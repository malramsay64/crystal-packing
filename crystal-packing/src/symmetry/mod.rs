//
// mod.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! The 17 plane (wallpaper) groups and their Wyckoff positions.
//!
//! The groups are defined by the coordinates of their general positions and the representative
//! coordinates of each Wyckoff position, taken from the International Tables for Crystallography
//! Volume A. Everything else (orbits, site symmetry, multiplicity, orientation constraints) is
//! derived from these definitions, so the only data which needs to be correct is what is in the
//! tables, and the tests check the derived properties match the tabulated multiplicities.

use std::f64::consts::PI;
use std::fmt;
use std::sync::OnceLock;

use anyhow::{anyhow, bail, Error};
use serde::{Deserialize, Serialize};

use crate::geometry::{angle_distance, wrap_angle, Mat2, Vec2};
use crate::hard_shape::ShapeSymmetry;

mod data;

/// The lattice systems of 2D crystals, which define the constraints on the unit cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LatticeFamily {
    /// No constraints on the cell
    Oblique,
    /// The cell angle is fixed at 90°
    Rectangular,
    /// The cell angle is fixed at 90° and both sides are the same length
    Square,
    /// The cell angle is fixed at 120° and both sides are the same length
    Hexagonal,
}

impl LatticeFamily {
    /// The cell matrix with unit length sides, used to convert symmetry operations from fractional
    /// to Cartesian coordinates.
    pub fn reference_cell(&self) -> Mat2 {
        match self {
            LatticeFamily::Hexagonal => Mat2::new(1., -0.5, 0., 3f64.sqrt() / 2.),
            _ => Mat2::IDENTITY,
        }
    }
}

/// An affine transformation of fractional coordinates, x' = R x + t
///
/// This represents both symmetry operations of a group, and the coordinates of a Wyckoff
/// position as a function of the free parameters (x, y). For the latter the matrix can be singular,
/// for example the position (x, 0).
///
/// ```
/// use crystal_packing::symmetry::SymOp;
/// use crystal_packing::geometry::Vec2;
/// let op = SymOp::parse("-y+1/2, x-y").unwrap();
/// let p = op.apply(Vec2::new(0.25, 0.5));
/// assert_eq!(p, Vec2::new(0., -0.25));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SymOp {
    pub rotation: [[i32; 2]; 2],
    pub translation: [f64; 2],
}

impl SymOp {
    pub const IDENTITY: SymOp = SymOp {
        rotation: [[1, 0], [0, 1]],
        translation: [0., 0.],
    };

    fn parse_component(s: &str) -> Result<([i32; 2], f64), Error> {
        let mut coeffs = [0i32; 2];
        let mut constant = 0.;
        let s: String = s.chars().filter(|c| !c.is_whitespace()).collect();
        if s.is_empty() {
            bail!("Empty coordinate in symmetry operation");
        }
        // Split into signed terms
        let mut terms: Vec<String> = vec![];
        let mut current = String::new();
        for c in s.chars() {
            if (c == '+' || c == '-') && !current.is_empty() {
                terms.push(current.clone());
                current.clear();
            }
            current.push(c);
        }
        terms.push(current);

        for term in terms {
            let (sign, body) = match term.strip_prefix('-') {
                Some(rest) => (-1, rest),
                None => (1, term.strip_prefix('+').unwrap_or(&term)),
            };
            let (multiplier, variable) = match body.chars().last() {
                Some(v @ ('x' | 'y')) => (&body[..body.len() - 1], Some(v)),
                _ => (body, None),
            };
            match variable {
                Some(v) => {
                    let m: i32 = if multiplier.is_empty() {
                        1
                    } else {
                        multiplier
                            .parse()
                            .map_err(|_| anyhow!("Invalid coefficient '{}'", multiplier))?
                    };
                    coeffs[if v == 'x' { 0 } else { 1 }] += sign * m;
                }
                None => {
                    let value = match multiplier.split_once('/') {
                        Some((n, d)) => {
                            let n: f64 = n.parse().map_err(|_| anyhow!("Invalid value '{}'", n))?;
                            let d: f64 = d.parse().map_err(|_| anyhow!("Invalid value '{}'", d))?;
                            n / d
                        }
                        None => multiplier
                            .parse()
                            .map_err(|_| anyhow!("Invalid value '{}'", multiplier))?,
                    };
                    constant += sign as f64 * value;
                }
            }
        }
        Ok((coeffs, constant))
    }

    /// Parse an operation in the notation of the International Tables, e.g. "-x+1/2,y"
    pub fn parse(s: &str) -> Result<SymOp, Error> {
        let s = s.trim().trim_matches(|c| c == '(' || c == ')');
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 2 {
            bail!("Expected two coordinates in '{}'", s);
        }
        let (rx, tx) = SymOp::parse_component(parts[0])?;
        let (ry, ty) = SymOp::parse_component(parts[1])?;
        Ok(SymOp {
            rotation: [rx, ry],
            translation: [tx, ty],
        })
    }

    pub fn apply(&self, p: Vec2) -> Vec2 {
        let r = &self.rotation;
        Vec2::new(
            r[0][0] as f64 * p.x + r[0][1] as f64 * p.y + self.translation[0],
            r[1][0] as f64 * p.x + r[1][1] as f64 * p.y + self.translation[1],
        )
    }

    /// The composition self ∘ other, applying other first
    pub fn compose(&self, other: &SymOp) -> SymOp {
        let a = &self.rotation;
        let b = &other.rotation;
        let rotation = [
            [
                a[0][0] * b[0][0] + a[0][1] * b[1][0],
                a[0][0] * b[0][1] + a[0][1] * b[1][1],
            ],
            [
                a[1][0] * b[0][0] + a[1][1] * b[1][0],
                a[1][0] * b[0][1] + a[1][1] * b[1][1],
            ],
        ];
        let t = self.apply(Vec2::new(other.translation[0], other.translation[1]));
        SymOp {
            rotation,
            translation: [t.x, t.y],
        }
    }

    pub fn linear(&self) -> Mat2 {
        let r = &self.rotation;
        Mat2::new(
            r[0][0] as f64,
            r[0][1] as f64,
            r[1][0] as f64,
            r[1][1] as f64,
        )
    }

    /// Which of the free parameters (x, y) affect the result
    pub fn uses(&self) -> [bool; 2] {
        let r = &self.rotation;
        [r[0][0] != 0 || r[1][0] != 0, r[0][1] != 0 || r[1][1] != 0]
    }
}

impl fmt::Display for SymOp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let component = |coeffs: [i32; 2], t: f64| {
            let mut s = String::new();
            for (c, v) in coeffs.iter().zip(["x", "y"]) {
                match c {
                    0 => {}
                    1 if s.is_empty() => s.push_str(v),
                    1 => s.push_str(&format!("+{}", v)),
                    -1 => s.push_str(&format!("-{}", v)),
                    c if *c > 0 && !s.is_empty() => s.push_str(&format!("+{}{}", c, v)),
                    c => s.push_str(&format!("{}{}", c, v)),
                }
            }
            if t.abs() > 1e-12 || s.is_empty() {
                if t >= 0. && !s.is_empty() {
                    s.push('+');
                }
                s.push_str(&format!("{}", t));
            }
            s
        };
        write!(
            f,
            "{},{}",
            component(self.rotation[0], self.translation[0]),
            component(self.rotation[1], self.translation[1])
        )
    }
}

/// Wrap fractional coordinates into the range [0, 1)
pub fn wrap_fractional(p: Vec2) -> Vec2 {
    Vec2::new(p.x - libm::floor(p.x), p.y - libm::floor(p.y))
}

/// Whether two fractional positions are the same point of the lattice
fn equivalent(a: Vec2, b: Vec2) -> bool {
    let d = a - b;
    let dx = d.x - libm::round(d.x);
    let dy = d.y - libm::round(d.y);
    dx.abs() < 1e-8 && dy.abs() < 1e-8
}

/// How the orientation of a shape on a site can vary
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum OrientationFreedom {
    /// Any orientation is allowed
    Free,
    /// Only a discrete set of orientations are allowed, due to mirror planes of the site
    Discrete(Vec<f64>),
}

/// A Wyckoff position of a plane group
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WyckoffPosition {
    pub letter: char,
    pub multiplicity: usize,
    /// The coordinates of the representative position as a function of the parameters (x, y)
    pub coordinates: SymOp,
    /// The operations generating every image of the site; image i is at `orbit[i](x, y)`.
    pub orbit: Vec<SymOp>,
    /// The Cartesian linear part of the group operation generating each image, which is applied
    /// to the orientation of the shape.
    pub orientations: Vec<Mat2>,
    /// The Cartesian linear parts of the operations leaving the site unchanged
    pub site_symmetry: Vec<Mat2>,
}

impl WyckoffPosition {
    /// Which of the parameters (x, y) are free to vary
    pub fn free_parameters(&self) -> [bool; 2] {
        self.coordinates.uses()
    }

    pub fn is_fixed(&self) -> bool {
        self.free_parameters() == [false, false]
    }

    /// The number of pure rotations (including the identity) in the site symmetry group
    pub fn rotation_order(&self) -> u32 {
        self.site_symmetry.iter().filter(|m| m.det() > 0.).count() as u32
    }

    /// The angles of the mirror lines through the site, in the range [0, π)
    pub fn mirror_angles(&self) -> Vec<f64> {
        self.site_symmetry
            .iter()
            .filter(|m| m.det() < 0.)
            .map(|m| wrap_angle(libm::atan2(m.m[1][0], m.m[0][0]) / 2., PI))
            .collect()
    }

    /// Determine whether a shape can occupy this site, and the orientations it can take.
    ///
    /// A shape on a special position needs to be invariant under the site symmetry. The rotations
    /// of the site must be a subgroup of the rotations of the shape, which is true for any
    /// orientation. Each mirror of the site needs to coincide with a mirror of the shape, which
    /// constrains the shape to a discrete set of orientations. Returns `None` when the shape can't
    /// occupy the site.
    pub fn orientation_freedom(&self, shape: &ShapeSymmetry) -> Option<OrientationFreedom> {
        if shape.continuous {
            return Some(OrientationFreedom::Discrete(vec![0.]));
        }
        if !shape.rotations.is_multiple_of(self.rotation_order()) {
            return None;
        }
        let site_mirrors = self.mirror_angles();
        let first = match site_mirrors.first() {
            None => return Some(OrientationFreedom::Free),
            Some(&m) => m,
        };
        // The orientation θ places the shape mirror β at θ + β, which needs to equal the site
        // mirror α modulo π.
        let period = 2. * PI / shape.rotations as f64;
        let mut allowed: Vec<f64> = vec![];
        for beta in shape.mirrors.iter() {
            for k in 0..2 {
                let theta = wrap_angle(first - beta + k as f64 * PI, 2. * PI);
                let valid = site_mirrors.iter().all(|alpha| {
                    shape.mirrors.iter().any(|b| {
                        let d = wrap_angle(theta + b - alpha, PI);
                        d < 1e-8 || PI - d < 1e-8
                    })
                });
                // Orientations related by a rotational symmetry of the shape are equivalent
                let duplicate = allowed.iter().any(|&a| {
                    angle_distance(wrap_angle(a, period), wrap_angle(theta, period), period) < 1e-8
                });
                if valid && !duplicate {
                    allowed.push(theta);
                }
            }
        }
        if allowed.is_empty() {
            None
        } else {
            allowed.sort_by(f64::total_cmp);
            Some(OrientationFreedom::Discrete(allowed))
        }
    }
}

/// One of the 17 plane groups
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlaneGroup {
    pub number: u8,
    pub name: String,
    pub family: LatticeFamily,
    /// All symmetry operations of the group modulo lattice translations, including centring.
    pub operations: Vec<SymOp>,
    pub wyckoff: Vec<WyckoffPosition>,
}

impl fmt::Display for PlaneGroup {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}

impl PlaneGroup {
    fn build(definition: &data::GroupDefinition) -> Result<PlaneGroup, Error> {
        let general: Vec<SymOp> = definition
            .operations
            .iter()
            .map(|s| SymOp::parse(s))
            .collect::<Result<_, _>>()?;
        let mut operations = general.clone();
        if definition.centred {
            let centring = SymOp {
                rotation: [[1, 0], [0, 1]],
                translation: [0.5, 0.5],
            };
            operations.extend(general.iter().map(|op| centring.compose(op)));
        }
        let cell = definition.family.reference_cell();
        let cell_inv = cell.inverse();
        let cartesian = |op: &SymOp| cell * op.linear() * cell_inv;

        // A generic point on a site, which is unlikely to coincide with any special position
        let generic = Vec2::new(0.1234567, 0.3141593);

        let mut wyckoff = vec![];
        for &(letter, multiplicity, coordinates) in definition.wyckoff.iter() {
            let coordinates = SymOp::parse(coordinates)?;
            let point = coordinates.apply(generic);
            let mut orbit: Vec<SymOp> = vec![];
            let mut orientations: Vec<Mat2> = vec![];
            let mut site_symmetry: Vec<Mat2> = vec![];
            for op in operations.iter() {
                let image = op.apply(point);
                if equivalent(image, point) {
                    site_symmetry.push(cartesian(op));
                }
                if !orbit
                    .iter()
                    .any(|o: &SymOp| equivalent(o.apply(generic), image))
                {
                    orbit.push(op.compose(&coordinates));
                    orientations.push(cartesian(op));
                }
            }
            if orbit.len() != multiplicity {
                bail!(
                    "Group {} site {} has multiplicity {} but {} images were generated",
                    definition.name,
                    letter,
                    multiplicity,
                    orbit.len()
                );
            }
            wyckoff.push(WyckoffPosition {
                letter,
                multiplicity,
                coordinates,
                orbit,
                orientations,
                site_symmetry,
            });
        }
        Ok(PlaneGroup {
            number: definition.number,
            name: definition.name.to_string(),
            family: definition.family,
            operations,
            wyckoff,
        })
    }

    /// All 17 plane groups, ordered by their number in the International Tables
    pub fn all() -> &'static [PlaneGroup] {
        static GROUPS: OnceLock<Vec<PlaneGroup>> = OnceLock::new();
        GROUPS.get_or_init(|| {
            data::GROUPS
                .iter()
                .map(|d| PlaneGroup::build(d).expect("Plane group data is invalid"))
                .collect()
        })
    }

    /// Find a plane group from its name, accepting both full and short symbols (e.g. p2mm, pmm)
    pub fn from_name(name: &str) -> Result<&'static PlaneGroup, Error> {
        let normalised = name.trim().to_lowercase();
        let lookup = data::ALIASES
            .iter()
            .find(|(alias, _)| *alias == normalised)
            .map(|(_, canonical)| *canonical)
            .unwrap_or(&normalised);
        PlaneGroup::all()
            .iter()
            .find(|g| g.name == lookup)
            .ok_or_else(|| anyhow!("Unknown plane group '{}'", name))
    }

    pub fn wyckoff_position(&self, letter: char) -> Result<&WyckoffPosition, Error> {
        self.wyckoff
            .iter()
            .find(|w| w.letter == letter)
            .ok_or_else(|| anyhow!("Group {} has no Wyckoff position {}", self.name, letter))
    }

    /// The number of symmetry operations, which is the multiplicity of the general position
    pub fn order(&self) -> usize {
        self.operations.len()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn parse_operations() {
        let op = SymOp::parse("x,y").unwrap();
        assert_eq!(op, SymOp::IDENTITY);
        let op = SymOp::parse("-x+1/2, y+1/2").unwrap();
        assert_eq!(op.rotation, [[-1, 0], [0, 1]]);
        assert_eq!(op.translation, [0.5, 0.5]);
        let op = SymOp::parse("x-y,-y").unwrap();
        assert_eq!(op.rotation, [[1, -1], [0, -1]]);
        let op = SymOp::parse("(2/3, 1/3)").unwrap();
        assert_eq!(op.rotation, [[0, 0], [0, 0]]);
        assert_abs_diff_eq!(op.translation[0], 2. / 3.);
        let op = SymOp::parse("x, x+1/2").unwrap();
        assert_eq!(op.rotation, [[1, 0], [1, 0]]);
        assert_eq!(op.uses(), [true, false]);
        assert!(SymOp::parse("x").is_err());
        assert!(SymOp::parse("x,z").is_err());
    }

    #[test]
    fn display_round_trip() {
        for s in ["x,y", "-y,x-y", "-x+0.5,y+0.5", "x,0", "0.5,0"] {
            let op = SymOp::parse(s).unwrap();
            assert_eq!(SymOp::parse(&op.to_string()).unwrap(), op);
        }
    }

    #[test]
    fn seventeen_groups() {
        let groups = PlaneGroup::all();
        assert_eq!(groups.len(), 17);
        for (i, g) in groups.iter().enumerate() {
            assert_eq!(g.number as usize, i + 1);
        }
    }

    #[test]
    fn group_orders() {
        let expected = [1, 2, 2, 2, 4, 4, 4, 4, 8, 4, 8, 8, 3, 6, 6, 6, 12];
        for (g, &order) in PlaneGroup::all().iter().zip(expected.iter()) {
            assert_eq!(g.order(), order, "{}", g.name);
        }
    }

    #[test]
    fn closure_under_composition() {
        // The product of any two operations is another operation, modulo lattice translations
        for g in PlaneGroup::all() {
            for a in g.operations.iter() {
                for b in g.operations.iter() {
                    let c = a.compose(b);
                    let found = g.operations.iter().any(|o| {
                        o.rotation == c.rotation
                            && equivalent(
                                Vec2::new(o.translation[0], o.translation[1]),
                                Vec2::new(c.translation[0], c.translation[1]),
                            )
                    });
                    assert!(found, "{}: {} ∘ {} = {} not in group", g.name, a, b, c);
                }
            }
        }
    }

    #[test]
    fn cartesian_operations_are_orthogonal() {
        for g in PlaneGroup::all() {
            for w in g.wyckoff.iter() {
                for m in w.orientations.iter().chain(w.site_symmetry.iter()) {
                    let product = *m * m.transpose();
                    assert!(
                        product.max_abs_diff(&Mat2::IDENTITY) < 1e-12,
                        "{} {}",
                        g.name,
                        w.letter
                    );
                }
            }
        }
    }

    #[test]
    fn multiplicity_matches_site_symmetry() {
        // |G| = multiplicity × |site symmetry group|
        for g in PlaneGroup::all() {
            for w in g.wyckoff.iter() {
                assert_eq!(
                    w.multiplicity * w.site_symmetry.len(),
                    g.order(),
                    "{} {}",
                    g.name,
                    w.letter
                );
            }
        }
    }

    #[test]
    fn general_position_first() {
        // The first listed position has the highest multiplicity and is the general position
        for g in PlaneGroup::all() {
            let general = &g.wyckoff[0];
            assert_eq!(general.multiplicity, g.order(), "{}", g.name);
            assert_eq!(general.free_parameters(), [true, true]);
            assert_eq!(general.site_symmetry.len(), 1);
        }
    }

    #[test]
    fn special_positions_invariant() {
        // Applying a site symmetry operation to any point on a site leaves it on the site
        for g in PlaneGroup::all() {
            for w in g.wyckoff.iter() {
                for params in [Vec2::new(0.2, 0.7), Vec2::new(0.9, 0.05)] {
                    let p = w.coordinates.apply(params);
                    let images: Vec<Vec2> = w.orbit.iter().map(|o| o.apply(params)).collect();
                    // Every group operation maps the orbit onto itself
                    for op in g.operations.iter() {
                        for image in images.iter() {
                            let q = op.apply(*image);
                            assert!(
                                images.iter().any(|i| equivalent(*i, q)),
                                "{} {}: orbit not closed",
                                g.name,
                                w.letter
                            );
                        }
                    }
                    assert!(equivalent(images[0], p));
                }
            }
        }
    }

    #[test]
    fn wyckoff_counts() {
        let expected = [1, 5, 3, 1, 2, 9, 4, 3, 6, 4, 7, 4, 4, 5, 4, 4, 6];
        for (g, &count) in PlaneGroup::all().iter().zip(expected.iter()) {
            assert_eq!(g.wyckoff.len(), count, "{}", g.name);
        }
    }

    #[test]
    fn lookup_names() {
        assert_eq!(PlaneGroup::from_name("p2mm").unwrap().number, 6);
        assert_eq!(PlaneGroup::from_name("pmm").unwrap().number, 6);
        assert_eq!(PlaneGroup::from_name("p1m1").unwrap().number, 3);
        assert_eq!(PlaneGroup::from_name("P6MM").unwrap().number, 17);
        assert!(PlaneGroup::from_name("p7").is_err());
    }

    #[test]
    fn orientation_square_on_p4mm() {
        let square = crate::hard_shape::ShapeSymmetry::detect(
            &(0..4)
                .map(|i| (Vec2::from_angle(PI / 2. * i as f64), 0.))
                .collect::<Vec<_>>(),
        );
        let group = PlaneGroup::from_name("p4mm").unwrap();
        // The 1a site has 4mm symmetry, a square can sit on it either aligned with the axes or
        // rotated by 45°, both of which keep the mirrors aligned.
        let site = group.wyckoff_position('a').unwrap();
        assert_eq!(site.rotation_order(), 4);
        match site.orientation_freedom(&square) {
            Some(OrientationFreedom::Discrete(angles)) => {
                assert_eq!(angles.len(), 2);
                assert_abs_diff_eq!(angles[1] - angles[0], PI / 4., epsilon = 1e-10);
            }
            other => panic!("Unexpected {:?}", other),
        }
        // The general position allows any orientation
        assert_eq!(
            group.wyckoff[0].orientation_freedom(&square),
            Some(OrientationFreedom::Free)
        );
    }

    #[test]
    fn orientation_incompatible() {
        let triangle = crate::hard_shape::ShapeSymmetry::detect(
            &(0..3)
                .map(|i| (Vec2::from_angle(2. * PI / 3. * i as f64), 0.))
                .collect::<Vec<_>>(),
        );
        let group = PlaneGroup::from_name("p4").unwrap();
        // A triangle doesn't have 4 fold or 2 fold symmetry
        assert_eq!(
            group
                .wyckoff_position('a')
                .unwrap()
                .orientation_freedom(&triangle),
            None
        );
        assert_eq!(
            group
                .wyckoff_position('c')
                .unwrap()
                .orientation_freedom(&triangle),
            None
        );
        let group = PlaneGroup::from_name("p3m1").unwrap();
        // A triangle fits on the 3m site, with two orientations (pointing either way)
        match group
            .wyckoff_position('a')
            .unwrap()
            .orientation_freedom(&triangle)
        {
            Some(OrientationFreedom::Discrete(angles)) => assert_eq!(angles.len(), 2),
            other => panic!("Unexpected {:?}", other),
        }
    }
}
