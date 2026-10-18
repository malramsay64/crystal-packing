//
// molecule.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use anyhow::{bail, Error};
use serde::{Deserialize, Serialize};

use super::{HardShape, Outline, ShapeSymmetry, OVERLAP_TOLERANCE};
use crate::geometry::{Mat2, Vec2};

/// A shape made from the union of disks
///
/// ```
/// use crystal_packing::hard_shape::{HardShape, Molecule};
/// let disk = Molecule::disk();
/// assert!((disk.area() - std::f64::consts::PI).abs() < 1e-12);
/// assert!(disk.symmetry().continuous);
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    name: String,
    /// The position and radius of each disk
    atoms: Vec<(Vec2, f64)>,
    area: f64,
    circumradius: f64,
    inradius: f64,
    symmetry: ShapeSymmetry,
}

/// A disk of unit radius
pub type Disk = Molecule;

/// The area of the union of a set of circles.
///
/// This uses Green's theorem, integrating ½(x dy - y dx) around the boundary of the union, which
/// is made of the arcs of each circle not inside any other circle. It is exact for any number of
/// overlapping circles.
fn union_area(circles: &[(Vec2, f64)]) -> f64 {
    let mut area = 0.;
    for (i, &(c, r)) in circles.iter().enumerate() {
        // A circle inside (or identical to, and listed after) another doesn't contribute
        let covered = circles.iter().enumerate().any(|(j, &(o, s))| {
            let d = (o - c).norm();
            j != i && (d + r < s || (d < 1e-12 && (r - s).abs() < 1e-12 && j < i))
        });
        if covered {
            continue;
        }
        // The angles where other circles cross this one split it into arcs
        let mut angles = vec![0., 2. * PI];
        for (j, &(o, s)) in circles.iter().enumerate() {
            let d = (o - c).norm();
            if j == i || d >= r + s || d + s <= r || d < 1e-12 {
                continue;
            }
            let base = (o - c).angle();
            let half = libm::acos(((d * d + r * r - s * s) / (2. * d * r)).clamp(-1., 1.));
            for angle in [base - half, base + half] {
                angles.push(crate::geometry::wrap_angle(angle, 2. * PI));
            }
        }
        angles.sort_by(f64::total_cmp);
        for pair in angles.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            if b - a < 1e-15 {
                continue;
            }
            // Keep the arc when its midpoint is outside every other circle
            let mid = c + Vec2::from_angle(0.5 * (a + b)) * r;
            let inside = circles
                .iter()
                .enumerate()
                .any(|(j, &(o, s))| j != i && (mid - o).norm_sq() < s * s);
            if !inside {
                let (sa, ca) = libm::sincos(a);
                let (sb, cb) = libm::sincos(b);
                area += 0.5 * (r * r * (b - a) + r * c.x * (sb - sa) - r * c.y * (cb - ca));
            }
        }
    }
    area
}

impl Molecule {
    pub fn new(name: &str, atoms: Vec<(Vec2, f64)>) -> Result<Molecule, Error> {
        if atoms.is_empty() {
            bail!("A molecule requires at least one atom");
        }
        if atoms.iter().any(|(_, r)| *r <= 0.) {
            bail!("Atom radii must be positive");
        }
        let circumradius = atoms.iter().map(|(p, r)| p.norm() + r).fold(0., f64::max);
        let inradius = atoms.iter().map(|(p, r)| r - p.norm()).fold(0., f64::max);
        let symmetry = if atoms.len() == 1 && atoms[0].0.norm() < 1e-12 {
            ShapeSymmetry::continuous()
        } else {
            ShapeSymmetry::detect(&atoms)
        };
        let area = union_area(&atoms);
        Ok(Molecule {
            name: name.to_string(),
            atoms,
            area,
            circumradius,
            inradius,
            symmetry,
        })
    }

    /// A disk with unit radius
    pub fn disk() -> Molecule {
        Molecule::new("Disk", vec![(Vec2::ZERO, 1.)]).expect("Valid disk")
    }

    /// A molecule of three atoms, a central atom of unit radius with two outer atoms
    ///
    /// The angle is between the two outer atoms in degrees, with the distance from the central
    /// atom to the outer atoms. The origin is placed at the centroid of the atom positions.
    pub fn trimer(radius: f64, angle: f64, distance: f64) -> Result<Molecule, Error> {
        let half = angle.to_radians() / 2.;
        let (s, c) = libm::sincos(half);
        Molecule::new(
            "Trimer",
            vec![
                (Vec2::new(0., -2. / 3. * distance * c), 1.),
                (Vec2::new(-distance * s, 1. / 3. * distance * c), radius),
                (Vec2::new(distance * s, 1. / 3. * distance * c), radius),
            ],
        )
    }

    pub fn atoms(&self) -> &[(Vec2, f64)] {
        &self.atoms
    }
}

#[derive(Debug, Clone, Default)]
pub struct PlacedMolecule {
    atoms: Vec<(Vec2, f64)>,
    scale: f64,
}

impl PlacedMolecule {
    /// The position and radius of each atom
    pub fn atoms(&self) -> &[(Vec2, f64)] {
        &self.atoms
    }
}

impl HardShape for Molecule {
    type Placed = PlacedMolecule;

    fn name(&self) -> &str {
        &self.name
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn circumradius(&self) -> f64 {
        self.circumradius
    }

    fn inradius(&self) -> f64 {
        self.inradius
    }

    fn symmetry(&self) -> &ShapeSymmetry {
        &self.symmetry
    }

    fn place(&self, orientation: &Mat2, position: Vec2, out: &mut PlacedMolecule) {
        out.atoms.clear();
        out.atoms.extend(
            self.atoms
                .iter()
                .map(|&(p, r)| (*orientation * p + position, r)),
        );
        out.scale = self.circumradius;
    }

    fn overlaps(a: &PlacedMolecule, b: &PlacedMolecule, shift: Vec2) -> bool {
        let tol = 1. - OVERLAP_TOLERANCE;
        a.atoms.iter().any(|&(pa, ra)| {
            b.atoms.iter().any(|&(pb, rb)| {
                let contact = (ra + rb) * tol;
                (pb + shift - pa).norm_sq() < contact * contact
            })
        })
    }

    fn support(placed: &PlacedMolecule, direction: Vec2) -> f64 {
        placed
            .atoms
            .iter()
            .map(|&(c, r)| c.dot(direction) + r)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    fn farthest(placed: &PlacedMolecule, point: Vec2) -> f64 {
        placed
            .atoms
            .iter()
            .map(|&(c, r)| (c - point).norm() + r)
            .fold(0., f64::max)
    }

    fn hull_circles(placed: &PlacedMolecule, out: &mut Vec<(Vec2, f64)>) {
        out.extend_from_slice(&placed.atoms);
    }

    fn outlines(&self, orientation: &Mat2, position: Vec2) -> Vec<Outline> {
        self.atoms
            .iter()
            .map(|&(p, r)| Outline::Circle(*orientation * p + position, r))
            .collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    #[test]
    fn disk_contact() {
        let disk = Molecule::disk();
        let mut a = PlacedMolecule::default();
        disk.place(&Mat2::IDENTITY, Vec2::ZERO, &mut a);
        assert!(!Molecule::overlaps(&a, &a, Vec2::new(2., 0.)));
        assert!(Molecule::overlaps(&a, &a, Vec2::new(1.999, 0.)));
    }

    #[test]
    fn support_and_farthest() {
        let disk = Molecule::disk();
        let mut placed = PlacedMolecule::default();
        disk.place(&Mat2::IDENTITY, Vec2::new(3., 1.), &mut placed);
        assert_abs_diff_eq!(
            Molecule::support(&placed, Vec2::new(1., 0.)),
            4.,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            Molecule::support(&placed, Vec2::new(0., -1.)),
            0.,
            epsilon = 1e-12
        );
        let (lo, hi) = Molecule::bounds(&placed);
        assert_abs_diff_eq!(lo.x, 2., epsilon = 1e-12);
        assert_abs_diff_eq!(lo.y, 0., epsilon = 1e-12);
        assert_abs_diff_eq!(hi.x, 4., epsilon = 1e-12);
        assert_abs_diff_eq!(hi.y, 2., epsilon = 1e-12);
        assert_abs_diff_eq!(
            Molecule::farthest(&placed, Vec2::ZERO),
            10f64.sqrt() + 1.,
            epsilon = 1e-12
        );
    }

    #[test]
    fn trimer_symmetry() {
        let trimer = Molecule::trimer(0.637556, 120., 1.).unwrap();
        assert_eq!(trimer.symmetry().rotations, 1);
        assert_eq!(trimer.symmetry().mirrors.len(), 1);
        assert_abs_diff_eq!(trimer.symmetry().mirrors[0], PI / 2., epsilon = 1e-10);
    }

    /// Estimate the area of a union of circles by sampling a fine grid
    fn grid_area(atoms: &[(Vec2, f64)]) -> f64 {
        let n = 1500;
        let extent = 4.;
        let h = 2. * extent / n as f64;
        let mut inside = 0;
        for i in 0..n {
            for j in 0..n {
                let p = Vec2::new(
                    -extent + (i as f64 + 0.5) * h,
                    -extent + (j as f64 + 0.5) * h,
                );
                if atoms.iter().any(|(c, r)| (p - *c).norm_sq() < r * r) {
                    inside += 1;
                }
            }
        }
        inside as f64 * h * h
    }

    #[test]
    fn union_area_three_overlapping() {
        // Every pair overlaps and all three share a common region, which the sum of pairwise
        // overlaps would count incorrectly
        let trimer = Molecule::trimer(0.9, 60., 0.3).unwrap();
        let expected = grid_area(trimer.atoms());
        assert_abs_diff_eq!(trimer.area(), expected, epsilon = 2e-3 * expected);
    }

    #[test]
    fn union_area_contained() {
        let molecule = Molecule::new(
            "contained",
            vec![
                (Vec2::ZERO, 1.),
                (Vec2::new(0.2, 0.1), 0.5),
                (Vec2::new(1.5, 0.), 0.7),
            ],
        )
        .unwrap();
        let expected = grid_area(molecule.atoms());
        assert_abs_diff_eq!(molecule.area(), expected, epsilon = 2e-3 * expected);
    }

    #[test]
    fn union_area_identical() {
        let molecule = Molecule::new("twice", vec![(Vec2::ZERO, 1.), (Vec2::ZERO, 1.)]).unwrap();
        assert_abs_diff_eq!(molecule.area(), PI, epsilon = 1e-12);
    }

    #[test]
    fn overlap_area() {
        // Two unit circles a unit distance apart overlap by a lens of area 2π/3 - √3/2
        let union = union_area(&[(Vec2::ZERO, 1.), (Vec2::new(1., 0.), 1.)]);
        let lens = 2. * PI / 3. - 3f64.sqrt() / 2.;
        assert_abs_diff_eq!(union, 2. * PI - lens, epsilon = 1e-12);
    }

    #[test]
    fn trimer_contact() {
        // A trimer and its mirror image, with their small atoms touching side by side
        let trimer = Molecule::trimer(0.5, 120., 1.).unwrap();
        let (mut a, mut b) = (PlacedMolecule::default(), PlacedMolecule::default());
        trimer.place(&Mat2::IDENTITY, Vec2::ZERO, &mut a);
        trimer.place(&Mat2::reflection(PI / 2.), Vec2::ZERO, &mut b);
        let reach = 2. * Molecule::support(&a, Vec2::new(1., 0.));
        assert!(!Molecule::overlaps(&a, &b, Vec2::new(reach, 0.)));
        assert!(Molecule::overlaps(&a, &b, Vec2::new(reach - 1e-6, 0.)));
        // Moved up or down they no longer touch
        assert!(!Molecule::overlaps(&a, &b, Vec2::new(reach - 1e-6, 0.1)));
    }

    proptest! {
        #[test]
        fn overlaps_match_brute_force(
            radius in 0.2f64..1.,
            angle in 60f64..180.,
            distance in 0.3f64..1.5,
            angle_a in 0. ..(2. * PI),
            angle_b in 0. ..(2. * PI),
            reflect in any::<bool>(),
            x in -3f64..3.,
            y in -3f64..3.,
        ) {
            let trimer = Molecule::trimer(radius, angle, distance);
            prop_assume!(trimer.is_ok());
            let trimer = trimer.unwrap();
            let a = Mat2::rotation(angle_a);
            let b = if reflect {
                Mat2::reflection(0.4) * Mat2::rotation(angle_b)
            } else {
                Mat2::rotation(angle_b)
            };
            let shift = Vec2::new(x, y);
            let (mut pa, mut pb) = (PlacedMolecule::default(), PlacedMolecule::default());
            trimer.place(&a, Vec2::ZERO, &mut pa);
            trimer.place(&b, Vec2::ZERO, &mut pb);
            let expected =
                crate::validate::brute_force_overlap(&trimer, (&a, Vec2::ZERO), (&b, shift));
            prop_assert_eq!(Molecule::overlaps(&pa, &pb, shift), expected);
            prop_assert_eq!(Molecule::overlaps(&pb, &pa, -shift), expected);
            // Shapes closer than the sum of their inradii always overlap, which the engine uses
            // to skip the full check, and further apart than their circumradii never do
            if shift.norm() < 2. * trimer.inradius() {
                prop_assert!(expected);
            }
            if shift.norm() >= 2. * trimer.circumradius() {
                prop_assert!(!expected);
            }
        }
    }
}
