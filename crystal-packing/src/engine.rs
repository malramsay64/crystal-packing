//
// engine.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! The state of a periodic packing of hard shapes and the overlap checks on that state.
//!
//! A [`Packing`] holds a plane group, a set of occupied Wyckoff sites, the lattice, and the free
//! parameters of each site. Every image of every site is placed into world coordinates and cached,
//! with the cache updated only for the parts of the structure which change.
//!
//! Overlaps are found using the symmetry of the structure. Every pair of shapes in the crystal is
//! related by a symmetry operation to a pair which includes the representative (first) image of
//! a site, so it is sufficient to check each representative against every image within the
//! contact distance. When a single site changes, only the pairs involving that site are checked.

use anyhow::{anyhow, bail, Error};
use serde::{Deserialize, Serialize};

use crate::geometry::{Mat2, Vec2};
use crate::hard_shape::HardShape;
use crate::lattice::Lattice;
use crate::random::Rng;
use crate::symmetry::{LatticeFamily, OrientationFreedom, PlaneGroup, WyckoffPosition};

/// The variable parameters of an occupied site
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SiteState {
    /// The free parameters (x, y) of the Wyckoff position in fractional coordinates
    pub position: Vec2,
    /// The orientation of the shape on the representative image, in radians
    pub angle: f64,
}

/// A Wyckoff position which is occupied by a shape, along with the constraints on the shape
#[derive(Debug, Clone)]
pub struct OccupiedSite {
    pub wyckoff: &'static WyckoffPosition,
    pub free: [bool; 2],
    pub orientation: OrientationFreedom,
}

#[derive(Debug, Clone, Default)]
struct Image<P> {
    centre: Vec2,
    placed: P,
}

/// A serialisable description of a packed structure
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Structure {
    pub shape: String,
    pub group: String,
    pub sites: Vec<SiteRecord>,
    pub lattice: Lattice,
    pub num_shapes: usize,
    pub packing_fraction: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SiteRecord {
    pub letter: char,
    pub x: f64,
    pub y: f64,
    pub angle: f64,
}

impl Structure {
    /// The Wyckoff letters of the occupied sites, e.g. "ab"
    pub fn sites_label(&self) -> String {
        self.sites.iter().map(|s| s.letter).collect()
    }
}

#[derive(Debug, Clone)]
pub struct Packing<S: HardShape> {
    shape: S,
    group: &'static PlaneGroup,
    sites: Vec<OccupiedSite>,
    lattice: Lattice,
    state: Vec<SiteState>,
    images: Vec<Vec<Image<S::Placed>>>,
    /// The reduced basis of the lattice, used to enumerate periodic images
    reduced: Mat2,
    reduced_inverse: Mat2,
    /// The range of the coefficients of each reduced basis vector within the contact distance
    range: Vec2,
    /// Images further apart than this distance can't overlap
    contact: f64,
    /// Images closer than this distance must overlap
    core: f64,
    num_shapes: usize,
}

impl<S: HardShape> Packing<S> {
    /// Create a packing with the shape occupying the given Wyckoff sites of the group.
    ///
    /// The sites are given by their Wyckoff letters. The parameters of each site are chosen at
    /// random and the lattice is made large enough that the shapes don't overlap.
    pub fn new(
        shape: S,
        group: &'static PlaneGroup,
        letters: &[char],
        rng: &mut Rng,
    ) -> Result<Packing<S>, Error> {
        if letters.is_empty() {
            bail!("At least one site needs to be occupied");
        }
        let mut sites = vec![];
        for &letter in letters {
            let wyckoff = group.wyckoff_position(letter)?;
            let orientation = wyckoff
                .orientation_freedom(shape.symmetry())
                .ok_or_else(|| {
                    anyhow!(
                        "The shape doesn't have the symmetry required for site {} of {}",
                        letter,
                        group.name
                    )
                })?;
            sites.push(OccupiedSite {
                wyckoff,
                free: wyckoff.free_parameters(),
                orientation,
            });
        }
        for (i, a) in sites.iter().enumerate() {
            for b in sites.iter().skip(i + 1) {
                if a.wyckoff.is_fixed() && a.wyckoff.letter == b.wyckoff.letter {
                    bail!(
                        "The fixed site {} can only be occupied once",
                        a.wyckoff.letter
                    );
                }
            }
        }
        let state: Vec<SiteState> = sites
            .iter()
            .map(|site| SiteState {
                position: Vec2::new(rng.uniform(), rng.uniform()),
                angle: match &site.orientation {
                    OrientationFreedom::Free => rng.uniform() * 2. * std::f64::consts::PI,
                    OrientationFreedom::Discrete(angles) => angles[rng.index(angles.len())],
                },
            })
            .collect();

        let num_shapes: usize = sites.iter().map(|s| s.wyckoff.multiplicity).sum();
        // Start from a dilute structure, expanding until there are no overlaps
        let target_area = num_shapes as f64 * shape.area() / 0.1;
        let lattice = Lattice::with_length(group.family, 1.);
        let lattice = lattice.scaled((target_area / lattice.area()).sqrt());
        let mut packing = Packing::build(shape, group, sites, lattice, state);
        for _ in 0..200 {
            if !packing.any_overlap() {
                return Ok(packing);
            }
            let expanded = packing.lattice.scaled(1.5);
            packing.set_lattice(expanded);
        }
        bail!("Unable to find a starting configuration without overlaps")
    }

    fn build(
        shape: S,
        group: &'static PlaneGroup,
        sites: Vec<OccupiedSite>,
        lattice: Lattice,
        state: Vec<SiteState>,
    ) -> Packing<S> {
        let images = sites
            .iter()
            .map(|s| vec![Image::default(); s.wyckoff.multiplicity])
            .collect();
        let num_shapes = sites.iter().map(|s| s.wyckoff.multiplicity).sum();
        let contact = 2. * shape.circumradius();
        // Reduced slightly so shapes in contact aren't considered overlapping
        let core = 2. * shape.inradius() * (1. - 1e-8);
        let mut packing = Packing {
            shape,
            group,
            sites,
            lattice: lattice.clone(),
            state,
            images,
            reduced: Mat2::IDENTITY,
            reduced_inverse: Mat2::IDENTITY,
            range: Vec2::ZERO,
            contact,
            core,
            num_shapes,
        };
        packing.set_lattice(lattice);
        packing
    }

    /// Recreate a packing from a saved structure
    pub fn from_structure(shape: S, structure: &Structure) -> Result<Packing<S>, Error> {
        let group = PlaneGroup::from_name(&structure.group)?;
        let mut sites = vec![];
        let mut state = vec![];
        for record in structure.sites.iter() {
            let wyckoff = group.wyckoff_position(record.letter)?;
            let orientation = wyckoff
                .orientation_freedom(shape.symmetry())
                .ok_or_else(|| anyhow!("Shape incompatible with site {}", record.letter))?;
            sites.push(OccupiedSite {
                wyckoff,
                free: wyckoff.free_parameters(),
                orientation,
            });
            state.push(SiteState {
                position: Vec2::new(record.x, record.y),
                angle: record.angle,
            });
        }
        if structure.lattice.family != group.family {
            bail!(
                "The lattice doesn't match the family of group {}",
                group.name
            );
        }
        Ok(Packing::build(
            shape,
            group,
            sites,
            structure.lattice.clone(),
            state,
        ))
    }

    pub fn structure(&self) -> Structure {
        Structure {
            shape: self.shape.name().to_string(),
            group: self.group.name.clone(),
            sites: self
                .sites
                .iter()
                .zip(self.state.iter())
                .map(|(site, state)| SiteRecord {
                    letter: site.wyckoff.letter,
                    x: state.position.x,
                    y: state.position.y,
                    angle: state.angle,
                })
                .collect(),
            lattice: self.lattice.clone(),
            num_shapes: self.num_shapes,
            packing_fraction: self.packing_fraction(),
        }
    }

    pub fn shape(&self) -> &S {
        &self.shape
    }

    pub fn group(&self) -> &'static PlaneGroup {
        self.group
    }

    pub fn sites(&self) -> &[OccupiedSite] {
        &self.sites
    }

    pub fn lattice(&self) -> &Lattice {
        &self.lattice
    }

    pub fn site_state(&self, index: usize) -> SiteState {
        self.state[index]
    }

    pub fn num_shapes(&self) -> usize {
        self.num_shapes
    }

    pub fn packing_fraction(&self) -> f64 {
        self.num_shapes as f64 * self.shape.area() / self.lattice.area()
    }

    /// The fractional coordinates of an image, wrapped into the unit cell
    fn image_fractional(&self, site: usize, image: usize) -> Vec2 {
        self.sites[site].wyckoff.orbit[image].apply(self.state[site].position)
    }

    fn update_site_images(&mut self, index: usize) {
        let site = &self.sites[index];
        let state = self.state[index];
        let rotation = Mat2::rotation(state.angle);
        for (i, image) in self.images[index].iter_mut().enumerate() {
            let fractional = site.wyckoff.orbit[i].apply(state.position);
            let centre = self.lattice.to_cartesian(fractional);
            let orientation = site.wyckoff.orientations[i] * rotation;
            self.shape.place(&orientation, centre, &mut image.placed);
            image.centre = centre;
        }
    }

    /// Update the parameters of a site, updating the cached images
    pub fn set_site(&mut self, index: usize, state: SiteState) {
        self.state[index] = state;
        self.update_site_images(index);
    }

    /// Update the lattice, updating all cached images
    pub fn set_lattice(&mut self, lattice: Lattice) {
        self.lattice = lattice;
        self.reduced = self.lattice.reduced_basis();
        self.reduced_inverse = self.reduced.inverse();
        // A vector t = n u + m v of length at most c has |m| ≤ c |u| / |u × v| (from
        // |u × t| = |m| |u × v| ≤ |u| |t|), and similarly for n.
        let (u, v) = (self.reduced.column(0), self.reduced.column(1));
        let det = self.reduced.det().abs();
        self.range = Vec2::new(self.contact * v.norm() / det, self.contact * u.norm() / det);
        for index in 0..self.sites.len() {
            self.update_site_images(index);
        }
    }

    /// Remove the shear of an oblique cell by changing the basis to b' = b - k a.
    ///
    /// This describes exactly the same structure, so the overlaps and packing fraction are
    /// unchanged, but keeps the cell parameters from drifting without bound during the
    /// optimisation. The fractional coordinates of sites transform as x' = x + k y. Fixed sites of
    /// p2 lie on half lattice points, which are only preserved when k is even.
    pub fn reduce_shear(&mut self) {
        if self.lattice.family != LatticeFamily::Oblique {
            return;
        }
        let p = self.lattice.parameters();
        let any_fixed = self.sites.iter().any(|s| s.wyckoff.is_fixed());
        let k = if any_fixed {
            2. * libm::round(p[1] / (2. * p[0]))
        } else {
            libm::round(p[1] / p[0])
        };
        if k == 0. {
            return;
        }
        for (site, state) in self.sites.iter().zip(self.state.iter_mut()) {
            if !site.wyckoff.is_fixed() {
                let x = state.position.x + k * state.position.y;
                state.position.x = x - libm::floor(x);
            }
        }
        self.set_lattice(Lattice::from_parameters(
            LatticeFamily::Oblique,
            &[p[0], p[1] - k * p[0], p[2]],
        ));
    }

    /// Replace an oblique cell with the reduced basis of the same lattice, rotated so the first
    /// vector lies along the x axis.
    ///
    /// The whole structure is rotated with the cell, so the site orientations change by the same
    /// angle, which is returned. This is only possible when every site is on a general position;
    /// fixed sites lie on half lattice points which can be relabelled by the change of basis, so
    /// only the shear is reduced and None is returned.
    pub fn reduce_lattice(&mut self) -> Option<f64> {
        if self.lattice.family != LatticeFamily::Oblique {
            return None;
        }
        if self.sites.iter().any(|s| s.wyckoff.is_fixed()) {
            self.reduce_shear();
            return None;
        }
        let (u, mut v) = crate::lattice::gauss_reduce(self.lattice.a(), self.lattice.b());
        if u.cross(v) < 0. {
            v = -v;
        }
        // The change of basis as an integer matrix, rounded to remove numerical noise
        let mut change = self.lattice.matrix.inverse() * Mat2::from_columns(u, v);
        for value in change.m.iter_mut().flatten() {
            *value = libm::round(*value);
        }
        if change.max_abs_diff(&Mat2::IDENTITY) == 0. {
            return Some(0.);
        }
        let inverse = change.inverse();
        let angle = u.angle();
        let rotated = Mat2::rotation(-angle) * self.lattice.matrix * change;
        for (site, state) in self.sites.iter().zip(self.state.iter_mut()) {
            let p = inverse * state.position;
            state.position = crate::symmetry::wrap_fractional(p);
            if site.orientation == OrientationFreedom::Free {
                state.angle =
                    crate::geometry::wrap_angle(state.angle - angle, 2. * std::f64::consts::PI);
            }
        }
        self.set_lattice(Lattice::from_parameters(
            LatticeFamily::Oblique,
            &[rotated.m[0][0], rotated.m[0][1], rotated.m[1][1]],
        ));
        Some(angle)
    }

    /// Check whether image `b` or any of its periodic copies overlaps image `a`
    ///
    /// When `same` is true a and b are the same image, so the zero translation is skipped.
    #[inline]
    fn images_overlap(&self, a: &Image<S::Placed>, b: &Image<S::Placed>, same: bool) -> bool {
        // Find the lattice translations t, as coefficients of the reduced basis, which bring b
        // within the contact distance of a, i.e. |d + t| < contact.
        let d = b.centre - a.centre;
        let f = self.reduced_inverse * d;
        let (u, v) = (self.reduced.column(0), self.reduced.column(1));
        let contact_sq = self.contact * self.contact;
        let core_sq = self.core * self.core;
        let m_min = libm::ceil(-f.y - self.range.y) as i64;
        let m_max = libm::floor(-f.y + self.range.y) as i64;
        let n_min = libm::ceil(-f.x - self.range.x) as i64;
        let n_max = libm::floor(-f.x + self.range.x) as i64;
        for m in m_min..=m_max {
            for n in n_min..=n_max {
                if same && n == 0 && m == 0 {
                    continue;
                }
                let t = u * n as f64 + v * m as f64;
                let distance_sq = (d + t).norm_sq();
                if distance_sq >= contact_sq {
                    continue;
                }
                if distance_sq < core_sq || S::overlaps(&a.placed, &b.placed, t) {
                    return true;
                }
            }
        }
        false
    }

    /// Check for overlaps involving any image of the site
    pub fn site_overlaps(&self, index: usize) -> bool {
        let representative = &self.images[index][0];
        for (j, images) in self.images.iter().enumerate() {
            for (l, image) in images.iter().enumerate() {
                if self.images_overlap(representative, image, j == index && l == 0) {
                    return true;
                }
            }
        }
        for (j, images) in self.images.iter().enumerate() {
            if j == index {
                continue;
            }
            for image in self.images[index].iter() {
                if self.images_overlap(&images[0], image, false) {
                    return true;
                }
            }
        }
        false
    }

    /// Check for any overlaps in the structure
    pub fn any_overlap(&self) -> bool {
        for (i, site_images) in self.images.iter().enumerate() {
            let representative = &site_images[0];
            for (j, images) in self.images.iter().enumerate() {
                for (l, image) in images.iter().enumerate() {
                    if self.images_overlap(representative, image, i == j && l == 0) {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Every placed shape within the given range of unit cells, as (site, centre, orientation,
    /// placed shape), used for rendering
    pub fn placed_shapes(&self, cells: i32) -> Vec<(usize, Vec2, Mat2, S::Placed)> {
        let mut shapes = vec![];
        for (index, site) in self.sites.iter().enumerate() {
            let rotation = Mat2::rotation(self.state[index].angle);
            for image in 0..site.wyckoff.multiplicity {
                let fractional =
                    crate::symmetry::wrap_fractional(self.image_fractional(index, image));
                let orientation = site.wyckoff.orientations[image] * rotation;
                for n in -cells..=cells {
                    for m in -cells..=cells {
                        let f = fractional + Vec2::new(n as f64, m as f64);
                        let centre = self.lattice.to_cartesian(f);
                        let mut placed = S::Placed::default();
                        self.shape.place(&orientation, centre, &mut placed);
                        shapes.push((index, centre, orientation, placed));
                    }
                }
            }
        }
        shapes
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::hard_shape::{Molecule, Polygon};
    use crate::random::Rng;
    use crate::symmetry::LatticeFamily;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    /// A brute force overlap check comparing every image with every image in a large supercell
    fn naive_overlap<S: HardShape>(packing: &Packing<S>) -> bool {
        let shapes = packing.placed_shapes(0);
        let lattice = packing.lattice();
        let k = 15;
        for (i, (_, ca, _, a)) in shapes.iter().enumerate() {
            for (j, (_, cb, _, b)) in shapes.iter().enumerate() {
                for n in -k..=k {
                    for m in -k..=k {
                        if i == j && n == 0 && m == 0 {
                            continue;
                        }
                        let shift = lattice.to_cartesian(Vec2::new(n as f64, m as f64));
                        if (*cb + shift - *ca).norm() < packing.contact && S::overlaps(a, b, shift)
                        {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    #[test]
    fn square_lattice() {
        let square = Polygon::regular(4).unwrap();
        let group = PlaneGroup::from_name("p1").unwrap();
        let mut rng = Rng::seed_from_u64(0);
        let mut packing = Packing::new(square, group, &['a'], &mut rng).unwrap();
        packing.set_site(
            0,
            SiteState {
                position: Vec2::ZERO,
                angle: std::f64::consts::FRAC_PI_4,
            },
        );
        let side = 2f64.sqrt();
        packing.set_lattice(Lattice::from_parameters(
            LatticeFamily::Oblique,
            &[side, 0., side],
        ));
        assert!(!packing.any_overlap());
        assert_abs_diff_eq!(packing.packing_fraction(), 1., epsilon = 1e-12);
        packing.set_lattice(Lattice::from_parameters(
            LatticeFamily::Oblique,
            &[side * 0.999, 0., side],
        ));
        assert!(packing.any_overlap());
    }

    #[test]
    fn hexagonal_disks() {
        let disk = Molecule::disk();
        let group = PlaneGroup::from_name("p6mm").unwrap();
        let mut rng = Rng::seed_from_u64(0);
        let mut packing = Packing::new(disk, group, &['a'], &mut rng).unwrap();
        packing.set_lattice(Lattice::from_parameters(LatticeFamily::Hexagonal, &[2.]));
        assert!(!packing.any_overlap());
        assert_abs_diff_eq!(
            packing.packing_fraction(),
            std::f64::consts::PI / (2. * 3f64.sqrt()),
            epsilon = 1e-12
        );
        packing.set_lattice(Lattice::from_parameters(LatticeFamily::Hexagonal, &[1.999]));
        assert!(packing.any_overlap());
    }

    #[test]
    fn incompatible_site() {
        let triangle = Polygon::regular(3).unwrap();
        let group = PlaneGroup::from_name("p4").unwrap();
        let mut rng = Rng::seed_from_u64(0);
        assert!(Packing::new(triangle, group, &['a'], &mut rng).is_err());
    }

    #[test]
    fn structure_round_trip() {
        let shape = Polygon::regular(5).unwrap();
        let group = PlaneGroup::from_name("p2gg").unwrap();
        let mut rng = Rng::seed_from_u64(4);
        let packing = Packing::new(shape.clone(), group, &['c'], &mut rng).unwrap();
        let structure = packing.structure();
        let json = serde_json::to_string(&structure).unwrap();
        let restored: Structure = serde_json::from_str(&json).unwrap();
        let packing2 = Packing::from_structure(shape, &restored).unwrap();
        assert_eq!(packing2.structure(), structure);
    }

    fn random_packing(group: &str, letters: &[char], seed: u64, scale: f64) -> Packing<Polygon> {
        random_polygon_packing(5, group, letters, seed, scale)
    }

    fn random_polygon_packing(
        sides: usize,
        group: &str,
        letters: &[char],
        seed: u64,
        scale: f64,
    ) -> Packing<Polygon> {
        let shape = Polygon::regular(sides).unwrap();
        let group = PlaneGroup::from_name(group).unwrap();
        let mut rng = Rng::seed_from_u64(seed);
        let mut packing = Packing::new(shape, group, letters, &mut rng).unwrap();
        // Compress so that overlaps are likely
        let lattice = packing.lattice().scaled(scale);
        packing.set_lattice(lattice);
        packing
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(300))]

        #[test]
        fn symmetry_reduced_check_matches_naive(
            seed in 0u64..10_000,
            scale in 0.2f64..1.,
            group in prop::sample::select(vec!["p1", "p2", "p1g1", "p2mg", "p2gg", "p4", "p3", "p6", "c1m1"]),
        ) {
            let letters: Vec<char> = vec![PlaneGroup::from_name(group).unwrap().wyckoff[0].letter];
            let packing = random_packing(group, &letters, seed, scale);
            prop_assert_eq!(packing.any_overlap(), naive_overlap(&packing));
        }

        #[test]
        fn multiple_sites_match_naive(seed in 0u64..10_000, scale in 0.2f64..1.) {
            let packing = random_polygon_packing(4, "p2", &['e', 'a'], seed, scale);
            prop_assert_eq!(packing.any_overlap(), naive_overlap(&packing));
            let packing = random_polygon_packing(4, "p4mm", &['g', 'b'], seed, scale);
            prop_assert_eq!(packing.any_overlap(), naive_overlap(&packing));
            let packing = random_polygon_packing(6, "p6mm", &['d', 'b'], seed, scale);
            prop_assert_eq!(packing.any_overlap(), naive_overlap(&packing));
            let packing = random_packing("p1", &['a', 'a'], seed, scale);
            prop_assert_eq!(packing.any_overlap(), naive_overlap(&packing));
        }

        #[test]
        fn sheared_cell_matches_naive(seed in 0u64..10_000, shear in -6f64..6., scale in 0.3f64..1.2) {
            let mut packing = random_packing("p2", &['e'], seed, 1.);
            let p = packing.lattice().parameters();
            let lattice = Lattice::from_parameters(
                LatticeFamily::Oblique,
                &[p[0] * scale, p[1] + shear * p[0] * scale, p[2] * scale],
            );
            packing.set_lattice(lattice);
            prop_assert_eq!(packing.any_overlap(), naive_overlap(&packing));
        }

        #[test]
        fn shear_reduction_preserves_structure(
            seed in 0u64..10_000,
            shear in -8f64..8.,
            scale in 0.3f64..1.,
            fixed in proptest::bool::ANY,
        ) {
            let sides = if fixed { 4 } else { 5 };
            let letters: &[char] = if fixed { &['e', 'b'] } else { &['e', 'e'] };
            let mut packing = random_polygon_packing(sides, "p2", letters, seed, 1.);
            let p = packing.lattice().parameters();
            packing.set_lattice(Lattice::from_parameters(
                LatticeFamily::Oblique,
                &[p[0] * scale, p[1] + shear * p[0] * scale, p[2] * scale],
            ));
            let overlap = packing.any_overlap();
            let phi = packing.packing_fraction();
            let lattice = packing.lattice().clone();
            let before: Vec<Vec2> = packing.placed_shapes(0).iter().map(|s| s.1).collect();
            packing.reduce_shear();
            let p = packing.lattice().parameters();
            let limit = if fixed { 1. } else { 0.5 };
            prop_assert!(p[1].abs() <= p[0] * limit + 1e-9);
            prop_assert_eq!(packing.any_overlap(), overlap);
            prop_assert!((packing.packing_fraction() - phi).abs() < 1e-12);
            // Every shape is at a position equivalent to a shape before the reduction
            let inverse = lattice.matrix.inverse();
            for s in packing.placed_shapes(0) {
                let found = before.iter().any(|&b| {
                    let f = inverse * (s.1 - b);
                    (f.x - f.x.round()).abs() < 1e-8 && (f.y - f.y.round()).abs() < 1e-8
                });
                prop_assert!(found);
            }
        }

        #[test]
        fn lattice_reduction_preserves_structure(
            seed in 0u64..10_000,
            ax in 0.3f64..3.,
            bx in -8f64..8.,
            by in 0.1f64..3.,
        ) {
            let mut packing = random_polygon_packing(5, "p2", &['e', 'e'], seed, 1.);
            let area = packing.lattice().area();
            // A cell of the same area with an arbitrary, possibly very poor, shape
            let raw = Lattice::from_parameters(LatticeFamily::Oblique, &[ax, bx, by]);
            packing.set_lattice(raw.scaled((area / raw.area()).sqrt() * 0.6));
            let overlap = packing.any_overlap();
            let phi = packing.packing_fraction();
            let lattice = packing.lattice().clone();
            let before: Vec<(Vec2, Mat2)> = packing
                .placed_shapes(0)
                .iter()
                .map(|s| (s.1, s.2))
                .collect();
            let angle = packing.reduce_lattice().unwrap();
            let p = packing.lattice().parameters();
            prop_assert!(p[1].abs() <= 0.5 * p[0] + 1e-9);
            prop_assert!(p[0] <= packing.lattice().b().norm() + 1e-9);
            prop_assert_eq!(packing.any_overlap(), overlap);
            prop_assert!((packing.packing_fraction() - phi).abs() < 1e-9);
            // Rotating back, every shape matches a shape from before the reduction
            let back = Mat2::rotation(angle);
            let inverse = lattice.matrix.inverse();
            for s in packing.placed_shapes(0) {
                let centre = back * s.1;
                let orientation = back * s.2;
                let found = before.iter().any(|(c, o)| {
                    let f = inverse * (centre - *c);
                    (f.x - f.x.round()).abs() < 1e-8
                        && (f.y - f.y.round()).abs() < 1e-8
                        && orientation.max_abs_diff(o) < 1e-8
                });
                prop_assert!(found);
            }
        }

        #[test]
        fn incremental_matches_full(seed in 0u64..10_000, scale in 0.3f64..1.) {
            let mut packing = random_packing("p2", &['e', 'e'], seed, scale);
            let mut rng = Rng::seed_from_u64(seed);
            // Start from a valid configuration
            while packing.any_overlap() {
                let lattice = packing.lattice().scaled(1.2);
                packing.set_lattice(lattice);
            }
            for _ in 0..20 {
                let index = rng.index(2);
                let mut state = packing.site_state(index);
                state.position.x += rng.centered() * 0.3;
                state.position.y += rng.centered() * 0.3;
                state.angle += rng.centered();
                packing.set_site(index, state);
                let incremental = packing.site_overlaps(index);
                prop_assert_eq!(incremental, packing.any_overlap());
                if incremental {
                    // Return to a valid state
                    while packing.any_overlap() {
                        let lattice = packing.lattice().scaled(1.2);
                        packing.set_lattice(lattice);
                    }
                }
            }
        }
    }
}
