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

fn lens_area(r: f64, d: f64) -> f64 {
    r * r * libm::acos(d / r) - d * libm::sqrt(r * r - d * d)
}

/// The area of the overlap of two circles
fn circle_overlap(a: (Vec2, f64), b: (Vec2, f64)) -> f64 {
    let d = (a.0 - b.0).norm();
    if d >= a.1 + b.1 {
        return 0.;
    }
    if d <= (a.1 - b.1).abs() {
        let r = a.1.min(b.1);
        return PI * r * r;
    }
    let d1 = (d * d + a.1 * a.1 - b.1 * b.1) / (2. * d);
    let d2 = (d * d + b.1 * b.1 - a.1 * a.1) / (2. * d);
    lens_area(a.1, d1) + lens_area(b.1, d2)
}

impl Molecule {
    pub fn new(name: &str, atoms: Vec<(Vec2, f64)>) -> Result<Molecule, Error> {
        if atoms.is_empty() {
            bail!("A molecule requires at least one atom");
        }
        if atoms.iter().any(|(_, r)| *r <= 0.) {
            bail!("Atom radii must be positive");
        }
        // The area is calculated as the sum of areas minus the pairwise overlaps, which is exact
        // provided no point is covered by three or more atoms.
        let total: f64 = atoms.iter().map(|(_, r)| PI * r * r).sum();
        let mut overlap = 0.;
        for (i, a) in atoms.iter().enumerate() {
            for b in atoms.iter().skip(i + 1) {
                overlap += circle_overlap(*a, *b);
            }
        }
        let circumradius = atoms.iter().map(|(p, r)| p.norm() + r).fold(0., f64::max);
        let inradius = atoms.iter().map(|(p, r)| r - p.norm()).fold(0., f64::max);
        let symmetry = if atoms.len() == 1 && atoms[0].0.norm() < 1e-12 {
            ShapeSymmetry::continuous()
        } else {
            ShapeSymmetry::detect(&atoms)
        };
        Ok(Molecule {
            name: name.to_string(),
            atoms,
            area: total - overlap,
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

    #[test]
    fn disk_contact() {
        let disk = Molecule::disk();
        let mut a = PlacedMolecule::default();
        disk.place(&Mat2::IDENTITY, Vec2::ZERO, &mut a);
        assert!(!Molecule::overlaps(&a, &a, Vec2::new(2., 0.)));
        assert!(Molecule::overlaps(&a, &a, Vec2::new(1.999, 0.)));
    }

    #[test]
    fn trimer_symmetry() {
        let trimer = Molecule::trimer(0.637556, 120., 1.).unwrap();
        assert_eq!(trimer.symmetry().rotations, 1);
        assert_eq!(trimer.symmetry().mirrors.len(), 1);
        assert_abs_diff_eq!(trimer.symmetry().mirrors[0], PI / 2., epsilon = 1e-10);
    }

    #[test]
    fn overlap_area() {
        // Two unit circles a unit distance apart
        let area = circle_overlap((Vec2::ZERO, 1.), (Vec2::new(1., 0.), 1.));
        assert_abs_diff_eq!(area, 2. * PI / 3. - 3f64.sqrt() / 2., epsilon = 1e-12);
    }
}
