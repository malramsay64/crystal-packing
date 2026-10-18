//
// lattice.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! The unit cell of a 2D crystal.
//!
//! The cell is stored as an upper triangular matrix whose columns are the lattice vectors, so
//! the first lattice vector always lies along the x axis. The free parameters of the matrix depend
//! on the lattice family. To find the periodic images which need to be checked for overlaps, the
//! lattice vectors are first reduced (Lagrange-Gauss reduction) to the shortest, most orthogonal
//! basis of the same lattice, which keeps the number of images small regardless of how sheared the
//! cell becomes.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::geometry::{Mat2, Vec2};
use crate::symmetry::LatticeFamily;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Lattice {
    pub family: LatticeFamily,
    /// The cell matrix, with columns being the lattice vectors a and b
    pub matrix: Mat2,
}

/// Reduce a basis to the shortest pair of vectors spanning the same lattice
///
/// ```
/// use crystal_packing::geometry::Vec2;
/// use crystal_packing::lattice::gauss_reduce;
/// let (u, v) = gauss_reduce(Vec2::new(1., 0.), Vec2::new(5.2, 1.));
/// assert!((u.norm() - 1.).abs() < 1e-12);
/// assert!(v.x.abs() <= 0.5);
/// ```
pub fn gauss_reduce(a: Vec2, b: Vec2) -> (Vec2, Vec2) {
    let (mut u, mut v) = (a, b);
    for _ in 0..100 {
        if u.norm_sq() > v.norm_sq() {
            std::mem::swap(&mut u, &mut v);
        }
        let mu = u.dot(v) / u.norm_sq();
        // The basis is reduced once |μ| ≤ 1/2; testing μ itself rather than its rounded value
        // stops the reduction cycling when μ = ±1/2 exactly, as in a hexagonal lattice.
        if mu.abs() <= 0.5 + 1e-12 {
            break;
        }
        v -= u * libm::round(mu);
    }
    (u, v)
}

impl Lattice {
    /// The number of free parameters for a lattice family
    pub fn num_parameters(family: LatticeFamily) -> usize {
        match family {
            LatticeFamily::Oblique => 3,
            LatticeFamily::Rectangular => 2,
            LatticeFamily::Square | LatticeFamily::Hexagonal => 1,
        }
    }

    /// Create a lattice from its free parameters
    ///
    /// - Oblique: [a_x, b_x, b_y]
    /// - Rectangular: [a, b]
    /// - Square and Hexagonal: [a]
    pub fn from_parameters(family: LatticeFamily, p: &[f64]) -> Lattice {
        let matrix = match family {
            LatticeFamily::Oblique => Mat2::new(p[0], p[1], 0., p[2]),
            LatticeFamily::Rectangular => Mat2::new(p[0], 0., 0., p[1]),
            LatticeFamily::Square => Mat2::new(p[0], 0., 0., p[0]),
            LatticeFamily::Hexagonal => family.reference_cell().scale(p[0]),
        };
        Lattice { family, matrix }
    }

    pub fn parameters(&self) -> Vec<f64> {
        let m = &self.matrix.m;
        match self.family {
            LatticeFamily::Oblique => vec![m[0][0], m[0][1], m[1][1]],
            LatticeFamily::Rectangular => vec![m[0][0], m[1][1]],
            LatticeFamily::Square | LatticeFamily::Hexagonal => vec![m[0][0]],
        }
    }

    /// A lattice with the side length `a` and the default shape for the family
    pub fn with_length(family: LatticeFamily, a: f64) -> Lattice {
        match family {
            LatticeFamily::Oblique => Lattice::from_parameters(family, &[a, 0., a]),
            LatticeFamily::Rectangular => Lattice::from_parameters(family, &[a, a]),
            _ => Lattice::from_parameters(family, &[a]),
        }
    }

    pub fn area(&self) -> f64 {
        self.matrix.det().abs()
    }

    pub fn a(&self) -> Vec2 {
        self.matrix.column(0)
    }

    pub fn b(&self) -> Vec2 {
        self.matrix.column(1)
    }

    /// The angle between the lattice vectors a and b
    pub fn angle(&self) -> f64 {
        let (a, b) = (self.a(), self.b());
        libm::acos((a.dot(b) / (a.norm() * b.norm())).clamp(-1., 1.))
    }

    pub fn to_cartesian(&self, fractional: Vec2) -> Vec2 {
        self.matrix * fractional
    }

    pub fn scaled(&self, factor: f64) -> Lattice {
        Lattice {
            family: self.family,
            matrix: self.matrix.scale(factor),
        }
    }

    /// The reduced basis of the lattice as the columns of a matrix
    pub fn reduced_basis(&self) -> Mat2 {
        let (u, v) = gauss_reduce(self.a(), self.b());
        Mat2::from_columns(u, v)
    }

    /// All lattice translations with a length at most `cutoff`, including the zero translation.
    pub fn translations(&self, cutoff: f64) -> Vec<Vec2> {
        let reduced = self.reduced_basis();
        let (u, v) = (reduced.column(0), reduced.column(1));
        let det = u.cross(v).abs();
        // From |u × t| = |m| |u × v| ≤ |u| |t| we get bounds on the coefficients of t = n u + m v
        let n_max = libm::ceil(cutoff * v.norm() / det) as i64;
        let m_max = libm::ceil(cutoff * u.norm() / det) as i64;
        let cutoff_sq = cutoff * cutoff;
        let mut translations = vec![];
        for n in -n_max..=n_max {
            for m in -m_max..=m_max {
                let t = u * n as f64 + v * m as f64;
                if t.norm_sq() <= cutoff_sq {
                    translations.push(t);
                }
            }
        }
        // Order by distance so the closest images, which are most likely to overlap, are first
        translations.sort_by(|a, b| a.norm_sq().total_cmp(&b.norm_sq()));
        translations
    }

    /// Check the lattice is a valid shape, with positive side lengths and a reasonable angle
    pub fn is_valid(&self) -> bool {
        let m = &self.matrix.m;
        m[0][0] > 0. && m[1][1] > 0. && m.iter().flatten().all(|v| v.is_finite())
    }

    /// The angle of the cell in degrees after reduction, useful for reporting
    pub fn reduced_angle_degrees(&self) -> f64 {
        let r = self.reduced_basis();
        let (u, v) = (r.column(0), r.column(1));
        libm::acos((u.dot(v) / (u.norm() * v.norm())).clamp(-1., 1.)) * 180. / PI
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    #[test]
    fn hexagonal() {
        let lattice = Lattice::from_parameters(LatticeFamily::Hexagonal, &[2.]);
        assert_abs_diff_eq!(lattice.angle(), 2. * PI / 3., epsilon = 1e-12);
        assert_abs_diff_eq!(lattice.b().norm(), 2., epsilon = 1e-12);
        assert_abs_diff_eq!(lattice.area(), 2. * 3f64.sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn hexagonal_reduction_terminates() {
        let lattice = Lattice::from_parameters(LatticeFamily::Hexagonal, &[1.]);
        let reduced = lattice.reduced_basis();
        assert_abs_diff_eq!(reduced.column(0).norm(), 1., epsilon = 1e-12);
        assert_abs_diff_eq!(reduced.column(1).norm(), 1., epsilon = 1e-12);
    }

    #[test]
    fn square_translations() {
        let lattice = Lattice::with_length(LatticeFamily::Square, 1.);
        // Within a radius of 1: the origin and the four nearest neighbours
        assert_eq!(lattice.translations(1.).len(), 5);
        // Within sqrt(2): the diagonals are included
        assert_eq!(lattice.translations(2f64.sqrt() + 1e-9).len(), 9);
    }

    proptest! {
        #[test]
        fn reduction_preserves_lattice(ax in 0.1f64..5., bx in -20f64..20., by in 0.1f64..5.) {
            let lattice = Lattice::from_parameters(LatticeFamily::Oblique, &[ax, bx, by]);
            let reduced = lattice.reduced_basis();
            // The area is unchanged
            prop_assert!((reduced.det().abs() - lattice.area()).abs() < 1e-9 * lattice.area().max(1.));
            // The reduced vectors are integer combinations of the original basis
            let coefficients = lattice.matrix.inverse() * reduced;
            for v in coefficients.m.iter().flatten() {
                prop_assert!((v - v.round()).abs() < 1e-6);
            }
            // The reduced basis is nearly orthogonal (angle between 60° and 120°)
            let (u, v) = (reduced.column(0), reduced.column(1));
            prop_assert!(u.dot(v).abs() <= 0.5 * u.norm_sq() + 1e-9);
        }

        #[test]
        fn translations_match_brute_force(
            ax in 0.5f64..3., bx in -3f64..3., by in 0.5f64..3., cutoff in 0.1f64..5.
        ) {
            let lattice = Lattice::from_parameters(LatticeFamily::Oblique, &[ax, bx, by]);
            let found = lattice.translations(cutoff);
            let mut expected = 0;
            for n in -80i64..=80 {
                for m in -80i64..=80 {
                    let t = lattice.a() * n as f64 + lattice.b() * m as f64;
                    if t.norm_sq() <= cutoff * cutoff {
                        expected += 1;
                        prop_assert!(found.iter().any(|f| (*f - t).norm() < 1e-9));
                    }
                }
            }
            prop_assert_eq!(found.len(), expected);
        }
    }
}
