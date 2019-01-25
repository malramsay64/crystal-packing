//
// cell.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

pub use crate::basis::{Basis, SharedValue, StandardBasis};
use nalgebra::{IsometryMatrix2, Point2};
use std::f64::consts::PI;

/// The different crystal families that can be represented
///
/// These are all the valid types of crystal symmetries which are valid in a 2D space.
///
#[derive(Debug, Clone, PartialEq)]
pub enum CrystalFamily {
    Monoclinic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
}

#[cfg(test)]
mod crystal_family_test {
    use super::*;

    #[test]
    fn equality() {
        assert_eq!(CrystalFamily::Monoclinic, CrystalFamily::Monoclinic);
        assert_eq!(CrystalFamily::Orthorhombic, CrystalFamily::Orthorhombic);
        assert_eq!(CrystalFamily::Hexagonal, CrystalFamily::Hexagonal);
        assert_eq!(CrystalFamily::Tetragonal, CrystalFamily::Tetragonal);
    }

    #[test]
    fn inequality() {
        assert_ne!(CrystalFamily::Orthorhombic, CrystalFamily::Monoclinic);
        assert_ne!(CrystalFamily::Hexagonal, CrystalFamily::Monoclinic);
        assert_ne!(CrystalFamily::Tetragonal, CrystalFamily::Monoclinic);
        assert_ne!(CrystalFamily::Hexagonal, CrystalFamily::Orthorhombic);
        assert_ne!(CrystalFamily::Tetragonal, CrystalFamily::Orthorhombic);
        assert_ne!(CrystalFamily::Tetragonal, CrystalFamily::Hexagonal);
    }
}

/// Representing the unit cell of a crystal packing
///
/// The unit cell holds the unit cell parameters, being the length of each side of the cell in
/// addition to the contained angles. Each cell belongs to one of the Crystal Families which
/// dictate the degrees of freedom the cell can take.
///
#[derive(Clone)]
pub struct Cell {
    x_len: SharedValue,
    y_len: SharedValue,
    angle: SharedValue,
    family: CrystalFamily,
}

impl std::fmt::Debug for Cell {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Cell {{ x: {}, y: {}, angle: {} }}",
            self.x_len.get_value(),
            self.y_len.get_value(),
            self.angle.get_value()
        )
    }
}

impl Default for Cell {
    fn default() -> Cell {
        Cell {
            x_len: SharedValue::new(1.),
            y_len: SharedValue::new(1.),
            angle: SharedValue::new(PI / 2.),
            family: CrystalFamily::Monoclinic,
        }
    }
}

impl Cell {
    /// Convert a transformation into cartesion coorsinates
    ///
    /// The positions of particles are stored in fractional coordinates, making changes to the
    /// unit cell simple. This function takes transformation to apply to a collection of points
    /// and converts the values of the fractional coordinates in the translation to real
    /// cartesian coordinates based on the current cell parameters.
    ///
    pub fn to_cartesian_isometry(&self, transform: &IsometryMatrix2<f64>) -> IsometryMatrix2<f64> {
        let (x, y) = self.to_cartesian(
            transform.translation.vector.x,
            transform.translation.vector.y,
        );
        IsometryMatrix2::from_parts(na::Translation2::new(x, y), transform.rotation)
    }

    pub fn to_cartesian_point(&self, point: Point2<f64>) -> Point2<f64> {
        let (x, y) = self.to_cartesian(point.x, point.y);
        Point2::new(x, y)
    }

    pub fn to_cartesian(&self, x: f64, y: f64) -> (f64, f64) {
        (
            x * self.x_len.get_value() + y * self.y_len.get_value() * self.angle.get_value().cos(),
            y * self.y_len.get_value() * self.angle.get_value().sin(),
        )
    }

    pub fn from_family(family: &CrystalFamily, length: f64) -> Cell {
        let (x_len, y_len, angle) = match family {
            // The Hexagonal Crystal has both sides equal with a fixed angle of 60 degrees.
            CrystalFamily::Hexagonal => {
                let len = SharedValue::new(length);
                (len.clone(), len.clone(), SharedValue::new(PI / 3.))
            }
            // The Tetragonal Crystal has both sides equal with a fixed angle of 90 degrees.
            CrystalFamily::Tetragonal => {
                let len = SharedValue::new(length);
                (len.clone(), len.clone(), SharedValue::new(PI / 2.))
            }
            // The Orthorhombic crystal has two variable sides with a fixed angle of 90 degrees.
            CrystalFamily::Orthorhombic => (
                SharedValue::new(length),
                SharedValue::new(length),
                SharedValue::new(PI / 2.),
            ),
            // The Monoclinic cell has two variable sides and a variable angle initialised to 90
            // degrees
            CrystalFamily::Monoclinic => (
                SharedValue::new(length),
                SharedValue::new(length),
                SharedValue::new(PI / 2.),
            ),
        };
        Cell {
            x_len,
            y_len,
            angle,
            family: family.clone(),
        }
    }

    pub fn get_basis(&self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];

        // All cells have at least a single variable cell length
        basis.push(StandardBasis::new(
            &self.x_len,
            0.01,
            self.x_len.get_value(),
        ));

        // Both the Orthorhombic and Monoclinic cells have a second variable cell length
        if (self.family == CrystalFamily::Orthorhombic) | (self.family == CrystalFamily::Monoclinic)
        {
            basis.push(StandardBasis::new(
                &self.y_len,
                0.01,
                self.y_len.get_value(),
            ));
        }

        // The Monoclinic family is the only one to have a variable cell angle.
        if self.family == CrystalFamily::Monoclinic {
            basis.push(StandardBasis::new(&self.angle, PI / 4., 3. * PI / 4.));
        }

        basis
    }

    pub fn center(&self) -> Point2<f64> {
        let (x, y) = self.to_cartesian(0.5, 0.5);
        Point2::new(x, y)
    }

    pub fn area(&self) -> f64 {
        self.angle.get_value().sin() * self.x_len.get_value() * self.y_len.get_value()
    }
}

#[cfg(test)]
mod cell_tests {
    use super::*;

    #[test]
    fn to_cartesian_test() {
        let cell = Cell::default();
        let trans = na::IsometryMatrix2::new(na::Vector2::new(0.5, 0.5), 0.);

        assert_eq!(cell.to_cartesian_isometry(&trans), trans);

        cell.angle.set_value(PI / 4.);
        let expected = na::IsometryMatrix2::new(
            na::Vector2::new(0.5 + 0.5 * 1. / f64::sqrt(2.), 0.5 * 1. / f64::sqrt(2.)),
            0.,
        );
        assert_abs_diff_eq!(cell.to_cartesian_isometry(&trans), expected);
    }
}
