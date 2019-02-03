//
// cell.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use nalgebra::{Point2, Vector2};

pub use crate::basis::{Basis, SharedValue, StandardBasis};
pub use crate::symmetry::SymmetryTransform;

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
    y_len: Option<SharedValue>,
    angle: SharedValue,
    family: CrystalFamily,
}

impl std::fmt::Debug for Cell {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Cell {{ x: {}, y: {}, angle: {} }}",
            self.x(),
            self.y(),
            self.angle.get_value()
        )
    }
}

impl Default for Cell {
    fn default() -> Cell {
        Cell {
            x_len: SharedValue::new(1.),
            y_len: Some(SharedValue::new(1.)),
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
    pub fn to_cartesian_isometry(&self, transform: &SymmetryTransform) -> SymmetryTransform {
        let (x, y) = self.to_cartesian(transform.translation.x, transform.translation.y);
        SymmetryTransform {
            rotation: transform.rotation,
            translation: Vector2::new(x, y),
        }
    }

    /// The $x$ component of the cell, also known as $a$
    ///
    /// This is the interface for accessing the length of the first crystallographic dimension of
    /// the unit cell.
    pub fn x(&self) -> f64 {
        self.x_len.get_value()
    }

    /// The $y$ component of the cell, also known as $a$
    ///
    /// This is the interface for accessing the length of the first crystallographic dimension of
    /// the unit cell.
    pub fn y(&self) -> f64 {
        match &self.y_len {
            Some(y) => y.get_value(),
            None => self.x_len.get_value(),
        }
    }

    /// The angle between the $x$ and $y$ components of the cell
    ///
    /// This is typically labelled $\theta$.
    pub fn angle(&self) -> f64 {
        self.angle.get_value()
    }

    /// Convert a point in relative coordinates to real coordinates
    pub fn to_cartesian_point(&self, point: Point2<f64>) -> Point2<f64> {
        let (x, y) = self.to_cartesian(point.x, point.y);
        Point2::new(x, y)
    }

    /// Convert two values in relative coordinates to real coordinates
    pub fn to_cartesian(&self, x: f64, y: f64) -> (f64, f64) {
        (
            x * self.x() + y * self.y() * self.angle().cos(),
            y * self.y() * self.angle().sin(),
        )
    }

    /// Initialise a Cell instance from the CrystalFamily the cell belongs to
    ///
    /// Initialising from the Crystal family configures the Cell to the restrictions that the
    /// crystal family impose upon the unit cell. This includes ensuring both sides of the unit
    /// cell are the same length, or restricting the angle to a specific value.
    ///
    pub fn from_family(family: &CrystalFamily, length: f64) -> Cell {
        let (x_len, y_len, angle) = match family {
            // The Hexagonal Crystal has both sides equal with a fixed angle of 60 degrees.
            CrystalFamily::Hexagonal => (SharedValue::new(length), None, SharedValue::new(PI / 3.)),
            // The Tetragonal Crystal has both sides equal with a fixed angle of 90 degrees.
            CrystalFamily::Tetragonal => {
                (SharedValue::new(length), None, SharedValue::new(PI / 2.))
            }
            // The Orthorhombic crystal has two variable sides with a fixed angle of 90 degrees.
            CrystalFamily::Orthorhombic => (
                SharedValue::new(length),
                Some(SharedValue::new(length)),
                SharedValue::new(PI / 2.),
            ),
            // The Monoclinic cell has two variable sides and a variable angle initialised to 90
            // degrees
            CrystalFamily::Monoclinic => (
                SharedValue::new(length),
                Some(SharedValue::new(length)),
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

    /// This finds the values of the unit cell which are allowed to be changed and how
    ///
    /// Each of the different crystal families impose different restrictions on the degrees of
    /// freedom of a unit cell. This compiles these degrees of freedom into a vector of Bases,
    /// which is the data structure used to modify the values.
    pub fn get_degrees_of_freedom(&self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];

        // All cells have at least a single variable cell length
        basis.push(StandardBasis::new(
            &self.x_len,
            0.01,
            self.x_len.get_value(),
        ));

        // Both the Orthorhombic and Monoclinic cells have a second variable cell length. This is
        // indicated by the presence of the optional value.
        if let Some(y) = &self.y_len {
            basis.push(StandardBasis::new(&y, 0.01, y.get_value()));
        }

        // The Monoclinic family is the only one to have a variable cell angle.
        if self.family == CrystalFamily::Monoclinic {
            basis.push(StandardBasis::new(&self.angle, PI / 4., 3. * PI / 4.));
        }

        basis
    }

    /// The center of the cell in real space
    ///
    /// This finds the center of the unit cell so it can be aligned when output.
    ///
    /// The cell is centered around the point (0, 0), although there are no instances within the
    /// calculations that this is required, when trying to plot the unit cell it should be plotted
    /// with the center at the appropriate position.
    ///
    pub fn center(&self) -> Point2<f64> {
        let (x, y) = self.to_cartesian(0.5, 0.5);
        Point2::new(x, y)
    }

    /// Calculates the area of the cell
    ///
    /// The general formula for the area of a rhombus.
    pub fn area(&self) -> f64 {
        self.angle().sin() * self.x() * self.y()
    }
}

#[cfg(test)]
mod cell_tests {
    use super::*;

    // TODO Cell area test

    // TODO center test

    #[test]
    fn to_cartesian_test() {
        let cell = Cell::default();
        let trans = SymmetryTransform::new(na::Vector2::new(0.5, 0.5), 0.);

        assert_eq!(cell.to_cartesian_isometry(&trans), trans);

        cell.angle.set_value(PI / 4.);
        let expected = SymmetryTransform::new(
            na::Vector2::new(0.5 + 0.5 * 1. / f64::sqrt(2.), 0.5 * 1. / f64::sqrt(2.)),
            0.,
        );
        assert_abs_diff_eq!(cell.to_cartesian_isometry(&trans), expected);
    }
}
