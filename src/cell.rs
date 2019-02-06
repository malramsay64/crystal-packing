//
// cell.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Point, Translation, VectorN};

pub use crate::basis::{Basis, SharedValue, StandardBasis};
use crate::shape::U2;
pub use crate::symmetry::Transform;

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
#[derive(Clone, Debug)]
pub struct Cell<D: U2>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    points: VectorN<f64, D>,
    angles: VectorN<f64, D>,
    family: CrystalFamily,
}

impl std::fmt::Display for Cell<na::U2> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Cell {{ x: {}, a: {}, angle: {} }}",
            self.x(),
            self.y(),
            self.angle()
        )
    }
}

impl Default for Cell<na::U2> {
    fn default() -> Self {
        Self {
            points: VectorN::<f64, na::U2>::new(1., 1.),
            angles: VectorN::<f64, na::U2>::new(PI / 2., 0.),
            family: CrystalFamily::Monoclinic,
        }
    }
}

impl Cell<na::U2> {
    /// Convert a transformation into Cartesian coordinates
    ///
    /// The positions of particles are stored in fractional coordinates, making changes to the
    /// unit cell simple. This function takes transformation to apply to a collection of points
    /// and converts the values of the fractional coordinates in the translation to real
    /// Cartesian coordinates based on the current cell parameters.
    ///
    pub fn to_cartesian_isometry(&self, transform: &Transform<na::U2>) -> Transform<na::U2> {
        let (x, y) = self.to_cartesian(
            transform.translation.vector.x,
            transform.translation.vector.y,
        );
        Transform::<na::U2>::from_parts(Translation::<f64, na::U2>::new(x, y), transform.rotation)
    }

    /// The $x$ component of the cell, also known as $a$
    ///
    /// This is the interface for accessing the length of the first crystallographic dimension of
    /// the unit cell.
    pub fn x(&self) -> f64 {
        self.points.x
    }

    /// The $y$ component of the cell, also known as $a$
    ///
    /// This is the interface for accessing the length of the first crystallographic dimension of
    /// the unit cell.
    pub fn y(&self) -> f64 {
        if self.points.y == 0. {
            return self.points.x;
        }
        self.points.y
    }

    /// The angle between the $x$ and $y$ components of the cell
    ///
    /// This is typically labelled $\theta$.
    pub fn angle(&self) -> f64 {
        self.angles.x
    }

    /// Convert a point in relative coordinates to real coordinates
    pub fn to_cartesian_point(&self, point: Point<f64, na::U2>) -> Point<f64, na::U2> {
        let (x, y) = self.to_cartesian(point.x, point.y);
        Point::<f64, na::U2>::new(x, y)
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
    pub fn from_family(family: &CrystalFamily, length: f64) -> Cell<na::U2> {
        let (x_len, y_len, angle) = match family {
            // The Hexagonal Crystal has both sides equal with a fixed angle of 60 degrees.
            CrystalFamily::Hexagonal => (length, 0., PI / 3.),
            // The Tetragonal Crystal has both sides equal with a fixed angle of 90 degrees.
            CrystalFamily::Tetragonal => (length, 0., PI / 2.),
            // The Orthorhombic crystal has two variable sides with a fixed angle of 90 degrees.
            CrystalFamily::Orthorhombic => (length, length, PI / 2.),
            // The Monoclinic cell has two variable sides and a variable angle initialised to 90
            // degrees
            CrystalFamily::Monoclinic => (length, length, PI / 2.),
        };
        Cell {
            points: VectorN::<f64, na::U2>::new(x_len, y_len),
            angles: VectorN::<f64, na::U2>::new(angle, 0.),
            family: family.clone(),
        }
    }

    /// This finds the values of the unit cell which are allowed to be changed and how
    ///
    /// Each of the different crystal families impose different restrictions on the degrees of
    /// freedom of a unit cell. This compiles these degrees of freedom into a vector of Bases,
    /// which is the data structure used to modify the values.
    pub fn get_degrees_of_freedom(&mut self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];

        // All cells have at least a single variable cell length
        basis.push(StandardBasis::new(
            SharedValue::new(&mut self.points.x),
            0.01,
            self.x(),
        ));

        // Both the Orthorhombic and Monoclinic cells have a second variable cell length. This is
        // indicated by the presence of the optional value.
        if self.points.y != 0. {
            basis.push(StandardBasis::new(
                SharedValue::new(&mut self.points.y),
                0.01,
                self.y(),
            ));
        }

        // The Monoclinic family is the only one to have a variable cell angle.
        if self.family == CrystalFamily::Monoclinic {
            basis.push(StandardBasis::new(
                SharedValue::new(&mut self.angles.x),
                PI / 4.,
                3. * PI / 4.,
            ));
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
    pub fn center(&self) -> Point<f64, na::U2> {
        let (x, y) = self.to_cartesian(0.5, 0.5);
        Point::<f64, na::U2>::new(x, y)
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
    use approx::assert_abs_diff_eq;

    use super::*;

    use crate::Transform2;

    // TODO Cell area test

    // TODO center test

    #[test]
    fn to_cartesian_test() {
        let mut cell = Cell::default();
        let trans = Transform2::new(na::Vector2::new(0.5, 0.5), 0.);

        assert_eq!(cell.to_cartesian_isometry(&trans), trans);

        cell.angles.x = PI / 4.;
        let expected = Transform2::new(
            na::Vector2::new(0.5 + 0.5 * 1. / f64::sqrt(2.), 0.5 * 1. / f64::sqrt(2.)),
            0.,
        );
        assert_abs_diff_eq!(cell.to_cartesian_isometry(&trans), expected);
    }
}
