//
// cell2.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use itertools::iproduct;
use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

use crate::traits::*;
use crate::{SharedValue, StandardBasis, Transform2};

/// The different crystal families that can be represented
///
/// These are all the valid types of crystal symmetries which are valid in a 2D space.
///
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Cell2 {
    points: Vector2<f64>,
    angles: Vector2<f64>,
    family: CrystalFamily,
}

impl std::fmt::Display for Cell2 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Cell2 {{ x: {}, y: {}, angle: {} }}",
            self.x(),
            self.y(),
            self.angle()
        )
    }
}

impl Default for Cell2 {
    fn default() -> Self {
        Self {
            points: Vector2::new(1., 1.),
            angles: Vector2::new(PI / 2., 0.),
            family: CrystalFamily::Monoclinic,
        }
    }
}

impl Cell for Cell2 {
    /// Convert a transformation into Cartesian coordinates
    ///
    /// The positions of particles are stored in fractional coordinates, making changes to the
    /// unit cell simple. This function takes transformation to apply to a collection of points
    /// and converts the values of the fractional coordinates in the translation to real
    /// Cartesian coordinates based on the current cell parameters.
    ///
    fn to_cartesian_isometry(&self, transform: &Transform2) -> Transform2 {
        transform.set_position(self.to_cartesian_point(transform.position()))
    }

    /// Convert a point in relative coordinates to real coordinates
    fn to_cartesian_point(&self, point: Vector2<f64>) -> Vector2<f64> {
        let (x, y) = self.to_cartesian(point.x, point.y);
        Vector2::new(x, y)
    }

    /// This finds the values of the unit cell which are allowed to be changed and how
    ///
    /// Each of the different crystal families impose different restrictions on the degrees of
    /// freedom of a unit cell. This compiles these degrees of freedom into a vector of Bases,
    /// which is the data structure used to modify the values.
    fn get_degrees_of_freedom(&mut self) -> Vec<StandardBasis> {
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
    fn center(&self) -> Vector2<f64> {
        let (x, y) = self.to_cartesian(0.5, 0.5);
        Vector2::new(x, y)
    }

    /// Calculates the area of the cell
    ///
    /// The general formula for the area of a rhombus.
    fn area(&self) -> f64 {
        self.angle().sin() * self.x() * self.y()
    }

    /// Initialise a Cell instance from the CrystalFamily the cell belongs to
    ///
    /// Initialising from the Crystal family configures the Cell to the restrictions that the
    /// crystal family impose upon the unit cell. This includes ensuring both sides of the unit
    /// cell are the same length, or restricting the angle to a specific value.
    ///
    fn from_family(family: &CrystalFamily, length: f64) -> Cell2 {
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
        Cell2 {
            points: Vector2::new(x_len, y_len),
            angles: Vector2::new(angle, 0.),
            family: family.clone(),
        }
    }

    fn periodic_images(&self, transform: &Transform2, zero: bool) -> Vec<Transform2> {
        // The periodic images to check. Checking the first and second shells i.e.
        // -2..=2, as this is necessary to ensure no intersections on tilted cells
        // and highly irregular cells.
        let iter_range = match (self.x() / self.y(), self.angle()) {
            (p, a) if 0.5 < p && p < 2. && f64::abs(a - PI / 2.) < 0.2 => -1..=1,
            (p, a) if 0.3 < p && p < 3. && f64::abs(a - PI / 2.) < 0.5 => -2..=2,
            _ => -3..=3,
        };

        if zero {
            iproduct!(iter_range.clone(), iter_range.clone())
                .map(|(x, y)| self.to_cartesian_translate(transform, x, y))
                .collect()
        } else {
            iproduct!(iter_range.clone(), iter_range.clone())
                .filter(|&(x, y)| !(x == 0 && y == 0))
                .map(|(x, y)| self.to_cartesian_translate(transform, x, y))
                .collect()
        }
    }

    fn get_corners(&self) -> Vec<Vector2<f64>> {
        let points = vec![
            Vector2::new(-0.5, -0.5),
            Vector2::new(-0.5, 0.5),
            Vector2::new(0.5, 0.5),
            Vector2::new(0.5, -0.5),
        ];

        points
            .into_iter()
            .map(|p| self.to_cartesian_point(p))
            .collect()
    }
}

impl Cell2 {
    fn to_cartesian_translate(&self, transform: &Transform2, x: i64, y: i64) -> Transform2 {
        let mut position = transform.position();
        position.x += x as f64;
        position.y += y as f64;
        transform.set_position(self.to_cartesian_point(position))
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

    /// Convert two values in relative coordinates to real coordinates
    pub fn to_cartesian(&self, x: f64, y: f64) -> (f64, f64) {
        (
            x * self.x() + y * self.y() * self.angle().cos(),
            y * self.y() * self.angle().sin(),
        )
    }
}

#[cfg(test)]
mod cell_tests {
    use approx::assert_abs_diff_eq;
    use itertools::izip;

    use super::*;

    use crate::LineShape;

    // TODO Cell area test

    // TODO center test

    #[test]
    fn to_cartesian_test() {
        let mut cell = Cell2::default();
        let trans = Transform2::new(0., (0.5, 0.5));

        assert_eq!(cell.to_cartesian_isometry(&trans), trans);

        cell.angles.x = PI / 4.;
        let expected = Transform2::new(
            0.,
            (0.5 + 0.5 * 1. / f64::sqrt(2.), 0.5 * 1. / f64::sqrt(2.)),
        );
        assert_abs_diff_eq!(cell.to_cartesian_isometry(&trans), expected);
    }

    #[test]
    fn periodic_intersection() {
        let shape = LineShape::from_radial("Square", vec![1.; 4]).unwrap();
        let cell = Cell2::default();
        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(&transform, false)
            .iter()
            .any(|t| shape.intersects(&shape.transform(t)));

        assert!(intersection)
    }

    #[test]
    fn periodic_edge_intersection() {
        let shape = LineShape::from_radial("Square", vec![0.5; 4]).unwrap();
        let cell = Cell2::default();
        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(&transform, false)
            .iter()
            .any(|t| shape.intersects(&shape.transform(t)));

        assert!(intersection)
    }

    #[test]
    fn no_periodic_intersection() {
        let shape = LineShape::from_radial("Square", vec![0.49; 4]).unwrap();
        let cell = Cell2::default();
        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(&transform, false)
            .iter()
            .any(|t| shape.intersects(&shape.transform(t)));

        assert!(!intersection)
    }

    #[test]
    fn periodic_images_nozero() {
        let translations = vec![
            Vector2::new(-1., -1.),
            Vector2::new(-1., 0.),
            Vector2::new(-1., 1.),
            Vector2::new(0., -1.),
            Vector2::new(0., 1.),
            Vector2::new(1., -1.),
            Vector2::new(1., -0.),
            Vector2::new(1., 1.),
        ];
        let cell = Cell2::default();
        let transform = Transform2::identity();
        for (calculated, expected) in izip!(cell.periodic_images(&transform, false), translations) {
            assert_abs_diff_eq!(calculated.position(), expected);
        }
    }

    #[test]
    fn periodic_images() {
        let translations = vec![
            Vector2::new(-1., -1.),
            Vector2::new(-1., 0.),
            Vector2::new(-1., 1.),
            Vector2::new(0., -1.),
            Vector2::new(0., 0.),
            Vector2::new(0., 1.),
            Vector2::new(1., -1.),
            Vector2::new(1., -0.),
            Vector2::new(1., 1.),
        ];
        let cell = Cell2::default();
        let transform = Transform2::identity();
        for (calculated, expected) in izip!(cell.periodic_images(&transform, true), translations) {
            assert_abs_diff_eq!(calculated.position(), expected);
        }
    }

    #[test]
    fn invalid_intersection() {
        let shape = LineShape::from_radial("Square", vec![1.; 4]).unwrap();
        let cell = Cell2 {
            points: Vector2::new(1.32, 1.59),
            angles: Vector2::new(1.21, 0.),
            family: CrystalFamily::Monoclinic,
        };

        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(&transform, false)
            .iter()
            .any(|t| shape.intersects(&shape.transform(t)));

        assert!(intersection)
    }
}
