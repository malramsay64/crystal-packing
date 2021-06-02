//
// cell2.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use itertools::iproduct;
use nalgebra::{Point2, Translation2};
use serde::{Deserialize, Serialize};

use crate::{Basis, SharedValue, Transform2};

/// The different crystal families that can be represented
///
/// These are all the valid types of crystal symmetries which are valid in a 2D space.
///
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
#[derive(Debug, Serialize, Deserialize)]
pub struct Cell2 {
    length: SharedValue,
    ratio: SharedValue,
    angle: SharedValue,
    family: CrystalFamily,
}

impl Clone for Cell2 {
    fn clone(&self) -> Self {
        Cell2 {
            length: SharedValue::new(self.length.get_value()),
            ratio: SharedValue::new(self.ratio.get_value()),
            angle: SharedValue::new(self.angle.get_value()),
            family: self.family,
        }
    }
}

impl std::fmt::Display for Cell2 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Cell2 {{ a: {}, b: {}, angle: {} }}",
            self.a(),
            self.b(),
            self.angle()
        )
    }
}

impl Default for Cell2 {
    fn default() -> Self {
        Self {
            length: SharedValue::new(1.),
            ratio: SharedValue::new(1.),
            angle: SharedValue::new(PI / 2.),
            family: CrystalFamily::Monoclinic,
        }
    }
}

impl Cell2 {
    pub fn a(&self) -> f64 {
        self.length.get_value()
    }

    pub fn b(&self) -> f64 {
        self.length.get_value() * self.ratio.get_value()
    }

    pub fn angle(&self) -> f64 {
        self.angle.get_value()
    }

    /// Convert a transformation into Cartesian coordinates
    ///
    /// The positions of particles are stored in fractional coordinates, making changes to the
    /// unit cell simple. This function takes transformation to apply to a collection of points
    /// and converts the values of the fractional coordinates in the translation to real
    /// Cartesian coordinates based on the current cell parameters.
    ///
    pub fn to_cartesian_isometry(&self, transform: Transform2) -> Transform2 {
        transform.set_position(self.to_cartesian_point(transform.position()))
    }

    /// Convert a point in relative coordinates to real coordinates
    ///
    /// # Example
    ///
    /// ```
    /// use crystal_packing::{Cell2, CrystalFamily};
    /// use nalgebra::Point2;
    /// let cell = Cell2::from_family(CrystalFamily::Monoclinic, 8.);
    /// let point = cell.to_cartesian_point(Point2::new(0.5, 0.5));
    /// assert_eq!(point, Point2::new(4., 4.));
    /// ```
    ///
    pub fn to_cartesian_point(&self, point: Point2<f64>) -> Point2<f64> {
        let (x, y) = self.to_cartesian(point.x, point.y);
        Point2::new(x, y)
    }

    /// This finds the values of the unit cell which are allowed to be changed and how
    ///
    /// Each of the different crystal families impose different restrictions on the degrees of
    /// freedom of a unit cell. This compiles these degrees of freedom into a vector of Bases,
    /// which is the data structure used to modify the values.
    pub fn get_degrees_of_freedom(&self) -> Vec<Basis> {
        let mut basis: Vec<Basis> = vec![
            // All cells have at least a single variable cell length
            Basis::StandardBasis {
            value: &self.length,
            min: 0.01,
            max: self.length.get_value(),
        }
        ];
        match self.family {
            // Monoclinic has both variable angle and varaible ratio of sides
            CrystalFamily::Monoclinic => {
                basis.push(Basis::StandardBasis {
                    value: &self.ratio,
                    min: 0.1,
                    max: self.ratio.get_value(),
                });
                basis.push(Basis::StandardBasis {
                    value: &self.angle,
                    min: PI / 4.,
                    max: PI / 2.,
                });
            }
            // The Orthorhombic have a second variable cell length in the ratio
            CrystalFamily::Orthorhombic => {
                basis.push(Basis::StandardBasis {
                    value: &self.ratio,
                    min: 0.1,
                    max: self.ratio.get_value(),
                });
            }
            _ => {}
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
    /// This uses the general formula for the area of a rhombus which is
    /// $ A = xy\sin(\theta) $
    ///
    pub fn area(&self) -> f64 {
        self.angle().sin() * self.a() * self.b()
    }

    /// Initialise a Cell instance from the CrystalFamily the cell belongs to
    ///
    /// Initialising from the Crystal family configures the Cell to the restrictions that the
    /// crystal family impose upon the unit cell. This includes ensuring both sides of the unit
    /// cell are the same length, or restricting the angle to a specific value.
    ///
    pub fn from_family(family: CrystalFamily, length: f64) -> Cell2 {
        let angle = match family {
            // The Hexagonal Crystal has both sides equal with a fixed angle of 60 degrees.
            CrystalFamily::Hexagonal => PI / 3.,
            // The Tetragonal, Orthorhombic, and Monoclinic all have an initial angle of 90 degrees
            _ => PI / 2.,
        };
        Cell2 {
            length: SharedValue::new(length),
            // The radio is initially always 1
            ratio: SharedValue::new(1.0),
            angle: SharedValue::new(angle),
            family,
        }
    }

    pub fn periodic_images(
        &self,
        transform: Transform2,
        shells: i64,
        zero: bool,
    ) -> impl Iterator<Item = Transform2> + '_ {
        iproduct!(-shells..=shells, -shells..=shells)
            .filter(move |&(x, y)| !(!zero && x == 0 && y == 0))
            .map(move |(x, y)| self.to_cartesian_translate(transform, x, y))
    }

    pub fn get_corners(&self) -> Vec<Point2<f64>> {
        let points = vec![
            Point2::new(-0.5, -0.5),
            Point2::new(-0.5, 0.5),
            Point2::new(0.5, 0.5),
            Point2::new(0.5, -0.5),
        ];

        points
            .into_iter()
            .map(|p| self.to_cartesian_point(p))
            .collect()
    }

    pub fn to_cartesian_translate(&self, transform: Transform2, x: i64, y: i64) -> Transform2 {
        let position = transform.position();
        transform
            .set_position(self.to_cartesian_point(Translation2::new(x as f64, y as f64) * position))
    }

    /// Convert two values in relative coordinates to real coordinates
    ///
    /// ```
    /// use crystal_packing::{Cell2, CrystalFamily};
    /// let cell = Cell2::from_family(CrystalFamily::Monoclinic, 8.);
    /// let point = cell.to_cartesian(0.25, 0.25);
    /// assert_eq!(point, (2., 2.));
    /// assert_eq!(cell.to_cartesian(0., 0.,), (0., 0.));
    /// assert_eq!(cell.to_cartesian(1., 1.,), (8., 8.));
    /// ```
    ///
    pub fn to_cartesian(&self, x: f64, y: f64) -> (f64, f64) {
        (
            x * self.a() + y * self.b() * self.angle().cos(),
            y * self.b() * self.angle().sin(),
        )
    }
}

#[cfg(test)]
mod cell_tests {
    use approx::assert_abs_diff_eq;
    use itertools::izip;

    use super::*;

    use crate::traits::{Intersect, Shape};
    use crate::LineShape;

    // TODO Cell area test

    // TODO center test

    #[test]
    fn to_cartesian_test() {
        let cell = Cell2::default();
        let trans = Transform2::new(0., (0.5, 0.5));

        assert_eq!(cell.to_cartesian_isometry(trans), trans);

        cell.angle.set_value(PI / 4.);
        let expected = Transform2::new(
            0.,
            (0.5 + 0.5 * 1. / f64::sqrt(2.), 0.5 * 1. / f64::sqrt(2.)),
        );
        assert_abs_diff_eq!(cell.to_cartesian_isometry(trans), expected);
    }

    #[test]
    fn periodic_intersection() {
        let shape = LineShape::from_radial("Square", vec![1.; 4]).unwrap();
        let cell = Cell2::default();
        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(transform, 1, false)
            .any(|t| shape.intersects(&shape.transform(&t)));

        assert!(intersection)
    }

    #[test]
    fn periodic_edge_intersection() {
        let shape = LineShape::from_radial("Square", vec![0.5; 4]).unwrap();
        let cell = Cell2::default();
        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(transform, 1, false)
            .any(|t| shape.intersects(&shape.transform(&t)));

        assert!(intersection)
    }

    #[test]
    fn no_periodic_intersection() {
        let shape = LineShape::from_radial("Square", vec![0.49; 4]).unwrap();
        let cell = Cell2::default();
        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(transform, 1, false)
            .any(|t| shape.intersects(&shape.transform(&t)));

        assert!(!intersection)
    }

    #[test]
    fn periodic_images_nozero() {
        let translations = vec![
            Point2::new(-1., -1.),
            Point2::new(-1., 0.),
            Point2::new(-1., 1.),
            Point2::new(0., -1.),
            Point2::new(0., 1.),
            Point2::new(1., -1.),
            Point2::new(1., -0.),
            Point2::new(1., 1.),
        ];
        let cell = Cell2::default();
        let transform = Transform2::identity();
        for (calculated, expected) in izip!(cell.periodic_images(transform, 1, false), translations)
        {
            assert_abs_diff_eq!(calculated.position(), expected);
        }
    }

    #[test]
    fn periodic_images() {
        let translations = vec![
            Point2::new(-1., -1.),
            Point2::new(-1., 0.),
            Point2::new(-1., 1.),
            Point2::new(0., -1.),
            Point2::new(0., 0.),
            Point2::new(0., 1.),
            Point2::new(1., -1.),
            Point2::new(1., -0.),
            Point2::new(1., 1.),
        ];
        let cell = Cell2::default();
        let transform = Transform2::identity();
        for (calculated, expected) in izip!(cell.periodic_images(transform, 1, true), translations)
        {
            assert_abs_diff_eq!(calculated.position(), expected);
        }
    }

    #[test]
    fn invalid_intersection() {
        let shape = LineShape::from_radial("Square", vec![1.; 4]).unwrap();
        let cell = Cell2 {
            length: SharedValue::new(1.59),
            ratio: SharedValue::new(0.83),
            angle: SharedValue::new(1.21),
            family: CrystalFamily::Monoclinic,
        };

        let transform = Transform2::new(0., (0., 0.));

        let intersection = cell
            .periodic_images(transform, 1, false)
            .any(|t| shape.intersects(&shape.transform(&t)));

        assert!(intersection)
    }
}
