//
// lj_shape.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::{fmt, slice, vec};

use itertools::iproduct;
use nalgebra::Point2;
use serde::{Deserialize, Serialize};

use super::{Transform2, LJ2};
use crate::traits::{Potential, Shape};

/// A shape defined by a collection of Atoms
///
/// This is a shape comprised of a series of circles which each have a position and radius.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LJShape2 {
    pub name: String,
    pub items: Vec<LJ2>,
}

impl<'a> IntoIterator for &'a LJShape2 {
    type Item = &'a LJ2;
    type IntoIter = slice::Iter<'a, LJ2>;

    fn into_iter(self) -> Self::IntoIter {
        self.items.iter()
    }
}

impl Potential for LJShape2 {
    fn energy(&self, other: &Self) -> f64 {
        iproduct!(self.items.iter(), other.items.iter())
            .map(|(s, o)| s.energy(o))
            .sum()
    }
}

impl Shape for LJShape2 {
    type Component = LJ2;

    fn score(&self, other: &Self) -> Option<f64> {
        Some(
            iproduct!(self.items.iter(), other.items.iter())
                .fold(0., |sum, (s, o)| sum + s.energy(o)),
        )
    }

    fn enclosing_radius(&self) -> f64 {
        self.items
            .iter()
            .map(|p| p.position.norm() + p.sigma / 2.)
            // The f64 type doesn't have complete ordering because of Nan and Inf, so the
            // standard min/max comparators don't work. Instead we use the f64::max which ignores
            // the NAN and max values.
            .fold(std::f64::MIN, f64::max)
    }

    fn get_items(&self) -> Vec<Self::Component> {
        self.items.clone()
    }

    fn iter(&self) -> slice::Iter<'_, Self::Component> {
        self.into_iter()
    }

    fn transform(&self, transform: &Transform2) -> Self {
        Self {
            name: self.name.clone(),
            items: self.into_iter().map(|i| i * transform).collect(),
        }
    }
}

impl fmt::Display for LJShape2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "LJShape2 {{ ")?;
        for item in self.items.iter() {
            write!(f, "{},", item)?;
        }
        write!(f, " }}")
    }
}

impl LJShape2 {
    /// Create a Trimer molecule instance
    ///
    /// A Trimer is a molecule consisting of three particles, a central particle of radius 1 and
    /// two smaller particles with a radius `radius`, separated by angle `angle` and at a distance
    /// `distance` from the center of the central particle. This is a class of particle I am
    /// studying in my research.
    ///
    /// # Arguments
    ///
    /// - `radius` - The radius of the smaller particle
    /// - `angle` - The angle between smaller particles in degrees
    /// - `distance` - The distance of the smaller particle from the center
    ///
    /// # Example
    ///
    /// ```
    /// # use packing::LJShape2;
    /// let shape = LJShape2::from_trimer(0.7, 120., 1.);
    /// # assert_eq!(shape.items.len(), 3);
    /// # assert_eq!(shape.name, "Trimer");
    /// ```
    ///
    pub fn from_trimer(radius: f64, angle: f64, distance: f64) -> Self {
        let x_base = distance * f64::sin(angle.to_radians() / 2.);
        let y_base = 1. / 3. * distance * f64::cos(angle.to_radians() / 2.);
        let positions = vec![
            (1., Point2::new(0., -2. * y_base)),
            (radius, Point2::new(-x_base, y_base)),
            (radius, Point2::new(x_base, y_base)),
        ];
        Self {
            name: String::from("Trimer"),
            items: positions
                .into_iter()
                .map(|(r, p)| LJ2 {
                    position: p.coords,
                    sigma: 2. * r,
                    cutoff: Some(3.5),
                    ..Default::default()
                })
                .collect(),
        }
    }

    /// Create an instance of a Circle
    ///
    /// This is the simplest molecular shape, a single circle at the origin with radius of 1.0.
    ///
    /// # Example
    ///
    /// ```
    /// # use packing::LJShape2;
    /// let shape = LJShape2::circle();
    /// # assert_eq!(shape.name, "circle");
    /// # assert_eq!(shape.items.len(), 1);
    /// ```
    ///
    pub fn circle() -> Self {
        Self {
            name: String::from("circle"),
            items: vec![LJ2::new(0., 0., 1.)],
        }
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;

    use super::*;
    use nalgebra::Vector2;

    #[test]
    fn from_trimer_test() {
        let shape = LJShape2::from_trimer(1., 180., 1.);
        assert_eq!(shape.items.len(), 3);

        assert_abs_diff_eq!(shape.items[0].position, Vector2::new(0., 0.));
        assert_abs_diff_eq!(shape.items[1].position, Vector2::new(-1., 0.));
        assert_abs_diff_eq!(shape.items[2].position, Vector2::new(1., 0.));

        let shape = LJShape2::from_trimer(0.637_556, 120., 1.);
        assert_abs_diff_eq!(shape.items[0].position, Vector2::new(0., -1. / 3.));
        assert_abs_diff_eq!(
            shape.items[1].position,
            Vector2::new(-0.866, 1. / 6.),
            epsilon = 1e-3,
        );
        assert_abs_diff_eq!(
            shape.items[2].position,
            Vector2::new(0.866, 1. / 6.),
            epsilon = 1e-3,
        );
    }

    #[test]
    fn trimer_not_hardcoded_sigma() {
        let shape = LJShape2::from_trimer(2., 180., 0.5);
        assert_abs_diff_eq!(shape.items[0].sigma, 2.);
        assert_abs_diff_eq!(shape.items[1].sigma, 4.);
        assert_abs_diff_eq!(shape.items[2].sigma, 4.);
    }

    #[test]
    fn trimer_cutoff() {
        let shape = LJShape2::from_trimer(2., 180., 0.5);
        assert_abs_diff_eq!(shape.items[0].cutoff.unwrap(), 3.5);
        assert_abs_diff_eq!(shape.items[1].cutoff.unwrap(), 3.5);
        assert_abs_diff_eq!(shape.items[2].cutoff.unwrap(), 3.5);
    }
}
