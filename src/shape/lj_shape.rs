//
// lj_shape.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;
use std::{fmt, slice, vec};

use itertools::{iproduct, Itertools};
use nalgebra::Point2;
use serde::{Deserialize, Serialize};

use super::{LJ2, Transform2};
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
    fn intersects(&self, other: &Self) -> f64 {
        iproduct!(self.items.iter(), other.items.iter()).map(|(s, o)| s.energy(o)).sum()
    }
}

impl Shape for LJShape2 {
    type Component = LJ2;
    type Transform = Transform2;

    fn area(&self) -> f64 {
        // TODO Implement an algorithm which takes into account multiple overlaps of circles, this
        // naive implementation is just a temporary measure.
        let total_area: f64 = self.items.iter().map(|a| PI * a.radius.powi(2)).sum();

        let naive_overlap: f64 = self
            .items
            .iter()
            .tuple_combinations()
            .map(|(a1, a2)| Self::circle_overlap(a1, a2))
            .sum();

        total_area - naive_overlap
    }

    fn enclosing_radius(&self) -> f64 {
        self.items
            .iter()
            .map(|p| na::distance(&Point2::origin(), &p.position) + p.radius)
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

    fn transform(&self, transform: &Self::Transform) -> Self {
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
    pub fn from_trimer(radius: f64, angle: f64, distance: f64) -> Self {
        Self {
            name: String::from("Trimer"),
            items: vec![
                Lj2::new(0., -2. / 3. * distance * f64::cos(angle / 2.), 1.),
                Lj2::new(
                    -distance * f64::sin(angle / 2.),
                    1. / 3. * distance * f64::cos(angle / 2.),
                    radius,
                ),
                Lj2::new(
                    distance * f64::sin(angle / 2.),
                    1. / 3. * distance * f64::cos(angle / 2.),
                    radius,
                ),
            ],
        }
    }

    /// Create an instance of a Circle
    ///
    /// This is the simplest molecular shape, a single circle at the origin with radius of 1.0.
    pub fn circle() -> Self {
        Self {
            name: String::from("circle"),
            items: vec![Lj2::new(0., 0., 1.)],
        }
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector2;

    use super::*;

    #[test]
    fn from_trimer_test() {
        let shape = LJShape2::from_trimer(1., PI, 1.);
        assert_eq!(shape.items.len(), 3);

        assert_abs_diff_eq!(shape.items[0].position, Point2::new(0., 0.));
        assert_abs_diff_eq!(shape.items[1].position, Point2::new(-1., 0.));
        assert_abs_diff_eq!(shape.items[2].position, Point2::new(1., 0.));

        let shape = LJShape2::from_trimer(0.637_556, 2. * PI / 3., 1.);
        assert_abs_diff_eq!(shape.items[0].position, Point2::new(0., -1. / 3.));
        assert_abs_diff_eq!(
            shape.items[1].position,
            Point2::new(-0.866, 1. / 6.),
            epsilon = 1e-3,
        );
        assert_abs_diff_eq!(
            shape.items[2].position,
            Point2::new(0.866, 1. / 6.),
            epsilon = 1e-3,
        );
    }

