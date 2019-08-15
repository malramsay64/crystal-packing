//
// molecular_shape2.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;
use std::{fmt, slice, vec};

use itertools::{iproduct, Itertools};
use serde::{Deserialize, Serialize};

use super::{Atom2, Transform2};
use crate::traits::{Intersect, Shape};

/// A shape defined by a collection of Atoms
///
/// This is a shape comprised of a series of circles which each have a position and radius.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MolecularShape2 {
    pub name: String,
    pub items: Vec<Atom2>,
}

impl<'a> IntoIterator for &'a MolecularShape2 {
    type Item = &'a Atom2;
    type IntoIter = slice::Iter<'a, Atom2>;

    fn into_iter(self) -> Self::IntoIter {
        self.items.iter()
    }
}

impl Intersect for MolecularShape2 {
    fn intersects(&self, other: &Self) -> bool {
        iproduct!(self.items.iter(), other.items.iter()).any(|(s, o)| s.intersects(o))
    }
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
}

impl Shape for MolecularShape2 {
    type Component = Atom2;

    fn score(&self, other: &Self) -> Result<f64, &'static str> {
        if self.intersects(other) {
            Err("Shape intersects")
        } else {
            Ok(self.area())
        }
    }

    fn enclosing_radius(&self) -> f64 {
        self.items
            .iter()
            .map(|p| p.position.norm() + p.radius)
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

impl fmt::Display for MolecularShape2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "MolShape {{ ")?;
        for item in self.items.iter() {
            write!(f, "{},", item)?;
        }
        write!(f, " }}")
    }
}

impl MolecularShape2 {
    fn overlap_area(r: f64, d: f64) -> f64 {
        r.powi(2) * f64::acos(d / r) - d * f64::sqrt(r.powi(2) - d.powi(2))
    }

    fn circle_overlap(a1: &Atom2, a2: &Atom2) -> f64 {
        let distance = (a1.position - a2.position).norm();
        // There is some overlap between the circles which needs to be calculated
        if distance < a1.radius + a2.radius {
            let d1 = (distance.powi(2) + a1.radius.powi(2) - a2.radius.powi(2)) / (2. * distance);
            let d2 = (distance.powi(2) + a2.radius.powi(2) - a1.radius.powi(2)) / (2. * distance);
            Self::overlap_area(a1.radius, d1) + Self::overlap_area(a2.radius, d2)
        } else {
            0.
        }
    }

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
                Atom2::new(0., -2. / 3. * distance * f64::cos(angle / 2.), 1.),
                Atom2::new(
                    -distance * f64::sin(angle / 2.),
                    1. / 3. * distance * f64::cos(angle / 2.),
                    radius,
                ),
                Atom2::new(
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
            items: vec![Atom2::new(0., 0., 1.)],
        }
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector2;

    use super::*;

    #[test]
    fn overlap_area_test() {
        assert_abs_diff_eq!(MolecularShape2::overlap_area(1., 1.), 0.);
    }

    #[test]
    fn circle_overlaps_test() {
        let a1 = Atom2::new(0., 0., 1.);
        let a2 = Atom2::new(2., 0., 1.);
        assert_abs_diff_eq!(MolecularShape2::circle_overlap(&a1, &a2), 0.);

        for i in 0..10 {
            let distance = f64::from(i + 1) / 10. * 2.;
            let a1 = Atom2::new(0., 0., 1.);
            let a2 = Atom2::new(distance, 0., 1.);
            // A known algorithm for confirming the area is calculated correctly, as found on
            // http://mathworld.wolfram.com/Circle-CircleIntersection.html
            let area = 2. * MolecularShape2::overlap_area(1., distance / 2.);
            assert_abs_diff_eq!(
                MolecularShape2::circle_overlap(&a1, &a2),
                area,
                epsilon = 1e-7
            );
        }
    }

    #[test]
    fn from_trimer_test() {
        let shape = MolecularShape2::from_trimer(1., PI, 1.);
        assert_eq!(shape.items.len(), 3);

        assert_abs_diff_eq!(shape.items[0].position, Vector2::new(0., 0.));
        assert_abs_diff_eq!(shape.items[1].position, Vector2::new(-1., 0.));
        assert_abs_diff_eq!(shape.items[2].position, Vector2::new(1., 0.));

        let shape = MolecularShape2::from_trimer(0.637_556, 2. * PI / 3., 1.);
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
    fn area_test() {
        let shape = MolecularShape2::from_trimer(1., PI, 2.);
        assert_abs_diff_eq!(shape.area(), 3. * PI);

        let shape = MolecularShape2::from_trimer(0.637_556, 2. * PI / 3., 1.);
        println!("{}", shape.area());
        assert!(shape.area() > 0.);
    }

    #[test]
    fn intersection() {
        let mol = MolecularShape2::circle();
        let transform = Transform2::new(1.0, 1.0, 0.);
        assert!(mol.intersects(&mol.transform(&transform)));
    }

    #[test]
    fn no_intersection_edge() {
        let mol = MolecularShape2::circle();
        let transform = Transform2::new(2., 0., 0.);
        assert!(!mol.intersects(&mol.transform(&transform)));
    }

    #[test]
    fn self_intersection() {
        let mol = MolecularShape2::circle();
        let transform = Transform2::new(0., 0., 0.);
        assert!(mol.intersects(&mol.transform(&transform)));
    }

    #[test]
    fn no_intersection() {
        let mol = MolecularShape2::circle();
        let transform = Transform2::new(2.01, 0., 0.);
        assert!(!mol.intersects(&mol.transform(&transform)));
    }

    #[test]
    fn no_intersection_corner() {
        let mol = MolecularShape2::circle();
        let transform = Transform2::new(2., 2., 0.);
        assert!(!mol.intersects(&mol.transform(&transform)));
    }

}
