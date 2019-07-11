//
// line_shape2.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;
use std::fmt;
use std::slice;
use std::vec;

use itertools::{iproduct, Itertools};
use nalgebra::Point2;
use serde::{Deserialize, Serialize};

use super::Line2;
use crate::traits::{Intersect, Shape};
use crate::Transform2;

/// A Shape constructed from a collection of Lines
///
/// This defines a collection of lines, from one point to another which define the area enclosed by
/// a shape. It is assumed that the lines completely enclose an area, and that the enclosed area is
/// close to the origin.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LineShape {
    pub name: String,
    pub items: Vec<Line2>,
}

impl<'a> IntoIterator for &'a LineShape {
    type Item = &'a Line2;
    type IntoIter = slice::Iter<'a, Line2>;

    fn into_iter(self) -> Self::IntoIter {
        self.items.iter()
    }
}

impl fmt::Display for LineShape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "LineShape {{ {} }}", self.items.iter().format(", "))
    }
}

impl Intersect for LineShape {
    /// Check whether this shape intersects with another shape
    ///
    /// A ShapeInstance is considered to intersect with another when one of it's components
    /// intersects with a component of the other shape. For a square, there is an intersection
    /// when a line from one square crosses the other. Each component item of `self` is
    /// checked against `other`.
    ///
    fn intersects(&self, other: &Self) -> bool {
        // We want to compare every item of the current shape with every item of the other shape.
        iproduct!(self.iter(), other.iter()).any(|(s, o)| s.intersects(o))
    }

    fn area(&self) -> f64 {
        // This is the sine of the angle between each point, this is used for every calculation
        // so pre-calculate here.
        let angle_term: f64 = f64::sin(2. * PI / self.items.len() as f64);
        let zero = Point2::origin();
        self.iter()
            // Calculate the area of the of triangle made by the line and the origin
            .map(|p| {
                0.5 * angle_term
                    * nalgebra::distance(&zero, &p.start)
                    * nalgebra::distance(&zero, &p.end)
            })
            .sum()
    }
}

impl Shape for LineShape {
    type Component = Line2;

    fn score(&self, other: &Self) -> Result<f64, &'static str> {
        if self.intersects(other) {
            Err("Shape intersects")
        } else {
            Ok(self.area())
        }
    }

    fn enclosing_radius(&self) -> f64 {
        self.iter()
            .map(|p| nalgebra::distance(&Point2::origin(), &p.start))
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

impl LineShape {
    /// Instantiate a LineShape from a collection of radial points
    ///
    /// The input is a Vector of points which are a radial distance from the origin, with the
    /// points separated by and equal angle. For example to create a Triangle, which is the shape
    /// with the fewest number of sides, we can run
    /// ```
    /// use packing::LineShape;
    /// let tri = LineShape::from_radial("Triangle", vec![1., 1., 1.]);
    /// ```
    /// More generally to create a regular polygon with an arbitrary number of sides
    /// ```
    /// use packing::LineShape;
    /// let sides = 10;
    /// let polygon = LineShape::from_radial("Polygon", vec![1.; sides]);
    /// ```
    ///
    pub fn from_radial(name: &str, points: Vec<f64>) -> Result<LineShape, &'static str> {
        if points.len() < 3 {
            return Err("The number of points provided is too few to create a 2D shape.");
        }
        let dtheta = 2. * PI / points.len() as f64;
        let mut items: Vec<Line2> = vec![];
        for (index, (r1, r2)) in points.iter().zip(points.iter().cycle().skip(1)).enumerate() {
            let angle = index as f64 * dtheta;
            items.push(Line2::new(
                (r1 * f64::sin(angle), r1 * f64::cos(angle)),
                (r2 * f64::sin(angle + dtheta), r2 * f64::cos(angle + dtheta)),
            ))
        }

        Ok(LineShape {
            name: String::from(name),
            items,
        })
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector2;

    use super::*;

    fn create_square() -> LineShape {
        LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap()
    }

    #[test]
    fn init() {
        let square = create_square();
        assert_eq!(square.name, "Square");
    }

    #[test]
    fn area() {
        let square = create_square();
        assert_abs_diff_eq!(square.area(), 2.);
    }

    #[test]
    fn max_radius() {
        let shape = LineShape::from_radial("iter_test", vec![1., 2., 3., 4.]).unwrap();
        assert_abs_diff_eq!(shape.enclosing_radius(), 4.);
        assert_abs_diff_eq!(shape.enclosing_radius(), 4.);
    }

    #[test]
    fn intersection() {
        let square = create_square();
        let transform = Transform2::new(Vector2::new(1., 1.), 0.);
        assert!(square.intersects(&square.transform(&transform)));
    }

    #[test]
    fn corner_no_intersection() {
        let square = create_square();
        let transform = Transform2::new(Vector2::new(2., 2.), 0.);
        assert!(!square.intersects(&square.transform(&transform)));
    }

    #[test]
    fn self_intersection() {
        let square = create_square();
        let transform = Transform2::new(Vector2::new(0., 0.), 0.);
        assert!(square.intersects(&square.transform(&transform)));
    }

    #[test]
    fn no_intersection() {
        let square = create_square();
        let transform = Transform2::new(Vector2::new(2.01, 2.01), 0.);
        assert!(!square.intersects(&square.transform(&transform)));
    }

}
