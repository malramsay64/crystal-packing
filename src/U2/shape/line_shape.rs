//
// line_shape.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;
use std::fmt;
use std::slice;
use std::vec;

use itertools::Itertools;
use nalgebra::Point2;

use super::{Line2, Transform2};
use crate::traits::Shape;

/// A Shape constructed from a collection of Lines
///
/// This defines a collection of lines, from one point to another which define the area enclosed by
/// a shape. It is assumed that the lines completely enclose an area, and that the enclosed area is
/// close to the origin.
#[derive(Debug, Clone, PartialEq)]
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

impl Shape for LineShape {
    type Component = Line2;
    type Transform = Transform2;

    fn area(&self) -> f64 {
        // This is the sine of the angle between each point, this is used for every calculation
        // so pre-calculate here.
        let angle_term: f64 = f64::sin(2. * PI / self.items.len() as f64);
        let zero = Point2::origin();
        self.iter()
            // Calculate the area of the of triangle made by the line and the origin
            .map(|p| 0.5 * angle_term * na::distance(&zero, &p.start) * na::distance(&zero, &p.end))
            .sum()
    }

    fn enclosing_radius(&self) -> f64 {
        self.iter()
            .map(|p| na::distance(&Point2::origin(), &p.start))
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

impl LineShape {
    /// Instantiate a LineShape from a collection of radial points
    ///
    /// The input is a Vector of points which are a radial distance from the origin, with the
    /// points separated by and equal angle. For example to create a Triangle, which is the shape
    /// with the fewest number of sides, we can run
    /// ```
    /// use packing::shape::LineShape;
    /// let tri = LineShape::from_radial("Triangle", vec![1., 1., 1.]);
    /// ```
    /// More generally to create a regular polygon with an arbitrary number of sides
    /// ```
    /// use packing::shape::LineShape;
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

}
