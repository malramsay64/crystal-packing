//
// shape.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Shape {
    pub name: String,
    pub radial_points: Vec<f64>,
    pub rotational_symmetries: u64,
    pub mirrors: u64,
}

impl<'a> IntoIterator for &'a Shape {
    type Item = (f64, f64);
    type IntoIter = ShapeIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        ShapeIter::new(&self)
    }
}

impl Shape {
    pub fn area(&self) -> f64 {
        // This is the sine of the angle between each point, this is used for every calculation
        // so pre-calculate here.
        let angle_term: f64 = f64::sin(2. * PI / self.radial_points.len() as f64);

        self.into_iter()
            .map(|(a, b)| 0.5 * angle_term * a * b)
            .sum()
    }

    pub fn max_radius(&self) -> f64 {
        return self
            .radial_points
            .iter()
            .cloned()
            // The f64 type doesn't have complete ordering because of Nan and Inf, so the
            // standard min/max comparators don't work. Instead we use the f64::max which ignores
            // the NAN and max values.
            .fold(std::f64::MIN, f64::max);
    }
}

/// Iterator for the Shape class
///
/// This iterates over all the line segments comprising a shape.
pub struct ShapeIter<'a> {
    shape: &'a Shape,
    index: usize,
    len: usize,
}

impl<'a> ShapeIter<'a> {
    fn new(shape: &'a Shape) -> Self {
        Self {
            shape,
            index: 0,
            len: shape.radial_points.len(),
        }
    }
}

impl<'a> Iterator for ShapeIter<'a> {
    type Item = (f64, f64);

    fn next(&mut self) -> Option<(f64, f64)> {
        if self.index >= self.len {
            return None;
        }
        let result = Some((
            self.shape.radial_points[self.index],
            self.shape.radial_points[(self.index + 1) % self.len],
        ));
        self.index += 1;
        result
    }
}

#[cfg(test)]
mod shape_tests {
    use super::*;

    fn create_square() -> Shape {
        Shape {
            name: String::from("Square"),
            radial_points: vec![1., 1., 1., 1.],
            rotational_symmetries: 4,
            mirrors: 4,
        }
    }

    #[test]
    fn init() {
        let square = create_square();
        assert_eq!(square.name, "Square");
        assert_eq!(square.radial_points, vec![1., 1., 1., 1.]);
        assert_eq!(square.rotational_symmetries, 4);
        assert_eq!(square.mirrors, 4);
    }

    #[test]
    fn area() {
        let square = create_square();
        assert_eq!(square.area(), 2.);
    }

    #[test]
    fn max_radius() {
        let shape = Shape {
            name: String::from("iter_test"),
            radial_points: vec![1., 2., 3., 4.],
            rotational_symmetries: 1,
            mirrors: 0,
        };
        assert_eq!(shape.max_radius(), 4.);
        assert_eq!(shape.max_radius(), 4.);
    }

    #[test]
    fn iter_values() {
        let shape = Shape {
            name: String::from("iter_test"),
            radial_points: vec![1., 2., 3., 4.],
            rotational_symmetries: 1,
            mirrors: 0,
        };
        let manual = vec![(1., 2.), (2., 3.), (3., 4.), (4., 1.)];
        assert_eq!(shape.into_iter().collect::<Vec<(f64, f64)>>(), manual);
    }
}
