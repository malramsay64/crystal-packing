//
// shape.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//
use approx::{AbsDiffEq, RelativeEq};
use nalgebra as na;
use nalgebra::{IsometryMatrix2, Point2};
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

    pub fn radial_step(&self) -> f64 {
        2. * PI / self.radial_points.len() as f64
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

#[derive(Debug, Clone, PartialEq)]
pub struct Line {
    pub start: Point2<f64>,
    pub end: Point2<f64>,
}

impl AbsDiffEq for Line {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        std::f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.start.abs_diff_eq(&other.start, epsilon) && self.end.abs_diff_eq(&other.end, epsilon)
    }
}

impl RelativeEq for Line {
    fn default_max_relative() -> Self::Epsilon {
        std::f64::EPSILON
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.start.relative_eq(&other.start, epsilon, max_relative)
            && self.end.relative_eq(&other.end, epsilon, max_relative)
    }
}

impl Line {
    pub fn new(start: Point2<f64>, end: Point2<f64>) -> Self {
        Self { start, end }
    }
}

pub struct ShapeInstance<'a> {
    shape: &'a Shape,
    isometry: IsometryMatrix2<f64>,
}

impl<'a> ShapeInstance<'a> {
    fn radii_to_line(&self, index: usize, radii: (f64, f64)) -> Line {
        let radial_step = self.shape.radial_step();
        let angle = index as f64 * radial_step;
        let (r1, r2) = radii;
        let start = na::Point2::new(r1 * f64::sin(angle), r1 * f64::cos(angle));
        let end = na::Point2::new(
            r2 * f64::sin(angle + radial_step),
            r2 * f64::cos(angle + radial_step),
        );

        Line {
            start: self.isometry * start,
            end: self.isometry * end,
        }
    }

    pub fn lines(&self) -> Vec<Line> {
        self.shape
            .into_iter()
            .enumerate()
            .map(|(index, r)| self.radii_to_line(index, r))
            .collect()
    }
}

#[cfg(test)]
mod shape_instance_tests {
    use super::*;

    #[test]
    fn lines() {
        let shape = Shape {
            name: String::from("Square"),
            radial_points: vec![1., 1., 1., 1.],
            rotational_symmetries: 4,
            mirrors: 4,
        };
        let shape_i = ShapeInstance {
            shape: &shape,
            isometry: IsometryMatrix2::identity(),
        };
        let expected_vec = vec![
            Line {
                start: Point2::new(0., 1.),
                end: Point2::new(1., 0.),
            },
            Line {
                start: Point2::new(1., 0.),
                end: Point2::new(0., -1.),
            },
            Line {
                start: Point2::new(0., -1.),
                end: Point2::new(-1., 0.),
            },
            Line {
                start: Point2::new(-1., 0.),
                end: Point2::new(0., 1.),
            },
        ];
        for (index, (expected, result)) in
            shape_i.lines().iter().zip(expected_vec.iter()).enumerate()
        {
            println!("{}", index);
            assert_abs_diff_eq!(*expected, *result, epsilon = 1e-8);
        }
    }

    #[test]
    fn lines_translate() {
        let shape = Shape {
            name: String::from("Square"),
            radial_points: vec![1., 1., 1., 1.],
            rotational_symmetries: 4,
            mirrors: 4,
        };
        let shape_i = ShapeInstance {
            shape: &shape,
            isometry: IsometryMatrix2::new(na::Vector2::new(-1., 0.), 0.),
        };
        let expected_vec = vec![
            Line {
                start: Point2::new(-1., 1.),
                end: Point2::new(0., 0.),
            },
            Line {
                start: Point2::new(0., 0.),
                end: Point2::new(-1., -1.),
            },
            Line {
                start: Point2::new(-1., -1.),
                end: Point2::new(-2., 0.),
            },
            Line {
                start: Point2::new(-2., 0.),
                end: Point2::new(-1., 1.),
            },
        ];
        for (index, (expected, result)) in
            shape_i.lines().iter().zip(expected_vec.iter()).enumerate()
        {
            println!("{}", index);
            assert_abs_diff_eq!(*expected, *result, epsilon = 1e-8);
        }
    }

    #[test]
    fn lines_translate_rotate() {
        let shape = Shape {
            name: String::from("Square"),
            radial_points: vec![1., 1., 1., 1.],
            rotational_symmetries: 4,
            mirrors: 4,
        };
        let shape_i = ShapeInstance {
            shape: &shape,
            isometry: IsometryMatrix2::new(na::Vector2::new(-1., 0.), PI),
        };
        let expected_vec = vec![
            Line {
                start: Point2::new(-1., -1.),
                end: Point2::new(-2., 0.),
            },
            Line {
                start: Point2::new(-2., 0.),
                end: Point2::new(-1., 1.),
            },
            Line {
                start: Point2::new(-1., 1.),
                end: Point2::new(0., 0.),
            },
            Line {
                start: Point2::new(0., 0.),
                end: Point2::new(-1., -1.),
            },
        ];
        for (index, (expected, result)) in
            shape_i.lines().iter().zip(expected_vec.iter()).enumerate()
        {
            println!("{}", index);
            assert_abs_diff_eq!(*expected, *result, epsilon = 1e-8);
        }
    }

}
