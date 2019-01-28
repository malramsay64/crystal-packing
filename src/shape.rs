// shape.rs
//
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//
use approx::{AbsDiffEq, RelativeEq};
use nalgebra as na;
use nalgebra::{IsometryMatrix2, Point2};
use std::f64::consts::PI;

trait Intersect {
    fn intersects(&self, other: &Self) -> bool;
}

#[derive(Clone, Copy, PartialEq)]
pub struct Line {
    pub start: Point2<f64>,
    pub end: Point2<f64>,
}

impl Intersect for Line {
    /// Determine whether two line segments intersect
    ///
    /// This calculates whether two lines intersect at a point contained within each line segment
    /// see [this](https://en.wikipedia.org/wiki/Intersection_%28Euclidean_geometry%29#Two_line_segments)
    /// Wikipedia article for more information on the algorithm used for this calculation.
    ///
    fn intersects(&self, other: &Self) -> bool {
        // Also see below links for other implementations of this algorithm
        // - https://github.com/georust/geo/blob/96c7846d703a74f59ba68e68929415cbce4a68d9/geo/src/algorithm/intersects.rs#L142
        // - https://github.com/brandonxiang/geojson-python-utils/blob/33b4c00c6cf27921fb296052d0c0341bd6ca1af2/geojson_utils.py
        // - http://www.kevlindev.com/gui/math/intersection/Intersection.js
        //
        let u_b = other.dy() * self.dx() - other.dx() * self.dy();
        // Where u_b == 0 the two lines are parallel. In this case we don't need any further checks
        // since we are only concerned with lines that cross, parralel is fine.
        if u_b == 0. {
            return false;
        }

        let ua_t = other.dx() * (self.start.y - other.start.y)
            - other.dy() * (self.start.x - other.start.x);
        let ub_t =
            self.dx() * (self.start.y - other.start.y) - self.dy() * (self.start.x - other.start.x);

        let ua = ua_t / u_b;
        let ub = ub_t / u_b;
        // Should the points ua, ub both lie on the interval [0, 1] the lines intersect.
        if 0. <= ua && ua <= 1. && 0. <= ub && ub <= 1. {
            return true;
        }
        false
    }
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

impl std::fmt::Debug for Line {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Line {{ ({:.5}, {:.5}), ({:.5}, {:.5}) }}",
            self.start[0], self.start[1], self.end[0], self.end[1]
        )
    }
}

impl std::ops::Mul<IsometryMatrix2<f64>> for Line {
    type Output = Self;

    fn mul(self, rhs: IsometryMatrix2<f64>) -> Self::Output {
        Self {
            start: rhs * self.start,
            end: rhs * self.end,
        }
    }
}

impl Line {
    pub fn new(start: (f64, f64), end: (f64, f64)) -> Self {
        Self {
            start: Point2::new(start.0, start.1),
            end: Point2::new(end.0, end.1),
        }
    }

    /// The diffence in the x values over the line.
    pub fn dx(&self) -> f64 {
        self.end.x - self.start.x
    }

    /// The diffence in the y values over the line.
    pub fn dy(&self) -> f64 {
        self.end.y - self.start.y
    }
}

#[cfg(test)]
mod line_tests {
    use super::*;

    #[test]
    fn new() {
        let line = Line::new((1., 0.), (0., 1.));
        assert_eq!(line.start, Point2::new(1., 0.));
        assert_eq!(line.end, Point2::new(0., 1.));
    }

    #[test]
    fn intersects_radial() -> Result<(), String> {
        // Testing lines to and from the same point don't intersect
        // Using values of -1, 0, 1 for the x and y axes
        let values: Vec<f64> = vec![-1., 0., 1.];
        let points: Vec<(f64, f64)> = values
            .iter()
            .zip(values.iter())
            .map(|(a, b)| (*a, *b))
            .collect();

        let mut result = Ok(());
        for start1 in points.iter() {
            for start2 in points.iter() {
                let l1 = Line::new(*start1, (0., 0.));
                let l2 = Line::new(*start2, (0., 0.));
                if l1.intersects(&l2) {
                    println!("{:?} {:?}", start1, start2);
                    result = Err(String::from(""));
                }
            }
        }
        result
    }

    #[test]
    fn isometry_matrix_mul() {
        let ident: IsometryMatrix2<f64> = IsometryMatrix2::identity();
        let line = Line::new((1., 1.), (0., 0.));
        assert_eq!(line * ident, line);

        let trans: IsometryMatrix2<f64> = IsometryMatrix2::new(na::Matrix2x1::new(1., 1.), 0.);
        assert_eq!(line * trans, Line::new((2., 2.), (1., 1.)));
    }

    //
    // +-------------------|-------------------+
    // |                   |                   |
    // |                   + 1                 |
    // |                   |                   |
    // |                   |                   |
    // |                   |                   |
    // |         p3        | p1      p2         |
    // |---------+---------+---------+---------|
    // |        -1         |         1         |
    // |                   |                   |
    // |                   |                   |
    // |          p5       |                   |
    // |         +      p4 + -1                |
    // |                   |                   |
    // |                   |                   |
    // +-------------------|-------------------+
    //
    #[test]
    fn intersects() {
        let line1 = Line::new((-1., 0.), (0., -1.));
        let line2 = Line::new((-1., -1.), (0., 0.));
        assert!(line1.intersects(&line2));
        assert!(line2.intersects(&line1));

        let line3 = Line::new((-2., -1.), (1., 0.));
        assert!(line2.intersects(&line3));
        assert!(line3.intersects(&line2));
        assert!(line1.intersects(&line3));
        assert!(line3.intersects(&line1));
    }

}

/// Puts an abstract shape object in a physical space
///
/// This matches a Shape to a transformation, placing it in 2D space.
#[derive(PartialEq)]
pub struct ShapeInstance<'a> {
    pub shape: &'a RadialShape,
    pub isometry: IsometryMatrix2<f64>,
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

    pub fn intersects(&self, other: &ShapeInstance) -> bool {
        // Same shape and same isometry => they intersect
        if self == other {
            return true;
        }
        for line_a in self.lines() {
            for line_b in other.lines() {
                if line_a.intersects(&line_b) {
                    return true;
                }
            }
        }
        false
    }
}

impl<'a> std::fmt::Debug for ShapeInstance<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "ShapeInstance{{ {:?} }}", self.lines())
    }
}

#[cfg(test)]
mod shape_instance_tests {
    use super::*;

    #[test]
    fn lines() {
        let shape = RadialShape {
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
        let shape = RadialShape {
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
        let shape = RadialShape {
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

    //
    // +-------------------|-------------------+
    // |                   |                   |
    // |          +        + 1       +         |
    // |           \       |       /           |
    // |             \     |     /             |
    // |               \   |   /               |
    // |                 \ | /                 |
    // |---------+---------+---------+---------|
    // |        -1       / | \       1         |
    // |               /   |   \               |
    // |             /     |     \             |
    // |           /       |       \           |
    // |         +         + -1      +         |
    // |                   |                   |
    // |                   |                   |
    // +-------------------|-------------------+
    //
    #[test]
    fn intersects() {
        let shape = RadialShape {
            name: String::from("Square"),
            radial_points: vec![1., 1., 1., 1.],
            rotational_symmetries: 4,
            mirrors: 4,
        };
        let shape_i1 = ShapeInstance {
            shape: &shape,
            isometry: IsometryMatrix2::new(na::Vector2::new(1., 0.), 0.),
        };
        assert!(shape_i1.intersects(&shape_i1));

        let shape_i2 = ShapeInstance {
            shape: &shape,
            isometry: IsometryMatrix2::new(na::Vector2::new(-1.001, 0.), 0.),
        };
        assert!(!shape_i1.intersects(&shape_i2));

        let shape_i3 = ShapeInstance {
            shape: &shape,
            isometry: IsometryMatrix2::new(na::Vector2::new(0., 0.), PI / 4.),
        };
        assert!(shape_i1.intersects(&shape_i3));
        assert!(shape_i2.intersects(&shape_i3));
    }
}
