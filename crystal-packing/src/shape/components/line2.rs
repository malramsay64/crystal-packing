//
// line2.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fmt;

#[cfg(test)]
use approx::AbsDiffEq;
use nalgebra::Point2;
use serde::{Deserialize, Serialize};

use crate::traits::Intersect;

#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
pub struct Line2 {
    pub start: Point2<f64>,
    pub end: Point2<f64>,
}

impl Intersect for Line2 {
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
        // since we are only concerned with lines that cross, parallel is fine.
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
        if (0. ..=1.).contains(&ua) && (0. ..=1.).contains(&ub) {
            return true;
        }
        false
    }

    fn area(&self) -> f64 {
        // TODO Implement some area calculation being the area to the origin or to the y axis.
        0.
    }
}

#[cfg(test)]
impl AbsDiffEq for Line2 {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        std::f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.start.abs_diff_eq(&other.start, epsilon) && self.end.abs_diff_eq(&other.end, epsilon)
    }
}

impl fmt::Display for Line2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Line2 {{ ({:.5}, {:.5}), ({:.5}, {:.5}) }}",
            self.start.x, self.start.y, self.end.x, self.end.y
        )
    }
}

impl Line2 {
    pub fn new(start: (f64, f64), end: (f64, f64)) -> Self {
        Self {
            start: Point2::new(start.0, start.1),
            end: Point2::new(end.0, end.1),
        }
    }

    /// The difference in the x values over the line.
    pub fn dx(&self) -> f64 {
        self.end.x - self.start.x
    }

    /// The difference in the y values over the line.
    pub fn dy(&self) -> f64 {
        self.end.y - self.start.y
    }
}

#[cfg(test)]
mod test {
    use itertools::iproduct;

    use super::*;
    use crate::Transform2;

    #[test]
    fn new() {
        let line = Line2::new((1., 0.), (0., 1.));
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
        for (start1, start2) in iproduct!(points.iter(), points.iter()) {
            let l1 = Line2::new(*start1, (0., 0.));
            let l2 = Line2::new(*start2, (0., 0.));
            if l1.intersects(&l2) {
                println!("{:?} {:?}", start1, start2);
                result = Err(String::from(""));
            }
        }
        result
    }

    #[test]
    fn isometry_matrix_mul() {
        let ident: Transform2 = Transform2::identity();
        let line = Line2::new((1., 1.), (0., 0.));
        assert_eq!(line * ident, line);

        let trans: Transform2 = Transform2::new(0., (1., 1.));
        assert_eq!(line * trans, Line2::new((2., 2.), (1., 1.)));
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
        let line1 = Line2::new((-1., 0.), (0., -1.));
        let line2 = Line2::new((-1., -1.), (0., 0.));
        assert!(line1.intersects(&line2));
        assert!(line2.intersects(&line1));

        let line3 = Line2::new((-2., -1.), (1., 0.));
        assert!(line2.intersects(&line3));
        assert!(line3.intersects(&line2));
        assert!(line1.intersects(&line3));
        assert!(line3.intersects(&line1));
    }
}
