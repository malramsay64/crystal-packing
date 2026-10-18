//
// geometry.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Small, copyable 2D vector and matrix types used throughout the hard shape packing engine.
//!
//! These are deliberately minimal; the hot loops of the overlap checks only need a handful of
//! operations, and keeping the types simple makes the generated code (native and wasm) easy for
//! the compiler to optimise. Trigonometric functions go through `libm` so results are bit-identical
//! across platforms.

use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Vec2 {
    pub const ZERO: Vec2 = Vec2 { x: 0., y: 0. };

    pub const fn new(x: f64, y: f64) -> Self {
        Vec2 { x, y }
    }

    /// A unit vector pointing in the direction `angle` measured anticlockwise from the x axis.
    pub fn from_angle(angle: f64) -> Self {
        let (s, c) = libm::sincos(angle);
        Vec2::new(c, s)
    }

    pub fn dot(self, other: Vec2) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// The z component of the 3D cross product
    pub fn cross(self, other: Vec2) -> f64 {
        self.x * other.y - self.y * other.x
    }

    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    pub fn norm(self) -> f64 {
        libm::sqrt(self.norm_sq())
    }

    /// The vector rotated by 90 degrees anticlockwise
    pub fn perp(self) -> Vec2 {
        Vec2::new(-self.y, self.x)
    }

    pub fn angle(self) -> f64 {
        libm::atan2(self.y, self.x)
    }

    pub fn normalised(self) -> Vec2 {
        let n = self.norm();
        Vec2::new(self.x / n, self.y / n)
    }
}

impl Add for Vec2 {
    type Output = Vec2;
    fn add(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl AddAssign for Vec2 {
    fn add_assign(&mut self, rhs: Vec2) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl Sub for Vec2 {
    type Output = Vec2;
    fn sub(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl SubAssign for Vec2 {
    fn sub_assign(&mut self, rhs: Vec2) {
        self.x -= rhs.x;
        self.y -= rhs.y;
    }
}

impl Neg for Vec2 {
    type Output = Vec2;
    fn neg(self) -> Vec2 {
        Vec2::new(-self.x, -self.y)
    }
}

impl Mul<f64> for Vec2 {
    type Output = Vec2;
    fn mul(self, rhs: f64) -> Vec2 {
        Vec2::new(self.x * rhs, self.y * rhs)
    }
}

/// A 2x2 matrix stored in row major order
///
/// ```
/// use crystal_packing::geometry::{Mat2, Vec2};
/// let r = Mat2::rotation(std::f64::consts::FRAC_PI_2);
/// let v = r * Vec2::new(1., 0.);
/// assert!((v.x).abs() < 1e-15 && (v.y - 1.).abs() < 1e-15);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Mat2 {
    pub m: [[f64; 2]; 2],
}

impl Default for Mat2 {
    fn default() -> Self {
        Mat2::IDENTITY
    }
}

impl Mat2 {
    pub const IDENTITY: Mat2 = Mat2 {
        m: [[1., 0.], [0., 1.]],
    };

    pub const fn new(m00: f64, m01: f64, m10: f64, m11: f64) -> Self {
        Mat2 {
            m: [[m00, m01], [m10, m11]],
        }
    }

    /// Construct a matrix from two column vectors
    pub fn from_columns(a: Vec2, b: Vec2) -> Self {
        Mat2::new(a.x, b.x, a.y, b.y)
    }

    /// Anticlockwise rotation by `angle` radians
    pub fn rotation(angle: f64) -> Self {
        let (s, c) = libm::sincos(angle);
        Mat2::new(c, -s, s, c)
    }

    /// Reflection across the line through the origin at `angle` radians from the x axis
    pub fn reflection(angle: f64) -> Self {
        let (s, c) = libm::sincos(2. * angle);
        Mat2::new(c, s, s, -c)
    }

    pub fn column(&self, index: usize) -> Vec2 {
        Vec2::new(self.m[0][index], self.m[1][index])
    }

    pub fn det(&self) -> f64 {
        self.m[0][0] * self.m[1][1] - self.m[0][1] * self.m[1][0]
    }

    pub fn inverse(&self) -> Mat2 {
        let d = self.det();
        Mat2::new(
            self.m[1][1] / d,
            -self.m[0][1] / d,
            -self.m[1][0] / d,
            self.m[0][0] / d,
        )
    }

    pub fn transpose(&self) -> Mat2 {
        Mat2::new(self.m[0][0], self.m[1][0], self.m[0][1], self.m[1][1])
    }

    pub fn scale(&self, factor: f64) -> Mat2 {
        Mat2::new(
            self.m[0][0] * factor,
            self.m[0][1] * factor,
            self.m[1][0] * factor,
            self.m[1][1] * factor,
        )
    }

    /// The largest absolute difference between the elements of two matrices
    pub fn max_abs_diff(&self, other: &Mat2) -> f64 {
        let mut max: f64 = 0.;
        for i in 0..2 {
            for j in 0..2 {
                max = max.max((self.m[i][j] - other.m[i][j]).abs());
            }
        }
        max
    }
}

impl Mul<Vec2> for Mat2 {
    type Output = Vec2;
    fn mul(self, v: Vec2) -> Vec2 {
        Vec2::new(
            self.m[0][0] * v.x + self.m[0][1] * v.y,
            self.m[1][0] * v.x + self.m[1][1] * v.y,
        )
    }
}

impl Mul<Mat2> for Mat2 {
    type Output = Mat2;
    fn mul(self, o: Mat2) -> Mat2 {
        Mat2::new(
            self.m[0][0] * o.m[0][0] + self.m[0][1] * o.m[1][0],
            self.m[0][0] * o.m[0][1] + self.m[0][1] * o.m[1][1],
            self.m[1][0] * o.m[0][0] + self.m[1][1] * o.m[1][0],
            self.m[1][0] * o.m[0][1] + self.m[1][1] * o.m[1][1],
        )
    }
}

/// Wrap an angle into the range [0, period)
pub fn wrap_angle(angle: f64, period: f64) -> f64 {
    let r = angle % period;
    if r < 0. {
        r + period
    } else {
        r
    }
}

/// The smallest distance between two angles on a circle with the given period
pub fn angle_distance(a: f64, b: f64, period: f64) -> f64 {
    let d = wrap_angle(a - b, period);
    d.min(period - d)
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn inverse() {
        let m = Mat2::new(2., 1., 0.5, 3.);
        let i = m * m.inverse();
        assert!(i.max_abs_diff(&Mat2::IDENTITY) < 1e-14);
    }

    #[test]
    fn reflection_is_involution() {
        for angle in [0., 0.3, PI / 4., 2.] {
            let r = Mat2::reflection(angle);
            assert!((r * r).max_abs_diff(&Mat2::IDENTITY) < 1e-14);
            assert_abs_diff_eq!(r.det(), -1., epsilon = 1e-14);
            // Points on the mirror line are unchanged
            let p = Vec2::from_angle(angle);
            let q = r * p;
            assert_abs_diff_eq!(q.x, p.x, epsilon = 1e-14);
            assert_abs_diff_eq!(q.y, p.y, epsilon = 1e-14);
        }
    }

    #[test]
    fn wrap() {
        assert_abs_diff_eq!(wrap_angle(-0.5, 2.), 1.5);
        assert_abs_diff_eq!(wrap_angle(4.5, 2.), 0.5);
        assert_abs_diff_eq!(angle_distance(0.1, 1.9, 2.), 0.2, epsilon = 1e-14);
    }
}
