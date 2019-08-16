//
// transform.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use log::warn;
use std::ops::Mul;

use approx::AbsDiffEq;
use nalgebra::{Matrix3, Vector2, Vector3};
use serde::{Deserialize, Serialize};

/// Perform coordinate tranforms on a point in space
///
/// This allows for defining a transformation of a point in space and allow for translations,
/// rotations, mirror planes, and any combination of the above.
///
/// Rotations and translations are the most common method with a simple construction, and when they
/// are the only transformations present this is known as an Isometric transform, or Isometry.
/// Which is a rotation, followed by a translation
///
/// ```
/// use packing::Transform2;
/// // Create a transform
/// let t = Transform2::new(0., (1., 1.))
/// ```
///
/// The order of rotation, followed by translation is followed in the initialisation, with the
/// angular rotation being the first argument, and the translation being the second argument.
///
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Transform2(Matrix3<f64>);

impl From<Matrix3<f64>> for Transform2 {
    fn from(matrix: Matrix3<f64>) -> Self {
        Self(matrix)
    }
}

impl Into<Matrix3<f64>> for Transform2 {
    fn into(self) -> Matrix3<f64> {
        self.0
    }
}

impl AbsDiffEq for Transform2 {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

binop_impl_all!(
    Mul, mul;
    self: Transform2, rhs: Vector2<f64>, Output = Vector2<f64>;
    [ref ref] => {
        let result = self.0 * Vector3::new(rhs.x, rhs.y, 1.);
        result.xy()
    };
);

binop_impl_all!(
    Mul, mul;
    self: Transform2, rhs: Transform2, Output = Transform2;
    [ref ref] => {
        Transform2::from(self.0 * rhs.0)
    };
);

impl Transform2 {
    pub fn new(rotation: f64, translation: (f64, f64)) -> Transform2 {
        Transform2::from(Matrix3::new(
            rotation.cos(),
            -rotation.sin(),
            translation.0,
            rotation.sin(),
            rotation.cos(),
            translation.1,
            0.,
            0.,
            1.,
        ))
    }

    pub fn identity() -> Self {
        Self(Matrix3::identity())
    }

    pub fn position(&self) -> Vector2<f64> {
        Vector2::new(self.0[(0, 2)], self.0[(1, 2)])
    }

    pub fn set_position(&self, position: Vector2<f64>) -> Transform2 {
        let mut transform = *self;
        transform.0[(0, 2)] = position.x;
        transform.0[(1, 2)] = position.y;
        transform
    }

    pub fn periodic(&self, period: f64, offset: f64) -> Transform2 {
        let mut position = self.position();
        position.x = (((position.x - offset) % period) + period) % period + offset;
        position.y = (((position.y - offset) % period) + period) % period + offset;
        self.set_position(position)
    }

    /// Convert the string representation of a symmetry operation to a vector.
    ///
    /// This converts the string representation of an operation to a Transform,
    /// extracting the rotation and translation components.
    ///
    /// ```
    /// use packing::{Transform2};
    /// let t2 = Transform2::from_operations("-x, y");
    /// ```
    ///
    pub fn from_operations(sym_ops: &str) -> Result<Transform2, &'static str> {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();

        match operations.len() {
            x if x < 2 => return Err("Not enough dimensions in input"),
            x if x > 2 => return Err("Too many dimensions in input"),
            _ => (),
        }

        let mut transform: Matrix3<f64> = Matrix3::zeros();

        for (index, op) in operations.iter().enumerate() {
            let mut sign = 1.;
            let mut constant = 0.;
            let mut operator: Option<char> = None;
            for c in op.chars() {
                match c {
                    'x' => {
                        transform[(index, 0)] = sign;
                        sign = 1.;
                    }
                    'y' => {
                        transform[(index, 1)] = sign;
                        sign = 1.;
                    }
                    '*' | '/' => {
                        operator = Some(c);
                    }
                    '-' => {
                        sign = -1.;
                    }
                    // This matches all digits from 0 to 9
                    '0'...'9' => {
                        let val = c.to_string().parse::<u64>().unwrap() as f64;
                        // Is there an operator defined, i.e. is this the first digit
                        constant = match operator {
                            Some(op) if op == '/' => sign * constant / val,
                            Some(op) if op == '*' => sign * constant / val,
                            Some(_) => 0.,
                            None => sign * val,
                        };
                        // Reset values
                        operator = None;
                        sign = 1.
                    }
                    ' ' | '+' => (),
                    // Default is do nothing (shouldn't encounter this at all)
                    x => {
                        warn!("Found '{}'", x);
                        return Err("Invalid value found");
                    }
                };
            }
            transform[(index, 2)] = constant;
        }
        Ok(Transform2::from(transform))
    }
}

#[cfg(test)]
mod test {
    use approx::{abs_diff_eq, assert_abs_diff_eq};
    use std::f64;

    use super::*;
    use quickcheck_macros::quickcheck;

    /// Init new with zeros should be the same as identity
    #[test]
    fn new_zeros() {
        assert_eq!(Transform2::new(0., (0., 0.)), Transform2::identity())
    }

    /// Ensure init with translation is same as returned
    #[quickcheck]
    fn new_translation(x: f64, y: f64) -> bool {
        let t = Transform2::new(0., (x, y));
        t.position() == Vector2::new(x, y)
    }

    /// Check get position is same as set
    #[quickcheck]
    fn set_get_translation(x: f64, y: f64) -> bool {
        let pos = Vector2::new(x, y);
        let t = Transform2::identity().set_position(pos);
        t.position() == pos
    }

    /// Transformation by identity matrix gives same result
    #[quickcheck]
    fn identity_transform(x: f64, y: f64) -> bool {
        let identity = Transform2::identity();
        let point = Vector2::new(x, y);
        identity * point == point
    }

    /// Rotation keeps point same distance from origin
    #[quickcheck]
    fn rotation_length(angle: f64) -> bool {
        let t = Transform2::new(angle, (0., 0.));
        let point = Vector2::new(1., 0.);
        abs_diff_eq!((t * point).norm(), 1.)
    }

    /// Translation from origin puts point in same location
    #[quickcheck]
    fn translation_from_origin(x: f64, y: f64) -> bool {
        let t = Transform2::new(0., (x, y));
        let origin = Vector2::zeros();
        t * origin == t.position()
    }

    /// Rotation followed by Translation is in circle about translation
    #[quickcheck]
    fn rotation_translation(rotation: f64, translation: (f64, f64)) -> bool {
        let t = Transform2::new(rotation, translation);
        let point = Vector2::new(1., 0.);
        abs_diff_eq!(
            (t * point - t.position()).norm(),
            1.,
            // Large transformations will have larger errors which this accounts for
            epsilon = t.position().norm() * Transform2::default_epsilon()
        )
    }

    #[test]
    fn transform() {
        let t = Transform2::new(f64::consts::PI / 2., (1., 1.));

        let point = Vector2::new(0.2, 0.2);
        assert_eq!(t * point, Vector2::new(0.8, 1.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Vector2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Vector2::new(0.1, 0.2));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Vector2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Vector2::new(-0.1, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Vector2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Vector2::new(0.6, -0.2));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Vector2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Vector2::new(-0.4, -0.2));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Vector2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Vector2::new(-0.2, 0.));
    }

    #[test]
    #[should_panic]
    fn parse_operation_z() {
        let input = String::from("(z, y)");
        Transform2::from_operations(&input).unwrap();
    }

    #[test]
    #[should_panic]
    fn parse_operation_3_inputs() {
        let input = String::from("(x, y, z)");
        Transform2::from_operations(&input).unwrap();
    }

    #[test]
    #[should_panic]
    fn parse_operation_1_input() {
        let input = String::from("(x)");
        Transform2::from_operations(&input).unwrap();
    }
}
