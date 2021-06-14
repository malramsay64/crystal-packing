//
// transform.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use anyhow::{bail, Error};
use std::ops::Mul;

#[cfg(test)]
use approx::AbsDiffEq;
use nalgebra::{Matrix3, Point2, Translation2};
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
/// use crystal_packing::Transform2;
/// // Create a transform
/// let t = Transform2::new(0., (1., 1.));
/// ```
///
/// The order of rotation, followed by translation is followed in the initialisation, with the
/// angular rotation being the first argument, and the translation being the second argument.
///
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Transform2(nalgebra::Transform2<f64>);

impl From<Matrix3<f64>> for Transform2 {
    fn from(matrix: Matrix3<f64>) -> Self {
        Self(nalgebra::Transform2::from_matrix_unchecked(matrix))
    }
}

// Since Matrix3 is a type defined in another crate, we run into the orphan rule where we can't
// define a trait for a type that we don't define.
#[allow(clippy::from_over_into)]
impl Into<Matrix3<f64>> for Transform2 {
    fn into(self) -> Matrix3<f64> {
        *self.0.matrix()
    }
}

#[cfg(test)]
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
    self: Transform2, rhs: Point2<f64>, Output = Point2<f64>;
    [ref ref] => {
        self.0 * rhs
    };
);

binop_impl_all!(
    Mul, mul;
    self: Transform2, rhs: Transform2, Output = Transform2;
    [ref ref] => {
        Transform2(self.0 * rhs.0)
    };
);

impl Transform2 {
    pub fn new(rotation: f64, translation: (f64, f64)) -> Transform2 {
        let translation = nalgebra::Translation2::new(translation.0, translation.1);
        let rotation = nalgebra::Rotation2::new(rotation);
        Transform2(nalgebra::Transform2::from_matrix_unchecked(
            nalgebra::IsometryMatrix2::from_parts(translation, rotation).to_homogeneous(),
        ))
    }

    pub fn identity() -> Self {
        Self(nalgebra::Transform2::identity())
    }

    pub fn position(&self) -> Point2<f64> {
        self.0 * Point2::origin()
    }

    pub fn get_translation(&self) -> Translation2<f64> {
        Translation2::new(self.0.matrix()[(0, 2)], self.0.matrix()[(1, 2)])
    }

    pub fn set_position(mut self, position: Point2<f64>) -> Transform2 {
        self.0[(0, 2)] = position.x;
        self.0[(1, 2)] = position.y;
        self
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
    /// use crystal_packing::{Transform2};
    /// let t2 = Transform2::from_operations("-x, y");
    /// ```
    ///
    pub fn from_operations(sym_ops: &str) -> Result<Transform2, Error> {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();

        match operations.len() {
            x if x < 2 => bail!("Not enough dimensions in input"),
            x if x > 2 => bail!("Too many dimensions in input"),
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
                    '0'..='9' => {
                        let val = c.to_string().parse::<u64>()? as f64;
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
                    x => bail!("Found invalid value: '{}'", x),
                };
            }
            transform[(index, 2)] = constant;
        }
        Ok(Transform2::from(transform))
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;
    use std::f64;

    use super::*;
    use nalgebra as na;
    use proptest::prelude::*;
    use proptest_attr_macro::proptest;

    /// Struct to aid in setting values for Proptest
    ///
    /// When testing the transformations, we want to be able to use property based testing,
    /// providing many examples around zero and throughout the space of numbers. However, some of
    /// these tests find the differece between two floating point values. This is fine close to
    /// zero, however the errors in this calculation get very large as we move towards the endges
    /// of float space. Rather than testing how the accuracy changes as a function of the size of
    /// the values, we instead clamp the values to a smaller range by implementing the Arbitrary
    /// trait for this struct.
    #[derive(Debug)]
    struct PropFloat {
        value: f64,
    }

    impl Arbitrary for PropFloat {
        type Parameters = ();
        type Strategy = BoxedStrategy<PropFloat>;

        fn arbitrary_with(_args: ()) -> Self::Strategy {
            (-100.0..=100.0)
                .prop_map(|x| PropFloat { value: x })
                .boxed()
        }
    }

    /// Init new with zeros should be the same as identity
    #[test]
    fn new_zeros() {
        assert_eq!(Transform2::new(0., (0., 0.)), Transform2::identity())
    }

    /// Ensure init with translation is same as returned
    #[proptest]
    fn new_translation(x: f64, y: f64) {
        let t = Transform2::new(0., (x, y));
        assert_eq!(t.position(), Point2::new(x, y))
    }

    /// Check get position is same as set
    #[proptest]
    fn set_get_translation(x: f64, y: f64) {
        let pos = Point2::new(x, y);
        let t = Transform2::identity().set_position(pos);
        assert_eq!(t.position(), pos);
    }

    /// Testing the multiplication of a Transform2 with a Vector2
    mod mul_transform_vec {
        use super::*;

        /// Transformation by identity matrix gives same result
        #[proptest]
        fn identity_transform(x: f64, y: f64) {
            let identity = Transform2::identity();
            let point = Point2::new(x, y);
            assert_eq!(identity * point, point);
        }

        /// Rotation keeps point same distance from origin
        #[proptest]
        fn rotation_length(angle: f64) {
            let t = Transform2::new(angle, (0., 0.));
            let point = Point2::new(1., 0.);
            assert_abs_diff_eq!(nalgebra::distance(&Point2::origin(), &(t * point)), 1.)
        }

        /// Translation from origin puts point in same location
        #[proptest]
        fn translation_from_origin(x: f64, y: f64) {
            let t = &Transform2::new(0., (x, y));
            assert_eq!(t * Point2::origin(), t.position())
        }

        /// Rotation followed by Translation is in circle about translation
        #[proptest]
        fn rotation_translation(rotation: PropFloat, translation: (PropFloat, PropFloat)) {
            let t = &Transform2::new(rotation.value, (translation.0.value, translation.1.value));
            let point = Point2::new(1., 0.);
            assert_abs_diff_eq!(
                nalgebra::distance(&Point2::origin(), &(t * point - t.position()).into()),
                1.,
                epsilon = Transform2::default_epsilon().sqrt()
            )
        }
    }

    /// Testing the multiplication of a Transform2 with a Transform2
    mod mul_transform_transform {
        use super::*;

        /// Translations with no rotations should be added
        #[proptest]
        fn mult_transform_translation(x: PropFloat, y: PropFloat) {
            let t = &Transform2::new(0., (x.value, y.value));
            assert_abs_diff_eq!((t * t).position(), 2. * Point2::new(x.value, y.value))
        }

        /// Translations should be added
        #[proptest]
        fn mult_transform_translations(t1: (f64, f64), t2: (f64, f64)) {
            let tf1 = &Transform2::new(0., t1);
            let tf2 = &Transform2::new(0., t2);
            assert_abs_diff_eq!(
                (tf1 * tf2).position(),
                tf1.get_translation() * tf2.position()
            )
        }

        /// Translations and rotations independent if translation multiplied first
        #[proptest]
        fn independent_trans_rot(rotation: f64, translation: (f64, f64)) {
            let tf1 = Transform2::new(0., translation);
            let tf2 = Transform2::new(rotation, (0., 0.));
            assert_abs_diff_eq!(tf1 * tf2, Transform2::new(rotation, translation))
        }

        /// Rotate a translation
        #[proptest]
        fn rotate_translation(rotation: f64, translation: (f64, f64)) {
            let tf1 = &Transform2::new(rotation, (0., 0.));
            let tf2 = &Transform2::new(0., translation);
            let rotated = tf1 * Point2::new(translation.0, translation.1);
            assert_abs_diff_eq!(tf1 * tf2, Transform2::new(rotation, (rotated.x, rotated.y)))
        }

        /// Rotations should be added
        #[proptest]
        fn combine_rotations(r1: PropFloat, r2: PropFloat) {
            let tf1 = &Transform2::new(r1.value, (0., 0.));
            let tf2 = &Transform2::new(r2.value, (0., 0.));
            assert_abs_diff_eq!(
                tf1 * tf2,
                Transform2::new(r1.value + r2.value, (0., 0.)),
                epsilon = 1e-10,
            )
        }

        #[proptest]
        // The value should be exact in this case
        #[allow(clippy::float_cmp)]
        fn rotation_and_trans_value(rotation: f64, translation: (f64, f64)) {
            let t = &Transform2::new(rotation, translation);
            dbg!(t * t);
            assert_eq!(t.0[(2, 2)], 1.)
        }

        #[proptest]
        fn rotation_and_trans(rotation: PropFloat, translation: (PropFloat, PropFloat)) {
            let t = &Transform2::new(rotation.value, (translation.0.value, translation.1.value));
            let position = t * Point2::new(translation.0.value, translation.1.value);
            // dbg!(position);
            println!("rotation: {:?}, translation: {:?}", rotation, translation);
            dbg!(
                (t * t).0,
                position,
                Transform2::new(rotation.value * 2., (position.x, position.y)).0
            );
            assert_abs_diff_eq!(
                t * t,
                Transform2::new(rotation.value + rotation.value, (position.x, position.y)),
                epsilon = 1e-10,
            );
        }

        #[proptest]
        fn rotation_and_trans_different(
            r1: PropFloat,
            r2: PropFloat,
            t1: (PropFloat, PropFloat),
            t2: (PropFloat, PropFloat),
        ) {
            let tf1 = &Transform2::new(r1.value, (t1.0.value, t1.1.value));
            let tf2 = &Transform2::new(r2.value, (t2.0.value, t2.1.value));
            let position = tf1 * Point2::new(t2.0.value, t2.1.value);
            println!("r1: {:?}, r2: {:?}, t1: {:?}, t2: {:?}", r1, r2, t1, t2);
            dbg!(
                (tf1 * tf2).0,
                Transform2::new(r1.value + r2.value, (position.x, position.y)).0
            );
            assert_abs_diff_eq!(
                tf1 * tf2,
                Transform2::new(r1.value + r2.value, (position.x, position.y)),
                epsilon = 1e-10,
            )
        }

        /// Compare with nalgebra implementation
        #[proptest]
        fn rotation_and_trans_nalgebra(r1: f64, r2: f64, t1: (f64, f64), t2: (f64, f64)) {
            let tf1 = Transform2::new(r1, t1);
            let tf2 = Transform2::new(r2, t2);
            let ntf1: na::Isometry2<f64> = na::Isometry2::from_parts(
                na::Translation2::new(t1.0, t1.1),
                na::UnitComplex::new(r1),
            );
            let ntf2 = na::Isometry2::from_parts(
                na::Translation2::new(t2.0, t2.1),
                na::UnitComplex::new(r2),
            );

            let mult = tf1 * tf2;
            assert_abs_diff_eq!(mult.0.matrix(), &(ntf1 * ntf2).to_homogeneous())
        }
    }

    #[test]
    fn transform() {
        let t = Transform2::new(f64::consts::PI / 2., (1., 1.));

        let point = Point2::new(0.2, 0.2);
        assert_eq!(t * point, Point2::new(0.8, 1.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(0.1, 0.2));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(-0.1, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(0.6, -0.2));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(-0.4, -0.2));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0)");
        let st = Transform2::from_operations(&input).unwrap();
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(-0.2, 0.));
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
