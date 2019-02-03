//
// symmetry.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use approx;
// Re-export these to allow importing along with the Transform struct
pub use nalgebra::{Matrix2, Point2, Vector2};

/// Define the transformations of particle positions
///
/// These transformations are divided into the 'rotation' component and the translation component.
/// The 'rotation' component is represented in matrix form, and is in quotes since matrix contains
/// both the rotation transformation in addition to the mirror transformations. The translation
/// component is represented as a vector, being applied after the rotation.
///
#[derive(Debug, Clone, PartialEq)]
pub struct Transform {
    pub rotation: Matrix2<f64>,
    pub translation: Vector2<f64>,
}

impl approx::AbsDiffEq for Transform {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.rotation.abs_diff_eq(&other.rotation, epsilon)
            && self.translation.abs_diff_eq(&other.translation, epsilon)
    }
}

impl Transform {
    /// Instantiate a symmetry transform from a translation vector and angle
    ///
    /// This is simpler than specifying an entire rotation matrix, although it doesn't allow for
    /// the specification of any mirror planes. It can be used as;
    /// ```
    /// use packing::symmetry::{Transform, Vector2};
    /// use std::f64::consts::PI;
    /// let t = Transform::new(Vector2::new(1., 1.), PI /2.);
    /// ```
    ///
    pub fn new(translation: Vector2<f64>, rotation: f64) -> Transform {
        Transform {
            // Convert a rotation angle in radians to a rotation matrix.
            rotation: Matrix2::new(
                rotation.cos(),
                -rotation.sin(),
                rotation.sin(),
                rotation.cos(),
            ),
            translation,
        }
    }

    /// Create the identity transform
    ///
    /// This is the transform which leaves the transformed vector unchanged.
    ///
    /// ```
    /// use packing::symmetry::{Transform, Vector2};
    /// let t = Transform::identity();
    /// let c = Vector2::new(-1., 1.);
    /// assert_eq!(t * &c, c);
    /// ```
    ///
    pub fn identity() -> Self {
        Self {
            translation: Vector2::zeros(),
            rotation: Matrix2::identity(),
        }
    }

    /// Convert the string representation of a symmetry operation to a vector.
    ///
    /// This converts the string representation of an operation to a Transform,
    /// extracting the rotation and translation components.
    ///
    /// ```
    /// use packing::symmetry::Transform;
    /// let t = Transform::from_operations("-x, y");
    /// ```
    ///
    pub fn from_operations(sym_ops: &str) -> Transform {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();
        let mut trans = Vector2::zeros();
        let mut rot: Matrix2<f64> = Matrix2::zeros();

        for (index, op) in operations.iter().enumerate() {
            let mut sign = 1.;
            let mut constant = 0.;
            let mut operator: Option<char> = None;
            for c in op.chars() {
                match c {
                    'x' => {
                        rot[(index, 0)] = sign;
                        sign = 1.;
                    }
                    'y' => {
                        rot[(index, 1)] = sign;
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
                            Some(op) => match op {
                                '/' => sign * constant / val,
                                '*' => sign * constant * val,
                                _ => 0.,
                            },
                            None => sign * val,
                        };
                        // Reset values
                        operator = None;
                        sign = 1.
                    }
                    // Default is do nothing (shouldn't encounter this at all)
                    _ => {}
                };
            }
            trans[index] = constant;
        }
        Transform {
            rotation: rot,
            translation: trans,
        }
    }

    pub fn transform(&self, position: &Point2<f64>) -> Point2<f64> {
        self * position
    }

    pub fn rotate(&self, vect: &Vector2<f64>) -> Vector2<f64> {
        self * vect
    }

    pub fn angle(&self) -> f64 {
        na::Rotation2::from_matrix_unchecked(self.rotation).angle()
    }
}

impl Default for Transform {
    fn default() -> Self {
        Self {
            rotation: Matrix2::identity(),
            translation: Vector2::zeros(),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::prelude::*;
    use rand::rngs::SmallRng;
    use rand::Rng;
    use std::f64::consts::PI;

    #[test]
    fn default() {
        let point = Point2::new(0.2, 0.2);
        let transform = Transform::default();
        assert_eq!(transform.transform(&point), point);
    }

    #[test]
    fn identity_transform() {
        let identity = Transform::default();
        let point = Point2::new(0.2, 0.2);
        assert_eq!(identity.transform(&point), point);

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(identity.rotate(&vec), vec);
    }

    #[test]
    fn transform() {
        let isometry = Transform::new(Vector2::new(1., 1.), PI / 2.);

        let point = Point2::new(0.2, 0.2);
        assert_eq!(isometry.transform(&point), Point2::new(0.8, 1.2));

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(isometry.rotate(&vec), Vector2::new(-0.2, 0.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y)");
        let st = Transform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.1, 0.2));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y)");
        let st = Transform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.1, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y)");
        let st = Transform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.6, -0.2));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y)");
        let st = Transform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.4, -0.2));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0)");
        let st = Transform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.2, 0.));
    }

    #[test]
    fn add_translation_test() {
        let ident = Transform::identity();
        let trans = Vector2::new(1., 1.);
        assert_abs_diff_eq!((&ident + trans).translation, trans);

        let trans = Vector2::new(-1., -1.);
        assert_abs_diff_eq!((&ident + trans).translation, trans);
    }

    #[test]
    fn mult_symmetry_transform_rotations() {
        for i in 0..10 {
            let ident = Transform::identity();
            let angle = f64::from(i) * PI / f64::from(10);
            let trans = Transform::new(Vector2::zeros(), angle);
            assert_abs_diff_eq!((ident * &trans), trans);
        }
    }

    #[test]
    fn mult_symmetry_transform_translations() {
        for i in 0..10 {
            let ident = Transform::identity();
            let trans = Transform::new(Vector2::new(f64::from(i), f64::from(i)), 0.);
            assert_abs_diff_eq!((ident * &trans), trans);
        }
    }

    fn rand_transform<R: Rng + ?Sized>(rng: &mut R) -> Transform {
        Transform::new(Vector2::new(rng.gen(), rng.gen()), rng.gen_range(0., PI))
    }

    #[test]
    fn mult_symmetry_transform_test() {
        let mut rng = SmallRng::seed_from_u64(0);

        for _ in 0..500 {
            let transform1 = rand_transform(&mut rng);
            let transform2 = rand_transform(&mut rng);

            let iso1 = na::IsometryMatrix2::from_parts(
                na::Translation2::from(transform1.translation),
                na::Rotation2::from_matrix_unchecked(transform1.rotation),
            );
            let iso2 = na::IsometryMatrix2::from_parts(
                na::Translation2::from(transform2.translation),
                na::Rotation2::from_matrix_unchecked(transform2.rotation),
            );

            let result: Transform = transform1 * transform2;
            let expected = iso1 * iso2;

            assert_abs_diff_eq!(result.translation, expected.translation.vector);
            assert_abs_diff_eq!(result.rotation, expected.rotation.matrix());
        }
    }

}
