//
// symmetry.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use std::ops;

use approx;
use nalgebra::{Matrix2, Point2, Vector2};

/// Define the transformations of particle positions
///
/// These transformations are divided into the 'rotation' component and the translation component.
/// The 'rotation' component is represented in matrix form, and is in quotes since matrix contains
/// both the rotation transformation in addition to the mirror transformations. The translation
/// component is represented as a vector, being applied after the rotation.
///
#[derive(Debug, Clone, PartialEq)]
pub struct SymmetryTransform {
    pub rotation: Matrix2<f64>,
    pub translation: Vector2<f64>,
}

impl approx::AbsDiffEq for SymmetryTransform {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.rotation.abs_diff_eq(&other.rotation, epsilon)
            && self.translation.abs_diff_eq(&other.translation, epsilon)
    }
}

impl ops::Add<Vector2<f64>> for SymmetryTransform {
    type Output = SymmetryTransform;

    fn add(self, other: Vector2<f64>) -> Self::Output {
        SymmetryTransform {
            rotation: self.rotation,
            translation: self.translation + other,
        }
    }
}

impl ops::Add<Vector2<f64>> for &SymmetryTransform {
    type Output = SymmetryTransform;

    fn add(self, other: Vector2<f64>) -> Self::Output {
        SymmetryTransform {
            rotation: self.rotation,
            translation: self.translation + other,
        }
    }
}

impl SymmetryTransform {
    pub fn new(translation: Vector2<f64>, rotation: f64) -> SymmetryTransform {
        SymmetryTransform {
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

    pub fn identity() -> Self {
        Self {
            translation: Vector2::zeros(),
            rotation: Matrix2::identity(),
        }
    }

    /// Convert the string representation of a symmetry operation to a vector.
    ///
    /// This converts the string representation of an operation to
    pub fn from_operations(sym_ops: &str) -> SymmetryTransform {
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
        SymmetryTransform {
            rotation: rot,
            translation: trans,
        }
    }

    pub fn transform(&self, position: &Point2<f64>) -> Point2<f64> {
        self.rotation * position + self.translation
    }

    pub fn rotate(&self, vect: &Vector2<f64>) -> Vector2<f64> {
        self.rotation * vect
    }

    pub fn angle(&self) -> f64 {
        na::Rotation2::from_matrix_unchecked(self.rotation).angle()
    }
}

impl Default for SymmetryTransform {
    fn default() -> Self {
        Self {
            rotation: Matrix2::identity(),
            translation: Vector2::zeros(),
        }
    }
}

impl ops::Mul<Vector2<f64>> for SymmetryTransform {
    type Output = Vector2<f64>;

    fn mul(self, other: Vector2<f64>) -> Self::Output {
        self.rotation * other
    }
}

impl ops::Mul<Vector2<f64>> for &SymmetryTransform {
    type Output = Vector2<f64>;

    fn mul(self, other: Vector2<f64>) -> Self::Output {
        self.rotation * other
    }
}

impl ops::Mul<Point2<f64>> for SymmetryTransform {
    type Output = Point2<f64>;

    fn mul(self, other: Point2<f64>) -> Self::Output {
        self.rotation * other + self.translation
    }
}

impl ops::Mul<Point2<f64>> for &SymmetryTransform {
    type Output = Point2<f64>;

    fn mul(self, other: Point2<f64>) -> Self::Output {
        self.rotation * other + self.translation
    }
}

impl ops::Mul<SymmetryTransform> for SymmetryTransform {
    type Output = SymmetryTransform;

    fn mul(self, other: SymmetryTransform) -> Self::Output {
        let shift = self.rotate(&other.translation);

        Self {
            translation: self.translation + shift,
            rotation: self.rotation * other.rotation,
        }
    }
}

impl ops::Mul<SymmetryTransform> for &SymmetryTransform {
    type Output = SymmetryTransform;

    fn mul(self, other: SymmetryTransform) -> Self::Output {
        let shift = self.rotate(&other.translation);

        SymmetryTransform {
            translation: self.translation + shift,
            rotation: self.rotation * other.rotation,
        }
    }
}

impl ops::Mul<&SymmetryTransform> for SymmetryTransform {
    type Output = SymmetryTransform;

    fn mul(self, other: &SymmetryTransform) -> Self::Output {
        let shift = self.rotate(&other.translation);

        Self {
            translation: self.translation + shift,
            rotation: self.rotation * other.rotation,
        }
    }
}

#[cfg(test)]
mod symmetry_transform_tests {
    use super::*;
    use rand::prelude::*;
    use rand::rngs::SmallRng;
    use rand::Rng;
    use std::f64::consts::PI;

    #[test]
    fn default() {
        let point = Point2::new(0.2, 0.2);
        let transform = SymmetryTransform::default();
        assert_eq!(transform.transform(&point), point);
    }

    #[test]
    fn identity_transform() {
        let identity = SymmetryTransform::default();
        let point = Point2::new(0.2, 0.2);
        assert_eq!(identity.transform(&point), point);

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(identity.rotate(&vec), vec);
    }

    #[test]
    fn transform() {
        let isometry = SymmetryTransform::new(Vector2::new(1., 1.), PI / 2.);

        let point = Point2::new(0.2, 0.2);
        assert_eq!(isometry.transform(&point), Point2::new(0.8, 1.2));

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(isometry.rotate(&vec), Vector2::new(-0.2, 0.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y)");
        let st = SymmetryTransform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.1, 0.2));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y)");
        let st = SymmetryTransform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.1, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y)");
        let st = SymmetryTransform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.6, -0.2));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y)");
        let st = SymmetryTransform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.4, -0.2));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0)");
        let st = SymmetryTransform::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.2, 0.));
    }

    #[test]
    fn add_translation_test() {
        let ident = SymmetryTransform::identity();
        let trans = Vector2::new(1., 1.);
        assert_abs_diff_eq!((&ident + trans).translation, trans);

        let trans = Vector2::new(-1., -1.);
        assert_abs_diff_eq!((&ident + trans).translation, trans);
    }

    #[test]
    fn mult_symmetry_transform_rotations() {
        for i in 0..10 {
            let ident = SymmetryTransform::identity();
            let angle = i as f64 * PI / 10 as f64;
            let trans = SymmetryTransform::new(Vector2::zeros(), angle);
            assert_abs_diff_eq!((ident * &trans), trans);
        }
    }

    #[test]
    fn mult_symmetry_transform_translations() {
        for i in 0..10 {
            let ident = SymmetryTransform::identity();
            let trans = SymmetryTransform::new(Vector2::new(i as f64, i as f64), 0.);
            assert_abs_diff_eq!((ident * &trans), trans);
        }
    }

    fn rand_transform<R: Rng + ?Sized>(rng: &mut R) -> SymmetryTransform {
        SymmetryTransform::new(Vector2::new(rng.gen(), rng.gen()), rng.gen_range(0., PI))
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

            let result = transform1 * transform2;
            let expected = iso1 * iso2;

            assert_abs_diff_eq!(result.translation, expected.translation.vector);
            assert_abs_diff_eq!(result.rotation, expected.rotation.matrix());
        }
    }

}
