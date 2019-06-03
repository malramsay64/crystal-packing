//
// symmetry.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use nalgebra::{IsometryMatrix3, Matrix3, Rotation3, Translation3, Vector3};

use crate::traits::FromSymmetry;

pub type Transform3 = IsometryMatrix3<f64>;

impl FromSymmetry for Transform3 {
    /// Convert the string representation of a symmetry operation to a vector.
    ///
    /// This converts the string representation of an operation to a Transform,
    /// extracting the rotation and translation components.
    ///
    /// ```
    /// use packing::symmetry::{Transform3, FromSymmetry};
    /// let t3 = Transform3::from_operations("-x, y, z+1/2");
    /// ```
    ///
    /// This should generate a warning, however will still work fine, ignoring the third dimension.
    fn from_operations(sym_ops: &str) -> Self {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();
        let mut trans = Vector3::zeros();
        let mut rot = Matrix3::zeros();

        // We are using the <D>::dim() to constrain the iterator to the right number of dimensions
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
                    'z' => {
                        rot[(index, 2)] = sign;
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
        Transform3::from_parts(
            Translation3::<f64>::from(trans),
            Rotation3::<f64>::from_matrix_unchecked(rot),
        )
    }
}

#[cfg(test)]
mod test {
    use std::f64::consts::PI;

    use approx::assert_abs_diff_eq;
    use nalgebra::{Point3, Vector3};

    use super::*;

    #[test]
    fn default() {
        let point = Point3::new(0.2, 0.2, 0.2);
        let transform = Transform3::identity();
        assert_eq!(transform * point, point);
    }

    #[test]
    fn identity_transform() {
        let identity = Transform3::identity();
        let point = Point3::new(0.2, 0.2, 0.2);
        assert_eq!(identity * point, point);

        let vec = Vector2::new(0.2, 0.2, 0.2);
        assert_eq!(identity * vec, vec);
    }

    #[test]
    fn transform() {
        let isometry = Transform3::new(Vector3::new(1., 1., 1.), PI / 2.);

        let point = Point3::new(0.2, 0.2, 0.2);
        assert_eq!(isometry * point, Point2::new(0.8, 1.2));

        let vec = Vector3::new(0.2, 0.2, 0.2);
        assert_eq!(isometry * vec, Vector3::new(-0.2, 0.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y, z)");
        let st = Transform3::from_operations(&input);
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.1, 0.2, 0.3));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y, z)");
        let st = Transform3::from_operations(&input);
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(-0.1, 0.3, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y, z)");
        let st = Transform3::from_operations(&input);
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.6, -0.2, 0.3));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y, z)");
        let st = Transform3::from_operations(&input);
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(-0.4, -0.2, 0.3));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0, 0)");
        let st = Transform3::from_operations(&input);
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point2::new(-0.2, 0., 0));
    }

    #[test]
    fn parse_operation_3d() {
        let input = String::from("(x, y, z)");
        let st = Transform3::from_operations(&input);
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.1, 0.2, 0.3));
    }

    #[test]
    fn parse_operation_from_3d_complex() {
        let input = String::from("(x+z, y-x, -z+y+1/2)");
        let st = Transform3::from_operations(&input);
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.4, 0.1, 0.4));
    }

    #[test]
    fn mult_symmetry_transform_rotations() {
        for i in 0..10 {
            let ident = Transform3::identity();
            let angle = f64::from(i) * PI / f64::from(10);
            let trans = Transform3::new(Vector3::zeros(), angle);
            assert_abs_diff_eq!((ident * trans), trans);
        }
    }

    #[test]
    fn mult_symmetry_transform_translations() {
        for i in 0..10 {
            let ident = Transform3::identity();
            let trans = Transform3::new(Vector3::new(f64::from(i), f64::from(i)), 0.);
            assert_abs_diff_eq!((ident * trans), trans);
        }
    }
}
