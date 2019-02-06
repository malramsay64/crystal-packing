//
// symmetry.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use log::warn;

// Re-export these to allow importing along with the Transform struct
use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, DimName, MatrixN, VectorN};
pub use nalgebra::{Isometry, Rotation, Translation, U2, U3};

pub type Transform<D> = Isometry<f64, D, Rotation<f64, D>>;

pub type Transform2 = Transform<U2>;
pub type Transform3 = Transform<U3>;

pub trait FromSymmetry {
    fn from_operations(ops: &str) -> Self;
}

impl<D: DimName> FromSymmetry for Transform<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    /// Convert the string representation of a symmetry operation to a vector.
    ///
    /// This converts the string representation of an operation to a Transform,
    /// extracting the rotation and translation components.
    ///
    /// ```
    /// use packing::symmetry::{Transform2, Transform3, FromSymmetry};
    /// let t2 = Transform2::from_operations("-x, y");
    ///
    /// let t3 = Transform3::from_operations("-x, y, z+1/2");
    /// ```
    ///
    /// Where a z dimension, or 3 dimensions worth of values are passed into this function, it
    /// handles the 'error' condition by ignoring those values. That is:
    /// ```
    /// use packing::symmetry::{Transform2, FromSymmetry};
    /// use nalgebra::Point2;
    /// let t = Transform2::from_operations("-x+z, y, z");
    /// assert_eq!(t * Point2::new(1., 1.), Point2::new(-1., 1.));
    /// ```
    /// This should generate a warning, however will still work fine, ignoring the third dimension.
    fn from_operations(sym_ops: &str) -> Transform<D> {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();
        let mut trans = VectorN::<f64, D>::zeros();
        let mut rot = MatrixN::<f64, D>::zeros();

        // We are using the <D>::dim() to constrain the iterator to the right number of dimensions
        for (index, op) in operations.iter().enumerate().take(<D>::dim()) {
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
                        if <D>::dim() < 3 {
                            warn!("Too many input dimensions, ignoring");
                            sign = 1.;
                        } else {
                            rot[(index, 2)] = sign;
                            sign = 1.;
                        }
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
        Transform::<D>::from_parts(
            Translation::<f64, D>::from(trans),
            Rotation::<f64, D>::from_matrix_unchecked(rot),
        )
    }
}

#[cfg(test)]
mod test {
    use std::f64::consts::PI;

    use approx::assert_abs_diff_eq;
    use nalgebra::{Point2, Point3, Vector2};

    use super::*;

    #[test]
    fn default() {
        let point = Point2::new(0.2, 0.2);
        let transform = Transform2::identity();
        assert_eq!(transform * point, point);
    }

    #[test]
    fn identity_transform() {
        let identity = Transform2::identity();
        let point = Point2::new(0.2, 0.2);
        assert_eq!(identity * point, point);

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(identity * vec, vec);
    }

    #[test]
    fn transform() {
        let isometry = Transform2::new(Vector2::new(1., 1.), PI / 2.);

        let point = Point2::new(0.2, 0.2);
        assert_eq!(isometry * point, Point2::new(0.8, 1.2));

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(isometry * vec, Vector2::new(-0.2, 0.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y)");
        let st = Transform2::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(0.1, 0.2));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y)");
        let st = Transform2::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(-0.1, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y)");
        let st = Transform2::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(0.6, -0.2));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y)");
        let st = Transform2::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(-0.4, -0.2));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0)");
        let st = Transform2::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(-0.2, 0.));
    }

    #[test]
    fn parse_operation_2d_from_3d() {
        let input = String::from("(x, y, z)");
        let st = Transform2::from_operations(&input);
        let point = Point2::new(0.1, 0.2);
        assert_abs_diff_eq!(st * point, Point2::new(0.1, 0.2));
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
            let ident = Transform2::identity();
            let angle = f64::from(i) * PI / f64::from(10);
            let trans = Transform2::new(Vector2::zeros(), angle);
            assert_abs_diff_eq!((ident * trans), trans);
        }
    }

    #[test]
    fn mult_symmetry_transform_translations() {
        for i in 0..10 {
            let ident = Transform2::identity();
            let trans = Transform2::new(Vector2::new(f64::from(i), f64::from(i)), 0.);
            assert_abs_diff_eq!((ident * trans), trans);
        }
    }
}
