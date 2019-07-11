//
// transform.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fmt::Write;

use log::warn;

use itertools::zip;
use nalgebra::allocator::Allocator;
use nalgebra::{
    DefaultAllocator, DimName, Isometry, MatrixN, Rotation, Translation, VectorN, U2, U3,
};

use crate::traits::*;

type Transform<D> = Isometry<f64, D, Rotation<f64, D>>;

impl Transformer for Transform2 {
    fn as_simple(&self) -> String {
        let mut s = String::new();
        write!(
            s,
            "x: {}, y: {}, angle: {}",
            self.translation.x,
            self.translation.y,
            self.rotation.angle()
        )
        .unwrap();
        s
    }
}

pub type Transform2 = Transform<U2>;
pub type Transform3 = Transform<U3>;

impl Periodic for f64 {
    type Output = f64;

    fn periodic(&self, rhs: Self) -> Self::Output {
        ((((self + 0.5) % rhs) + rhs) % rhs) - 0.5
    }
}

impl PeriodicAssign for f64 {
    fn periodic_assign(&mut self, rhs: Self) {
        *self = ((((*self + 0.5) % rhs) + rhs) % rhs) - 0.5
    }
}

impl<D: DimName> PeriodicAssign<f64> for VectorN<f64, D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    fn periodic_assign(&mut self, rhs: f64) {
        for x in self.iter_mut() {
            x.periodic_assign(rhs)
        }
    }
}

impl<D: DimName> Periodic<f64> for VectorN<f64, D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = VectorN<f64, D>;

    fn periodic(&self, rhs: f64) -> Self::Output {
        let mut tmp = self.clone();
        tmp.periodic_assign(rhs);
        tmp
    }
}

impl<D: DimName> PeriodicAssign<f64> for Translation<f64, D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    fn periodic_assign(&mut self, rhs: f64) {
        self.vector.periodic_assign(rhs)
    }
}

impl<D: DimName> Periodic<f64> for Translation<f64, D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Translation<f64, D>;

    fn periodic(&self, rhs: f64) -> Self::Output {
        let mut tmp = self.clone();
        tmp.periodic_assign(rhs);
        tmp
    }
}

impl<D: DimName> PeriodicAssign<f64> for Transform<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    fn periodic_assign(&mut self, rhs: f64) {
        self.translation.periodic_assign(rhs)
    }
}

impl<D: DimName> Periodic<f64> for Transform<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Transform<D>;

    fn periodic(&self, rhs: f64) -> Self::Output {
        let mut tmp = self.clone();
        tmp.periodic_assign(rhs);
        tmp
    }
}

impl<D: DimName> AdjustPeriod<D> for Translation<f64, D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Translation<f64, D>;

    fn adjust_period(&self, adjustment: VectorN<f64, D>) -> Self::Output {
        let mut tmp = self.clone();
        for (v, a) in zip(tmp.vector.iter_mut(), adjustment.iter()) {
            *v += a;
        }
        tmp
    }
}

impl<D: DimName> AdjustPeriod<D> for Transform<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Transform<D>;

    fn adjust_period(&self, adjustment: VectorN<f64, D>) -> Self::Output {
        let mut tmp = self.clone();
        tmp.translation = tmp.translation.adjust_period(adjustment);
        tmp
    }
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
    /// use packing::{Transform2, FromSymmetry};
    /// let t2 = Transform2::from_operations("-x, y");
    /// ```
    ///
    fn from_operations(sym_ops: &str) -> Result<Transform<D>, &'static str> {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();

        match operations.len() {
            x if x < <D>::dim() => return Err("Not enough dimensions in input"),
            x if x > <D>::dim() => return Err("Too many dimensions in input"),
            _ => (),
        }

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
                            return Err("No Z dimension in a 2D transform");
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
            trans[index] = constant;
        }
        Ok(Transform::<D>::from_parts(
            Translation::<f64, D>::from(trans),
            Rotation::<f64, D>::from_matrix_unchecked(rot),
        ))
    }
}

#[cfg(test)]
mod test_2d {
    use std::f64::consts::PI;

    use approx::assert_abs_diff_eq;
    use nalgebra::{Point2, Vector2};

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

#[cfg(test)]
mod test_3d {
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

        let vec = Vector3::new(0.2, 0.2, 0.2);
        assert_eq!(identity * vec, vec);
    }

    #[test]
    fn transform() {
        let isometry = Transform3::new(Vector3::new(1., 1., 1.), Vector3::new(0., 0., PI / 2.));

        let point = Point3::new(0.2, 0.2, 0.2);
        assert_eq!(isometry * point, Point3::new(0.8, 1.2, 1.2));

        let vec = Vector3::new(0.2, 0.2, 0.2);
        assert_eq!(isometry * vec, Vector3::new(-0.2, 0.2, 0.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y, z)");
        let st = Transform3::from_operations(&input).unwrap();
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.1, 0.2, 0.3));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y, z)");
        let st = Transform3::from_operations(&input).unwrap();
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(-0.1, 0.3, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y, z)");
        let st = Transform3::from_operations(&input).unwrap();
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.6, -0.2, 0.3));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y, z)");
        let st = Transform3::from_operations(&input).unwrap();
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(-0.4, -0.2, 0.3));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0, 0)");
        let st = Transform3::from_operations(&input).unwrap();
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(-0.2, 0., 0.));
    }

    #[test]
    fn parse_operation_3d() {
        let input = String::from("(x, y, z)");
        let st = Transform3::from_operations(&input).unwrap();
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.1, 0.2, 0.3));
    }

    #[test]
    fn parse_operation_from_3d_complex() {
        let input = String::from("(x+z, y-x, -z+y+1/2)");
        let st = Transform3::from_operations(&input).unwrap();
        let point = Point3::new(0.1, 0.2, 0.3);
        assert_abs_diff_eq!(st * point, Point3::new(0.4, 0.1, 0.4));
    }

    #[test]
    fn mult_symmetry_transform_rotations() {
        for i in 0..10 {
            let ident = Transform3::identity();
            let angle = f64::from(i) * PI / f64::from(10);
            let trans = Transform3::new(Vector3::zeros(), Vector3::new(angle, angle, angle));
            assert_abs_diff_eq!((ident * trans), trans);
        }
    }

    #[test]
    fn mult_symmetry_transform_translations() {
        for i in 0..10 {
            let ident = Transform3::identity();
            let trans = Transform3::new(
                Vector3::new(f64::from(i), f64::from(i), f64::from(i)),
                Vector3::zeros(),
            );
            assert_abs_diff_eq!((ident * trans), trans);
        }
    }

    #[test]
    #[should_panic]
    fn parse_operation_2_inputs() {
        let input = String::from("(x, y)");
        Transform3::from_operations(&input).unwrap();
    }

    #[test]
    #[should_panic]
    fn parse_operation_invalid() {
        let input = String::from("(x, y, w)");
        Transform3::from_operations(&input).unwrap();
    }
}
