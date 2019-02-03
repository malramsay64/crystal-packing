//
// symmetry_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
#![allow(clippy::suspicious_arithmetic_impl)]
use std::ops::{Add, Mul};

use nalgebra::{Point2, Vector2};

use crate::symmetry::SymmetryTransform;

/// Implement a binary operation
///
/// https://github.com/rustsim/nalgebra/blob/911ddca588d2c89e4f7ec9b45aa0aa7f787209c4/src/geometry/isometry_ops.rs
macro_rules! symmetry_binop_impl(
    ($Op: ident, $op: ident;
     $lhs: ident: $Lhs: ty, $rhs: ident: $Rhs: ty, Output = $Output: ty;
     $action: expr; $($lives: tt),*) => {
        impl<$($lives ,)*> $Op<$Rhs> for $Lhs {
            type Output = $Output;

            #[inline]
            fn $op($lhs, $rhs: $Rhs) -> Self::Output {
                $action
            }
        }
    }
);

/// Implement all combinations of operations
///
/// When implementing an operation there are four different operations that need to be
/// implemented to cover all use cases, since an operation of a reference is considered different
/// to on a value. This is a macro which makes writing the implementation of that operation much
/// simpler, only needing to write out the implementation for the reference * reference values.
///
/// This macro is adapted from the
/// [nalgebra](https://github.com/rustsim/nalgebra/blob/911ddca588d2c89e4f7ec9b45aa0aa7f787209c4/src/geometry/isometry_ops.rs)
/// source code.
///
macro_rules! symmetry_binop_impl_all(
    // The trait we implementing, along with the function required to implement the trait
    ($Op: ident, $op: ident;
     // This is types and variable names of each side of the operand, and then the type of the
     // output.
     $lhs: ident: $Lhs: ty, $rhs: ident: $Rhs: ty, Output = $Output: ty;
     [val val] => $action_val_val: expr;
     [ref val] => $action_ref_val: expr;
     [val ref] => $action_val_ref: expr;
     [ref ref] => $action_ref_ref: expr;) => {
        symmetry_binop_impl!(
            $Op, $op;
            $lhs: $Lhs, $rhs: $Rhs, Output = $Output;
            $action_val_val; );

        symmetry_binop_impl!(
            $Op, $op;
            $lhs: &'a $Lhs, $rhs: $Rhs, Output = $Output;
            $action_ref_val; 'a);

        symmetry_binop_impl!(
            $Op, $op;
            $lhs: $Lhs, $rhs: &'b $Rhs, Output = $Output;
            $action_val_ref; 'b);

        symmetry_binop_impl!(
            $Op, $op;
            $lhs: &'a $Lhs, $rhs: &'b $Rhs, Output = $Output;
            $action_ref_ref; 'a, 'b);
    }
);

// SymmetryTransform + Vector
symmetry_binop_impl_all!(
    Add, add;
    self: SymmetryTransform, rhs: Vector2<f64>, Output = SymmetryTransform;
    // It is always safe to take a reference, so reduce all instances to a reference * reference
    // then have a single implementation.
    [val val] => &self + &rhs;
    [ref val] => self + &rhs;
    [val ref] => &self + rhs;
    [ref ref] => {
        SymmetryTransform {
            rotation: self.rotation,
            translation: self.translation + rhs,
        }
    };
);

// SymmetryTransform * Vector
// Multiplying by a vector doesn't apply the translational component, only the rotation.
symmetry_binop_impl_all!(
    Mul, mul;
    self: SymmetryTransform, rhs: Vector2<f64>, Output = Vector2<f64>;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        self.rotation * rhs
    };
);

// SymmetryTransform * Point
symmetry_binop_impl_all!(
    Mul, mul;
    self: SymmetryTransform, rhs: Point2<f64>, Output = Point2<f64>;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        self.rotation * rhs + self.translation
    };
);

// SymmetryTransform * SymmetryTransform
symmetry_binop_impl_all!(
    Mul, mul;
    self: SymmetryTransform, rhs: SymmetryTransform, Output = SymmetryTransform;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        let shift = self.rotate(&rhs.translation);

        SymmetryTransform {
            translation: self.translation + shift,
            rotation: self.rotation * rhs.rotation,
        }
    };
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {}
}
