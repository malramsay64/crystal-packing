//
// macros_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

#![macro_use]

/// Implement a binary operation
///
/// This converts a set of values and types to an impl block for a specific operation. In doing
/// this it is possible to implement all the operations at once with the binop_impl_all macro.
///
macro_rules! _binop_impl(
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
/// # Example
///
/// ```
/// # use crate::ops_macros::binop_impl_all;
/// struct Wrapper(f64);
///
/// binop_impl_all!(
///     Mul, mul;
///     self: Wrapper, rhs: f64, Output = f64;
///     [ref ref] => self * rhs
/// )
///
///
/// let w = Wrapper(8.);
/// w * 2.
/// # assert!(w * 2. == 16.)
///
/// ```
macro_rules! binop_impl_all(
    // The trait we implementing, along with the function required to implement the trait
    ($Op: ident, $op: ident;
     // This is types and variable names of each side of the operand, and then the type of the
     // output.
     $lhs: ident: $Lhs: ty, $rhs: ident: $Rhs: ty, Output = $Output: ty;
     // These are the expressions used to perform each of the operations on the different
     // combinations of value and reference. When implementing, by converting all values to the
     // form ref * ref, there is only a single implementation required.
     [ref ref] => $action_ref_ref: expr;) => {
        // This is the block where the transformation actually occurs.
        _binop_impl!(
            $Op, $op;
            $lhs: $Lhs, $rhs: $Rhs, Output = $Output;
            &$lhs * &$rhs; );

        _binop_impl!(
            $Op, $op;
            $lhs: &'a $Lhs, $rhs: $Rhs, Output = $Output;
            $lhs * &$rhs; 'a);

        _binop_impl!(
            $Op, $op;
            $lhs: $Lhs, $rhs: &'b $Rhs, Output = $Output;
            &$lhs * $rhs; 'b);

        _binop_impl!(
            $Op, $op;
            $lhs: &'a $Lhs, $rhs: &'b $Rhs, Output = $Output;
            $action_ref_ref; 'a, 'b);
    }
);
