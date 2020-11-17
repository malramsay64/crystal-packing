//
// lj2_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]
use std::ops::Mul;

use super::LJ2;
use crate::Transform2;

binop_impl_all!(
    Mul, mul;
    self: Transform2, rhs: LJ2, Output = LJ2;
    [ref ref] => {
        LJ2 {
            position: self * rhs.position,
            sigma: rhs.sigma,
            epsilon: rhs.epsilon,
            cutoff: rhs.cutoff
        }
    };
);

binop_impl_all!(
    Mul, mul;
    self: LJ2, rhs: Transform2, Output = LJ2;
    [ref ref] => {
        LJ2 {
            position: rhs * self.position,
            sigma: self.sigma,
            epsilon: self.epsilon,
            cutoff: self.cutoff
        }
    };
);
