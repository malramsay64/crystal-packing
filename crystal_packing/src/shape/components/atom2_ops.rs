//
// atom2_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]

use std::ops::Mul;

use super::Atom2;
use crate::Transform2;

binop_impl_all!(
    Mul, mul;
    self: Transform2, rhs: Atom2, Output = Atom2;
    [ref ref] => {
        Atom2 {
            position: self * rhs.position,
            radius: rhs.radius,
        }
    };
);

binop_impl_all!(
    Mul, mul;
    self: Atom2, rhs: Transform2, Output = Atom2;
    [ref ref] => {
        Atom2 {
            position: rhs * self.position,
            radius: self.radius,
        }
    };
);
