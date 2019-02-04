//
// atom_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]

use std::ops::Mul;

use crate::shape::Atom;
use crate::symmetry::Transform;

binop_impl_all!(
    Mul, mul;
    self: Transform, rhs: Atom, Output = Atom;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Atom {
            position: self * rhs.position,
            radius: rhs.radius,
        }
    };
);

binop_impl_all!(
    Mul, mul;
    self: Atom, rhs: Transform, Output = Atom;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Atom {
            position: rhs * self.position,
            radius: self.radius,
        }
    };
);
