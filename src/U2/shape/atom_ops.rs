//
// atom_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]

use std::ops::Mul;

use super::{Atom2, Transform2};

impl<'a, 'b> Mul<&'b Transform2> for &'a Atom2 {
    type Output = Atom2;

    fn mul(self, rhs: &Transform2) -> Self::Output {
        Atom2 {
            position: rhs * self.position,
            radius: self.radius,
        }
    }
}

impl<'a> Mul<Transform2> for &'a Atom2 {
    type Output = Atom2;

    fn mul(self, rhs: Transform2) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a Transform2> for Atom2 {
    type Output = Atom2;

    fn mul(self, rhs: &Transform2) -> Self::Output {
        &self * rhs
    }
}

impl Mul<Transform2> for Atom2 {
    type Output = Atom2;

    fn mul(self, rhs: Transform2) -> Self::Output {
        &self * &rhs
    }
}
