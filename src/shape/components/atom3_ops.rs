//
// atom3_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::ops::Mul;

use super::Atom3;
use crate::Transform3;

impl<'a, 'b> Mul<&'b Transform3> for &'a Atom3 {
    type Output = Atom3;

    fn mul(self, rhs: &Transform3) -> Self::Output {
        Atom3 {
            position: rhs * self.position,
            radius: self.radius,
        }
    }
}

impl<'a> Mul<Transform3> for &'a Atom3 {
    type Output = Atom3;

    fn mul(self, rhs: Transform3) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a Transform3> for Atom3 {
    type Output = Atom3;

    fn mul(self, rhs: &Transform3) -> Self::Output {
        &self * rhs
    }
}

impl Mul<Transform3> for Atom3 {
    type Output = Atom3;

    fn mul(self, rhs: Transform3) -> Self::Output {
        &self * &rhs
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn it_works() {}
}
