//
// atom3_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::ops::Mul;

use super::Atom3;
use crate::Transform3;

binop_impl_all!(
    Mul, mul;
    self: Transform3, rhs: Atom3, Output = Atom3;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Atom3 {
            position: self * rhs.position,
            radius: rhs.radius,
        }
    };
);

binop_impl_all!(
    Mul, mul;
    self: Atom3, rhs: Transform3, Output = Atom3;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Atom3 {
            position: rhs * self.position,
            radius: self.radius,
        }
    };
);

#[cfg(test)]
mod tests {

    #[test]
    fn it_works() {}
}
