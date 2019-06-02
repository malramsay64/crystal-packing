//
// line_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]
use std::ops::Mul;

use super::{Line2, Transform2};

binop_impl_all!(
    Mul, mul;
    self: Transform2, rhs: Line2, Output = Line2;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Line2 {
            start: self * rhs.start,
            end: self * rhs.end,
        }
    };
);

binop_impl_all!(
    Mul, mul;
    self: Line2, rhs: Transform2, Output = Line2;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Line2 {
            start: rhs * self.start,
            end: rhs * self.end,
        }
    };
);
