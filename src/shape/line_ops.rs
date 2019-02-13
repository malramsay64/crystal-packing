//
// line_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]
use std::ops::Mul;

use crate::{Line, Transform2};

binop_impl_all!(
    Mul, mul;
    self: Transform2, rhs: Line, Output = Line;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Line {
            start: self * rhs.start,
            end: self * rhs.end,
        }
    };
);

binop_impl_all!(
    Mul, mul;
    self: Line, rhs: Transform2, Output = Line;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        Line {
            start: rhs * self.start,
            end: rhs * self.end,
        }
    };
);
