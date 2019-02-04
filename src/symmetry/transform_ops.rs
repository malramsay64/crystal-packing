//
// symmetry_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
#![allow(clippy::suspicious_arithmetic_impl)]
#![allow(clippy::op_ref)]
use std::ops::{Add, Mul};

use nalgebra::{Point2, Vector2};

use crate::symmetry::Transform;

// Transform + Vector
binop_impl_all!(
    Add, add;
    self: Transform, rhs: Vector2<f64>, Output = Transform;
    // It is always safe to take a reference, so reduce all instances to a reference * reference
    // then have a single implementation.
    [val val] => &self + &rhs;
    [ref val] => self + &rhs;
    [val ref] => &self + rhs;
    [ref ref] => {
        Transform {
            rotation: self.rotation,
            translation: self.translation + rhs,
        }
    };
);

// Transform * Vector
// Multiplying by a vector doesn't apply the translational component, only the rotation.
binop_impl_all!(
    Mul, mul;
    self: Transform, rhs: Vector2<f64>, Output = Vector2<f64>;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        self.rotation * rhs
    };
);

// Transform * Point
binop_impl_all!(
    Mul, mul;
    self: Transform, rhs: Point2<f64>, Output = Point2<f64>;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        self.rotation * rhs + self.translation
    };
);

// Transform * Transform
binop_impl_all!(
    Mul, mul;
    self: Transform, rhs: Transform, Output = Transform;
    [val val] => &self * &rhs;
    [ref val] => self * &rhs;
    [val ref] => &self * rhs;
    [ref ref] => {
        let shift = self.rotate(&rhs.translation);

        Transform {
            translation: self.translation + shift,
            rotation: self.rotation * rhs.rotation,
        }
    };
);
