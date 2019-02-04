// shape.rs
//
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use std::fmt;
use std::ops;
use std::slice;

use crate::symmetry::Transform;

pub trait Intersect: ops::Mul<Transform> {
    fn intersects(&self, other: &Self) -> bool;
}

pub trait Shape: PartialEq + fmt::Debug + Clone + fmt::Display {
    type Component: Intersect
        + fmt::Debug
        + fmt::Display
        + ops::Mul<Transform, Output = Self::Component>;

    fn area(&self) -> f64;
    fn enclosing_radius(&self) -> f64;
    fn get_items(&self) -> Vec<Self::Component>;
    fn rotational_symmetries(&self) -> u64 {
        1
    }
    fn iter(&self) -> slice::Iter<'_, Self::Component>;
}

impl fmt::Display for LineShape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Shape({}) {{ ", self.name)?;
        for item in self.items.iter() {
            write!(f, "{}, ", item)?;
        }
        write!(f, "}}")
    }
}

mod atom;
mod atom_ops;

mod line;
mod line_ops;

mod line_shape;
mod molecular_shape;
mod shape_instance;

pub use atom::*;
pub use line::*;
pub use line_shape::*;
pub use molecular_shape::*;
pub use shape_instance::*;
