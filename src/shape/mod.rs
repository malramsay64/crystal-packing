//
// shape.rs
//
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fmt;
use std::ops::Mul;
use std::slice;

use nalgebra::allocator::Allocator;
use nalgebra::{DefaultAllocator, DimName};

use crate::Transform;

pub trait U2: DimName {}

impl U2 for na::U2 {}

pub trait Intersect<D: DimName>: Sized + Mul<Transform<D>, Output = Self>
where
    for<'a> Self: Mul<&'a Transform<D>, Output = Self>,
    for<'a, 'b> &'a Self: Mul<&'b Transform<D>, Output = Self>,
    for<'a> &'a Self: Mul<Transform<D>, Output = Self>,
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    fn intersects(&self, other: &Self) -> bool;
}

pub trait Shape<D: DimName>: PartialEq + fmt::Debug + Clone + fmt::Display
where
    for<'a> Self::Component: Mul<&'a Transform<D>, Output = Self::Component>,
    for<'a, 'b> &'a Self::Component: Mul<&'b Transform<D>, Output = Self::Component>,
    for<'a> &'a Self::Component: Mul<Transform<D>, Output = Self::Component>,
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Component: Intersect<D> + fmt::Debug + fmt::Display;

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
