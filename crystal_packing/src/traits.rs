//
// traits.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::{fmt, ops, slice};

use anyhow::Error;
use nalgebra::SVector;
use serde::Serialize;
use svg::node::element::Group;
use svg::Document;

use crate::{Basis, Transform2};

pub trait Transformer {
    fn as_simple(&self) -> String;
}

pub trait Periodic<Rhs = Self> {
    type Output;

    fn periodic(&self, rhs: Rhs) -> Self::Output;
}

pub trait PeriodicAssign<Rhs = Self> {
    fn periodic_assign(&mut self, rhs: Rhs);
}

pub trait AdjustPeriod<const D: usize> {
    type Output;
    fn adjust_period(&self, adjustment: SVector<f64, D>) -> Self::Output;
}

pub trait Intersect {
    fn intersects(&self, other: &Self) -> bool;
    fn area(&self) -> f64;
}

pub trait Potential {
    fn energy(&self, other: &Self) -> f64;
}

pub trait Shape:
    Clone + Send + Sync + Serialize + fmt::Debug + fmt::Display + ToSVG<Value = Group>
{
    type Component: Clone
        + Send
        + Sync
        + Serialize
        + fmt::Debug
        + fmt::Display
        + ops::Mul<Transform2, Output = Self::Component>
        + ToSVG;

    fn score(&self, other: &Self) -> Option<f64>;
    fn enclosing_radius(&self) -> f64;
    fn get_items(&self) -> Vec<Self::Component>;
    fn rotational_symmetries(&self) -> u64 {
        1
    }
    fn iter(&self) -> slice::Iter<'_, Self::Component>;
    fn transform(&self, transform: &Transform2) -> Self;
}

pub trait FromSymmetry: Sized {
    fn from_operations(ops: &str) -> Result<Self, Error>;
}

pub trait State:
    Eq
    + PartialEq
    + PartialOrd
    + Ord
    + Clone
    + Send
    + Sync
    + Serialize
    + fmt::Debug
    + ToSVG<Value = Document>
{
    fn score(&self) -> Option<f64>;
    fn generate_basis(&self) -> Vec<Basis>;
    fn total_shapes(&self) -> usize;
    fn as_positions(&self) -> Result<String, Error>;
}

pub trait ToSVG {
    type Value: svg::Node;
    fn as_svg(&self) -> Self::Value;
}
