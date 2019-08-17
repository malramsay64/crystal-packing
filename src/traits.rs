//
// traits.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::{fmt, ops, slice};

use nalgebra::allocator::Allocator;
use nalgebra::{DefaultAllocator, DimName, Vector2, VectorN};
use rand::Rng;
use serde::Serialize;
use svg::node::element::Group;
use svg::Document;

use crate::wallpaper::WyckoffSite;
use crate::{CrystalFamily, StandardBasis, Transform2};

pub trait Transformer {
    fn as_simple(&self) -> String;
}

pub trait Basis {
    fn set_value(&mut self, new_value: f64);
    fn get_value(&self) -> f64;
    fn reset_value(&self);
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R, step_size: f64) -> f64;
}

pub trait Periodic<Rhs = Self> {
    type Output;

    fn periodic(&self, rhs: Rhs) -> Self::Output;
}

pub trait PeriodicAssign<Rhs = Self> {
    fn periodic_assign(&mut self, rhs: Rhs);
}

pub trait AdjustPeriod<D: DimName>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output;
    fn adjust_period(&self, adjustment: VectorN<f64, D>) -> Self::Output;
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

    fn score(&self, other: &Self) -> Result<f64, &'static str>;
    fn enclosing_radius(&self) -> f64;
    fn get_items(&self) -> Vec<Self::Component>;
    fn rotational_symmetries(&self) -> u64 {
        1
    }
    fn iter(&self) -> slice::Iter<'_, Self::Component>;
    fn transform(&self, transform: &Transform2) -> Self;
}

pub trait FromSymmetry: Sized {
    fn from_operations(ops: &str) -> Result<Self, &'static str>;
}

pub trait Cell:
    Clone + Send + Sync + Serialize + fmt::Debug + fmt::Display + ToSVG<Value = Group>
{
    fn periodic_images<'a>(
        &'a self,
        transform: Transform2,
        zero: bool,
    ) -> Box<dyn Iterator<Item = Transform2> + 'a>;
    fn from_family(group: CrystalFamily, max_size: f64) -> Self;
    fn to_cartesian_isometry(&self, transform: &Transform2) -> Transform2;
    fn to_cartesian_point(&self, point: Vector2<f64>) -> Vector2<f64>;
    fn get_degrees_of_freedom(&mut self) -> Vec<StandardBasis>;
    fn center(&self) -> Vector2<f64>;
    fn area(&self) -> f64;
    fn get_corners(&self) -> Vec<Vector2<f64>>;
}

pub trait Site: Clone + Send + Sync + Serialize + fmt::Debug {
    fn transform(&self) -> Transform2;
    fn positions<'a>(&'a self) -> Box<dyn Iterator<Item = Transform2> + 'a>;
    fn multiplicity(&self) -> usize;
    fn from_wyckoff(wyckoff: &WyckoffSite) -> Self;
    fn get_basis(&mut self, rot_symmetry: u64) -> Vec<StandardBasis>;
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
    fn score(&self) -> Result<f64, &'static str>;
    fn generate_basis(&mut self) -> Vec<StandardBasis>;
    fn total_shapes(&self) -> usize;
    fn as_positions(&self) -> Result<String, fmt::Error>;
}

pub trait ToSVG {
    type Value: svg::Node;
    fn as_svg(&self) -> Self::Value;
}
