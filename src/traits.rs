//
// traits.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use crate::wallpaper::WyckoffSite;
use crate::StandardBasis;
use crate::U2::CrystalFamily;
use std::{fmt, ops, slice};

pub trait Intersect {
    fn intersects(&self, other: &Self) -> bool;
}

pub trait Shape: Clone + Send + Sync + fmt::Debug + fmt::Display + Intersect {
    type Transform: Clone + Send + Sync + fmt::Debug;
    type Component: Intersect
        + Clone
        + Send
        + Sync
        + fmt::Debug
        + fmt::Display
        + ops::Mul<Self::Transform, Output = Self::Component>;

    fn area(&self) -> f64;
    fn enclosing_radius(&self) -> f64;
    fn get_items(&self) -> Vec<Self::Component>;
    fn rotational_symmetries(&self) -> u64 {
        1
    }
    fn iter(&self) -> slice::Iter<'_, Self::Component>;
    fn transform(&self, transform: &Self::Transform) -> Self;
}

pub trait FromSymmetry {
    fn from_operations(ops: &str) -> Self;
}

pub trait Cell: Clone + Send + Sync + fmt::Debug {
    type Transform: Clone + Send + Sync + fmt::Debug;
    type Point: Clone + Send + Sync + fmt::Display;

    fn periodic_images(&self, transform: &Self::Transform, zero: bool) -> Vec<Self::Transform>;
    fn from_family(group: &CrystalFamily, max_size: f64) -> Self;
    fn to_cartesian_isometry(&self, transform: &Self::Transform) -> Self::Transform;
    fn to_cartesian_point(&self, point: Self::Point) -> Self::Point;
    fn get_degrees_of_freedom(&mut self) -> Vec<StandardBasis>;
    fn center(&self) -> Self::Point;
    fn area(&self) -> f64;
    fn get_corners(&self) -> Vec<Self::Point>;
}

pub trait Site: Clone + Send + Sync + fmt::Debug {
    type Transform: Clone + Send + Sync + fmt::Debug;

    fn transform(&self) -> Self::Transform;
    fn positions(&self) -> Vec<Self::Transform>;
    fn multiplicity(&self) -> usize;
    fn from_wyckoff(wyckoff: &WyckoffSite) -> Self;
    fn get_basis(&mut self, rot_symmetry: u64) -> Vec<StandardBasis>;
}
