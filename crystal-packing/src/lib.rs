//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
// This file is primarily for the definition of the public API of the crate, with most of the
// functions and data structures being defined in submodules. This should additionally mean that
// most of the imports throughout the rest of the crate can just be from the top level and nicely
// grouped together.

pub mod basis;
pub mod cell;
pub mod engine;
pub mod geometry;
pub mod hard_shape;
pub mod lattice;
pub mod monte_carlo;
pub mod ops_macros;
pub mod optimisation;
pub mod random;
pub mod shape;
pub mod site;
pub mod state;
pub mod symmetry;
pub mod to_svg;
pub mod traits;
pub mod transform;
pub mod wallpaper;

pub use crate::basis::*;
pub use crate::cell::*;
pub use crate::optimisation::*;
pub use crate::random::Rng;
pub use crate::shape::*;
pub use crate::site::*;
pub use crate::state::*;
pub use crate::traits::{FromSymmetry, Intersect, Shape};
pub use crate::transform::Transform2;
pub use crate::wallpaper::WallpaperGroup;
