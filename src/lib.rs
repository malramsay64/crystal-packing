//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
// This file is primarily for the definition of the public API of the crate, with most of the
// functions and data structures being defined in submodules. This should additionally mean that
// most of the imports throughout the rest of the crate can just be from the top level and nicely
// grouped together.

extern crate approx;
extern crate clap;
extern crate itertools;
extern crate log;
extern crate nalgebra as na;
extern crate rand;

pub mod ops_macros;

pub mod basis;
pub mod cell;
pub mod optimisation;
pub mod packing;
pub mod shape;
pub mod site;
pub mod traits;
pub mod transform;
pub mod wallpaper;

pub use crate::basis::*;
pub use crate::cell::*;
pub use crate::optimisation::{monte_carlo_best_packing, MCVars};
pub use crate::packing::PackedState;
pub use crate::shape::*;
pub use crate::site::*;
pub use crate::traits::{FromSymmetry, Intersect, Shape};
pub use crate::transform::*;
pub use crate::wallpaper::WallpaperGroup;
