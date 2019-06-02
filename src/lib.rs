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

mod ops_macros;

pub mod U2;
pub mod U3;
pub mod basis;
pub mod packing;
pub mod traits;
pub mod wallpaper;

pub use crate::basis::StandardBasis;
pub use crate::cell::{Cell, CrystalFamily};
pub use crate::packing::{monte_carlo_best_packing, MCVars, PackedState};
pub use crate::shape::*;
pub use crate::site::OccupiedSite;
pub use crate::symmetry::*;
pub use crate::traits::*;
pub use crate::wallpaper::WallpaperGroup;
