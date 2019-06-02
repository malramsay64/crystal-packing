//
// shape.rs
//
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fmt;
use std::slice;

pub use super::Transform2;

mod atom;
mod atom_ops;

mod line;
mod line_ops;

mod line_shape;
mod molecular_shape;

pub use atom::*;
pub use line::*;
pub use line_shape::*;
pub use molecular_shape::*;
