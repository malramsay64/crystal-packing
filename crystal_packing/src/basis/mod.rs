//
// basis.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use anyhow::{anyhow, Error};
pub mod shared_value;

pub use shared_value::SharedValue;

#[non_exhaustive]
pub enum Basis<'a> {
    CellBasis {
        value: &'a SharedValue,
    },
    StandardBasis {
        value: &'a SharedValue,
        min: f64,
        max: f64,
    },
}

impl<'a> Basis<'a> {
    pub fn set_value(&self, new_value: f64) -> Result<(), Error> {
        match self {
            Basis::CellBasis { value, .. } => {
                value.set_value(new_value);
                Ok(())
            }
            Basis::StandardBasis { value, min, max } if min <= &new_value && &new_value <= max => {
                value.set_value(new_value);
                Ok(())
            }
            Basis::StandardBasis { .. } => Err(anyhow!("Out of Bounds")),
        }
    }

    pub fn get_value(&self) -> f64 {
        match self {
            Basis::CellBasis { value, .. } => value.get_value(),
            Basis::StandardBasis { value, .. } => value.get_value(),
        }
    }

    pub fn scale(&self) -> f64 {
        1.
    }
}
