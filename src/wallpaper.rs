//
// wallpaper.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use crate::cell::CrystalFamily;
use crate::symmetry::SymmetryTransform;

/// Defining one of the Crystallographic wallpaper groups.
///
/// This is the highest level description of the symmetry operations of a crystal structure.
///
#[derive(Debug, Clone)]
pub struct Wallpaper {
    pub name: String,
    pub family: CrystalFamily,
}

#[derive(Debug, Clone)]
pub struct WyckoffSite {
    pub letter: char,
    pub symmetries: Vec<SymmetryTransform>,
    pub num_rotations: u64,
    pub mirror_primary: bool,
    pub mirror_secondary: bool,
}

impl WyckoffSite {
    pub fn multiplicity(&self) -> usize {
        self.symmetries.len()
    }

    pub fn degrees_of_freedom(&self) -> &[bool] {
        // TODO implement -> This is only required for the non-general wyckoff sites since all the
        // general sites have 3 degrees-of-freedom.
        //
        // This will be checked as a method of the SymmetryTransform struct.
        &[true, true, true]
    }
}

#[cfg(test)]
mod wyckoff_site_tests {
    use super::*;

    pub fn create_wyckoff() -> WyckoffSite {
        WyckoffSite {
            letter: 'a',
            symmetries: vec![SymmetryTransform::default()],
            num_rotations: 1,
            mirror_primary: false,
            mirror_secondary: false,
        }
    }

    #[test]
    fn multiplicity() {
        let wyckoff = create_wyckoff();
        assert_eq!(wyckoff.multiplicity(), 1);
    }

}
