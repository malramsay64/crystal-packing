//
// site.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::basis::{SharedValue, StandardBasis};
use crate::wallpaper::WyckoffSite;
use crate::Transform2;

#[derive(Debug, Serialize, Deserialize)]
pub struct OccupiedSite {
    wyckoff: WyckoffSite,
    x: SharedValue,
    y: SharedValue,
    angle: SharedValue,
}

impl Clone for OccupiedSite {
    fn clone(&self) -> Self {
        OccupiedSite {
            wyckoff: self.wyckoff.clone(),
            x: SharedValue::new(self.x.get_value()),
            y: SharedValue::new(self.y.get_value()),
            angle: SharedValue::new(self.angle.get_value()),
        }
    }
}

impl OccupiedSite {
    pub fn transform(&self) -> Transform2 {
        Transform2::new(
            self.angle.get_value(),
            (self.x.get_value(), self.y.get_value()),
        )
    }

    pub fn positions<'a>(&'a self) -> impl Iterator<Item = Transform2> + 'a {
        let transform = self.transform();
        self.symmetries()
            .map(move |sym| sym * &transform)
            .map(|sym| sym.periodic(1., -0.5))
    }

    pub fn multiplicity(&self) -> usize {
        self.wyckoff.symmetries.len() as usize
    }

    pub fn from_wyckoff(wyckoff: &WyckoffSite) -> Self {
        let position = -0.5 + 0.5 / wyckoff.multiplicity() as f64;
        let x = SharedValue::new(position);
        let y = SharedValue::new(position);
        let angle = SharedValue::new(0.);

        OccupiedSite {
            wyckoff: wyckoff.clone(),
            x,
            y,
            angle,
        }
    }

    pub fn get_basis(&self, rot_symmetry: u64) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];
        let dof = self.wyckoff.degrees_of_freedom();

        if dof[0] {
            basis.push(StandardBasis::new(&self.x, -0.5, 0.5));
        }
        if dof[1] {
            basis.push(StandardBasis::new(&self.y, -0.5, 0.5));
        }
        if dof[2] {
            basis.push(StandardBasis::new(
                &self.angle,
                0.,
                2. * PI / rot_symmetry as f64,
            ));
        }
        basis
    }
    pub fn symmetries<'a>(&'a self) -> impl Iterator<Item = &Transform2> + 'a {
        self.wyckoff.symmetries.iter()
    }
}
