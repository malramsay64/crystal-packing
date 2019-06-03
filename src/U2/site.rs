//
// site.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use std::f64::consts::PI;

use nalgebra::Vector2;

use super::Transform2;
use crate::basis::{SharedValue, StandardBasis};
use crate::traits::Site;
use crate::wallpaper::WyckoffSite;

#[derive(Clone, Debug)]
pub struct OccupiedSite {
    wyckoff: WyckoffSite,
    x: f64,
    y: f64,
    angle: f64,
}

impl Site for OccupiedSite {
    type Transform = Transform2;

    fn transform(&self) -> Self::Transform {
        Transform2::new(Vector2::new(self.x, self.y), self.angle)
    }

    fn positions(&self) -> Vec<Self::Transform> {
        let transform = self.transform();
        self.symmetries()
            .iter()
            .map(|sym| sym * transform)
            .collect()
    }

    fn multiplicity(&self) -> usize {
        self.wyckoff.symmetries.len() as usize
    }

    fn from_wyckoff(wyckoff: &WyckoffSite) -> Self {
        let position = -0.5 + 0.5 / wyckoff.multiplicity() as f64;
        let x = position;
        let y = position;
        let angle = 0.;

        OccupiedSite {
            wyckoff: wyckoff.clone(),
            x,
            y,
            angle,
        }
    }

    fn get_basis(&mut self, rot_symmetry: u64) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];
        let dof = self.wyckoff.degrees_of_freedom();

        if dof[0] {
            basis.push(StandardBasis::new(SharedValue::new(&mut self.x), -0.5, 0.5));
        }
        if dof[1] {
            basis.push(StandardBasis::new(SharedValue::new(&mut self.y), -0.5, 0.5));
        }
        if dof[2] {
            basis.push(StandardBasis::new(
                SharedValue::new(&mut self.angle),
                0.,
                2. * PI / rot_symmetry as f64,
            ));
        }
        basis
    }
}

impl OccupiedSite {
    fn symmetries(&self) -> &[Transform2] {
        &self.wyckoff.symmetries
    }
}
