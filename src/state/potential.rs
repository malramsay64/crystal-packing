//
// potential.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

// This is an issue with the Derived traits
#![allow(clippy::type_repetition_in_bounds)]

use std::cmp::Ordering;
use std::fmt::Write;

use anyhow::Error;
use log::debug;
use serde::{Deserialize, Serialize};

use crate::traits::{Potential, Shape, State};
use crate::wallpaper::{Wallpaper, WallpaperGroup, WyckoffSite};
use crate::{Cell2, OccupiedSite, StandardBasis, Transform2};

pub type PotentialState2<S> = PotentialState<S>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PotentialState<S>
where
    S: Shape + Potential,
{
    pub wallpaper: Wallpaper,
    pub shape: S,
    pub cell: Cell2,
    occupied_sites: Vec<OccupiedSite>,
}

impl<S> Eq for PotentialState<S> where S: Shape + Potential {}

impl<S> PartialEq for PotentialState<S>
where
    S: Shape + Potential,
{
    fn eq(&self, other: &Self) -> bool {
        match (self.score(), other.score()) {
            (Some(s), Some(o)) => s.eq(&o),
            (_, _) => false,
        }
    }
}

impl<S> PartialOrd for PotentialState<S>
where
    S: Shape + Potential,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.score(), other.score()) {
            (Some(s), Some(o)) => s.partial_cmp(&o),
            (_, _) => None,
        }
    }
}

impl<S> Ord for PotentialState<S>
where
    S: Shape + Potential,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(&other).unwrap()
    }
}

impl<S> State for PotentialState<S>
where
    S: Shape + Potential,
{
    fn generate_basis(&self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];
        basis.append(&mut self.cell.get_degrees_of_freedom());
        for site in self.occupied_sites.iter() {
            basis.append(&mut site.get_basis(1));
        }
        basis
    }

    fn score(&self) -> Option<f64> {
        let mut sum = 0.;

        // Compare within the current cell
        for (index, shape1) in self
            .cartesian_positions()
            .map(|p| self.shape.transform(&p))
            .enumerate()
        {
            for shape2 in self
                .cartesian_positions()
                .map(|p| self.shape.transform(&p))
                .skip(index + 1)
            {
                sum += shape1.energy(&shape2);
            }
        }

        // Compare in periodic cells
        for shape1 in self.cartesian_positions().map(|p| self.shape.transform(&p)) {
            for position in self.relative_positions() {
                for shape2 in self
                    .cell
                    .periodic_images(position, 3, false)
                    .map(|p| self.shape.transform(&p))
                {
                    sum += shape1.energy(&shape2);
                }
            }
        }
        // We want to minimize the potential energy, so the score we want to maximize is the
        // negation of the potential energy.
        Some(-sum / self.total_shapes() as f64)
    }

    fn total_shapes(&self) -> usize {
        self.occupied_sites
            .iter()
            .fold(0, |sum, site| sum + site.multiplicity())
    }

    fn as_positions(&self) -> Result<String, Error> {
        let mut output = String::new();
        writeln!(&mut output, "{}", self.cell)?;
        writeln!(&mut output, "Positions")?;

        for transform in self.cartesian_positions() {
            writeln!(&mut output, "{:?}", transform)?;
        }
        Ok(output)
    }
}

impl<S> PotentialState<S>
where
    S: Shape + Potential,
{
    pub fn cartesian_positions<'a>(&'a self) -> impl Iterator<Item = Transform2> + 'a {
        self.relative_positions()
            .map(move |position| self.cell.to_cartesian_isometry(&position))
    }

    pub fn relative_positions<'a>(&'a self) -> impl Iterator<Item = Transform2> + 'a {
        self.occupied_sites.iter().flat_map(OccupiedSite::positions)
    }

    pub fn from_group(shape: S, group: &WallpaperGroup) -> Result<Self, Error> {
        let wallpaper = Wallpaper::new(&group);
        let isopointal = &[WyckoffSite::new(group)?];
        Ok(Self::initialise(
            shape.clone(),
            wallpaper.clone(),
            isopointal,
        ))
    }

    pub fn initialise(
        shape: S,
        wallpaper: Wallpaper,
        isopointal: &[WyckoffSite],
    ) -> PotentialState<S> {
        let num_shapes = isopointal.iter().fold(0, |acc, x| acc + x.multiplicity());
        let max_cell_size = 2. * shape.enclosing_radius() * num_shapes as f64;

        let cell = Cell2::from_family(wallpaper.family, max_cell_size);

        debug!("Cell: {:?}", cell);

        let occupied_sites: Vec<_> = isopointal.iter().map(OccupiedSite::from_wyckoff).collect();

        PotentialState {
            wallpaper,
            shape,
            cell,
            occupied_sites,
        }
    }
}

#[cfg(test)]
mod packed_state_tests {
    use super::*;
    use crate::{CrystalFamily, LJShape2, Transform2};

    fn create_wallpaper_p1() -> (Wallpaper, Vec<WyckoffSite>) {
        let wallpaper = Wallpaper {
            name: String::from("p1"),
            family: CrystalFamily::Monoclinic,
        };
        let isopointal = vec![WyckoffSite {
            letter: 'a',
            symmetries: vec![Transform2::from_operations("x,y").unwrap()],
            num_rotations: 1,
            mirror_primary: false,
            mirror_secondary: false,
        }];

        (wallpaper, isopointal)
    }

    fn create_wallpaper_p2mg() -> (Wallpaper, Vec<WyckoffSite>) {
        let wallpaper = Wallpaper {
            name: String::from("p2mg"),
            family: CrystalFamily::Monoclinic,
        };
        let isopointal = vec![WyckoffSite {
            letter: 'd',
            symmetries: vec![
                Transform2::from_operations("x,y").unwrap(),
                Transform2::from_operations("-x,-y").unwrap(),
                Transform2::from_operations("-x+1/2,y").unwrap(),
                Transform2::from_operations("x+1/2,-y").unwrap(),
            ],
            num_rotations: 1,
            mirror_primary: false,
            mirror_secondary: false,
        }];

        (wallpaper, isopointal)
    }

    fn init_state(group: &str) -> PotentialState<LJShape2> {
        let circle = LJShape2::circle();

        let (wallpaper, isopointal) = (match group {
            "p1" => Some(create_wallpaper_p1()),
            "p2mg" => Some(create_wallpaper_p2mg()),
            _ => None,
        })
        .unwrap();
        PotentialState::initialise(circle, wallpaper, &isopointal)
    }

    #[test]
    fn total_shapes_p1() {
        let state = init_state("p1");
        assert_eq!(state.total_shapes(), 1);
    }

    #[test]
    fn total_shapes_p2mg() {
        let state = init_state("p2mg");
        assert_eq!(state.total_shapes(), 4);
    }
}
