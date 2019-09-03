//
// potential.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::cmp::Ordering;
use std::fmt;
use std::fmt::Write;

use log::debug;
use serde::{Deserialize, Serialize};

use crate::traits::*;
use crate::wallpaper::{Wallpaper, WallpaperGroup, WyckoffSite};
use crate::{Cell2, OccupiedSite, StandardBasis, Transform2};

pub type PotentialState2<S> = PotentialState<S, Cell2, OccupiedSite>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
    pub wallpaper: Wallpaper,
    pub shape: S,
    pub cell: C,
    occupied_sites: Vec<T>,
}

impl<S, C, T> Eq for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
}

impl<S, C, T> PartialEq for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
    fn eq(&self, other: &Self) -> bool {
        self.score() == other.score()
    }
}

impl<S, C, T> PartialOrd for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.score().partial_cmp(&other.score())
    }
}

impl<S, C, T> Ord for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.score().partial_cmp(&other.score()).unwrap()
    }
}

impl<S, C, T> State for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
    fn generate_basis(&self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];
        basis.append(&mut self.cell.get_degrees_of_freedom());
        for site in self.occupied_sites.iter() {
            basis.append(&mut site.get_basis(1));
        }
        basis
    }

    fn score(&self) -> Result<f64, &'static str> {
        let mut sum = 0.;
        for (index1, position1) in self.relative_positions().enumerate() {
            let shape1 = self
                .shape
                .transform(&self.cell.to_cartesian_isometry(&position1));
            for (index2, position2) in self.relative_positions().enumerate().skip(index1) {
                for transform in self.cell.periodic_images(position2, index1 != index2) {
                    let shape2 = self.shape.transform(&transform);
                    sum += shape1.energy(&shape2)
                }
            }
        }
        // We want to minimize the potential energy, so the score we want to maximize is the
        // negation of the potential energy.
        Ok(-sum / 15. / self.total_shapes() as f64)
    }

    fn total_shapes(&self) -> usize {
        self.occupied_sites
            .iter()
            .fold(0, |sum, site| sum + site.multiplicity())
    }

    fn as_positions(&self) -> Result<String, fmt::Error> {
        let mut output = String::new();
        writeln!(&mut output, "{}", self.cell)?;
        writeln!(&mut output, "Positions")?;

        for transform in self.cartesian_positions() {
            writeln!(&mut output, "{:?}", transform)?;
        }
        Ok(output)
    }
}

impl<S, C, T> PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
    pub fn cartesian_positions<'a>(&'a self) -> impl Iterator<Item = Transform2> + 'a {
        self.relative_positions()
            .map(move |position| self.cell.to_cartesian_isometry(&position))
    }

    pub fn relative_positions<'a>(&'a self) -> impl Iterator<Item = Transform2> + 'a {
        self.occupied_sites.iter().flat_map(Site::positions)
    }

    pub fn from_group(shape: S, group: &WallpaperGroup) -> Self {
        let wallpaper = Wallpaper::new(&group);
        let isopointal = &[WyckoffSite::new(group.clone())];
        Self::initialise(shape.clone(), wallpaper.clone(), isopointal)
    }

    pub fn initialise(
        shape: S,
        wallpaper: Wallpaper,
        isopointal: &[WyckoffSite],
    ) -> PotentialState<S, C, T> {
        let num_shapes = isopointal.iter().fold(0, |acc, x| acc + x.multiplicity());
        let max_cell_size = 2. * shape.enclosing_radius() * num_shapes as f64;

        let cell = C::from_family(wallpaper.family, max_cell_size);

        debug!("Cell: {:?}", cell);

        let mut occupied_sites: Vec<T> = Vec::new();
        for wyckoff in isopointal.iter() {
            let site = T::from_wyckoff(wyckoff);
            occupied_sites.push(site);
        }

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
    use crate::{Cell2, CrystalFamily, LJShape2, OccupiedSite, Transform2};

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

    fn init_state(group: &str) -> PotentialState<LJShape2, Cell2, OccupiedSite> {
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
