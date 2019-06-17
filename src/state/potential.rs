//
// potential.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::cmp::Ordering;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

use log::debug;
use serde::{Deserialize, Serialize};

use crate::traits::*;
use crate::wallpaper::{Wallpaper, WyckoffSite};
use crate::StandardBasis;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    pub wallpaper: Wallpaper,
    pub shape: S,
    pub cell: C,
    occupied_sites: Vec<T>,
}

impl<S, C, T> Eq for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
}

impl<S, C, T> PartialEq for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn eq(&self, other: &Self) -> bool {
        self.score() == other.score()
    }
}

impl<S, C, T> PartialOrd for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.score().partial_cmp(&other.score())
    }
}

impl<S, C, T> Ord for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.score().partial_cmp(&other.score()).unwrap()
    }
}

impl<S, C, T> State for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn generate_basis(&mut self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];
        basis.append(&mut self.cell.get_degrees_of_freedom());
        for site in self.occupied_sites.iter_mut() {
            basis.append(&mut site.get_basis(1));
        }
        basis
    }

    fn score(&self) -> Result<f64, &'static str> {
        let mut sum = 0.;
        for (index1, position1) in self.relative_positions().iter().enumerate() {
            let shape1 = self
                .shape
                .transform(&self.cell.to_cartesian_isometry(position1));
            for (index2, position2) in self.relative_positions().iter().enumerate().skip(index1) {
                for transform in self.cell.periodic_images(position2, index1 != index2) {
                    let shape2 = self
                        .shape
                        .transform(&self.cell.to_cartesian_isometry(&transform));
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
}

impl<S, C, T> PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn cartesian_positions(&self) -> Vec<S> {
        self.relative_positions()
            .iter()
            .map(|position| {
                self.shape
                    .transform(&self.cell.to_cartesian_isometry(position))
            })
            .collect()
    }

    fn relative_positions(&self) -> Vec<T::Transform> {
        self.occupied_sites
            .iter()
            .flat_map(Site::positions)
            .collect()
    }

    pub fn initialise(
        shape: S,
        wallpaper: Wallpaper,
        isopointal: &[WyckoffSite],
    ) -> PotentialState<S, C, T> {
        let num_shapes = isopointal.iter().fold(0, |acc, x| acc + x.multiplicity());
        let max_cell_size = 2. * shape.enclosing_radius() * num_shapes as f64;

        let cell = C::from_family(&wallpaper.family, max_cell_size);

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

    pub fn to_figure(&self, filename: &str) {
        let path = Path::new(filename);
        let display = path.display();
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why.description()),
            Ok(file) => file,
        };

        let points = self.cell.get_corners();

        for (p0, p1) in points.iter().zip(points.iter().cycle().skip(1)) {
            let colour = 'k';
            writeln!(file, "{}, {}, {}", p0, p1, colour).unwrap();
        }

        for position in self.relative_positions().iter() {
            let shape_i = self
                .shape
                .transform(&self.cell.to_cartesian_isometry(&position));

            for item in shape_i.iter() {
                writeln!(file, "{}, b", item).unwrap();
            }

            for transform in self.cell.periodic_images(position, false) {
                let shape_i = self
                    .shape
                    .transform(&self.cell.to_cartesian_isometry(&transform));

                for item in shape_i.iter() {
                    writeln!(file, "{}, g", item).unwrap();
                }
            }
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
