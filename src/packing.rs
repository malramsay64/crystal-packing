//
// packing.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::cmp::Ordering;
use std::error::Error;
use std::fmt::{Debug, Display};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

use log::{debug, trace};
use serde::{Deserialize, Serialize};

use crate::basis::StandardBasis;
use crate::traits::*;
use crate::wallpaper::{Wallpaper, WyckoffSite};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PackedState<S, C, T>
where
    S: Shape,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    pub wallpaper: Wallpaper,
    pub shape: S,
    pub cell: C,
    occupied_sites: Vec<T>,
}

impl<S, C, T> Eq for PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
}

impl<S, C, T> PartialEq for PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn eq(&self, other: &Self) -> bool {
        self.score() == other.score()
    }
}

impl<S, C, T> PartialOrd for PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.score().partial_cmp(&other.score())
    }
}

impl<S, C, T> Ord for PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.score().partial_cmp(&other.score()).unwrap()
    }
}

impl<S, C, T> PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    pub fn cartesian_positions(&self) -> Vec<S> {
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

    /// Check for intersections of shapes in the current state.
    ///
    /// This checks for intersections between any shapes, checking all occupied sites and their
    /// symmetry defined copies for the current cell and the neighbouring cells. Checking the
    /// neighbouring cells ensures there are no intersections of when tiling space.
    ///
    fn check_intersection(&self) -> bool {
        for (index1, position1) in self.relative_positions().iter().enumerate() {
            let shape_i1 = self
                .shape
                .transform(&self.cell.to_cartesian_isometry(position1));

            // We only need to check the positions after that of index, since the previous ones
            // have already been checked, hence `.skip(index)`
            for (index2, position2) in self.relative_positions().iter().enumerate().skip(index1) {
                trace!("Checking {} against {}", index1, index2);
                for transform in self.cell.periodic_images(position2, index1 != index2) {
                    let shape_i2 = self
                        .shape
                        .transform(&self.cell.to_cartesian_isometry(&transform));

                    if shape_i1.intersects(&shape_i2) {
                        return true;
                    }
                }
            }
        }
        false
    }

    pub fn total_shapes(&self) -> usize {
        self.occupied_sites
            .iter()
            .fold(0, |sum, site| sum + site.multiplicity())
    }

    pub fn score(&self) -> Result<f64, &'static str> {
        if self.check_intersection() {
            Err("Intersection in packing")
        } else {
            Ok((self.shape.area() * self.total_shapes() as f64) / self.cell.area())
        }
    }

    pub fn initialise(
        shape: S,
        wallpaper: Wallpaper,
        isopointal: &[WyckoffSite],
    ) -> PackedState<S, C, T> {
        let num_shapes = isopointal.iter().fold(0, |acc, x| acc + x.multiplicity());
        let max_cell_size = 4. * shape.enclosing_radius() * num_shapes as f64;

        let cell = C::from_family(&wallpaper.family, max_cell_size);

        debug!("Cell: {:?}", cell);

        let mut occupied_sites: Vec<T> = Vec::new();
        for wyckoff in isopointal.iter() {
            let site = T::from_wyckoff(wyckoff);
            occupied_sites.push(site);
        }

        PackedState {
            wallpaper,
            shape,
            cell,
            occupied_sites,
        }
    }

    pub fn generate_basis(&mut self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];
        basis.append(&mut self.cell.get_degrees_of_freedom());
        for site in self.occupied_sites.iter_mut() {
            basis.append(&mut site.get_basis(1));
        }
        basis
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
    use crate::U2::{Cell2, CrystalFamily, LineShape, OccupiedSite, Transform2};
    use approx::assert_abs_diff_eq;

    fn create_square() -> LineShape {
        LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap()
    }

    fn create_wallpaper_p1() -> (Wallpaper, Vec<WyckoffSite>) {
        let wallpaper = Wallpaper {
            name: String::from("p1"),
            family: CrystalFamily::Monoclinic,
        };
        let isopointal = vec![WyckoffSite {
            letter: 'a',
            symmetries: vec![Transform2::from_operations("x,y")],
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
                Transform2::from_operations("x,y"),
                Transform2::from_operations("-x,-y"),
                Transform2::from_operations("-x+1/2,y"),
                Transform2::from_operations("x+1/2,-y"),
            ],
            num_rotations: 1,
            mirror_primary: false,
            mirror_secondary: false,
        }];

        (wallpaper, isopointal)
    }

    fn init_packed_state(group: &str) -> PackedState<LineShape, Cell2, OccupiedSite> {
        let square: LineShape = create_square();

        let (wallpaper, isopointal) = (match group {
            "p1" => Some(create_wallpaper_p1()),
            "p2mg" => Some(create_wallpaper_p2mg()),
            _ => None,
        })
        .unwrap();
        PackedState::initialise(square, wallpaper, &isopointal)
    }

    #[test]
    fn total_shapes_p1() {
        let state = init_packed_state("p1");
        assert_eq!(state.total_shapes(), 1);
    }

    #[test]
    fn packing_fraction_p1() {
        let state = init_packed_state("p1");
        assert_abs_diff_eq!(state.score().unwrap(), 1. / 8.);
    }

    #[test]
    fn total_shapes_p2mg() {
        let state = init_packed_state("p2mg");
        assert_eq!(state.total_shapes(), 4);
    }

    #[test]
    fn packing_fraction_p2mg() {
        let state = init_packed_state("p2mg");
        assert_abs_diff_eq!(state.score().unwrap(), 1. / 32.);
    }

}
