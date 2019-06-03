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

use log::{debug, trace, warn};
use rand::distributions::Uniform;
use rand::prelude::*;
use rand::rngs::SmallRng;

use crate::basis::{Basis, StandardBasis};
use crate::traits::*;
use crate::wallpaper::{Wallpaper, WyckoffSite};

#[derive(Clone, Debug)]
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
        self.packing_fraction() == other.packing_fraction()
    }
}

impl<S, C, T> PartialOrd for PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.packing_fraction()
            .partial_cmp(&other.packing_fraction())
    }
}

impl<S, C, T> Ord for PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.packing_fraction()
            .partial_cmp(&other.packing_fraction())
            .unwrap()
    }
}

impl<S, C, T> PackedState<S, C, T>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    pub fn cartesian_positions(&self) -> Vec<S> {
        let mut positions: Vec<S> = vec![];
        for position in self.relative_positions().iter() {
            let shape_i = self
                .shape
                .transform(&self.cell.to_cartesian_isometry(position));
            positions.push(shape_i);
        }
        positions
    }

    fn relative_positions(&self) -> Vec<T::Transform> {
        self.occupied_sites
            .iter()
            .flat_map(|site| site.positions())
            .collect()
    }

    /// Check for intersections of shapes in the current state.
    ///
    /// This checks for intersections between any shapes, checking all occupied sites and their
    /// symmetry defined copies for the current cell and the neighbouring cells. Checking the
    /// neighbouring cells ensures there are no intersections of when tiling space.
    ///
    pub fn check_intersection(&self) -> bool {
        for (index1, position1) in self.relative_positions().iter().enumerate() {
            // We only need to check the positions after that of index, since the previous ones
            // have already been checked, hence `.skip(index)`
            trace!(
                "Creating shape from: {:?}, results in {:?}",
                position1,
                self.cell.to_cartesian_isometry(position1)
            );
            let shape_i1 = self
                .shape
                .transform(&self.cell.to_cartesian_isometry(position1));

            for (index2, position2) in self.relative_positions().iter().enumerate().skip(index1) {
                debug!("Checking {} against {}", index1, index2);
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

    pub fn packing_fraction(&self) -> Result<f64, &'static str> {
        match (self.shape.area() * self.total_shapes() as f64) / self.cell.area() {
            x => Ok(x),
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

    fn generate_basis(&mut self) -> Vec<StandardBasis> {
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
            // The list of periodic images to check. Currently only checking the first shell,
            // i.e. -1, 0, 1. For highly tilted cells checking the second shell may also be
            // necessary, although this is currently not an issue due to the limiting of the
            // value of the cell angle.

            let shape_i = self.shape.clone().transform(&position);

            for item in shape_i.iter() {
                writeln!(file, "{}, b", item).unwrap();
            }

            for transform in self.cell.periodic_images(position, false) {
                let shape_i = self
                    .shape
                    .clone()
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
    use crate::U2::LineShape;
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
            symmetries: vec![na::Transform2::from_operations("x,y")],
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

    fn init_packed_state(group: &str) -> PackedState<LineShape> {
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
        assert_abs_diff_eq!(state.packing_fraction().unwrap(), 1. / 8.);
    }

    #[test]
    fn total_shapes_p2mg() {
        let state = init_packed_state("p2mg");
        assert_eq!(state.total_shapes(), 4);
    }

    #[test]
    fn packing_fraction_p2mg() {
        let state = init_packed_state("p2mg");
        assert_abs_diff_eq!(state.packing_fraction().unwrap(), 1. / 32.);
    }

}

pub struct MCVars {
    pub kt_start: f64,
    pub kt_finish: f64,
    pub max_step_size: f64,
    pub num_start_configs: u64,
    pub steps: u64,
    pub seed: Option<u64>,
}

impl Default for MCVars {
    fn default() -> MCVars {
        MCVars {
            kt_start: 0.1,
            kt_finish: 0.0005,
            max_step_size: 0.01,
            num_start_configs: 32,
            steps: 100,
            seed: None,
        }
    }
}

impl MCVars {
    fn kt_ratio(&self) -> f64 {
        f64::powf(self.kt_finish / self.kt_start, 1.0 / self.steps as f64)
    }
}

fn mc_temperature(old: f64, new: f64, kt: f64, n: u64) -> f64 {
    f64::exp((1. / old - 1. / new) / kt) * (old / new).powi(n as i32)
}

pub fn monte_carlo_best_packing<S, C, T>(
    vars: &MCVars,
    mut state: PackedState<S, C, T>,
) -> Result<PackedState<S, C, T>, &'static str>
where
    S: Shape + Debug + Display,
    C: Cell<Transform = S::Transform>,
    T: Site<Transform = S::Transform>,
{
    // When a random seed is provided, use it, otherwise seed the random number generator from the
    // system entropy.
    let mut rng = match vars.seed {
        Some(x) => SmallRng::seed_from_u64(x),
        None => SmallRng::from_entropy(),
    };
    let mut rejections: u64 = 0;

    let mut kt: f64 = vars.kt_start;
    let kt_ratio: f64 = vars.kt_ratio();
    let total_shapes: u64 = state.total_shapes() as u64;

    let mut basis = state.generate_basis();
    let basis_distribution = Uniform::new(0, basis.len() as u64);

    let mut packing: f64 = state.packing_fraction()?;

    let mut packing_prev: f64 = 0.;
    let mut packing_max: f64 = 0.;

    let mut best_state = state.clone();

    for _ in 0..vars.steps {
        let basis_index: usize = basis_distribution.sample(&mut rng) as usize;
        if let Some(basis_current) = basis.get_mut(basis_index) {
            basis_current.set_value(basis_current.sample(&mut rng, vars.max_step_size));
        }

        if state.check_intersection() {
            debug!("Rejected for intersection.");
            rejections += 1;
            basis[basis_index].reset_value();
        } else {
            packing = match state.packing_fraction() {
                Err(_) => {
                    warn!("Rejected for invalid packing fraction.");
                    rejections += 1;
                    basis[basis_index].reset_value();

                    // Set packing to it's previous value
                    packing_prev
                }
                Ok(new_packing) => {
                    if rng.gen::<f64>()
                        > mc_temperature(packing_prev, new_packing, kt, total_shapes)
                    {
                        // Packing fraction was increased too much so reject the step
                        debug!("Rejected for Increasing packing fraction.");
                        rejections += 1;
                        basis[basis_index].reset_value();

                        // Set packing to it's previous value
                        packing_prev
                    } else {
                        // This is where we update the packing fraction cause the test was
                        // successful
                        packing_prev = new_packing;
                        new_packing
                    }
                }
            };
        }
        if packing > packing_max {
            best_state = state.clone();
            packing_max = packing;
        }
        kt *= kt_ratio;
    }
    println!(
        "Packing Fraction: {:.4}, Rejections: {:.2} %",
        packing_max,
        100. * rejections as f64 / vars.steps as f64,
    );
    Ok(best_state)
}
