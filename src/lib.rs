//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use std::cmp::Ordering;
use std::error::Error;
use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

extern crate approx;
extern crate clap;
extern crate itertools;
extern crate log;
extern crate nalgebra as na;
extern crate rand;

use itertools::iproduct;
use log::{debug, trace, warn};
use nalgebra::{Point2, Vector2};
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand::rngs::SmallRng;
use rand::Rng;

pub mod basis;
pub mod cell;
pub mod ops_macros;
pub mod shape;
pub mod symmetry;
pub mod wallpaper;

pub use crate::basis::{Basis, SharedValue, StandardBasis};
pub use crate::cell::{Cell, CrystalFamily};
pub use crate::shape::{LineShape, Shape, ShapeInstance};
pub use crate::symmetry::Transform;
pub use crate::wallpaper::{Wallpaper, WyckoffSite};

#[derive(Clone, Debug)]
struct OccupiedSite {
    wyckoff: WyckoffSite,
    x: f64,
    y: f64,
    angle: f64,
}

impl OccupiedSite {
    fn multiplicity(&self) -> usize {
        self.wyckoff.symmetries.len()
    }

    pub fn transform(&self) -> Transform {
        Transform::new(Vector2::new(self.x, self.y), self.angle)
    }

    fn from_wyckoff(wyckoff: &WyckoffSite) -> OccupiedSite {
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

#[derive(Clone, Debug)]
pub struct PackedState<T: shape::Shape> {
    pub wallpaper: Wallpaper,
    pub shape: T,
    pub cell: Cell,
    occupied_sites: Vec<OccupiedSite>,
}

impl<T: shape::Shape> Eq for PackedState<T> {}

impl<T: shape::Shape> PartialEq for PackedState<T> {
    fn eq(&self, other: &Self) -> bool {
        self.packing_fraction() == other.packing_fraction()
    }
}

impl<T: shape::Shape> PartialOrd for PackedState<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.packing_fraction()
            .partial_cmp(&other.packing_fraction())
    }
}

impl<T: shape::Shape> Ord for PackedState<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.packing_fraction()
            .partial_cmp(&other.packing_fraction())
            .unwrap()
    }
}

impl<T: shape::Shape> PackedState<T> {
    pub fn cartesian_positions(&self) -> Vec<Vec<T::Component>> {
        let mut positions: Vec<Vec<T::Component>> = vec![];
        for position in self.relative_positions().iter() {
            let shape_i: ShapeInstance<T::Component> =
                ShapeInstance::from(&self.shape, self.cell.to_cartesian_isometry(position));
            positions.push(shape_i.items);
        }
        positions
    }

    fn relative_positions(&self) -> Vec<Transform> {
        let mut transforms: Vec<Transform> = vec![];
        for site in self.occupied_sites.iter() {
            for symmetry in site.wyckoff.symmetries.iter() {
                transforms.push(symmetry * site.transform());
            }
        }
        transforms
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
            let shape_i1 =
                ShapeInstance::from(&self.shape, self.cell.to_cartesian_isometry(position1));
            for (index2, position2) in self.relative_positions().iter().enumerate().skip(index1) {
                debug!("Checking {} against {}", index1, index2);
                trace!(
                    "Creating shape from: {:?}, results in {:?}",
                    position2,
                    self.cell.to_cartesian_isometry(position2)
                );

                // The periodic images to check. Checking the first and second shells i.e.
                // -2..=2, as this is necessary to ensure no intersections on tilted cells.
                for (x_periodic, y_periodic) in iproduct!(-2..=2, -2..=2) {
                    // A shape is always going to intersect with itself. This skips the check for
                    // a shape intersecting with itself, while still checking the periodic
                    // copies.
                    if index1 == index2 && x_periodic == 0 && y_periodic == 0 {
                        continue;
                    }

                    let periodic_position2 =
                        position2 + Vector2::new(f64::from(x_periodic), f64::from(y_periodic));
                    let shape_i2 = ShapeInstance::from(
                        &self.shape,
                        self.cell.to_cartesian_isometry(&periodic_position2),
                    );
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
            x if 0. < x && x <= 1. => Ok(x),
            _ => Err("Invalid packing fraction"),
        }
    }

    pub fn initialise(
        shape: T,
        wallpaper: Wallpaper,
        isopointal: &[WyckoffSite],
    ) -> PackedState<T> {
        let num_shapes = isopointal.iter().fold(0, |acc, x| acc + x.multiplicity());
        let max_cell_size = 4. * shape.enclosing_radius() * num_shapes as f64;

        let cell = Cell::from_family(&wallpaper.family, max_cell_size);

        debug!("Cell: {:?}", cell);

        let mut occupied_sites: Vec<OccupiedSite> = Vec::new();
        for wyckoff in isopointal.iter() {
            let site = OccupiedSite::from_wyckoff(wyckoff);
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

        // Write cell lines to file
        let pc = na::Translation2::from(-self.cell.center().coords);

        let mut points = vec![
            Point2::new(0., 0.),
            Point2::new(1., 0.),
            Point2::new(1., 1.),
            Point2::new(0., 1.),
        ];

        points = points
            .into_iter()
            .map(|p| pc * self.cell.to_cartesian_point(p))
            .collect();

        for (p0, p1) in points.iter().zip(points.iter().cycle().skip(1)) {
            let colour = 'k';

            writeln!(file, "{}, {}, {}", p0, p1, colour).unwrap();
        }

        for position in self.relative_positions().iter() {
            // The list of periodic images to check. Currently only checking the first shell,
            // i.e. -1, 0, 1. For highly tilted cells checking the second shell may also be
            // necessary, although this is currently not an issue due to the limiting of the
            // value of the cell angle.
            let (x_r, y_r) = (position.translation.x, position.translation.y);
            debug!(
                "Relative Position: {} {} {}",
                x_r,
                y_r,
                position.angle() * 180. / PI
            );
            let (x_c, y_c) = self.cell.to_cartesian(x_r, y_r);
            debug!("Cartesian Position: {} {}", x_c, y_c);
            for (x_periodic, y_periodic) in iproduct!(-1..=1, -1..=1) {
                let p_position =
                    position + Vector2::new(f64::from(x_periodic), f64::from(y_periodic));

                let shape_i =
                    ShapeInstance::from(&self.shape, self.cell.to_cartesian_isometry(&p_position));

                let colour = match () {
                    // These are the 'original' coordinates so differentiate
                    _ if x_periodic == 0 && y_periodic == 0 => 'b',
                    _ => 'g',
                };
                for item in shape_i.items.iter() {
                    writeln!(file, "{}, {}", item, colour).unwrap();
                }
            }
        }
    }
}

#[cfg(test)]
mod packed_state_tests {
    use super::*;
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
            symmetries: vec![Transform::from_operations("x,y")],
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
                Transform::from_operations("x,y"),
                Transform::from_operations("-x,-y"),
                Transform::from_operations("-x+1/2,y"),
                Transform::from_operations("x+1/2,-y"),
            ],
            num_rotations: 1,
            mirror_primary: false,
            mirror_secondary: false,
        }];

        (wallpaper, isopointal)
    }

    fn init_packed_state(group: &str) -> PackedState<shape::LineShape> {
        let square: shape::LineShape = create_square();

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

pub fn monte_carlo_best_packing<T: shape::Shape>(
    vars: &MCVars,
    mut state: PackedState<T>,
) -> Result<PackedState<T>, &'static str> {
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
