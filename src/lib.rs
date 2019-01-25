//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

#[allow(unused_imports)]
#[macro_use]
extern crate approx;
#[macro_use]
extern crate log;
#[macro_use]
extern crate clap;
extern crate nalgebra as na;
extern crate rand;

use nalgebra::{IsometryMatrix2, Point2};
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand::rngs::SmallRng;
use rand::Rng;
use std::error::Error;
use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

pub mod basis;
pub mod cell;
pub mod shape;
pub mod symmetry;
pub mod wallpaper;

pub use crate::basis::{Basis, SharedValue, StandardBasis};
pub use crate::cell::{Cell, CrystalFamily};
pub use crate::shape::{RadialShape, ShapeInstance};
pub use crate::symmetry::SymmetryTransform;
pub use crate::wallpaper::{Wallpaper, WyckoffSite};

#[derive(Clone, Debug)]
struct OccupiedSite {
    wyckoff: WyckoffSite,
    x: SharedValue,
    y: SharedValue,
    angle: SharedValue,
}

impl OccupiedSite {
    fn multiplicity(&self) -> usize {
        self.wyckoff.symmetries.len()
    }

    fn from_wyckoff(wyckoff: &WyckoffSite) -> OccupiedSite {
        let x = SharedValue::new(0.25);
        let y = SharedValue::new(0.25);
        let angle = SharedValue::new(0.);

        OccupiedSite {
            wyckoff: wyckoff.clone(),
            x,
            y,
            angle,
        }
    }

    fn get_basis(&self, rot_symmetry: u64) -> Vec<StandardBasis> {
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
}

#[derive(Clone, Debug)]
pub struct PackedState {
    pub wallpaper: Wallpaper,
    pub shape: RadialShape,
    pub cell: Cell,
    occupied_sites: Vec<OccupiedSite>,
    basis: Vec<StandardBasis>,
}

impl PackedState {
    pub fn cartesian_positions(&self) -> Vec<Vec<shape::Line>> {
        let mut positions: Vec<Vec<shape::Line>> = vec![];
        for position in self.relative_positions().iter() {
            let shape_i = ShapeInstance {
                shape: &self.shape,
                isometry: self.cell.to_cartesian_isometry(position),
            };
            positions.push(shape_i.lines());
        }
        positions
    }

    fn relative_positions(&self) -> Vec<IsometryMatrix2<f64>> {
        let mut transforms: Vec<IsometryMatrix2<f64>> = vec![];
        for site in self.occupied_sites.iter() {
            let position_transform = na::IsometryMatrix2::new(
                na::Vector2::new(site.x.get_value(), site.y.get_value()),
                site.angle.get_value(),
            );
            for symmetry in site.wyckoff.symmetries.iter() {
                transforms.push(symmetry.isometry * position_transform);
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
            let shape_i1 = ShapeInstance {
                shape: &self.shape,
                isometry: self.cell.to_cartesian_isometry(position1),
            };
            // We only need to check the positions after that of index1, since the previous ones
            // have already been checked, hence `.skip(index1)`
            for (index2, position2) in self.relative_positions().iter().enumerate().skip(index1) {
                // The list of periodic images to check. Currently only checking the first shell,
                // i.e. -1, 0, 1. For highly tilted cells checking the second shell may also be
                // nessecary, although this is currently not an issue due to the limiting of the
                // value of the cell angle.
                debug!("Checking {} against {}", index1, index2);
                let periodic_images: &[f64] = &[-1., 0., 1.];
                for x_periodic in periodic_images {
                    for y_periodic in periodic_images {
                        // A shape is always going to intersect with itself. This skips the check
                        // for a shape intersecting with itself, while still checking the periodic
                        // copies.
                        if index1 == index2 && *x_periodic == 0. && *y_periodic == 0. {
                            continue;
                        }
                        let iso = IsometryMatrix2::from_parts(
                            position2.translation * na::Translation2::new(*x_periodic, *y_periodic),
                            position2.rotation,
                        );

                        let shape_i2 = ShapeInstance {
                            shape: &self.shape,
                            isometry: self.cell.to_cartesian_isometry(&iso),
                        };
                        if shape_i1.intersects(&shape_i2) {
                            return true;
                        }
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

    pub fn packing_fraction(&self) -> f64 {
        (self.shape.area() * self.total_shapes() as f64) / self.cell.area()
    }

    pub fn initialise(
        shape: RadialShape,
        wallpaper: Wallpaper,
        isopointal: &[WyckoffSite],
    ) -> PackedState {
        let mut basis: Vec<StandardBasis> = Vec::new();

        let num_shapes = isopointal.iter().fold(0, |acc, x| acc + x.multiplicity());
        let max_cell_size = 4. * shape.max_radius() * num_shapes as f64;

        let cell = Cell::from_family(&wallpaper.family, max_cell_size);
        basis.append(&mut cell.get_basis());

        debug!("Cell: {:?}", cell);

        let mut occupied_sites: Vec<OccupiedSite> = Vec::new();
        for wyckoff in isopointal.iter() {
            let site = OccupiedSite::from_wyckoff(wyckoff);
            basis.append(&mut site.get_basis(shape.rotational_symmetries));
            occupied_sites.push(site);
        }

        PackedState {
            wallpaper,
            shape,
            cell,
            occupied_sites,
            basis,
        }
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

            writeln!(file, "{} {} {}", p0, p1, colour).unwrap();
        }

        for position in self.relative_positions().iter() {
            // The list of periodic images to check. Currently only checking the first shell,
            // i.e. -1, 0, 1. For highly tilted cells checking the second shell may also be
            // nessecary, although this is currently not an issue due to the limiting of the
            // value of the cell angle.
            let (x_r, y_r) = (position.translation.vector.x, position.translation.vector.y);
            let rotation = position.rotation;
            debug!(
                "Relative Position: {} {} {}",
                x_r,
                y_r,
                rotation.angle() * 180. / PI
            );
            let (x_c, y_c) = self.cell.to_cartesian(x_r, y_r);
            debug!("Cartesian Position: {} {}", x_c, y_c);
            let periodic_images: &[f64] = &[-1., 0., 1.];
            for x_periodic in periodic_images {
                for y_periodic in periodic_images {
                    let iso = IsometryMatrix2::from_parts(
                        position.translation * na::Translation2::new(*x_periodic, *y_periodic),
                        position.rotation,
                    );

                    let shape_i = ShapeInstance {
                        shape: &self.shape,
                        isometry: self.cell.to_cartesian_isometry(&iso),
                    };

                    let colour = match () {
                        // This is the 'original' coordinates so differentiate
                        _ if *x_periodic == 0. && *y_periodic == 0. => 'b',
                        _ => 'g',
                    };
                    for line in shape_i.lines() {
                        writeln!(file, "{} {} {}", line.start, line.end, colour).unwrap();
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod packed_state_tests {
    use super::*;

    fn create_square() -> RadialShape {
        RadialShape {
            name: String::from("Square"),
            radial_points: vec![1., 1., 1., 1.],
            rotational_symmetries: 4,
            mirrors: 4,
        }
    }

    fn create_wallpaper_p1() -> (Wallpaper, Vec<WyckoffSite>) {
        let wallpaper = Wallpaper {
            name: String::from("p1"),
            family: CrystalFamily::Monoclinic,
        };
        let isopointal = vec![WyckoffSite {
            letter: 'a',
            symmetries: vec![SymmetryTransform::new("x,y")],
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
                SymmetryTransform::new("x,y"),
                SymmetryTransform::new("-x,-y"),
                SymmetryTransform::new("-x+1/2,y"),
                SymmetryTransform::new("x+1/2,-y"),
            ],
            num_rotations: 1,
            mirror_primary: false,
            mirror_secondary: false,
        }];

        (wallpaper, isopointal)
    }

    fn init_packed_state(group: &str) -> PackedState {
        let square = create_square();

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
        assert_relative_eq!(state.packing_fraction(), 1. / 8.);
    }

    #[test]
    fn total_shapes_p2mg() {
        let state = init_packed_state("p2mg");
        assert_eq!(state.total_shapes(), 4);
    }

    #[test]
    fn packing_fraction_p2mg() {
        let state = init_packed_state("p2mg");
        assert_relative_eq!(state.packing_fraction(), 1. / 32.);
    }

}

pub struct MCVars {
    pub kt_start: f64,
    pub kt_finish: f64,
    pub max_step_size: f64,
    pub num_start_configs: u64,
    pub steps: u64,
}

impl Default for MCVars {
    fn default() -> MCVars {
        MCVars {
            kt_start: 0.1,
            kt_finish: 0.0005,
            max_step_size: 0.01,
            num_start_configs: 32,
            steps: 100,
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

pub fn monte_carlo_best_packing(vars: &MCVars, state: &mut PackedState) -> PackedState {
    let mut rng = SmallRng::seed_from_u64(0);
    let mut rejections: u64 = 0;

    let mut kt: f64 = vars.kt_start;
    let kt_ratio: f64 = vars.kt_ratio();
    let total_shapes: u64 = state.total_shapes() as u64;
    let basis_distribution = Uniform::new(0, state.basis.len() as u64);

    let mut packing: f64 = state.packing_fraction();
    let mut packing_prev: f64 = 0.;
    let mut packing_max: f64 = 0.;

    let mut best_state = state.clone();

    for _ in 0..vars.steps {
        let basis_index: usize = basis_distribution.sample(&mut rng) as usize;
        if let Some(basis_current) = state.basis.get_mut(basis_index) {
            basis_current.set_value(basis_current.sample(&mut rng, vars.max_step_size));
        }

        if state.check_intersection() {
            debug!("Rejected for intersection.");
            rejections += 1;
            state.basis[basis_index].reset_value();
        } else {
            packing = state.packing_fraction();
            if rng.gen::<f64>() > mc_temperature(packing_prev, packing, kt, total_shapes) {
                debug!("Rejected for Increasing packing fraction.");
                rejections += 1;
                state.basis[basis_index].reset_value();
                packing = packing_prev;
            } else {
                // Keep current state, so update previous packing
                packing_prev = packing;
            }
        }
        if packing > packing_max {
            best_state = state.clone();
            packing_max = packing;
        }
        kt *= kt_ratio;
    }
    println!(
        "Packing Fraction: {}, Rejection Percentage: {:.2}%",
        packing_max,
        100. * rejections as f64 / vars.steps as f64,
    );
    best_state
}
