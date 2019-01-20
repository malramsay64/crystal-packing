//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//
extern crate nalgebra;
extern crate rand;

pub mod basis;

use nalgebra::{Matrix2, Vector2};
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::f64::consts::PI;

pub use crate::basis::{Basis, SharedValue, StandardBasis};

/// The different crystal families that can be represented
///
/// These are all the valid types of crystal symmetries which are valid in a 2D space.
///
#[derive(Debug, Clone, PartialEq)]
pub enum CrystalFamily {
    Monoclinic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
}

#[cfg(test)]
mod crystal_family_test {
    use super::*;

    #[test]
    fn equality() {
        assert_eq!(CrystalFamily::Monoclinic, CrystalFamily::Monoclinic);
        assert_ne!(CrystalFamily::Hexagonal, CrystalFamily::Monoclinic);
    }
}

/// Defining one of the Crystallographic wallpaper groups.
///
/// This is the highest level description of the symmetry operations of a crystal structure.
///
#[derive(Debug, Clone)]
pub struct Wallpaper {
    pub name: String,
    pub family: CrystalFamily,
}

/// Define the transformations of particle positions
///
/// These
#[derive(Debug, Clone)]
pub struct SymmetryTransform {
    pub rotation: Matrix2<f64>,
    pub translation: Vector2<f64>,
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
    fn multiplicity(&self) -> usize {
        return self.symmetries.len();
    }

    fn degrees_of_freedom(&self) -> &[bool] {
        // TODO implement
        // Check x
        // Check y
        // Check rotations
        return &[true, true, true];
    }
}

#[cfg(test)]
mod wyckoff_site_tests {
    use super::*;

    pub fn create_wyckoff() -> WyckoffSite {
        WyckoffSite {
            letter: 'a',
            symmetries: vec![SymmetryTransform {
                rotation: Matrix2::new(1., 0., 1., 0.),
                translation: Vector2::new(0., 0.),
            }],
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

#[derive(Debug, Clone)]
pub struct Shape {
    pub name: String,
    pub radial_points: Vec<f64>,
    pub rotational_symmetries: u64,
    pub mirrors: u64,
}

pub struct ShapeIter<'a> {
    shape: &'a Shape,
    index: usize,
    len: usize,
}

impl<'a> ShapeIter<'a> {
    fn new(shape: &'a Shape) -> Self {
        Self {
            shape,
            index: 0,
            len: shape.radial_points.len(),
        }
    }
}

impl<'a> Iterator for ShapeIter<'a> {
    type Item = (f64, f64);

    fn next(&mut self) -> Option<(f64, f64)> {
        if self.index >= self.len {
            return None;
        }
        let result = Some((
            self.shape.radial_points[self.index],
            self.shape.radial_points[(self.index + 1) % self.len],
        ));
        self.index += 1;
        result
    }
}

impl<'a> IntoIterator for &'a Shape {
    type Item = (f64, f64);
    type IntoIter = ShapeIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        ShapeIter::new(&self)
    }
}

impl Shape {
    fn area(&self) -> f64 {
        // This is the sine of the angle between each point, this is used for every calculation
        // so pre-calculate here.
        let angle_term: f64 = f64::sin(2. * PI / self.radial_points.len() as f64);

        self.into_iter()
            .map(|(a, b)| 0.5 * angle_term * a * b)
            .sum()
    }

    fn max_radius(&self) -> f64 {
        return self
            .radial_points
            .iter()
            .cloned()
            // The f64 type doesn't have complete ordering because of Nan and Inf, so the
            // standard min/max comparators don't work. Instead we use the f64::max which ignores
            // the NAN and max values.
            .fold(std::f64::MIN, f64::max);
    }
}

#[cfg(test)]
mod shape_tests {
    use super::*;

    pub fn create_square() -> Shape {
        Shape {
            name: String::from("Square"),
            radial_points: vec![1., 1., 1., 1.],
            rotational_symmetries: 4,
            mirrors: 4,
        }
    }

    #[test]
    fn init() {
        let square = create_square();
        assert_eq!(square.name, "Square");
        assert_eq!(square.radial_points, vec![1., 1., 1., 1.]);
        assert_eq!(square.rotational_symmetries, 4);
        assert_eq!(square.mirrors, 4);
    }

    #[test]
    fn area() {
        let square = create_square();
        assert_eq!(square.area(), 2.);
    }

    #[test]
    fn max_radius() {
        let shape = Shape {
            name: String::from("iter_test"),
            radial_points: vec![1., 2., 3., 4.],
            rotational_symmetries: 1,
            mirrors: 0,
        };
        assert_eq!(shape.max_radius(), 4.);
        assert_eq!(shape.max_radius(), 4.);
    }

    #[test]
    fn iter_values() {
        let shape = Shape {
            name: String::from("iter_test"),
            radial_points: vec![1., 2., 3., 4.],
            rotational_symmetries: 1,
            mirrors: 0,
        };
        let manual = vec![(1., 2.), (2., 3.), (3., 4.), (4., 1.)];
        assert_eq!(shape.into_iter().collect::<Vec<(f64, f64)>>(), manual);
    }
}

#[derive(Clone)]
struct OccupiedSite {
    wyckoff: WyckoffSite,
    x: SharedValue,
    y: SharedValue,
    angle: SharedValue,
}

impl OccupiedSite {
    fn multiplicity(&self) -> usize {
        return self.wyckoff.symmetries.len();
    }

    fn from_wyckoff(wyckoff: &WyckoffSite) -> OccupiedSite {
        let x = SharedValue::new(0.);
        let y = SharedValue::new(0.);
        let angle = SharedValue::new(0.);

        return OccupiedSite {
            wyckoff: wyckoff.clone(),
            x,
            y,
            angle,
        };
    }

    fn get_basis(&self, rot_symmetry: u64) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];
        let dof = self.wyckoff.degrees_of_freedom();

        if dof[0] {
            basis.push(StandardBasis::new(&self.x, 0., 1.));
        }
        if dof[1] {
            basis.push(StandardBasis::new(&self.y, 0., 1.));
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

#[derive(Clone)]
struct Cell {
    x_len: SharedValue,
    y_len: SharedValue,
    angle: SharedValue,
    family: CrystalFamily,
}
impl Cell {
    fn from_family(family: &CrystalFamily, length: f64) -> Cell {
        let (x_len, y_len, angle) = match family {
            // The Hexagonal Crystal has both sides equal with a fixed angle of 60 degrees.
            CrystalFamily::Hexagonal => {
                let len = SharedValue::new(length);
                (len.clone(), len.clone(), SharedValue::new(PI / 3.))
            }
            // The Tetragonal Crystal has both sides equal with a fixed angle of 90 degrees.
            CrystalFamily::Tetragonal => {
                let len = SharedValue::new(length);
                (len.clone(), len.clone(), SharedValue::new(PI / 2.))
            }
            // The Orthorhombic crystal has two variable sides with a fixed angle of 90 degrees.
            CrystalFamily::Orthorhombic => (
                SharedValue::new(length),
                SharedValue::new(length),
                SharedValue::new(PI / 2.),
            ),
            // The Monoclinic cell has two variable sides and a variable angle initialised to 90
            // degrees
            CrystalFamily::Monoclinic => (
                SharedValue::new(length),
                SharedValue::new(length),
                SharedValue::new(PI / 2.),
            ),
        };
        return Cell {
            x_len,
            y_len,
            angle,
            family: family.clone(),
        };
    }

    fn get_basis(&self) -> Vec<StandardBasis> {
        let mut basis: Vec<StandardBasis> = vec![];

        // All cells have at least a single variable cell length
        basis.push(StandardBasis::new(
            &self.x_len,
            0.01,
            self.x_len.get_value(),
        ));

        // Both the Orthorhombic and Monoclinic cells have a second variable cell length
        if (self.family == CrystalFamily::Orthorhombic) | (self.family == CrystalFamily::Monoclinic)
        {
            basis.push(StandardBasis::new(
                &self.y_len,
                0.01,
                self.y_len.get_value(),
            ));
        }

        // The Monoclinic family is the only one to have a variable cell angle.
        if self.family == CrystalFamily::Monoclinic {
            basis.push(StandardBasis::new(&self.angle, PI / 4., 3. * PI / 4.));
        }

        basis
    }

    fn area(&self) -> f64 {
        self.angle.get_value().sin() * self.x_len.get_value() * self.y_len.get_value()
    }
}

#[derive(Clone)]
pub struct PackedState {
    wallpaper: Wallpaper,
    shape: Shape,
    cell: Cell,
    occupied_sites: Vec<OccupiedSite>,
    basis: Vec<StandardBasis>,
}

impl PackedState {
    pub fn check_intersection(&self) -> bool {
        // TODO Implement
        return true;
    }

    pub fn total_shapes(&self) -> usize {
        let mut sum: usize = 0;
        for site in self.occupied_sites.iter() {
            sum += site.multiplicity();
        }
        return sum;
    }

    pub fn packing_fraction(&self) -> f64 {
        self.cell.area() / (self.shape.area() * self.total_shapes() as f64)
    }

    pub fn initialise(
        shape: Shape,
        wallpaper: Wallpaper,
        isopointal: Vec<WyckoffSite>,
        step_size: f64,
    ) -> PackedState {
        let mut basis: Vec<StandardBasis> = Vec::new();

        let num_shapes = isopointal.iter().fold(0, |acc, x| acc + x.multiplicity());
        let max_cell_size = 4. * shape.max_radius() * num_shapes as f64;

        let cell = Cell::from_family(&wallpaper.family, max_cell_size);
        basis.append(&mut cell.get_basis());

        let mut occupied_sites: Vec<OccupiedSite> = Vec::new();
        for wyckoff in isopointal.iter() {
            let site = OccupiedSite::from_wyckoff(wyckoff);
            basis.append(&mut site.get_basis(shape.rotational_symmetries));
            occupied_sites.push(site);
        }

        return PackedState {
            wallpaper,
            shape,
            cell,
            occupied_sites,
            basis,
        };
    }
}

#[cfg(test)]
mod packed_state_tests {
    use super::*;

    fn init_packed_state() -> PackedState {
        let square = shape_tests::create_square();

        let wallpaper = Wallpaper {
            name: String::from("p1"),
            family: CrystalFamily::Monoclinic,
        };

        let isopointal = vec![WyckoffSite {
            letter: 'a',
            symmetries: vec![SymmetryTransform {
                rotation: Matrix2::new(1., 0., 1., 0.),
                translation: Vector2::new(0., 0.),
            }],
            num_rotations: 1,
            mirror_primary: false,
            mirror_secondary: false,
        }];

        PackedState::initialise(square, wallpaper, isopointal, 0.1)
    }

    #[test]
    fn packed_state_total_shapes() {
        let state = init_packed_state();
        assert_eq!(state.total_shapes(), 1);
    }

}

struct MCVars {
    kt_start: f64,
    kt_finish: f64,
    max_step_size: f64,
    num_start_configs: u64,
    steps: u64,
}

impl MCVars {
    fn kt_ratio(&self) -> f64 {
        return f64::powf(self.kt_finish / self.kt_start, 1.0 / self.steps as f64);
    }
}

fn mc_temperature(old: f64, new: f64, kt: f64, n: u64) -> f64 {
    return f64::exp((1. / old - 1. / new) / kt) * (old / new).powi(n as i32);
}

fn monte_carlo_best_packing<'a, 'b>(vars: &'a MCVars, state: &'b mut PackedState) -> PackedState {
    let mut rng = rand::thread_rng();
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
            basis_current.set_value(basis_current.sample(&mut rng));
        }

        if state.check_intersection() {
            rejections += 1;
            state.basis[basis_index].reset_value();
        } else {
            packing = state.packing_fraction();
            if rng.gen::<f64>() > mc_temperature(packing_prev, packing, kt, total_shapes) {
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
        "Packing Fraction: {}, Rejection Percentage {}",
        packing_max,
        rejections as f64 / vars.steps as f64,
    );
    return best_state;
}
