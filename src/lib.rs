//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//
extern crate nalgebra as na;
extern crate rand;

use nalgebra::{Matrix2, Vector2};
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::f64::consts::PI;

/// A Value which can be modified in many places
///
/// Rust's ownership rules preclude multiple mutable borrows from taking place. The SharedValue
/// struct is a wrapper around some of Rust's built in methods of handling this, providing a
/// nicer interface for that tasks I require.
///
pub struct SharedValue {
    value: std::rc::Rc<std::cell::RefCell<f64>>,
}

impl Clone for SharedValue {
    fn clone(&self) -> Self {
        return Self {
            value: std::rc::Rc::clone(&self.value),
        };
    }
}

impl SharedValue {
    pub fn get_value(&self) -> f64 {
        return *self.value.borrow();
    }
    pub fn set_value(&self, value: f64) {
        *self.value.borrow_mut() = value;
    }

    pub fn new(val: f64) -> Self {
        return Self {
            value: std::rc::Rc::new(std::cell::RefCell::new(val)),
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_shared_value() {
        let value = SharedValue::new(1.);
        assert_eq!(value.get_value(), 1.);
    }

    #[test]
    fn update_shared_value() {
        let value = SharedValue::new(1.);
        assert_eq!(value.get_value(), 1.);
        value.set_value(0.5);
        assert_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn cloned_shared_value() {
        let value1 = SharedValue::new(1.);
        let value2 = value1.clone();
        assert_eq!(value1.get_value(), 1.);
        assert_eq!(value2.get_value(), 1.);

        value2.set_value(0.5);
        assert_eq!(value1.get_value(), 0.5);
        assert_eq!(value2.get_value(), 0.5);
    }
}

pub trait Basis {
    fn set_value(&mut self, new_value: f64);
    fn get_value(&self) -> f64;
    fn reset_value(&self);
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64;
}

#[derive(Clone)]
pub struct StandardBasis {
    value: SharedValue,
    old: f64,
    min: f64,
    max: f64,
}

impl StandardBasis {
    pub fn new(value: f64, min: f64, max: f64) -> Self {
        Self {
            value: SharedValue::new(value),
            old: value,
            min,
            max,
        }
    }

    fn value_range(&self) -> f64 {
        return self.max - self.min;
    }
}

impl Basis for StandardBasis {
    fn get_value(&self) -> f64 {
        return self.value.get_value();
    }

    fn set_value(&mut self, new_value: f64) {
        self.old = self.get_value();
        self.value.set_value(match new_value {
            x if x < self.min => self.min,
            x if x > self.max => self.max,
            x => x,
        })
    }

    fn reset_value(&self) {
        self.value.set_value(self.old);
    }

    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        return self.get_value() + self.value_range() * rng.gen_range(-0.5, 0.5);
    }
}

/// The different crystal families that can be represented
///
/// These are all the valid types of crystal symmetries which are valid in a 2D space.
///
#[derive(Debug, Clone)]
enum CrystalFamily {
    Monoclinic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
}

/// Defining one of the Crystallographic wallpaper groups.
///
/// This is the highest level description of the symmetry operations of a crystal structure.
///
#[derive(Debug, Clone)]
struct Wallpaper {
    name: String,
    family: CrystalFamily,
}

/// Define the transformations of particle positions
///
/// These
#[derive(Debug, Clone)]
struct SymmetryTransform {
    rotation: Matrix2<f64>,
    translation: Vector2<f64>,
}

#[derive(Debug, Clone)]
struct WyckoffSite {
    letter: char,
    symmetries: Vec<SymmetryTransform>,
    num_rotations: u64,
    mirror_primary: bool,
    mirror_secondary: bool,
}

impl WyckoffSite {
    fn multiplicity(&self) -> usize {
        return self.symmetries.len();
    }

    fn degrees_of_freedom(&self) -> &[bool] {
        // Check x
        // Check y
        // Check rotations
        return &[true, true, true];
    }
}

#[derive(Debug, Clone)]
struct Shape {
    name: String,
    radial_points: Vec<f64>,
    rotational_symmetries: u64,
    mirrors: u64,
}

impl Shape {
    fn area(&self) -> f64 {
        // This is the sine of the angle between each point
        let angle_term: f64 = f64::sin(2. * PI / self.radial_points.len() as f64);
        return self
            .radial_points
            .iter()
            // Get the next value of the iterator to compare, cycling back to the start
            .zip(&self.radial_points.iter().cycle().nth(1))
            // Calculate the area of each triangle
            .map(|(a, b)| a * *b * angle_term)
            // Sum everything together
            .sum();
    }

    fn max_radius(&self) -> f64 {
        return self.radial_points.iter().cloned().fold(0. / 0., f64::max);
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
        let basis: &[StandardBasis] = &[
            StandardBasis {
                value: self.x.clone(),
                old: self.x.get_value(),
                min: 0.,
                max: 1.,
            },
            StandardBasis {
                value: self.y.clone(),
                old: self.y.get_value(),
                min: 0.,
                max: 1.,
            },
            StandardBasis {
                value: self.angle.clone(),
                old: self.angle.get_value(),
                min: 0.,
                max: 2. * PI / rot_symmetry as f64,
            },
        ];

        return self
            .wyckoff
            .degrees_of_freedom()
            .iter()
            .zip(basis)
            .filter(|&(b, _)| *b)
            .map(|(_, val)| val.clone())
            .collect::<Vec<StandardBasis>>();
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
        let mut basis: Vec<StandardBasis> = vec![
            StandardBasis {
                value: self.x_len.clone(),
                old: self.x_len.get_value(),
                min: 0.01,
                max: self.x_len.get_value(),
            },
            StandardBasis {
                value: self.y_len.clone(),
                old: self.y_len.get_value(),
                min: 0.01,
                max: self.y_len.get_value(),
            },
            StandardBasis {
                value: self.angle.clone(),
                old: self.angle.get_value(),
                min: PI / 4.,
                max: 3. * PI / 4.,
            },
        ];

        match self.family {
            CrystalFamily::Hexagonal | CrystalFamily::Tetragonal => basis.truncate(1),
            CrystalFamily::Orthorhombic => basis.truncate(2),
            CrystalFamily::Monoclinic => basis.truncate(3),
        };

        return basis;
    }
}

#[derive(Clone)]
struct PackedState {
    wallpaper: Wallpaper,
    shape: Shape,
    cell: Cell,
    occupied_sites: Vec<OccupiedSite>,
    basis: Vec<StandardBasis>,
}

impl PackedState {
    fn check_intersection(&self) -> bool {
        // TODO Implement
        return true;
    }

    fn total_shapes(&self) -> usize {
        let mut sum: usize = 0;
        for site in self.occupied_sites.iter() {
            sum += site.multiplicity();
        }
        return sum;
    }

    fn packing_fraction(&self) -> f64 {
        // TODO Implement
        return 0.;
    }

    fn initialise(
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
