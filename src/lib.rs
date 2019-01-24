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

extern crate nalgebra as na;
extern crate rand;

pub mod basis;
pub mod shape;

use nalgebra::{IsometryMatrix2, Matrix2, Point2, Vector2};
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::f64::consts::PI;

pub use crate::basis::{Basis, SharedValue, StandardBasis};
pub use crate::shape::RadialShape;

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
        assert_eq!(CrystalFamily::Orthorhombic, CrystalFamily::Orthorhombic);
        assert_eq!(CrystalFamily::Hexagonal, CrystalFamily::Hexagonal);
        assert_eq!(CrystalFamily::Tetragonal, CrystalFamily::Tetragonal);
    }

    #[test]
    fn inequality() {
        assert_ne!(CrystalFamily::Orthorhombic, CrystalFamily::Monoclinic);
        assert_ne!(CrystalFamily::Hexagonal, CrystalFamily::Monoclinic);
        assert_ne!(CrystalFamily::Tetragonal, CrystalFamily::Monoclinic);
        assert_ne!(CrystalFamily::Hexagonal, CrystalFamily::Orthorhombic);
        assert_ne!(CrystalFamily::Tetragonal, CrystalFamily::Orthorhombic);
        assert_ne!(CrystalFamily::Tetragonal, CrystalFamily::Hexagonal);
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
    isometry: IsometryMatrix2<f64>,
}

impl SymmetryTransform {
    fn parse_ops(ops: &str) -> (Vector2<f64>, f64) {
        let mut vec = Vector2::zeros();
        let mut sign = 1.;
        let mut constant = 0.;
        let mut operator: Option<char> = None;
        for c in ops.chars() {
            match c {
                'x' => {
                    vec[0] = sign;
                    sign = 1.;
                }
                'y' => {
                    vec[1] = sign;
                    sign = 1.;
                }
                '*' | '/' => {
                    operator = Some(c);
                }
                '-' => {
                    sign = -1.;
                }
                // This matches all digits from 0 to 9
                '0'...'9' => {
                    let val = c.to_string().parse::<u64>().unwrap() as f64;
                    // Is there an operator defined, i.e. is this the first digit
                    if let Some(op) = operator {
                        constant = match op {
                            '/' => sign * constant / val,
                            '*' => sign * constant * val,
                            _ => 0.,
                        };
                        operator = None;
                    } else {
                        constant = sign * val;
                    }
                    sign = 1.
                }
                // Default is do nothing (shouldn't encounter this at all)
                _ => {}
            };
        }

        (vec, constant)
    }

    pub fn new(sym_ops: &str) -> SymmetryTransform {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();
        let mut trans = Vector2::new(0., 0.);
        let mut rot: Matrix2<f64> = Matrix2::new(1., 0., 0., 1.);

        for (index, op) in operations.iter().enumerate() {
            let (r, t) = SymmetryTransform::parse_ops(op);
            rot.set_row(index, &r.transpose());
            trans[index] = t;
        }
        SymmetryTransform {
            isometry: IsometryMatrix2::from_parts(
                na::Translation2::from(trans),
                na::Rotation2::from_matrix_unchecked(rot),
            ),
        }
    }

    pub fn transform(&self, position: &Point2<f64>) -> Point2<f64> {
        self.isometry * position
    }

    pub fn rotate(&self, vect: &Vector2<f64>) -> Vector2<f64> {
        self.isometry * vect
    }
}

impl Default for SymmetryTransform {
    fn default() -> Self {
        Self {
            isometry: IsometryMatrix2::identity(),
        }
    }
}

#[cfg(test)]
mod symmetry_transform_tests {
    use super::*;

    fn create_identity() -> SymmetryTransform {
        SymmetryTransform {
            isometry: IsometryMatrix2::identity(),
        }
    }

    #[test]
    fn default() {
        let point = Point2::new(0.2, 0.2);
        let transform = SymmetryTransform::default();
        assert_eq!(transform.transform(&point), point);
    }

    #[test]
    fn identity_transform() {
        let identity = create_identity();
        let point = Point2::new(0.2, 0.2);
        assert_eq!(identity.transform(&point), point);

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(identity.rotate(&vec), vec);
    }

    #[test]
    fn transform() {
        let isometry = SymmetryTransform {
            isometry: IsometryMatrix2::new(Vector2::new(1., 1.), PI / 2.),
        };

        let point = Point2::new(0.2, 0.2);
        assert_eq!(isometry.transform(&point), Point2::new(0.8, 1.2));

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(isometry.rotate(&vec), Vector2::new(-0.2, 0.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.1, 0.2));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.1, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.6, -0.2));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.4, -0.2));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.2, 0.));
    }
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
        self.symmetries.len()
    }

    fn degrees_of_freedom(&self) -> &[bool] {
        // TODO implement -> This is only required for the non-general wyckoff sites since all the
        // general sites have 3 degrees-of-freedom.
        //
        // This will be checked as a method of the SymmetryTransform struct.
        &[true, true, true]
    }
}

#[cfg(test)]
mod wyckoff_site_tests {
    use super::*;

    pub fn create_wyckoff() -> WyckoffSite {
        WyckoffSite {
            letter: 'a',
            symmetries: vec![SymmetryTransform::default()],
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

#[derive(Clone)]
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

/// Representing the unit cell of a crystal packing
///
/// The unit cell holds the unit cell parameters, being the length of each side of the cell in
/// addition to the contained angles. Each cell belongs to one of the Crystal Families which
/// dictate the degrees of freedom the cell can take.
///
#[derive(Clone)]
pub struct Cell {
    x_len: SharedValue,
    y_len: SharedValue,
    angle: SharedValue,
    family: CrystalFamily,
}

impl Default for Cell {
    fn default() -> Cell {
        Cell {
            x_len: SharedValue::new(1.),
            y_len: SharedValue::new(1.),
            angle: SharedValue::new(PI / 2.),
            family: CrystalFamily::Monoclinic,
        }
    }
}

impl Cell {
    /// Convert a transformation into cartesion coorsinates
    ///
    /// The positions of particles are stored in fractional coordinates, making changes to the
    /// unit cell simple. This function takes transformation to apply to a collection of points
    /// and converts the values of the fractional coordinates in the translation to real
    /// cartesian coordinates based on the current cell parameters.
    ///
    pub fn to_cartesian(&self, transform: IsometryMatrix2<f64>) -> IsometryMatrix2<f64> {
        let x = transform.translation.vector.x * self.x_len.get_value()
            + transform.translation.vector.y
                * self.y_len.get_value()
                * self.angle.get_value().cos();
        let y =
            transform.translation.vector.y * self.y_len.get_value() * self.angle.get_value().sin();

        IsometryMatrix2::from_parts(na::Translation2::new(x, y), transform.rotation)
    }

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
        Cell {
            x_len,
            y_len,
            angle,
            family: family.clone(),
        }
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

    pub fn center(&self) -> Point2<f64> {
        let (x, y) = self.to_cartesian(0.5, 0.5);
        Point2::new(x, y)
    }

    pub fn area(&self) -> f64 {
        self.angle.get_value().sin() * self.x_len.get_value() * self.y_len.get_value()
    }
}

#[cfg(test)]
mod cell_tests {
    use super::*;

    #[test]
    fn to_cartesian_test() {
        let cell = Cell::default();
        let trans = na::IsometryMatrix2::new(na::Vector2::new(0.5, 0.5), 0.);

        assert_eq!(cell.to_cartesian(trans), trans);

        cell.angle.set_value(PI / 4.);
        let expected = na::IsometryMatrix2::new(
            na::Vector2::new(0.5 + 0.5 * 1. / f64::sqrt(2.), 0.5 * 1. / f64::sqrt(2.)),
            0.,
        );
        assert_abs_diff_eq!(cell.to_cartesian(trans), expected);
    }
}

#[derive(Clone)]
pub struct PackedState {
    pub wallpaper: Wallpaper,
    pub shape: RadialShape,
    pub cell: Cell,
    occupied_sites: Vec<OccupiedSite>,
    basis: Vec<StandardBasis>,
}

impl PackedState {
    fn all_positions(&self) -> Vec<SymmetryTransform> {
        let mut transforms: Vec<SymmetryTransform> = vec![];
        for site in self.occupied_sites.iter() {
            for symmetry in site.wyckoff.symmetries.iter() {
                transforms.push(symmetry.clone());
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
        for (index1, position1) in self.all_positions().iter().enumerate() {
            let shape_i1 = shape::ShapeInstance {
                shape: &self.shape,
                isometry: self.cell.to_cartesian(position1.isometry),
            };
            // We only need to check the positions after that of index1, since the previous ones
            // have already been checked, hence `.skip(index1)`
            for (index2, position2) in self.all_positions().iter().enumerate().skip(index1) {
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
                        let translation = na::Translation2::new(*x_periodic, *y_periodic);
                        let shape_i2 = shape::ShapeInstance {
                            shape: &self.shape,
                            isometry: self.cell.to_cartesian(position2.isometry * translation),
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

impl MCVars {
    fn kt_ratio(&self) -> f64 {
        f64::powf(self.kt_finish / self.kt_start, 1.0 / self.steps as f64)
    }
}

fn mc_temperature(old: f64, new: f64, kt: f64, n: u64) -> f64 {
    f64::exp((1. / old - 1. / new) / kt) * (old / new).powi(n as i32)
}

pub fn monte_carlo_best_packing(vars: &MCVars, state: &mut PackedState) -> PackedState {
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
        "Packing Fraction: {}, Rejection Percentage {:.2}",
        packing_max,
        100. * rejections as f64 / vars.steps as f64,
    );
    best_state
}
