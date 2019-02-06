//
// atom.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fmt;

use nalgebra as na;
use nalgebra::allocator::Allocator;
use nalgebra::{DefaultAllocator, DimName, Point, U2, U3};

use crate::shape::Intersect;

pub type Atom2 = Atom<U2>;
pub type Atom3 = Atom<U3>;

#[derive(Clone, PartialEq, Debug)]
pub struct Atom<D: DimName>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    pub position: Point<f64, D>,
    pub radius: f64,
}

impl<D: DimName> Intersect<D> for Atom<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    fn intersects(&self, other: &Self) -> bool {
        let r_squared = (self.radius + other.radius).powi(2);
        // We have an intersection when the distance between the particles is less than the
        // combined radius of the two particles.
        na::distance_squared(&self.position, &other.position) < r_squared
    }
}

impl fmt::Display for Atom<U2> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Atom2 {{ {}, {}, {} }}",
            self.position.x, self.position.y, self.radius
        )
    }
}

impl fmt::Display for Atom<U3> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Atom3 {{ {}, {}, {}, {} }}",
            self.position.x, self.position.y, self.position.z, self.radius
        )
    }
}

impl Atom<U2> {
    pub fn new(x: f64, y: f64, radius: f64) -> Self {
        Atom {
            position: Point::<f64, U2>::new(x, y),
            radius,
        }
    }
}

impl Atom<U3> {
    pub fn new(x: f64, y: f64, z: f64, radius: f64) -> Self {
        Atom {
            position: Point::<f64, U3>::new(x, y, z),
            radius,
        }
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn init_test() {
        let a = Atom2::new(0., 0., 1.);
        assert_abs_diff_eq!(a.position.x, 0.);
        assert_abs_diff_eq!(a.position.y, 0.);
        assert_abs_diff_eq!(a.radius, 1.);
    }

    #[test]
    fn distance_squared_test() {
        let a0 = Atom2::new(0., 0., 1.);
        let a1 = Atom2::new(0.5, 0., 1.);
        assert_abs_diff_eq!(na::distance_squared(&a0.position, &a1.position), 0.25);
        assert!(a0.intersects(&a1));
    }

    #[test]
    fn intersection_test() {
        let a0 = Atom2::new(0., 0., 1.);
        let a1 = Atom2::new(1., 0., 1.);
        let a2 = Atom2::new(0.5, 0.5, 1.);
        let a3 = Atom2::new(1.5, 1.5, 1.);

        assert!(a0.intersects(&a1));
        assert!(a1.intersects(&a2));
        assert!(a3.intersects(&a2));

        assert!(!a0.intersects(&a3));
    }

    #[test]
    fn intersection_calculation_test() {
        let a0 = Atom2::new(0., 0., f64::sqrt(2.) / 2.);
        let a1 = Atom2::new(1., 1., f64::sqrt(2.) / 2.);
        let a2 = Atom2::new(1., 1., f64::sqrt(2.) / 2. - 2. * std::f64::EPSILON);
        println!("Radii: {}", a0.radius * a0.radius + a1.radius * a1.radius);
        assert!(a0.intersects(&a1));
        assert!(a1.intersects(&a2));

        assert!(!a0.intersects(&a2));
    }

}
