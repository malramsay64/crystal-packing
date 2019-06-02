//
// atom.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use std::fmt;

use nalgebra::Point3;

use crate::Intersect;

#[derive(Clone, PartialEq, Debug)]
pub struct Atom3 {
    pub position: Point3<f64>,
    pub radius: f64,
}

impl Intersect for Atom3 {
    fn intersects(&self, other: &Self) -> bool {
        let r_squared = (self.radius + other.radius).powi(2);
        // We have an intersection when the distance between the particles is less than the
        // combined radius of the two particles.
        na::distance_squared(&self.position, &other.position) < r_squared
    }
}

impl fmt::Display for Atom3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Atom3 {{ {}, {}, {}, {} }}",
            self.position.x, self.position.y, self.position.z, self.radius
        )
    }
}

impl Atom3 {
    pub fn new(x: f64, y: f64, z: f64, radius: f64) -> Self {
        Atom3 {
            position: Point3::<f64>::new(x, y, z),
            radius,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {}
}
