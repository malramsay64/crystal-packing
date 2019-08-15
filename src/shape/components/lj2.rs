//
// lj2.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fmt;

use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

use crate::traits::Potential;

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct LJ2 {
    pub position: Vector2<f64>,
    pub sigma: f64,
    pub epsilon: f64,
}

impl Potential for LJ2 {
    fn energy(&self, other: &Self) -> f64 {
        let sigma_squared = self.sigma.powi(2);
        let r_squared = (&self.position - &other.position).norm_squared();
        let sigma2_r2_cubed = (sigma_squared / r_squared).powi(3);

        4. * self.epsilon * (sigma2_r2_cubed.powi(2) - sigma2_r2_cubed)
    }
}

impl fmt::Display for LJ2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "LJ2 {{ {}, {}, {}, {}}}",
            self.position.x, self.position.y, self.sigma, self.epsilon
        )
    }
}

impl LJ2 {
    pub fn new(x: f64, y: f64, sigma: f64) -> Self {
        LJ2 {
            position: Vector2::new(x, y),
            sigma,
            epsilon: 1.,
        }
    }
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn init_test() {
        let a = LJ2::new(0., 0., 1.);
        assert_abs_diff_eq!(a.position.x, 0.);
        assert_abs_diff_eq!(a.position.y, 0.);
        assert_abs_diff_eq!(a.sigma, 1.);
    }

    #[test]
    fn distance_squared_test() {
        let a0 = LJ2::new(0., 0., 1.);
        let a1 = LJ2::new(1., 0., 1.);
        assert_abs_diff_eq!((a0.position - a1.position).norm_squared(), 1.);
        assert_abs_diff_eq!(a0.energy(&a1), 0.);
    }

}
