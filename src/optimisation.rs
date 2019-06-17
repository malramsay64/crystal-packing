//
// optimisation.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use log::{info, trace};
use rand::distributions::Uniform;
use rand::prelude::*;
use rand::rngs::SmallRng;

use crate::traits::*;

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

pub fn monte_carlo_best_packing(
    vars: &MCVars,
    mut state: impl State,
) -> Result<impl State, &'static str> {
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

    let mut score_prev: f64 = state.score()?;
    let mut score: f64 = 0.;
    let mut score_max: f64 = 0.;

    let mut best_state = state.clone();

    for _ in 0..vars.steps {
        let basis_index: usize = basis_distribution.sample(&mut rng) as usize;
        if let Some(basis_current) = basis.get_mut(basis_index) {
            basis_current.set_value(basis_current.sample(&mut rng, vars.max_step_size));
        }

        score = match state.score() {
            Err(_) => {
                trace!("Trace rejected for increasing packing fraction");
                score_prev
            }
            Ok(score_new) => {
                if rng.gen::<f64>() > mc_temperature(score_prev, score_new, kt, total_shapes) {
                    // Packing fraction was increased too much so reject the step
                    trace!("Rejected for Increasing packing fraction.");
                    rejections += 1;
                    basis[basis_index].reset_value();

                    // Set packing to it's previous value
                    score_prev
                } else {
                    // This is where we update the score cause the test was successful
                    score_prev = score_new;
                    score_new
                }
            }
        };
        if score > score_max {
            best_state = state.clone();
            score_max = score;
        }
        kt *= kt_ratio;
    }
    info!(
        "Score: {:.4}, Rejections: {:.2} %",
        score,
        100. * rejections as f64 / vars.steps as f64,
    );
    Ok(best_state)
}
