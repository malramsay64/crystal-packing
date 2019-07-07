//
// optimisation.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use log::{info, trace};
use rand::distributions::Uniform;
use rand::prelude::*;
use rand_pcg::Pcg64Mcg;

use crate::traits::*;

#[derive(Debug, Clone, Copy)]
pub struct BuildOptimiser {
    kt_start: f64,
    kt_finish: f64,
    max_step_size: f64,
    steps: u64,
    seed: Option<u64>,
}

impl Default for BuildOptimiser {
    fn default() -> Self {
        Self {
            kt_start: 0.1,
            kt_finish: 0.0005,
            max_step_size: 0.01,
            steps: 100,
            seed: None,
        }
    }
}

impl BuildOptimiser {
    pub fn kt_start(&mut self, kt_start: f64) -> &mut Self {
        self.kt_start = kt_start;
        self
    }

    pub fn kt_finish(&mut self, kt_finish: f64) -> &mut Self {
        self.kt_finish = kt_finish;
        self
    }

    pub fn max_step_size(&mut self, max_step_size: f64) -> &mut Self {
        self.max_step_size = max_step_size;
        self
    }

    pub fn steps(&mut self, steps: u64) -> &mut Self {
        self.steps = steps;
        self
    }

    pub fn seed(&mut self, seed: u64) -> &mut Self {
        self.seed = Some(seed);
        self
    }

    pub fn build(&self) -> MCOptimiser {
        let kt_ratio = f64::powf(self.kt_finish / self.kt_start, 1.0 / self.steps as f64);
        let seed = match self.seed {
            None => Pcg64Mcg::from_entropy().gen(),
            Some(x) => x,
        };

        MCOptimiser {
            kt_start: self.kt_start,
            kt_ratio,
            max_step_size: self.max_step_size,
            steps: self.steps,
            seed,
        }
    }
}

pub struct MCOptimiser {
    kt_start: f64,
    kt_ratio: f64,
    max_step_size: f64,
    steps: u64,
    seed: u64,
}

impl MCOptimiser {
    fn temperature(&self, old: f64, new: f64, kt: f64, num_shapes: u64) -> f64 {
        f64::exp((new - old) / kt)
    }

    pub fn optimise_state(&self, mut state: impl State) -> Result<impl State, &'static str> {
        let mut rng = Pcg64Mcg::seed_from_u64(self.seed);
        let mut rejections: u64 = 0;

        let mut kt: f64 = self.kt_start;
        let total_shapes: u64 = state.total_shapes() as u64;

        let mut basis = state.generate_basis();
        let basis_distribution = Uniform::new(0, basis.len() as u64);

        let mut score_prev: f64 = state.score()?;
        let mut score: f64 = 0.;
        let mut score_max: f64 = 0.;

        let mut best_state = state.clone();

        for _ in 0..self.steps {
            let basis_index: usize = basis_distribution.sample(&mut rng) as usize;
            if let Some(basis_current) = basis.get_mut(basis_index) {
                basis_current.set_value(basis_current.sample(&mut rng, self.max_step_size));
            }

            score = match state.score() {
                Err(_) => {
                    trace!("Rejected for invalid score.");
                    rejections += 1;
                    score_prev
                }
                Ok(score_new) => {
                    if rng.gen::<f64>() > self.temperature(score_prev, score_new, kt, total_shapes)
                    {
                        trace!(
                            "Rejected for lowering score at t={:.4}\nscore_prev: {:.4}, score_new: {:.4}, kt: {:.4}",
                            self.temperature(score_prev, score_new, kt, total_shapes),
                            score_prev,
                            score_new,
                            kt
                        );
                        rejections += 1;
                        basis[basis_index].reset_value();

                        // Set packing to it's previous value
                        score_prev
                    } else {
                        // This is where we update the score cause the test was successful
                        score_prev = score;
                        score_new
                    }
                }
            };
            if score > score_max {
                best_state = state.clone();
                score_max = score;
            }
            kt *= self.kt_ratio;
        }
        info!(
            "Score: {:.4}, Rejections: {:.2} %",
            score_max,
            100. * rejections as f64 / self.steps as f64,
        );
        Ok(best_state)
    }
}
