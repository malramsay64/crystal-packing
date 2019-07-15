//
// optimisation.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use log::{debug, info, trace};
use rand::distributions::Uniform;
use rand::prelude::*;
use rand_pcg::Pcg64Mcg;

use crate::traits::*;

#[derive(Debug, Clone, Copy)]
pub struct BuildOptimiser {
    kt_start: f64,
    kt_finish: Option<f64>,
    kt_ratio: Option<f64>,
    max_step_size: f64,
    steps: u64,
    inner_steps: u64,
    seed: Option<u64>,
}

impl Default for BuildOptimiser {
    fn default() -> Self {
        Self {
            kt_start: 0.1,
            kt_finish: Some(0.001),
            kt_ratio: None,
            max_step_size: 0.01,
            steps: 1000,
            inner_steps: 1000,
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
        self.kt_finish = Some(kt_finish);
        self
    }

    pub fn kt_ratio(&mut self, kt_ratio: Option<f64>) -> &mut Self {
        self.kt_ratio = kt_ratio;
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

    pub fn inner_steps(&mut self, inner_steps: u64) -> &mut Self {
        self.inner_steps = inner_steps;
        self
    }

    pub fn seed(&mut self, seed: u64) -> &mut Self {
        self.seed = Some(seed);
        self
    }

    pub fn build(&self) -> MCOptimiser {
        let kt_ratio = match (self.kt_ratio, self.kt_finish) {
            (Some(ratio), _) => 1. - ratio,
            (None, Some(finish)) => f64::powf(finish / self.kt_start, 1. / self.steps as f64),
            (None, None) => 0.1,
        };
        debug!("Setting kt_ratio to: {}", kt_ratio);
        let seed = match self.seed {
            None => Pcg64Mcg::from_entropy().gen(),
            Some(x) => x,
        };

        MCOptimiser {
            kt_start: self.kt_start,
            kt_ratio,
            max_step_size: self.max_step_size,
            steps: self.steps,
            inner_steps: self.inner_steps,
            seed,
        }
    }
}

pub struct MCOptimiser {
    kt_start: f64,
    kt_ratio: f64,
    max_step_size: f64,
    steps: u64,
    inner_steps: u64,
    seed: u64,
}

impl MCOptimiser {
    fn accept_probability(&self, new: f64, old: f64, kt: f64) -> f64 {
        f64::min(f64::exp((new - old) / kt), 1.)
    }

    pub fn optimise_state(&self, mut state: impl State) -> Result<impl State, &'static str> {
        let mut rng = Pcg64Mcg::seed_from_u64(self.seed);
        let mut rejections: u64 = 0;

        let mut kt: f64 = self.kt_start;

        let mut basis = state.generate_basis();
        let basis_distribution = Uniform::new(0, basis.len() as u64);

        let mut score_prev: f64 = state.score()?;
        let mut score: f64 = 0.;
        let mut score_max: f64 = 0.;
        let mut step_ratio = 1.;

        let mut best_state = state.clone();

        let mut total_steps = 0;
        for loops in 1..=self.steps / self.inner_steps {
            let mut loop_rejections: u64 = 0;
            for _ in 0..self.inner_steps {
                let basis_index: usize = basis_distribution.sample(&mut rng) as usize;
                if let Some(basis_current) = basis.get_mut(basis_index) {
                    basis_current
                        .set_value(basis_current.sample(&mut rng, self.max_step_size * step_ratio));
                }

                score = match state.score() {
                    Ok(score_new) => {
                        // When the new score is better we accept the new state or
                        // When the score is higher, there is a probability of accepting the new
                        // score, which is evaluated based on the difference between the new score,
                        // and the old score along with the temperature value (kt).
                        // The acceptance occurs when a random number smaller than the acceptance
                        // probability is drawn.
                        if score_new > score_prev
                            || rng.gen::<f64>() < self.accept_probability(score_new, score_prev, kt)
                        {
                            trace!("Accepted new score: {}, previous: {}", score_new, score);
                            score_prev = score;
                            score_new
                        }
                        // If the first two tests fail, then the score is rejected.
                        else {
                            trace!(
                                "Rejected for lowering score at t={:.4}\nscore_prev: {:.4}, score_new: {:.4}, kt: {:.4}",
                                self.accept_probability(score_prev, score_new, kt),
                                score_prev,
                                score_new,
                                kt
                            );
                            loop_rejections += 1;
                            basis[basis_index].reset_value();

                            // Set packing to it's previous value
                            score_prev
                        }
                    }
                    Err(_) => {
                        trace!("Rejected for invalid score.");
                        loop_rejections += 1;
                        score_prev
                    }
                };
                if score > score_max {
                    best_state = state.clone();
                    score_max = score;
                }
            }
            rejections += loop_rejections;
            kt *= self.kt_ratio;
            total_steps = loops * self.inner_steps;

            // Scale step ratio with goal of 50% rejections
            if step_ratio > 1e-3 {
                step_ratio *= self.inner_steps as f64 / (2. * loop_rejections as f64 + 1.);
            }
        }
        info!(
            "Score: {:.4}, Rejected Fraction: {:.2}%",
            score_max,
            100. * rejections as f64 / total_steps as f64,
        );
        Ok(best_state)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use quickcheck_macros::quickcheck;
    use approx::abs_diff_eq;

    static OPT: MCOptimiser = MCOptimiser {
        kt_start: 0.,
        kt_ratio: 0.,
        max_step_size: 0.,
        steps: 0,
        inner_steps: 0,
        seed: 0,
    };

    #[quickcheck]
    fn test_accept_probability(new: f64, old: f64) -> bool {
        let result = OPT.accept_probability(new, old, 0.);
        if new < old {
            abs_diff_eq!(result, 0.)
        } else if new >= old {
            abs_diff_eq!(result, 1.)
        } else {
            false
        }
    }

    #[quickcheck]
    fn test_accept_probability_temperature(new: f64, old: f64) -> bool {
        let result = OPT.accept_probability(new, old, 0.5);
        if new < old {
            0. < result && result < 1.
        } else if new >= old {
            abs_diff_eq!(result, 1.)
        } else {
            false
        }
    }
}
