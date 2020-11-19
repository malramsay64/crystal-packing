//
// optimiser.rs
// Copyright (C) 2020 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use rand::distributions::Uniform;
use rand::prelude::*;
use rand_pcg::Pcg64Mcg;
use wasm_bindgen::prelude::*;

use crate::state::JSState;
use packing::traits::Basis;

#[wasm_bindgen]
extern "C" {
    // Use `js_namespace` here to bind `console.log(..)` instead of just
    // `log(..)`
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

#[wasm_bindgen]
pub struct Optimiser {
    pub kt: f64,
    pub step_size: f64,
    pub steps: u64,
    pub seed: u64,
}

// These are the function that we want to keep internal
impl Optimiser {
    #[inline]
    fn energy_surface(&self, new: f64, old: f64) -> f64 {
        f64::min(f64::exp((new - old) / self.kt), 1.)
    }

    #[inline]
    fn test_acceptance(&self, threshold: f64, new: f64, old: f64) -> bool {
        threshold < self.energy_surface(new, old)
    }

    fn accept_score<R: Rng + ?Sized>(
        &self,
        new: Option<f64>,
        old: f64,
        rng: &mut R,
    ) -> Option<f64> {
        let threshold: f64 = rng.gen();

        match new {
            // New score is better, keep updated state
            Some(new_score) if new_score > old => Some(new_score),
            // When the score increases, there is a probability of accepting the new
            // score, which is evaluated based on the difference between the new score,
            // and the old score along with the temperature value (kt).
            // The acceptance occurs when a random number smaller than the acceptance
            // probability is drawn.
            Some(new_score) if self.test_acceptance(threshold, new_score, old) => Some(new_score),
            // If the first two tests fail, then the score is rejected.
            _ => None,
        }
    }
}

#[wasm_bindgen]
impl Optimiser {
    pub fn new(kt: f64, step_size: f64, steps: u64, seed: u64) -> Optimiser {
        Optimiser {
            kt,
            step_size,
            steps,
            seed,
        }
    }

    pub fn optimise_state(&mut self, state: JSState) -> JSState {
        let mut score_current = match state.score() {
            Some(score) => score,
            _ => panic!("Invalid configuration passed to function, exiting."),
        };

        let mut rng = Pcg64Mcg::seed_from_u64(self.seed);

        let mut basis = state.generate_basis();
        let basis_distribution = Uniform::new(0, basis.len() as usize);
        let mut rejections: u64 = 0;

        for _ in 0..self.steps {
            // Choose a basis at random to modify
            // This is needed later if we need to undo the change
            let basis_index: usize = basis_distribution.sample(&mut rng);

            // Make a random modification to the selected basis
            basis
                .get_mut(basis_index)
                // There was some error in accessing the basis,
                // This should never occur in normal operation so panic and exit
                .expect("Trying to access basis which doesn't exist")
                .set_sampled(&mut rng, self.step_size);

            // Check if modification was good
            score_current = match self.accept_score(state.score(), score_current, &mut rng) {
                Some(score) => score,
                // Score was rejected so we have to undo the change
                None => {
                    basis
                        .get(basis_index)
                        // There was some error in accessing the basis,
                        // This should never occur in normal operation so panic and exit
                        .expect("Trying to access basis which doesn't exist.")
                        .reset_value();
                    // Increment counter of rejections
                    rejections += 1;
                    score_current
                }
            };
        }
        // Scale step size with goal of 75% rejections
        // Taking shinking the cell as an example, 50% of steps will  increase the cell, so
        // we want 50% of the steps which can improve the performance to be accepted.
        // There is a limit to the usefulness though and 1e-4 has been good.
        self.step_size *= self.steps as f64 / (2 * rejections + 1) as f64;
        assert!(
            state.score().is_some(),
            "Final score is invalid, this shouldn't occur in normal operation"
        );
        state
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::abs_diff_eq;
    use quickcheck_macros::quickcheck;

    static OPT: Optimiser = Optimiser {
        kt: 0.,
        step_size: 0.,
        steps: 0,
        seed: 0,
    };

    #[quickcheck]
    fn test_energy_surface(new: f64, old: f64) -> bool {
        let result = OPT.energy_surface(new, old, 0.);
        if new < old {
            abs_diff_eq!(result, 0.)
        } else if new >= old {
            abs_diff_eq!(result, 1.)
        } else {
            false
        }
    }

    #[quickcheck]
    fn test_energy_surface_temperature(new: f64, old: f64) -> bool {
        let result = OPT.energy_surface(new, old, 0.5);
        if new < old {
            0. < result && result < 1.
        } else if new >= old {
            abs_diff_eq!(result, 1.)
        } else {
            false
        }
    }
}
