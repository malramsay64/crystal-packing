//
// random.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! A small, dependency free, deterministic random number generator.
//!
//! The optimisation needs to give bit-identical results on every platform for a given seed,
//! including when compiled to WebAssembly. Implementing the generator here, rather than relying on
//! an external crate, ensures the stream of numbers never changes underneath us.

use serde::{Deserialize, Serialize};

/// The SplitMix64 generator, used to expand a single `u64` seed into the full generator state.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The xoshiro256++ generator by Blackman and Vigna.
///
/// ```
/// use crystal_packing::Rng;
/// let mut a = Rng::seed_from_u64(42);
/// let mut b = Rng::seed_from_u64(42);
/// assert_eq!(a.next_u64(), b.next_u64());
/// let x = a.uniform();
/// assert!((0. ..1.).contains(&x));
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn seed_from_u64(seed: u64) -> Self {
        let mut state = seed;
        let s = [
            splitmix64(&mut state),
            splitmix64(&mut state),
            splitmix64(&mut state),
            splitmix64(&mut state),
        ];
        Rng { s }
    }

    /// Combine a seed with a sequence of values to create an independent stream.
    ///
    /// This is used to give every unit of work (group, sites, replica) its own stream of random
    /// numbers which doesn't depend on the order in which the work is scheduled.
    pub fn derive_seed(seed: u64, values: &[u64]) -> u64 {
        let mut state = seed;
        let mut out = splitmix64(&mut state);
        for v in values {
            state ^= v.wrapping_mul(0xD6E8_FEB8_6659_FD93);
            out ^= splitmix64(&mut state);
        }
        out
    }

    pub fn next_u64(&mut self) -> u64 {
        let result = (self.s[0].wrapping_add(self.s[3]))
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// A uniformly distributed value in the range [0, 1)
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// A uniformly distributed value in the range [-0.5, 0.5)
    pub fn centered(&mut self) -> f64 {
        self.uniform() - 0.5
    }

    /// A uniformly distributed index in the range [0, n)
    pub fn index(&mut self, n: usize) -> usize {
        assert!(n > 0, "Cannot choose from an empty range");
        // Lemire's multiply-shift method; the bias is negligible for the small n used here.
        ((self.next_u64() as u128 * n as u128) >> 64) as usize
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn deterministic() {
        let mut a = Rng::seed_from_u64(1);
        let mut b = Rng::seed_from_u64(1);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds() {
        let mut a = Rng::seed_from_u64(1);
        let mut b = Rng::seed_from_u64(2);
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn uniform_moments() {
        let mut rng = Rng::seed_from_u64(7);
        let n = 100_000;
        let mean: f64 = (0..n).map(|_| rng.uniform()).sum::<f64>() / n as f64;
        assert!((mean - 0.5).abs() < 0.01);
    }

    #[test]
    fn index_in_range() {
        let mut rng = Rng::seed_from_u64(3);
        let mut counts = [0usize; 5];
        for _ in 0..10_000 {
            counts[rng.index(5)] += 1;
        }
        assert!(counts.iter().all(|&c| c > 1800 && c < 2200));
    }

    #[test]
    fn derived_seeds_differ() {
        let a = Rng::derive_seed(0, &[1, 2]);
        let b = Rng::derive_seed(0, &[2, 1]);
        let c = Rng::derive_seed(0, &[1, 2]);
        assert_ne!(a, b);
        assert_eq!(a, c);
    }
}
