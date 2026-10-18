//
// random.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! The deterministic random number generator used by the optimisation.
//!
//! The generator is xoshiro256++ from `rand_xoshiro`, seeded with SplitMix64. It's a specified
//! algorithm, so a seed gives the same numbers on every platform (including WebAssembly) and with
//! every version of the crate. This wraps it with the few ways the optimisation draws numbers.

use rand_core::{Rng as _, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// One step of SplitMix64, used to mix seeds together in [`Rng::derive_seed`]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A xoshiro256++ generator
///
/// ```
/// use crystal_packing::Rng;
/// let mut a = Rng::seed_from_u64(42);
/// let mut b = Rng::seed_from_u64(42);
/// assert_eq!(a.next_u64(), b.next_u64());
/// let x = a.uniform();
/// assert!((0. ..1.).contains(&x));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rng(Xoshiro256PlusPlus);

impl Rng {
    /// A generator with its state expanded from a seed by SplitMix64
    pub fn seed_from_u64(seed: u64) -> Self {
        Rng(Xoshiro256PlusPlus::seed_from_u64(seed))
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
        self.0.next_u64()
    }

    /// A uniformly distributed value in the range [0, 1), from the top 53 bits
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
        // Lemire's multiply-shift method; the bias of at most n / 2^64 is negligible.
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

    #[test]
    fn matches_reference_values() {
        // The xoshiro256++ reference implementation from a state of [1, 2, 3, 4]
        let mut state = [0u8; 32];
        for (i, word) in [1u64, 2, 3, 4].iter().enumerate() {
            state[8 * i..8 * i + 8].copy_from_slice(&word.to_le_bytes());
        }
        let mut rng = Rng(Xoshiro256PlusPlus::from_seed(state));
        let expected = [
            41943041,
            58720359,
            3588806011781223,
            3591011842654386,
            9228616714210784205,
        ];
        for value in expected {
            assert_eq!(rng.next_u64(), value);
        }
        // Seeding with SplitMix64, as recorded before using rand_xoshiro, so a seed gives the
        // same structures it always has
        assert_eq!(Rng::seed_from_u64(0).next_u64(), 0x5317_5d61_490b_23df);
    }
}
