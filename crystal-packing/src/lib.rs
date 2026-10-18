//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Finding dense packings of 2D shapes in periodic crystal structures.
//!
//! The crystal structures are described by one of the 17 plane groups ([`symmetry`]) with a set of
//! occupied Wyckoff positions, a unit cell ([`lattice`]) and the free parameters of each occupied
//! site. A [`engine::Packing`] holds such a structure for a [`hard_shape::HardShape`], checking for
//! overlaps using the symmetry of the structure. The structures are compressed using isobaric
//! Monte Carlo ([`monte_carlo`]), or for molecules interacting through a Lennard-Jones potential,
//! the lattice energy is minimised with simulated annealing ([`potential`]).
//!
//! A [`search::Search`] finds the best structure over every compatible plane group and choice of
//! Wyckoff sites, splitting the work into independent jobs which can run in parallel.
//!
//! ```
//! use crystal_packing::search::{Search, SearchConfig, ShapeSpec};
//!
//! let config = SearchConfig {
//!     groups: vec!["p2".into()],
//!     replicas: 2,
//!     steps: 5_000,
//!     ..Default::default()
//! };
//! let mut search = Search::new(ShapeSpec::Polygon { sides: 4 }, config).unwrap();
//! let best = search.run().unwrap();
//! assert!(best.packing_fraction > 0.99);
//! ```

#![forbid(unsafe_code)]

pub mod compare;
pub mod container;
pub mod engine;
pub mod geometry;
pub mod hard_shape;
pub mod lattice;
pub mod monte_carlo;
pub mod potential;
pub mod random;
pub mod render;
pub mod search;
pub mod symmetry;
pub mod validate;

pub use crate::random::Rng;
