//
// potential.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Minimising the lattice energy of molecules interacting with a Lennard-Jones potential.
//!
//! The molecules are made of atoms, where the radius of each atom is half its Lennard-Jones σ, so
//! the σ of a pair of atoms is the sum of their radii (Lorentz mixing). The potential is truncated
//! and shifted at a cutoff given in units of σ. The energy per molecule (plus a pressure term for
//! the area of the cell) is minimised using simulated annealing with the same moves as the hard
//! shape compression.

use serde::{Deserialize, Serialize};

use crate::engine::{Packing, SiteState, Structure};
use crate::geometry::Vec2;
use crate::hard_shape::{HardShape, Molecule, PlacedMolecule};
use crate::lattice::Lattice;
use crate::monte_carlo::{
    keep_reduced, propose_lattice, propose_reorientation, propose_rotation, propose_translation,
    site_moves, Progress, StepSize,
};
use crate::random::Rng;

/// The Lennard-Jones potential and the parameters of the annealing schedule
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LennardJones {
    /// The depth of the potential well
    pub epsilon: f64,
    /// The cutoff of the potential in units of σ
    pub cutoff: f64,
    /// The temperature at the start of the annealing
    pub kt_start: f64,
    /// The temperature at the end of the annealing
    pub kt_end: f64,
    /// The pressure applied to the cell, in units of ε / σ² of an atom with unit radius
    pub pressure: f64,
}

impl Default for LennardJones {
    fn default() -> Self {
        LennardJones {
            epsilon: 1.,
            cutoff: 2.5,
            kt_start: 0.1,
            kt_end: 1e-3,
            pressure: 0.,
        }
    }
}

impl LennardJones {
    /// The energy of an atom pair separated by a squared distance, with a given σ
    #[inline]
    fn atom_energy(&self, distance_sq: f64, sigma: f64) -> f64 {
        let cutoff = self.cutoff * sigma;
        if distance_sq >= cutoff * cutoff {
            return 0.;
        }
        let s6 = libm::pow(sigma * sigma / distance_sq, 3.);
        let c6 = libm::pow(1. / self.cutoff, 6.);
        4. * self.epsilon * ((s6 * s6 - s6) - (c6 * c6 - c6))
    }

    /// The interaction energy of two placed molecules, translating `b` by `shift`
    pub fn pair_energy(&self, a: &PlacedMolecule, b: &PlacedMolecule, shift: Vec2) -> f64 {
        let mut energy = 0.;
        for &(pa, ra) in a.atoms() {
            for &(pb, rb) in b.atoms() {
                energy += self.atom_energy((pb + shift - pa).norm_sq(), ra + rb);
            }
        }
        energy
    }

    /// The largest distance between the centres of two interacting molecules
    pub fn range(&self, shape: &Molecule) -> f64 {
        let reach = shape
            .atoms()
            .iter()
            .map(|(p, _)| p.norm())
            .fold(0., f64::max);
        let sigma = 2. * shape.atoms().iter().map(|(_, r)| *r).fold(0., f64::max);
        2. * reach + self.cutoff * sigma
    }

    /// The energy per molecule of a packing
    pub fn energy(&self, packing: &Packing<Molecule>) -> f64 {
        packing.energy(self.range(packing.shape()), |a, b, t| {
            self.pair_energy(a, b, t)
        })
    }
}

/// Simulated annealing of the lattice energy of a molecular crystal
pub struct Annealer {
    packing: Packing<Molecule>,
    potential: LennardJones,
    steps: u64,
    rng: Rng,
    translate: StepSize,
    rotate: StepSize,
    cell: StepSize,
    step: u64,
    /// The enthalpy per molecule of the current state
    enthalpy: f64,
    best: (Vec<SiteState>, Lattice),
    best_enthalpy: f64,
}

const TARGET_ACCEPTANCE: f64 = 0.3;
const ADAPT_INTERVAL: u64 = 100;

impl Annealer {
    pub fn new(packing: Packing<Molecule>, potential: LennardJones, steps: u64, seed: u64) -> Self {
        let mut annealer = Annealer {
            packing,
            potential,
            steps,
            rng: Rng::seed_from_u64(seed),
            translate: StepSize::new(0.1, 1e-12, 1.),
            rotate: StepSize::new(0.1, 1e-12, std::f64::consts::PI),
            cell: StepSize::new(0.01, 1e-12, 0.2),
            step: 0,
            enthalpy: 0.,
            best: (
                vec![],
                Lattice::with_length(crate::symmetry::LatticeFamily::Square, 1.),
            ),
            best_enthalpy: f64::INFINITY,
        };
        annealer.enthalpy = annealer.compute_enthalpy();
        annealer.update_best();
        annealer
    }

    fn compute_enthalpy(&self) -> f64 {
        let area_per_molecule = self.packing.lattice().area() / self.packing.num_shapes() as f64;
        self.potential.energy(&self.packing) + self.potential.pressure * area_per_molecule
    }

    fn snapshot(&self) -> (Vec<SiteState>, Lattice) {
        (
            (0..self.packing.sites().len())
                .map(|i| self.packing.site_state(i))
                .collect(),
            self.packing.lattice().clone(),
        )
    }

    fn update_best(&mut self) {
        if self.enthalpy < self.best_enthalpy {
            self.best_enthalpy = self.enthalpy;
            self.best = self.snapshot();
        }
    }

    fn temperature(&self) -> f64 {
        let fraction = (self.step as f64 / self.steps.max(1) as f64).min(1.);
        self.potential.kt_start
            * libm::pow(self.potential.kt_end / self.potential.kt_start, fraction)
    }

    /// Metropolis acceptance of a change of the enthalpy per molecule
    fn accept(&mut self, new_enthalpy: f64, kt: f64) -> bool {
        let delta = (new_enthalpy - self.enthalpy) * self.packing.num_shapes() as f64;
        if !new_enthalpy.is_finite() {
            return false;
        }
        delta <= 0. || (kt > 0. && self.rng.uniform() < libm::exp(-delta / kt))
    }

    fn try_site(&mut self, index: usize, new: SiteState, kt: f64) -> bool {
        self.packing.set_site(index, new);
        let enthalpy = self.compute_enthalpy();
        if self.accept(enthalpy, kt) {
            self.enthalpy = enthalpy;
            true
        } else {
            self.packing.undo();
            false
        }
    }

    fn try_cell(&mut self, kt: f64) -> bool {
        let new = propose_lattice(self.packing.lattice(), self.cell.size, &mut self.rng);
        if !new.is_valid() {
            return false;
        }
        self.packing.set_lattice(new);
        let enthalpy = self.compute_enthalpy();
        if self.accept(enthalpy, kt) {
            self.enthalpy = enthalpy;
            keep_reduced(&mut self.packing);
            true
        } else {
            self.packing.undo();
            false
        }
    }

    fn step_once(&mut self, kt: f64) {
        let num_sites = self.packing.sites().len();
        let choice = self.rng.index(num_sites + 1);
        if choice == num_sites {
            let accepted = self.try_cell(kt);
            self.cell
                .record(accepted, ADAPT_INTERVAL, TARGET_ACCEPTANCE);
        } else {
            let (translate, rotate, reorient) = site_moves(&self.packing, choice);
            if reorient && self.rng.uniform() < 0.1 {
                if let Some(new) = propose_reorientation(&self.packing, choice, &mut self.rng) {
                    self.try_site(choice, new, kt);
                }
            } else if translate && (!rotate || self.rng.index(2) == 0) {
                let new =
                    propose_translation(&self.packing, choice, self.translate.size, &mut self.rng);
                let accepted = self.try_site(choice, new, kt);
                self.translate
                    .record(accepted, ADAPT_INTERVAL, TARGET_ACCEPTANCE);
            } else if rotate {
                let state = self.packing.site_state(choice);
                let new = propose_rotation(state, self.rotate.size, &mut self.rng);
                let accepted = self.try_site(choice, new, kt);
                self.rotate
                    .record(accepted, ADAPT_INTERVAL, TARGET_ACCEPTANCE);
            }
        }
        self.update_best();
    }

    /// Run up to `steps` steps of the annealing, returning the progress
    pub fn run(&mut self, steps: u64) -> Progress {
        let end = (self.step + steps).min(self.steps);
        while self.step < end {
            let kt = self.temperature();
            self.step_once(kt);
            self.step += 1;
        }
        self.progress()
    }

    pub fn progress(&self) -> Progress {
        Progress {
            steps: self.step,
            total_steps: self.steps,
            pressure: self.potential.pressure,
            packing_fraction: self.packing.packing_fraction(),
            best_packing_fraction: self.best_packing_fraction(),
            best_energy: Some(self.best_enthalpy),
        }
    }

    pub fn is_complete(&self) -> bool {
        self.step >= self.steps
    }

    fn best_packing_fraction(&self) -> f64 {
        self.packing.num_shapes() as f64 * self.packing.shape().area() / self.best.1.area()
    }

    /// The current structure, with its energy
    pub fn current_structure(&self) -> Structure {
        let mut structure = self.packing.structure();
        structure.energy = Some(self.potential.energy(&self.packing));
        structure
    }

    pub fn packing(&self) -> &Packing<Molecule> {
        &self.packing
    }

    /// Quench the best structure found, returning it with its energy
    pub fn finish(&mut self) -> Structure {
        let best = self.best.clone();
        self.packing.set_lattice(best.1.clone());
        for (i, state) in best.0.iter().enumerate() {
            self.packing.set_site(i, *state);
        }
        self.enthalpy = self.compute_enthalpy();
        // A short zero temperature minimisation with small steps
        self.translate.size = self.translate.size.min(1e-3);
        self.rotate.size = self.rotate.size.min(1e-3);
        self.cell.size = self.cell.size.min(1e-4);
        for _ in 0..(self.steps / 10).clamp(500, 5000) {
            self.step_once(0.);
        }
        let best = self.best.clone();
        self.packing.set_lattice(best.1);
        for (i, state) in best.0.iter().enumerate() {
            self.packing.set_site(i, *state);
        }
        let mut structure = self.packing.structure();
        structure.energy = Some(self.potential.energy(&self.packing));
        structure
    }

    pub fn shape(&self) -> &Molecule {
        self.packing.shape()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::symmetry::{LatticeFamily, PlaneGroup};
    use approx::assert_abs_diff_eq;

    #[test]
    fn pair_minimum() {
        let lj = LennardJones::default();
        // For unit radius atoms σ = 2, the minimum is at 2^(1/6) σ with depth ε (before shifting)
        let r_min = 2f64.powf(1. / 6.) * 2.;
        let c6 = (1. / 2.5f64).powi(6);
        let shift = 4. * (c6 * c6 - c6);
        assert_abs_diff_eq!(
            lj.atom_energy(r_min * r_min, 2.),
            -1. - shift,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(lj.atom_energy(25. - 1e-9, 2.), 0., epsilon = 1e-9);
        assert_eq!(lj.atom_energy(25.1, 2.), 0.);
    }

    #[test]
    fn hexagonal_lattice_energy() {
        // Disks on a hexagonal lattice with spacing r have six neighbours at r, six at √3 r, ...
        let lj = LennardJones {
            cutoff: 1.8,
            ..Default::default()
        };
        let group = PlaneGroup::from_name("p6mm").unwrap();
        let mut rng = Rng::seed_from_u64(0);
        let mut packing = Packing::new(Molecule::disk(), group, &['a'], &mut rng).unwrap();
        let spacing = 2.3;
        packing.set_lattice(Lattice::from_parameters(
            LatticeFamily::Hexagonal,
            &[spacing],
        ));
        // With a cutoff of 1.8σ = 3.6 only the nearest (2.3) and second nearest (3.98 > 3.6)
        // neighbours are within range, so just the six nearest contribute.
        let expected = 6. * lj.atom_energy(spacing * spacing, 2.) / 2.;
        assert_abs_diff_eq!(lj.energy(&packing), expected, epsilon = 1e-12);
    }

    #[test]
    fn annealing_finds_hexagonal_lattice() {
        // The lowest energy structure of LJ disks is the hexagonal lattice near the minimum
        let lj = LennardJones::default();
        let group = PlaneGroup::from_name("p2").unwrap();
        let mut rng = Rng::seed_from_u64(1);
        let packing = Packing::new(Molecule::disk(), group, &['e'], &mut rng).unwrap();
        let mut annealer = Annealer::new(packing, lj, 20_000, 1);
        annealer.run(20_000);
        let structure = annealer.finish();
        let energy = structure.energy.unwrap();
        // The energy of the hexagonal lattice at its optimal spacing is about -3.3 ε per molecule
        // for a cutoff of 2.5σ (six nearest neighbours contribute most of it).
        assert!(energy < -3.0, "{}", energy);
        // A dense structure close to hexagonal packing of the σ = 2 disks
        let phi = structure.packing_fraction;
        assert!(phi > 0.7 && phi < 0.85, "{}", phi);
    }

    #[test]
    fn deterministic() {
        let run = || {
            let group = PlaneGroup::from_name("p2").unwrap();
            let mut rng = Rng::seed_from_u64(5);
            let trimer = Molecule::trimer(0.637556, 120., 1.).unwrap();
            let packing = Packing::new(trimer, group, &['e'], &mut rng).unwrap();
            let mut annealer = Annealer::new(packing, LennardJones::default(), 2000, 5);
            annealer.run(2000);
            annealer.finish()
        };
        assert_eq!(run(), run());
    }
}
