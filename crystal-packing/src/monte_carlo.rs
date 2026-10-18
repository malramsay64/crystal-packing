//
// monte_carlo.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Isobaric ("floppy box") Monte Carlo compression of a periodic packing.
//!
//! The optimisation samples the isothermal-isobaric ensemble of hard shapes, slowly increasing the
//! pressure so the structure compresses towards its densest packing. Moves either change the
//! parameters of a single site (translation, rotation, or a discrete re-orientation) or change the
//! unit cell. Moves creating an overlap are always rejected, and cell moves are accepted with the
//! probability
//!
//! ```text
//! min(1, exp(-P* ΔA / A_shape + N ln(A' / A)))
//! ```
//!
//! where P* is the reduced pressure and N the number of shapes in the cell. The step size of each
//! type of move adapts to keep the acceptance rate near a target, so moves stay effective as the
//! structure becomes dense. The best structure seen is tracked, and at the end a deterministic
//! squeeze compresses it until the shapes are in contact.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::engine::{Packing, SiteState, Structure};
use crate::geometry::Vec2;
use crate::hard_shape::HardShape;
use crate::lattice::Lattice;
use crate::random::Rng;
use crate::symmetry::{LatticeFamily, OrientationFreedom};

/// Parameters controlling a Monte Carlo compression
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MonteCarloConfig {
    /// The number of attempted moves
    pub steps: u64,
    /// The reduced pressure at the start of the compression
    pub pressure_start: f64,
    /// The reduced pressure at the end of the compression
    pub pressure_end: f64,
    /// The acceptance rate the step sizes are adjusted towards
    pub target_acceptance: f64,
    /// The number of attempts of a move type between step size adjustments
    pub adapt_interval: u64,
    /// The number of moves in the final squeeze to contact
    pub squeeze_steps: u64,
}

/// The pressure through a compression, increasing geometrically from the start to the end
#[derive(Debug, Clone, Copy)]
pub(crate) struct PressureSchedule {
    start: f64,
    /// The logarithm of the ratio of the final to the initial pressure
    log_ratio: f64,
    steps: u64,
}

impl PressureSchedule {
    pub(crate) fn new(config: &MonteCarloConfig) -> PressureSchedule {
        PressureSchedule {
            start: config.pressure_start,
            log_ratio: libm::log(config.pressure_end / config.pressure_start),
            steps: config.steps.max(1),
        }
    }

    /// The pressure after a number of steps, staying at the final pressure after the end
    pub(crate) fn at(&self, step: u64) -> f64 {
        let fraction = (step as f64 / self.steps as f64).min(1.);
        self.start * libm::exp(self.log_ratio * fraction)
    }
}

impl Default for MonteCarloConfig {
    fn default() -> Self {
        MonteCarloConfig {
            steps: 20_000,
            pressure_start: 0.1,
            pressure_end: 1e6,
            target_acceptance: 0.3,
            adapt_interval: 100,
            squeeze_steps: 2_000,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MoveKind {
    Translate,
    Rotate,
    Cell,
}

/// The adaptive step size for one type of move
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StepSize {
    pub size: f64,
    min: f64,
    max: f64,
    attempts: u64,
    accepted: u64,
    pub total_attempts: u64,
    pub total_accepted: u64,
}

impl StepSize {
    pub(crate) fn new(size: f64, min: f64, max: f64) -> Self {
        StepSize {
            size,
            min,
            max,
            attempts: 0,
            accepted: 0,
            total_attempts: 0,
            total_accepted: 0,
        }
    }

    pub(crate) fn record(&mut self, accepted: bool, interval: u64, target: f64) {
        self.attempts += 1;
        self.total_attempts += 1;
        if accepted {
            self.accepted += 1;
            self.total_accepted += 1;
        }
        if self.attempts >= interval {
            let rate = self.accepted as f64 / self.attempts as f64;
            // Scale the step size towards the target acceptance, limiting the change per update
            let factor = ((rate + 0.01) / target).clamp(0.5, 1.5);
            self.size = (self.size * factor).clamp(self.min, self.max);
            self.attempts = 0;
            self.accepted = 0;
        }
    }

    pub fn acceptance(&self) -> f64 {
        if self.total_attempts == 0 {
            0.
        } else {
            self.total_accepted as f64 / self.total_attempts as f64
        }
    }
}

/// The current progress of an optimisation
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Progress {
    pub steps: u64,
    pub total_steps: u64,
    pub pressure: f64,
    pub packing_fraction: f64,
    pub best_packing_fraction: f64,
    /// The lowest enthalpy per molecule, when minimising an energy
    pub best_energy: Option<f64>,
}

/// A resumable Monte Carlo compression of a packing
pub struct Optimiser<S: HardShape> {
    packing: Packing<S>,
    config: MonteCarloConfig,
    schedule: PressureSchedule,
    rng: Rng,
    translate: StepSize,
    rotate: StepSize,
    cell: StepSize,
    step: u64,
    best: (Vec<SiteState>, Lattice),
    best_packing_fraction: f64,
}

/// Displace a site by a random amount in Cartesian space, in units of the shape's circumradius
pub(crate) fn propose_translation<S: HardShape>(
    packing: &Packing<S>,
    index: usize,
    step: f64,
    rng: &mut Rng,
) -> SiteState {
    let mut new = packing.site_state(index);
    let free = packing.sites()[index].free;
    let scale = step * packing.shape().circumradius();
    let displacement = Vec2::new(rng.centered(), rng.centered()) * (2. * scale);
    let fractional = packing.lattice().matrix.inverse() * displacement;
    if free[0] {
        new.position.x += fractional.x;
    }
    if free[1] {
        new.position.y += fractional.y;
    }
    new.position.x -= libm::floor(new.position.x);
    new.position.y -= libm::floor(new.position.y);
    new
}

/// Rotate a site by a random angle
pub(crate) fn propose_rotation(state: SiteState, step: f64, rng: &mut Rng) -> SiteState {
    let mut new = state;
    new.angle = crate::geometry::wrap_angle(new.angle + 2. * step * rng.centered(), 2. * PI);
    new
}

/// Choose one of the allowed orientations of a site constrained by mirror symmetry
pub(crate) fn propose_reorientation<S: HardShape>(
    packing: &Packing<S>,
    index: usize,
    rng: &mut Rng,
) -> Option<SiteState> {
    match &packing.sites()[index].orientation {
        OrientationFreedom::Discrete(angles) => {
            let mut new = packing.site_state(index);
            new.angle = angles[rng.index(angles.len())];
            Some(new)
        }
        OrientationFreedom::Free => None,
    }
}

/// Propose a random change to the lattice, with the step size being the logarithm of the
/// relative change of a cell parameter.
pub(crate) fn propose_lattice(lattice: &Lattice, step: f64, rng: &mut Rng) -> Lattice {
    let mut p = lattice.parameters();
    let delta = 2. * step * rng.centered();
    match lattice.family {
        LatticeFamily::Oblique => match rng.index(4) {
            0 => p[0] *= libm::exp(delta),
            1 => p[2] *= libm::exp(delta),
            // A shear of the cell, which doesn't change the area
            2 => p[1] += delta * p[0],
            _ => p.iter_mut().for_each(|v| *v *= libm::exp(delta)),
        },
        LatticeFamily::Rectangular => match rng.index(3) {
            0 => p[0] *= libm::exp(delta),
            1 => p[1] *= libm::exp(delta),
            _ => p.iter_mut().for_each(|v| *v *= libm::exp(delta)),
        },
        LatticeFamily::Square | LatticeFamily::Hexagonal => p[0] *= libm::exp(delta),
    }
    Lattice::from_parameters(lattice.family, &p)
}

/// Keep an oblique cell close to its reduced form, so the cell moves remain effective
pub(crate) fn keep_reduced<S: HardShape>(packing: &mut Packing<S>) {
    let lattice = packing.lattice();
    let p = lattice.parameters();
    if lattice.family == LatticeFamily::Oblique
        && (p[1].abs() > p[0] || p[0] > 1.5 * lattice.b().norm())
    {
        packing.reduce_lattice();
    }
}

/// Whether a site has any free parameters, translational or rotational
pub(crate) fn site_moves(packing: &Packing<impl HardShape>, index: usize) -> (bool, bool, bool) {
    let site = &packing.sites()[index];
    let translate = site.free[0] || site.free[1];
    let (rotate, reorient) = match &site.orientation {
        OrientationFreedom::Free => (true, false),
        OrientationFreedom::Discrete(angles) => (false, angles.len() > 1),
    };
    (translate, rotate, reorient)
}

impl<S: HardShape> Optimiser<S> {
    pub fn new(packing: Packing<S>, config: MonteCarloConfig, seed: u64) -> Self {
        let best = Optimiser::snapshot(&packing);
        let best_packing_fraction = packing.packing_fraction();
        Optimiser {
            packing,
            schedule: PressureSchedule::new(&config),
            config,
            rng: Rng::seed_from_u64(seed),
            // Translations are in units of the shape's circumradius
            translate: StepSize::new(0.2, 1e-12, 1.),
            rotate: StepSize::new(0.2, 1e-12, PI),
            // Cell moves are the logarithm of the relative change
            cell: StepSize::new(0.02, 1e-12, 0.5),
            step: 0,
            best,
            best_packing_fraction,
        }
    }

    fn snapshot(packing: &Packing<S>) -> (Vec<SiteState>, Lattice) {
        (
            (0..packing.sites().len())
                .map(|i| packing.site_state(i))
                .collect(),
            packing.lattice().clone(),
        )
    }

    fn restore(&mut self, snapshot: &(Vec<SiteState>, Lattice)) {
        self.packing.set_lattice(snapshot.1.clone());
        for (i, state) in snapshot.0.iter().enumerate() {
            self.packing.set_site(i, *state);
        }
    }

    pub fn packing(&self) -> &Packing<S> {
        &self.packing
    }

    pub fn step_sizes(&self) -> (StepSize, StepSize, StepSize) {
        (self.translate, self.rotate, self.cell)
    }

    /// The reduced pressure at the current step, increasing geometrically through the run
    fn pressure(&self) -> f64 {
        self.schedule.at(self.step)
    }

    fn update_best(&mut self) {
        let phi = self.packing.packing_fraction();
        if phi > self.best_packing_fraction {
            self.best_packing_fraction = phi;
            self.best = Optimiser::snapshot(&self.packing);
        }
    }

    /// Apply a change to a site, keeping it when there are no overlaps
    fn try_site(&mut self, index: usize, new: SiteState) -> bool {
        self.packing.set_site(index, new);
        if self.packing.site_overlaps(index) {
            self.packing.undo();
            false
        } else {
            true
        }
    }

    fn translate_site(&mut self, index: usize) -> bool {
        let new = propose_translation(&self.packing, index, self.translate.size, &mut self.rng);
        self.try_site(index, new)
    }

    fn rotate_site(&mut self, index: usize) -> bool {
        let new = propose_rotation(
            self.packing.site_state(index),
            self.rotate.size,
            &mut self.rng,
        );
        self.try_site(index, new)
    }

    fn reorient_site(&mut self, index: usize) -> bool {
        match propose_reorientation(&self.packing, index, &mut self.rng) {
            Some(new) => self.try_site(index, new),
            None => false,
        }
    }

    fn propose_lattice(&mut self) -> Lattice {
        propose_lattice(self.packing.lattice(), self.cell.size, &mut self.rng)
    }

    fn cell_move(&mut self, pressure: Option<f64>) -> bool {
        let new = self.propose_lattice();
        if !new.is_valid() {
            return false;
        }
        let pressure = pressure.unwrap_or_else(|| self.pressure());
        let (area_old, area_new) = (self.packing.lattice().area(), new.area());
        let n = self.packing.num_shapes() as f64;
        let shape_area = self.packing.shape().area();
        let log_acceptance = if pressure.is_infinite() {
            if area_new <= area_old {
                0.
            } else {
                f64::NEG_INFINITY
            }
        } else {
            -pressure * (area_new - area_old) / shape_area + n * libm::log(area_new / area_old)
        };
        // Check the acceptance before overlaps since it is much cheaper
        if log_acceptance < 0. && self.rng.uniform() >= libm::exp(log_acceptance) {
            return false;
        }
        self.packing.set_lattice(new);
        if self.packing.any_overlap() {
            self.packing.undo();
            return false;
        }
        keep_reduced(&mut self.packing);
        true
    }

    /// Perform a single Monte Carlo move at a fixed pressure, or None for the pressure of the
    /// schedule, which is only found for the moves which need it
    fn step_once(&mut self, pressure: Option<f64>) {
        let num_sites = self.packing.sites().len();
        let interval = self.config.adapt_interval;
        let target = self.config.target_acceptance;
        let choice = self.rng.index(num_sites + 1);
        if choice == num_sites {
            let accepted = self.cell_move(pressure);
            self.cell.record(accepted, interval, target);
            if accepted {
                self.update_best();
            }
            return;
        }
        let (translate, rotate, reorient) = site_moves(&self.packing, choice);
        let options = [translate, rotate];
        let available = options.iter().filter(|&&o| o).count();
        if reorient && self.rng.uniform() < 0.1 {
            self.reorient_site(choice);
            return;
        }
        if available == 0 {
            return;
        }
        let pick = if available == 2 {
            self.rng.index(2)
        } else if translate {
            0
        } else {
            1
        };
        if pick == 0 {
            let accepted = self.translate_site(choice);
            self.translate.record(accepted, interval, target);
        } else {
            let accepted = self.rotate_site(choice);
            self.rotate.record(accepted, interval, target);
        }
    }

    /// Run up to `steps` moves of the compression, returning the progress.
    ///
    /// This can be called repeatedly to run the optimisation in slices.
    pub fn run(&mut self, steps: u64) -> Progress {
        let end = (self.step + steps).min(self.config.steps);
        while self.step < end {
            self.step_once(None);
            self.step += 1;
        }
        self.progress()
    }

    pub fn progress(&self) -> Progress {
        Progress {
            steps: self.step,
            total_steps: self.config.steps,
            pressure: self.pressure(),
            packing_fraction: self.packing.packing_fraction(),
            best_packing_fraction: self.best_packing_fraction,
            best_energy: None,
        }
    }

    pub fn is_complete(&self) -> bool {
        self.step >= self.config.steps
    }

    /// Uniformly scale the cell down until the shapes come into contact
    fn squeeze_cell(&mut self) {
        let lattice = self.packing.lattice().clone();
        let valid = |optimiser: &mut Self, scale: f64| {
            optimiser.packing.set_lattice(lattice.scaled(scale));
            !optimiser.packing.any_overlap()
        };
        let mut high = 1.;
        let mut low = 0.99;
        while valid(self, low) && low > 1e-3 {
            high = low;
            low *= 0.99;
        }
        for _ in 0..50 {
            let mid = 0.5 * (low + high);
            if valid(self, mid) {
                high = mid;
            } else {
                low = mid;
            }
        }
        self.packing.set_lattice(lattice.scaled(high));
        self.update_best();
    }

    /// Finish the optimisation, compressing the best structure until it is in contact.
    ///
    /// This alternates uniform compression of the cell with moves at infinite pressure, where the
    /// cell can only shrink. The resulting structure is returned.
    pub fn finish(&mut self) -> Structure {
        let best = self.best.clone();
        self.restore(&best);
        self.squeeze_cell();
        let rounds = 4;
        for _ in 0..rounds {
            // Small steps are needed to make progress once the shapes are in contact
            self.translate.size = self.translate.size.min(1e-3);
            self.rotate.size = self.rotate.size.min(1e-3);
            self.cell.size = self.cell.size.min(1e-4);
            for _ in 0..self.config.squeeze_steps / rounds {
                self.step_once(Some(f64::INFINITY));
            }
            self.squeeze_cell();
        }
        let best = self.best.clone();
        self.restore(&best);
        debug_assert!(!self.packing.any_overlap());
        self.packing.structure()
    }

    pub fn best_packing_fraction(&self) -> f64 {
        self.best_packing_fraction
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::hard_shape::{Molecule, Polygon};
    use crate::symmetry::PlaneGroup;

    fn optimise<S: HardShape>(
        shape: S,
        group: &str,
        sites: &[char],
        steps: u64,
        seed: u64,
    ) -> (Structure, Optimiser<S>) {
        let group = PlaneGroup::from_name(group).unwrap();
        let mut rng = Rng::seed_from_u64(seed);
        let packing = Packing::new(shape, group, sites, &mut rng).unwrap();
        let config = MonteCarloConfig {
            steps,
            ..Default::default()
        };
        let mut optimiser = Optimiser::new(packing, config, seed);
        optimiser.run(steps);
        let structure = optimiser.finish();
        (structure, optimiser)
    }

    #[test]
    fn step_size_converges_to_target() {
        let mut step = StepSize::new(1., 1e-9, 10.);
        // A move where the acceptance probability falls as the step size grows
        let mut rng = Rng::seed_from_u64(1);
        for _ in 0..20_000 {
            let accepted = rng.uniform() < (0.1 / step.size).min(1.);
            step.record(accepted, 100, 0.3);
        }
        assert!((step.size - 0.33).abs() < 0.15, "{}", step.size);
    }

    #[test]
    fn deterministic() {
        let (a, _) = optimise(Polygon::regular(5).unwrap(), "p2", &['e'], 3000, 7);
        let (b, _) = optimise(Polygon::regular(5).unwrap(), "p2", &['e'], 3000, 7);
        assert_eq!(a, b);
    }

    #[test]
    fn best_is_never_worse_than_final() {
        let group = PlaneGroup::from_name("p2").unwrap();
        let mut rng = Rng::seed_from_u64(3);
        let packing = Packing::new(Polygon::regular(4).unwrap(), group, &['e'], &mut rng).unwrap();
        let mut optimiser = Optimiser::new(packing, MonteCarloConfig::default(), 3);
        let progress = optimiser.run(5000);
        assert!(progress.best_packing_fraction >= progress.packing_fraction);
        let structure = optimiser.finish();
        assert!(structure.packing_fraction >= progress.best_packing_fraction);
    }

    #[test]
    fn disks_hexagonal() {
        // A single disk in a hexagonal cell can only reach the hexagonal packing
        let (structure, optimiser) = optimise(Molecule::disk(), "p6mm", &['a'], 2000, 1);
        let expected = PI / (2. * 3f64.sqrt());
        assert!(
            (structure.packing_fraction - expected).abs() < 1e-6,
            "{}",
            structure.packing_fraction
        );
        assert!(!optimiser.packing().any_overlap());
    }

    #[test]
    fn squares_fill_space() {
        let (structure, _) = optimise(Polygon::regular(4).unwrap(), "p1", &['a'], 20_000, 2);
        assert!(
            structure.packing_fraction > 0.99,
            "{}",
            structure.packing_fraction
        );
    }
}
