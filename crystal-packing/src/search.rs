//
// search.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Searching for the densest packing over every plane group and choice of Wyckoff sites.
//!
//! A [`Candidate`] is a plane group with a set of occupied Wyckoff sites (an isopointal set). The
//! candidates compatible with a shape are enumerated, removing choices of sites which are related
//! by a shift of the origin and therefore describe the same structures.
//!
//! The [`Search`] spends the computational budget using successive halving: every candidate is
//! optimised for a short time, then the best fraction are optimised for longer, repeating until a
//! single candidate remains. The search is split into independent [`Job`]s which can be run in any
//! order or in parallel, then submitted back to the search. This allows the same search to be
//! driven by threads on the command line or by web workers in a browser.

use std::cmp::Ordering;
use std::collections::BTreeSet;

use anyhow::{anyhow, bail, Error};
use serde::{Deserialize, Serialize};

use crate::engine::{Packing, Structure};
use crate::geometry::Vec2;
use crate::hard_shape::{HardShape, Molecule, Polygon, ShapeSymmetry};
use crate::monte_carlo::{MonteCarloConfig, Optimiser, Progress};
use crate::potential::{Annealer, LennardJones};
use crate::random::Rng;
use crate::symmetry::{PlaneGroup, SymOp};
use crate::validate::validate;

/// A serialisable description of a shape
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum ShapeSpec {
    /// A regular polygon with vertices a unit distance from the centre
    Polygon { sides: usize },
    /// A polygon with vertices equally spaced in angle at the given distances from the centre
    Radial { name: String, radii: Vec<f64> },
    /// A polygon with the given vertices
    Vertices {
        name: String,
        vertices: Vec<[f64; 2]>,
    },
    /// A disk with unit radius
    Disk,
    /// A molecule of three disks
    Trimer {
        radius: f64,
        angle: f64,
        distance: f64,
    },
}

/// A shape built from a [`ShapeSpec`]
#[derive(Debug, Clone)]
pub enum AnyShape {
    Polygon(Polygon),
    Molecule(Molecule),
}

impl ShapeSpec {
    pub fn build(&self) -> Result<AnyShape, Error> {
        Ok(match self {
            ShapeSpec::Polygon { sides } => AnyShape::Polygon(Polygon::regular(*sides)?),
            ShapeSpec::Radial { name, radii } => {
                AnyShape::Polygon(Polygon::from_radial(name, radii)?)
            }
            ShapeSpec::Vertices { name, vertices } => AnyShape::Polygon(Polygon::from_vertices(
                name,
                vertices.iter().map(|v| Vec2::new(v[0], v[1])).collect(),
            )?),
            ShapeSpec::Disk => AnyShape::Molecule(Molecule::disk()),
            ShapeSpec::Trimer {
                radius,
                angle,
                distance,
            } => AnyShape::Molecule(Molecule::trimer(*radius, *angle, *distance)?),
        })
    }
}

impl AnyShape {
    pub fn symmetry(&self) -> &ShapeSymmetry {
        match self {
            AnyShape::Polygon(s) => s.symmetry(),
            AnyShape::Molecule(s) => s.symmetry(),
        }
    }

    pub fn area(&self) -> f64 {
        match self {
            AnyShape::Polygon(s) => s.area(),
            AnyShape::Molecule(s) => s.area(),
        }
    }

    pub fn name(&self) -> &str {
        match self {
            AnyShape::Polygon(s) => s.name(),
            AnyShape::Molecule(s) => s.name(),
        }
    }
}

/// A plane group with a set of occupied Wyckoff sites
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Candidate {
    pub group: String,
    pub sites: Vec<char>,
}

impl Candidate {
    pub fn label(&self) -> String {
        format!("{}:{}", self.group, self.sites.iter().collect::<String>())
    }

    fn seed_values(&self) -> Vec<u64> {
        let group = PlaneGroup::from_name(&self.group)
            .map(|g| g.number as u64)
            .unwrap_or(0);
        std::iter::once(group)
            .chain(self.sites.iter().map(|&c| c as u64))
            .collect()
    }
}

/// Find the most special Wyckoff position containing a point.
///
/// Any point lies on the general position, so the position with the lowest multiplicity which
/// contains the point is the one the point belongs to.
fn locate(group: &PlaneGroup, point: Vec2) -> Option<char> {
    let mut best: Option<(usize, char)> = None;
    for site in group.wyckoff.iter() {
        let contains = site.orbit.iter().any(|op| solve_on(op, point));
        if contains && best.is_none_or(|(m, _)| site.multiplicity < m) {
            best = Some((site.multiplicity, site.letter));
        }
    }
    best.map(|(_, letter)| letter)
}

/// Check whether there are parameters p such that op(p) is equivalent to the point
fn solve_on(op: &SymOp, point: Vec2) -> bool {
    let r = op.linear();
    let used = op.uses();
    for nx in -2..=2 {
        for ny in -2..=2 {
            let target = point - Vec2::new(op.translation[0], op.translation[1])
                + Vec2::new(nx as f64, ny as f64);
            let params = match used {
                [true, true] => r.inverse() * target,
                [true, false] | [false, true] => {
                    let column = if used[0] { r.column(0) } else { r.column(1) };
                    let t = target.dot(column) / column.norm_sq();
                    if used[0] {
                        Vec2::new(t, 0.)
                    } else {
                        Vec2::new(0., t)
                    }
                }
                [false, false] => Vec2::ZERO,
            };
            if (r * params - target).norm() < 1e-8 {
                return true;
            }
        }
    }
    false
}

/// The shifts of the origin which map a group onto itself (the translations of its normaliser),
/// as a mapping of Wyckoff letters for each shift.
fn origin_shifts(group: &PlaneGroup) -> Vec<Vec<(char, char)>> {
    let values = [0., 0.5, 1. / 3., 2. / 3.];
    let mut mappings = vec![];
    for &tx in values.iter() {
        for &ty in values.iter() {
            let t = Vec2::new(tx, ty);
            // Conjugating x → Rx + τ by the shift gives x → Rx + τ + (I - R)t
            let normalises = group.operations.iter().all(|op| {
                let shift = t - op.linear() * t;
                let translation = Vec2::new(op.translation[0], op.translation[1]) + shift;
                group.operations.iter().any(|o| {
                    let d = Vec2::new(o.translation[0], o.translation[1]) - translation;
                    o.rotation == op.rotation
                        && (d.x - d.x.round()).abs() < 1e-8
                        && (d.y - d.y.round()).abs() < 1e-8
                })
            });
            if !normalises {
                continue;
            }
            let generic = Vec2::new(0.1234567, 0.3141593);
            let mapping: Option<Vec<(char, char)>> = group
                .wyckoff
                .iter()
                .map(|w| locate(group, w.coordinates.apply(generic) + t).map(|l| (w.letter, l)))
                .collect();
            if let Some(mapping) = mapping {
                mappings.push(mapping);
            }
        }
    }
    mappings
}

/// Enumerate the candidates compatible with a shape.
///
/// Each candidate occupies between 1 and `max_sites` Wyckoff sites, with at most `max_shapes`
/// shapes in the unit cell. Sites with free parameters can be occupied multiple times, while fixed
/// sites can be occupied at most once. Candidates related by a shift of the origin are removed.
pub fn enumerate_candidates(
    symmetry: &ShapeSymmetry,
    groups: &[&PlaneGroup],
    max_sites: usize,
    max_shapes: usize,
) -> Vec<Candidate> {
    let mut candidates = vec![];
    for group in groups {
        let compatible: Vec<_> = group
            .wyckoff
            .iter()
            .filter(|w| w.orientation_freedom(symmetry).is_some())
            .collect();
        let shifts = origin_shifts(group);
        let mut seen: BTreeSet<Vec<char>> = BTreeSet::new();

        // Multisets of compatible sites, as non-decreasing sequences of indices
        let mut stack: Vec<Vec<usize>> = (0..compatible.len()).map(|i| vec![i]).collect();
        while let Some(choice) = stack.pop() {
            let shapes: usize = choice.iter().map(|&i| compatible[i].multiplicity).sum();
            let repeated_fixed = choice
                .windows(2)
                .any(|w| w[0] == w[1] && compatible[w[0]].is_fixed());
            if shapes > max_shapes || repeated_fixed {
                continue;
            }
            let mut letters: Vec<char> = choice.iter().map(|&i| compatible[i].letter).collect();
            letters.sort();
            // The canonical form is the smallest relabelling under the origin shifts
            let canonical = shifts
                .iter()
                .map(|mapping| {
                    let mut mapped: Vec<char> = letters
                        .iter()
                        .map(|l| {
                            mapping
                                .iter()
                                .find(|(from, _)| from == l)
                                .map_or(*l, |(_, to)| *to)
                        })
                        .collect();
                    mapped.sort();
                    mapped
                })
                .min()
                .unwrap_or_else(|| letters.clone());
            if seen.insert(canonical.clone()) {
                candidates.push(Candidate {
                    group: group.name.clone(),
                    sites: canonical,
                });
            }
            if choice.len() < max_sites {
                let last = *choice.last().unwrap();
                for next in last..compatible.len() {
                    let mut extended = choice.clone();
                    extended.push(next);
                    stack.push(extended);
                }
            }
        }
    }
    candidates.sort_by(|a, b| {
        let ga = PlaneGroup::from_name(&a.group)
            .map(|g| g.number)
            .unwrap_or(0);
        let gb = PlaneGroup::from_name(&b.group)
            .map(|g| g.number)
            .unwrap_or(0);
        ga.cmp(&gb).then(a.sites.cmp(&b.sites))
    });
    candidates
}

/// A unit of work: optimise a single candidate from a random or given starting structure
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Job {
    pub shape: ShapeSpec,
    pub candidate: Candidate,
    pub seed: u64,
    pub config: MonteCarloConfig,
    pub start: Option<Structure>,
    /// Minimise the energy of this potential rather than maximising the packing fraction
    #[serde(default)]
    pub potential: Option<LennardJones>,
}

/// The result of running a [`Job`]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Outcome {
    pub candidate: Candidate,
    pub structure: Option<Structure>,
    pub error: Option<String>,
}

impl Outcome {
    pub fn packing_fraction(&self) -> f64 {
        self.structure.as_ref().map_or(0., |s| s.packing_fraction)
    }

    /// The value being maximised, see [`Structure::score`]
    pub fn score(&self) -> f64 {
        self.structure
            .as_ref()
            .map_or(f64::NEG_INFINITY, |s| s.score())
    }
}

fn start_packing<S: HardShape>(shape: S, job: &Job) -> Result<(Packing<S>, u64), Error> {
    let group = PlaneGroup::from_name(&job.candidate.group)?;
    let mut rng = Rng::seed_from_u64(job.seed);
    let packing = match &job.start {
        Some(structure) => Packing::from_structure(shape, structure)?,
        None => Packing::new(shape, group, &job.candidate.sites, &mut rng)?,
    };
    Ok((packing, rng.next_u64()))
}

fn start_optimiser<S: HardShape>(shape: S, job: &Job) -> Result<Optimiser<S>, Error> {
    let (packing, seed) = start_packing(shape, job)?;
    Ok(Optimiser::new(packing, job.config, seed))
}

fn start_annealer(shape: AnyShape, potential: LennardJones, job: &Job) -> Result<Annealer, Error> {
    let AnyShape::Molecule(shape) = shape else {
        bail!("A potential can only be used with molecules made of disks");
    };
    let (packing, seed) = start_packing(shape, job)?;
    Ok(Annealer::new(packing, potential, job.config.steps, seed))
}

fn finish_annealer(annealer: &mut Annealer) -> Result<Structure, Error> {
    let structure = annealer.finish();
    match structure.energy {
        Some(e) if e.is_finite() && structure.lattice.is_valid() => Ok(structure),
        _ => bail!("Annealing produced an invalid structure"),
    }
}

fn finish_optimiser<S: HardShape>(optimiser: &mut Optimiser<S>) -> Result<Structure, Error> {
    let structure = optimiser.finish();
    validate(optimiser.packing().shape(), &structure)
        .map_err(|e| anyhow!("Optimised structure failed validation: {}", e))?;
    Ok(structure)
}

enum AnyOptimiser {
    Polygon(Box<Optimiser<Polygon>>),
    Molecule(Box<Optimiser<Molecule>>),
    Annealer(Box<Annealer>),
}

/// Runs a [`Job`] in slices, so the progress can be reported and the work interleaved with
/// other tasks, such as in a web worker.
pub struct JobRunner {
    candidate: Candidate,
    optimiser: Result<AnyOptimiser, String>,
}

impl JobRunner {
    pub fn new(job: &Job) -> JobRunner {
        let optimiser = job
            .shape
            .build()
            .and_then(|shape| {
                Ok(match (shape, job.potential) {
                    (shape, Some(potential)) => {
                        AnyOptimiser::Annealer(Box::new(start_annealer(shape, potential, job)?))
                    }
                    (AnyShape::Polygon(s), None) => {
                        AnyOptimiser::Polygon(Box::new(start_optimiser(s, job)?))
                    }
                    (AnyShape::Molecule(s), None) => {
                        AnyOptimiser::Molecule(Box::new(start_optimiser(s, job)?))
                    }
                })
            })
            .map_err(|e| e.to_string());
        JobRunner {
            candidate: job.candidate.clone(),
            optimiser,
        }
    }

    /// Run up to `steps` Monte Carlo steps, returning the progress
    pub fn run(&mut self, steps: u64) -> Option<Progress> {
        match self.optimiser.as_mut().ok()? {
            AnyOptimiser::Polygon(o) => Some(o.run(steps)),
            AnyOptimiser::Molecule(o) => Some(o.run(steps)),
            AnyOptimiser::Annealer(o) => Some(o.run(steps)),
        }
    }

    pub fn is_complete(&self) -> bool {
        match &self.optimiser {
            Ok(AnyOptimiser::Polygon(o)) => o.is_complete(),
            Ok(AnyOptimiser::Molecule(o)) => o.is_complete(),
            Ok(AnyOptimiser::Annealer(o)) => o.is_complete(),
            Err(_) => true,
        }
    }

    /// Finish the job, compressing and validating the best structure found
    pub fn finish(mut self) -> Outcome {
        let result = match self.optimiser.as_mut() {
            Ok(AnyOptimiser::Polygon(o)) => finish_optimiser(o).map_err(|e| e.to_string()),
            Ok(AnyOptimiser::Molecule(o)) => finish_optimiser(o).map_err(|e| e.to_string()),
            Ok(AnyOptimiser::Annealer(o)) => finish_annealer(o).map_err(|e| e.to_string()),
            Err(e) => Err(e.clone()),
        };
        match result {
            Ok(structure) => Outcome {
                candidate: self.candidate,
                structure: Some(structure),
                error: None,
            },
            Err(error) => Outcome {
                candidate: self.candidate,
                structure: None,
                error: Some(error),
            },
        }
    }
}

/// Run a job to completion
pub fn run_job(job: &Job) -> Outcome {
    let mut runner = JobRunner::new(job);
    runner.run(job.config.steps);
    runner.finish()
}

/// How the computational budget is divided between the candidates
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Strategy {
    /// Successive halving, focusing the budget on the most promising candidates
    #[default]
    Halving,
    /// Every candidate is optimised for every round, to survey all the groups
    Exhaustive,
}

/// Configuration of a search over candidates
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct SearchConfig {
    /// The groups to search, or all groups when empty
    pub groups: Vec<String>,
    /// The maximum number of occupied Wyckoff sites
    pub max_sites: usize,
    /// The maximum number of shapes in the unit cell
    pub max_shapes: usize,
    /// The number of independent optimisations of each candidate per round
    pub replicas: usize,
    /// The number of Monte Carlo steps per optimisation in the first round
    pub steps: u64,
    /// The fraction of candidates kept after each round is 1 / eta, while the number of steps
    /// increases by a factor of eta.
    pub eta: usize,
    /// The number of steps grows by at most eta to this power
    pub max_step_growth: u32,
    pub strategy: Strategy,
    /// Once the candidates are narrowed down, the search finishes after this many rounds
    /// without the best score improving
    pub patience: usize,
    /// The relative change of the best score counted as an improvement
    pub improvement: f64,
    /// The maximum number of rounds, excluding the refinement round
    pub max_rounds: usize,
    /// Candidates within this fraction of the best score are refined with the final budget
    pub within: f64,
    pub seed: u64,
    pub monte_carlo: MonteCarloConfig,
    /// Minimise the energy of this potential rather than maximising the packing fraction
    pub potential: Option<LennardJones>,
}

impl Default for SearchConfig {
    fn default() -> Self {
        SearchConfig {
            groups: vec![],
            max_sites: 1,
            max_shapes: 12,
            replicas: 4,
            steps: 10_000,
            eta: 3,
            max_step_growth: 3,
            strategy: Strategy::Halving,
            patience: 2,
            improvement: 1e-5,
            max_rounds: 8,
            within: 0.02,
            seed: 0,
            monte_carlo: MonteCarloConfig::default(),
            potential: None,
        }
    }
}

/// The best result found for a candidate
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CandidateResult {
    pub candidate: Candidate,
    pub best: Option<Structure>,
    /// The last round of the search the candidate took part in
    pub rounds: usize,
    /// The largest number of steps of any optimisation of the candidate
    #[serde(default)]
    pub steps: u64,
    pub errors: Vec<String>,
}

impl CandidateResult {
    pub fn score(&self) -> f64 {
        self.best.as_ref().map_or(f64::NEG_INFINITY, |s| s.score())
    }
}

/// The stages of a search
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Phase {
    /// Narrowing down the candidates
    Explore,
    /// Optimising the remaining candidates until they stop improving
    Converge,
    /// Giving every candidate close to the best the same budget as the best
    Refine,
    Finished,
}

/// A search over every candidate for a shape
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Search {
    pub config: SearchConfig,
    pub shape: ShapeSpec,
    pub results: Vec<CandidateResult>,
    active: Vec<usize>,
    round: usize,
    phase: Phase,
    /// The number of consecutive rounds without improvement in the converge phase
    stale: usize,
    /// The number of steps of the latest round before refinement
    final_steps: u64,
}

impl Search {
    pub fn new(shape: ShapeSpec, config: SearchConfig) -> Result<Search, Error> {
        let built = shape.build()?;
        let groups: Vec<&PlaneGroup> = if config.groups.is_empty() {
            PlaneGroup::all().iter().collect()
        } else {
            config
                .groups
                .iter()
                .map(|g| PlaneGroup::from_name(g))
                .collect::<Result<_, _>>()?
        };
        if config.eta < 2 {
            bail!("eta needs to be at least 2");
        }
        let candidates = enumerate_candidates(
            built.symmetry(),
            &groups,
            config.max_sites,
            config.max_shapes,
        );
        if candidates.is_empty() {
            bail!("No candidate structures are compatible with the shape");
        }
        let results: Vec<CandidateResult> = candidates
            .into_iter()
            .map(|candidate| CandidateResult {
                candidate,
                best: None,
                rounds: 0,
                steps: 0,
                errors: vec![],
            })
            .collect();
        let phase = match (config.strategy, results.len()) {
            (Strategy::Exhaustive, _) | (_, 1) => Phase::Converge,
            _ => Phase::Explore,
        };
        Ok(Search {
            active: (0..results.len()).collect(),
            config,
            shape,
            results,
            round: 0,
            phase,
            stale: 0,
            final_steps: 0,
        })
    }

    pub fn is_finished(&self) -> bool {
        self.phase == Phase::Finished
    }

    pub fn phase(&self) -> Phase {
        self.phase
    }

    /// The number of steps for each optimisation in the current round
    fn round_steps(&self) -> u64 {
        if self.phase == Phase::Refine {
            return self.final_steps;
        }
        let growth = (self.round as u32).min(self.config.max_step_growth);
        self.config.steps * (self.config.eta as u64).pow(growth)
    }

    pub fn round(&self) -> usize {
        self.round
    }

    pub fn num_active(&self) -> usize {
        self.active.len()
    }

    /// The jobs for the current round, which is empty once the search is finished
    pub fn jobs(&self) -> Vec<Job> {
        if self.is_finished() {
            return vec![];
        }
        let steps = self.round_steps();
        let mut jobs = vec![];
        for &index in self.active.iter() {
            let result = &self.results[index];
            for replica in 0..self.config.replicas {
                let mut values = result.candidate.seed_values();
                values.push(self.round as u64);
                values.push(replica as u64);
                let mut config = self.config.monte_carlo;
                config.steps = steps;
                // After the first round, half the replicas continue from the best structure,
                // starting from a higher pressure so they explore around it.
                let start = if replica % 2 == 1 {
                    result.best.clone()
                } else {
                    None
                };
                if start.is_some() {
                    config.pressure_start = config.pressure_start.max(10.);
                }
                jobs.push(Job {
                    shape: self.shape.clone(),
                    candidate: result.candidate.clone(),
                    seed: Rng::derive_seed(self.config.seed, &values),
                    config,
                    start,
                    potential: self.config.potential,
                });
            }
        }
        jobs
    }

    /// Submit the outcomes of the current round's jobs, advancing to the next round
    pub fn submit(&mut self, outcomes: Vec<Outcome>) {
        let steps = self.round_steps();
        let previous_best = self.best_score();
        let mut improved = false;
        for outcome in outcomes {
            let Some(result) = self
                .results
                .iter_mut()
                .find(|r| r.candidate == outcome.candidate)
            else {
                continue;
            };
            result.rounds = self.round + 1;
            result.steps = result.steps.max(steps);
            if let Some(error) = outcome.error {
                result.errors.push(error);
            }
            if let Some(structure) = outcome.structure {
                let old = result.score();
                if structure.score() > old {
                    if structure.score() > old + self.config.improvement * old.abs() {
                        improved = true;
                    }
                    result.best = Some(structure);
                }
            }
        }
        // In the halving phase only the change of the overall best matters
        let best_improved =
            self.best_score() > previous_best + self.config.improvement * previous_best.abs();

        match self.phase {
            Phase::Explore => {
                let mut ranked = self.active.clone();
                ranked.sort_by(|a, b| {
                    self.results[*b]
                        .score()
                        .partial_cmp(&self.results[*a].score())
                        .unwrap_or(Ordering::Equal)
                });
                ranked.truncate(ranked.len().div_ceil(self.config.eta).max(1));
                self.active = ranked;
                if self.active.len() == 1 {
                    self.phase = Phase::Converge;
                }
                self.final_steps = steps;
                self.round += 1;
            }
            Phase::Converge => {
                let progress = match self.config.strategy {
                    Strategy::Halving => best_improved,
                    Strategy::Exhaustive => improved,
                };
                self.stale = if progress { 0 } else { self.stale + 1 };
                self.final_steps = steps;
                self.round += 1;
                if self.stale >= self.config.patience || self.round >= self.config.max_rounds {
                    self.start_refinement();
                }
            }
            Phase::Refine => {
                self.round += 1;
                self.phase = Phase::Finished;
            }
            Phase::Finished => {}
        }
        if self.phase != Phase::Finished && self.round >= self.config.max_rounds + 1 {
            self.phase = Phase::Finished;
        }
    }

    /// Give every candidate near the best which had a smaller budget a round with the final budget
    fn start_refinement(&mut self) {
        let threshold = self.threshold(self.config.within);
        let refine: Vec<usize> = (0..self.results.len())
            .filter(|&i| {
                let result = &self.results[i];
                result.score() >= threshold && result.steps < self.final_steps
            })
            .collect();
        if refine.is_empty() {
            self.phase = Phase::Finished;
        } else {
            self.active = refine;
            self.phase = Phase::Refine;
        }
    }

    fn best_score(&self) -> f64 {
        self.results
            .iter()
            .map(|r| r.score())
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// The lowest score within a relative tolerance of the best score
    fn threshold(&self, within: f64) -> f64 {
        let best = self.best_score();
        best - within * best.abs()
    }

    /// The candidates with a score within a fraction `within` of the best, best first
    pub fn near_best(&self, within: f64) -> Vec<&CandidateResult> {
        let threshold = self.threshold(within);
        self.ranking()
            .into_iter()
            .filter(|r| r.best.is_some() && r.score() >= threshold)
            .collect()
    }

    /// The best candidate of every group searched, in the order of the International Tables
    pub fn best_per_group(&self) -> Vec<(String, Option<&CandidateResult>)> {
        let mut groups: Vec<&PlaneGroup> = if self.config.groups.is_empty() {
            PlaneGroup::all().iter().collect()
        } else {
            self.config
                .groups
                .iter()
                .filter_map(|g| PlaneGroup::from_name(g).ok())
                .collect()
        };
        groups.sort_by_key(|g| g.number);
        groups
            .into_iter()
            .map(|group| {
                let best = self
                    .results
                    .iter()
                    .filter(|r| r.candidate.group == group.name && r.best.is_some())
                    .max_by(|a, b| a.score().total_cmp(&b.score()));
                (group.name.clone(), best)
            })
            .collect()
    }

    /// All candidates ranked by their best packing fraction
    pub fn ranking(&self) -> Vec<&CandidateResult> {
        let mut ranked: Vec<&CandidateResult> = self.results.iter().collect();
        ranked.sort_by(|a, b| b.score().total_cmp(&a.score()));
        ranked
    }

    pub fn best(&self) -> Option<&Structure> {
        self.ranking().first().and_then(|r| r.best.as_ref())
    }

    /// Run the whole search sequentially
    pub fn run(&mut self) -> Option<&Structure> {
        while !self.is_finished() {
            let outcomes = self.jobs().iter().map(run_job).collect();
            self.submit(outcomes);
        }
        self.best()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn candidates(spec: ShapeSpec, groups: &[&str], max_sites: usize) -> Vec<Candidate> {
        let shape = spec.build().unwrap();
        let groups: Vec<_> = groups
            .iter()
            .map(|g| PlaneGroup::from_name(g).unwrap())
            .collect();
        enumerate_candidates(shape.symmetry(), &groups, max_sites, 24)
    }

    #[test]
    fn p2_special_sites_are_equivalent() {
        // The four fixed sites of p2 are related by origin shifts, leaving a and e
        let found = candidates(ShapeSpec::Polygon { sides: 4 }, &["p2"], 1);
        let labels: Vec<String> = found.iter().map(|c| c.label()).collect();
        assert_eq!(labels, vec!["p2:a", "p2:e"]);
    }

    #[test]
    fn p2_pairs() {
        let found = candidates(ShapeSpec::Polygon { sides: 4 }, &["p2"], 2);
        let labels: Vec<String> = found.iter().map(|c| c.label()).collect();
        // Origin shifts reduce the pairs of distinct fixed sites to ab, ac and ad, and pairs
        // including the general position to ae or ee. Only translations of the normaliser are
        // used, so ab, ac and ad (related by a change of basis) are kept separate.
        assert_eq!(
            labels,
            vec!["p2:a", "p2:ab", "p2:ac", "p2:ad", "p2:ae", "p2:e", "p2:ee"]
        );
    }

    #[test]
    fn incompatible_sites_removed() {
        let found = candidates(ShapeSpec::Polygon { sides: 5 }, &["p2", "p4", "p6mm"], 1);
        let labels: Vec<String> = found.iter().map(|c| c.label()).collect();
        // A pentagon has no 2-fold or 4-fold symmetry, though it can sit on a mirror line
        assert_eq!(labels, vec!["p2:e", "p4:d", "p6mm:d", "p6mm:e", "p6mm:f"]);
    }

    #[test]
    fn mirror_lines_for_triangles() {
        // Triangles have mirrors so can sit on the mirror lines of p1m1
        let found = candidates(ShapeSpec::Polygon { sides: 3 }, &["p1m1"], 1);
        let labels: Vec<String> = found.iter().map(|c| c.label()).collect();
        assert_eq!(labels, vec!["p1m1:a", "p1m1:c"]);
    }

    #[test]
    fn all_groups_disk() {
        let found = candidates(ShapeSpec::Disk, &[], 1);
        assert!(found.is_empty());
        let shape = ShapeSpec::Disk.build().unwrap();
        let groups: Vec<_> = PlaneGroup::all().iter().collect();
        let found = enumerate_candidates(shape.symmetry(), &groups, 1, 12);
        // Every group has at least one candidate for a disk
        for g in PlaneGroup::all() {
            assert!(found.iter().any(|c| c.group == g.name), "{}", g.name);
        }
    }

    #[test]
    fn search_finds_hexagonal_disks() {
        let config = SearchConfig {
            groups: vec!["p1".into(), "p2".into(), "p6mm".into()],
            replicas: 2,
            steps: 5000,
            ..Default::default()
        };
        let mut search = Search::new(ShapeSpec::Disk, config).unwrap();
        let best = search.run().unwrap().clone();
        let expected = std::f64::consts::PI / (2. * 3f64.sqrt());
        assert!(
            (best.packing_fraction - expected).abs() < 1e-4,
            "{}",
            best.packing_fraction
        );
        assert!(search.is_finished());
    }

    #[test]
    fn jobs_are_deterministic() {
        let config = SearchConfig {
            groups: vec!["p2".into()],
            replicas: 2,
            steps: 1000,
            ..Default::default()
        };
        let search = Search::new(ShapeSpec::Polygon { sides: 5 }, config).unwrap();
        let a: Vec<Outcome> = search.jobs().iter().map(run_job).collect();
        let b: Vec<Outcome> = search.jobs().iter().rev().map(run_job).rev().collect();
        assert_eq!(a, b);
    }

    #[test]
    fn sliced_job_matches_complete_job() {
        let job = Job {
            shape: ShapeSpec::Polygon { sides: 5 },
            candidate: Candidate {
                group: "p2".into(),
                sites: vec!['e'],
            },
            seed: 3,
            config: MonteCarloConfig {
                steps: 3000,
                ..Default::default()
            },
            start: None,
            potential: None,
        };
        let mut runner = JobRunner::new(&job);
        while !runner.is_complete() {
            runner.run(250);
        }
        assert_eq!(runner.finish(), run_job(&job));
    }

    #[test]
    fn invalid_job_reports_error() {
        let job = Job {
            shape: ShapeSpec::Polygon { sides: 3 },
            candidate: Candidate {
                group: "p4".into(),
                sites: vec!['a'],
            },
            seed: 0,
            config: MonteCarloConfig::default(),
            start: None,
            potential: None,
        };
        let outcome = run_job(&job);
        assert!(outcome.structure.is_none());
        assert!(outcome.error.unwrap().contains("symmetry"));
    }

    /// A structure with a given packing fraction, standing in for the result of a job
    fn fake_structure(candidate: &Candidate, phi: f64) -> Structure {
        let group = PlaneGroup::from_name(&candidate.group).unwrap();
        Structure {
            shape: "fake".into(),
            group: candidate.group.clone(),
            sites: vec![],
            lattice: crate::lattice::Lattice::with_length(group.family, 1.),
            num_shapes: 1,
            packing_fraction: phi,
            energy: None,
        }
    }

    /// Submit an outcome for every job, with the packing fraction given by `score`
    fn submit_round(search: &mut Search, score: impl Fn(&Job) -> f64) -> Vec<Job> {
        let jobs = search.jobs();
        let outcomes = jobs
            .iter()
            .map(|job| Outcome {
                candidate: job.candidate.clone(),
                structure: Some(fake_structure(&job.candidate, score(job))),
                error: None,
            })
            .collect();
        search.submit(outcomes);
        jobs
    }

    fn disk_search(config: SearchConfig) -> Search {
        Search::new(ShapeSpec::Disk, config).unwrap()
    }

    fn all_groups() -> Vec<String> {
        PlaneGroup::all().iter().map(|g| g.name.clone()).collect()
    }

    /// A score for each candidate which doesn't change between rounds
    fn fixed_score(job: &Job) -> f64 {
        let label = job.candidate.label();
        let hash = label
            .bytes()
            .fold(7u64, |h, b| h.wrapping_mul(31).wrapping_add(b as u64));
        0.5 + (hash % 1000) as f64 / 4000.
    }

    #[test]
    fn halving_keeps_a_fraction_of_candidates() {
        let mut search = disk_search(SearchConfig {
            groups: all_groups(),
            ..Default::default()
        });
        let n = search.num_active();
        assert!(n > 9, "{}", n);
        submit_round(&mut search, fixed_score);
        assert_eq!(search.num_active(), n.div_ceil(3));
        assert!(!search.is_finished());
    }

    #[test]
    fn continues_while_improving() {
        let mut search = disk_search(SearchConfig {
            groups: vec!["p1".into()],
            patience: 2,
            max_rounds: 20,
            ..Default::default()
        });
        assert_eq!(search.num_active(), 1);
        // Keep improving for 5 rounds
        for round in 0..5 {
            submit_round(&mut search, |_| 0.5 + 0.01 * round as f64);
            assert!(
                !search.is_finished(),
                "Finished while improving in round {}",
                round
            );
        }
        // Then stop improving, finishing after `patience` rounds
        submit_round(&mut search, |_| 0.5);
        assert!(!search.is_finished());
        submit_round(&mut search, |_| 0.5);
        assert!(search.is_finished());
    }

    #[test]
    fn stops_at_max_rounds() {
        let mut search = disk_search(SearchConfig {
            groups: vec!["p1".into()],
            max_rounds: 3,
            ..Default::default()
        });
        for round in 0..3 {
            assert!(!search.is_finished());
            submit_round(&mut search, |_| 0.5 + 0.01 * round as f64);
        }
        assert!(search.is_finished());
    }

    #[test]
    fn step_growth_is_capped() {
        let mut search = disk_search(SearchConfig {
            groups: vec!["p1".into()],
            steps: 100,
            eta: 3,
            max_step_growth: 2,
            max_rounds: 10,
            ..Default::default()
        });
        let mut steps = vec![];
        for round in 0..5 {
            let jobs = submit_round(&mut search, |_| 0.5 + 0.01 * round as f64);
            steps.push(jobs[0].config.steps);
        }
        assert_eq!(steps, vec![100, 300, 900, 900, 900]);
    }

    #[test]
    fn refines_candidates_near_the_best() {
        let mut search = disk_search(SearchConfig {
            groups: all_groups(),
            within: 0.02,
            ..Default::default()
        });
        // Two candidates are close to the best, and the rest far below
        let score = |job: &Job| match job.candidate.label().as_str() {
            "p6mm:a" => 0.90,
            "p3:a" => 0.89,
            "p1:a" => 0.885,
            _ => 0.5,
        };
        let mut final_steps = 0;
        while !search.is_finished() {
            let jobs = submit_round(&mut search, score);
            final_steps = final_steps.max(jobs[0].config.steps);
        }
        let near: Vec<String> = search
            .near_best(0.02)
            .iter()
            .map(|r| r.candidate.label())
            .collect();
        assert_eq!(near, vec!["p6mm:a", "p3:a", "p1:a"]);
        // Every candidate within the tolerance had the full budget
        for result in search.near_best(0.02) {
            assert_eq!(result.steps, final_steps, "{}", result.candidate.label());
        }
        // Candidates far from the best were not refined
        let far = search
            .results
            .iter()
            .find(|r| r.candidate.label() == "p2mm:a")
            .unwrap();
        assert!(far.steps < final_steps);
    }

    #[test]
    fn exhaustive_runs_every_candidate() {
        let mut search = disk_search(SearchConfig {
            groups: all_groups(),
            strategy: Strategy::Exhaustive,
            patience: 1,
            ..Default::default()
        });
        let n = search.num_active();
        submit_round(&mut search, fixed_score);
        assert_eq!(search.num_active(), n);
        submit_round(&mut search, fixed_score);
        assert!(search.is_finished());
        assert!(search.results.iter().all(|r| r.best.is_some()));
    }

    #[test]
    fn best_per_group_lists_every_group() {
        let mut search = disk_search(SearchConfig {
            groups: all_groups(),
            strategy: Strategy::Exhaustive,
            patience: 1,
            ..Default::default()
        });
        while !search.is_finished() {
            submit_round(&mut search, fixed_score);
        }
        let per_group = search.best_per_group();
        assert_eq!(per_group.len(), 17);
        for (group, result) in per_group {
            let result = result.expect("Every group has a candidate for disks");
            assert_eq!(result.candidate.group, group);
            // It is the best of the group's candidates
            for other in search.results.iter().filter(|r| r.candidate.group == group) {
                assert!(
                    other.best.as_ref().unwrap().score() <= result.best.as_ref().unwrap().score()
                );
            }
        }
    }

    #[test]
    fn near_best_uses_relative_tolerance() {
        let mut search = disk_search(SearchConfig {
            groups: all_groups(),
            strategy: Strategy::Exhaustive,
            patience: 1,
            ..Default::default()
        });
        while !search.is_finished() {
            submit_round(&mut search, fixed_score);
        }
        let best = search.best().unwrap().score();
        let near = search.near_best(0.05);
        assert!(!near.is_empty());
        let scores: Vec<f64> = near
            .iter()
            .map(|r| r.best.as_ref().unwrap().score())
            .collect();
        assert!(scores.windows(2).all(|w| w[0] >= w[1]), "Not sorted");
        assert!(scores.iter().all(|&s| s >= best * 0.95));
        let excluded = search
            .results
            .iter()
            .filter(|r| r.best.as_ref().unwrap().score() < best * 0.95)
            .count();
        assert_eq!(near.len() + excluded, search.results.len());
    }

    #[test]
    fn shape_spec_serialisation() {
        let spec = ShapeSpec::Trimer {
            radius: 0.6,
            angle: 120.,
            distance: 1.,
        };
        let json = serde_json::to_string(&spec).unwrap();
        assert_eq!(
            json,
            r#"{"type":"trimer","radius":0.6,"angle":120.0,"distance":1.0}"#
        );
        assert_eq!(serde_json::from_str::<ShapeSpec>(&json).unwrap(), spec);
    }
}
