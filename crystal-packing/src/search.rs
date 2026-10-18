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
}

fn start_optimiser<S: HardShape>(shape: S, job: &Job) -> Result<Optimiser<S>, Error> {
    let group = PlaneGroup::from_name(&job.candidate.group)?;
    let mut rng = Rng::seed_from_u64(job.seed);
    let packing = match &job.start {
        Some(structure) => Packing::from_structure(shape, structure)?,
        None => Packing::new(shape, group, &job.candidate.sites, &mut rng)?,
    };
    Ok(Optimiser::new(packing, job.config, rng.next_u64()))
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
                Ok(match shape {
                    AnyShape::Polygon(s) => {
                        AnyOptimiser::Polygon(Box::new(start_optimiser(s, job)?))
                    }
                    AnyShape::Molecule(s) => {
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
        }
    }

    pub fn is_complete(&self) -> bool {
        match &self.optimiser {
            Ok(AnyOptimiser::Polygon(o)) => o.is_complete(),
            Ok(AnyOptimiser::Molecule(o)) => o.is_complete(),
            Err(_) => true,
        }
    }

    /// Finish the job, compressing and validating the best structure found
    pub fn finish(mut self) -> Outcome {
        let result = match self.optimiser.as_mut() {
            Ok(AnyOptimiser::Polygon(o)) => finish_optimiser(o).map_err(|e| e.to_string()),
            Ok(AnyOptimiser::Molecule(o)) => finish_optimiser(o).map_err(|e| e.to_string()),
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

/// Configuration of a search over candidates
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
    pub seed: u64,
    pub monte_carlo: MonteCarloConfig,
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
            seed: 0,
            monte_carlo: MonteCarloConfig::default(),
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
    pub errors: Vec<String>,
}

/// A successive halving search over every candidate for a shape
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Search {
    pub config: SearchConfig,
    pub shape: ShapeSpec,
    pub results: Vec<CandidateResult>,
    active: Vec<usize>,
    round: usize,
    finished: bool,
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
                errors: vec![],
            })
            .collect();
        Ok(Search {
            active: (0..results.len()).collect(),
            config,
            shape,
            results,
            round: 0,
            finished: false,
        })
    }

    pub fn is_finished(&self) -> bool {
        self.finished
    }

    pub fn round(&self) -> usize {
        self.round
    }

    pub fn num_active(&self) -> usize {
        self.active.len()
    }

    /// The jobs for the current round, which is empty once the search is finished
    pub fn jobs(&self) -> Vec<Job> {
        if self.finished {
            return vec![];
        }
        let steps = self.config.steps * (self.config.eta as u64).pow(self.round as u32);
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
                });
            }
        }
        jobs
    }

    /// Submit the outcomes of the current round's jobs, advancing to the next round
    pub fn submit(&mut self, outcomes: Vec<Outcome>) {
        for outcome in outcomes {
            let Some(result) = self
                .results
                .iter_mut()
                .find(|r| r.candidate == outcome.candidate)
            else {
                continue;
            };
            result.rounds = self.round + 1;
            if let Some(error) = outcome.error {
                result.errors.push(error);
            }
            if let Some(structure) = outcome.structure {
                let better = result
                    .best
                    .as_ref()
                    .is_none_or(|b| structure.packing_fraction > b.packing_fraction);
                if better {
                    result.best = Some(structure);
                }
            }
        }
        if self.active.len() <= 1 {
            self.finished = true;
            return;
        }
        let score = |i: &usize| {
            self.results[*i]
                .best
                .as_ref()
                .map_or(0., |s| s.packing_fraction)
        };
        let mut ranked = self.active.clone();
        ranked.sort_by(|a, b| score(b).partial_cmp(&score(a)).unwrap_or(Ordering::Equal));
        let keep = ranked.len().div_ceil(self.config.eta).max(1);
        ranked.truncate(keep);
        self.active = ranked;
        self.round += 1;
    }

    /// All candidates ranked by their best packing fraction
    pub fn ranking(&self) -> Vec<&CandidateResult> {
        let mut ranked: Vec<&CandidateResult> = self.results.iter().collect();
        let score = |r: &CandidateResult| r.best.as_ref().map_or(0., |s| s.packing_fraction);
        ranked.sort_by(|a, b| score(b).total_cmp(&score(a)));
        ranked
    }

    pub fn best(&self) -> Option<&Structure> {
        self.ranking().first().and_then(|r| r.best.as_ref())
    }

    /// Run the whole search sequentially
    pub fn run(&mut self) -> Option<&Structure> {
        while !self.finished {
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
        };
        let outcome = run_job(&job);
        assert!(outcome.structure.is_none());
        assert!(outcome.error.unwrap().contains("symmetry"));
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
