//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use anyhow::{bail, Context, Error};
use clap::{Args, Parser, Subcommand, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, LevelFilter};
use rayon::prelude::*;

use crystal_packing::engine::Structure;
use crystal_packing::monte_carlo::MonteCarloConfig;
use crystal_packing::potential::LennardJones;
use crystal_packing::render::render_structure;
use crystal_packing::search::{
    run_job, Candidate, CandidateResult, Job, Search, SearchConfig, ShapeSpec, Strategy,
};
use crystal_packing::Rng;

#[derive(Debug, Parser)]
#[command(
    name = "packing",
    version,
    about = "Find dense packings of 2D shapes in crystals"
)]
struct Cli {
    /// Pass many times for more log output
    #[arg(long, short, action = clap::ArgAction::Count, global = true)]
    verbosity: u8,

    /// The number of threads to use, defaulting to the number of cores
    #[arg(long, global = true)]
    threads: Option<usize>,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Search the plane groups for the densest packing of a shape
    Search(SearchArgs),
    /// Optimise every candidate of every plane group, reporting the best of each group
    Survey(SurveyArgs),
    /// Pack independent shapes in a periodic cell, without any symmetry
    Unconstrained(UnconstrainedArgs),
    /// Optimise a single plane group and choice of Wyckoff sites
    Optimise(OptimiseArgs),
    /// Render a saved structure as an SVG image
    Render(RenderArgs),
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum PotentialKind {
    /// Hard shapes, maximising the packing fraction
    Hard,
    /// Lennard-Jones atoms (σ is twice the atom radius), minimising the lattice energy
    Lj,
}

#[derive(Debug, Args)]
struct PotentialArgs {
    /// The interaction between shapes
    #[arg(long, value_enum, default_value = "hard")]
    potential: PotentialKind,
    /// The depth of the Lennard-Jones potential
    #[arg(long, default_value = "1")]
    epsilon: f64,
    /// The cutoff of the Lennard-Jones potential in units of σ
    #[arg(long, default_value = "2.5")]
    cutoff: f64,
    /// The temperature at the start of the annealing
    #[arg(long, default_value = "0.1")]
    kt_start: f64,
    /// The temperature at the end of the annealing
    #[arg(long, default_value = "0.001")]
    kt_end: f64,
    /// The pressure applied to the cell when annealing
    #[arg(long, default_value = "0")]
    pressure: f64,
}

impl PotentialArgs {
    fn config(&self) -> Option<LennardJones> {
        match self.potential {
            PotentialKind::Hard => None,
            PotentialKind::Lj => Some(LennardJones {
                epsilon: self.epsilon,
                cutoff: self.cutoff,
                kt_start: self.kt_start,
                kt_end: self.kt_end,
                pressure: self.pressure,
            }),
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ShapeKind {
    /// A regular polygon, with the number of sides given by --sides
    Polygon,
    /// A unit disk
    Disk,
    /// A molecule of three disks, see --radius, --angle and --distance
    Trimer,
    /// A shape defined by distances from the centre, given by --radii
    Radial,
    /// A polygon with vertices read from a JSON file of [x, y] pairs, given by --vertices
    Vertices,
}

#[derive(Debug, Args)]
struct ShapeArgs {
    /// The shape to pack
    #[arg(long, value_enum, default_value = "polygon")]
    shape: ShapeKind,
    /// The number of sides of a regular polygon
    #[arg(long, default_value = "4")]
    sides: usize,
    /// The radius of the outer disks of a trimer
    #[arg(long, default_value = "0.637556")]
    radius: f64,
    /// The angle between the outer disks of a trimer in degrees
    #[arg(long, default_value = "120")]
    angle: f64,
    /// The distance from the central disk to the outer disks of a trimer
    #[arg(long, default_value = "1")]
    distance: f64,
    /// Comma separated distances from the centre for a radial shape
    #[arg(long, value_delimiter = ',')]
    radii: Vec<f64>,
    /// A JSON file with a list of [x, y] vertices
    #[arg(long)]
    vertices: Option<PathBuf>,
    /// The name of a radial or vertex shape
    #[arg(long, default_value = "Shape")]
    name: String,
}

impl ShapeArgs {
    fn spec(&self) -> Result<ShapeSpec, Error> {
        Ok(match self.shape {
            ShapeKind::Polygon => ShapeSpec::Polygon { sides: self.sides },
            ShapeKind::Disk => ShapeSpec::Disk,
            ShapeKind::Trimer => ShapeSpec::Trimer {
                radius: self.radius,
                angle: self.angle,
                distance: self.distance,
            },
            ShapeKind::Radial => {
                if self.radii.is_empty() {
                    bail!("--radii is required for a radial shape");
                }
                ShapeSpec::Radial {
                    name: self.name.clone(),
                    radii: self.radii.clone(),
                }
            }
            ShapeKind::Vertices => {
                let path = self
                    .vertices
                    .as_ref()
                    .context("--vertices is required for a vertex shape")?;
                let vertices: Vec<[f64; 2]> = serde_json::from_str(&fs::read_to_string(path)?)
                    .context("Expected a JSON list of [x, y] pairs")?;
                ShapeSpec::Vertices {
                    name: self.name.clone(),
                    vertices,
                }
            }
        })
    }
}

#[derive(Debug, Args)]
struct MonteCarloArgs {
    /// The number of Monte Carlo steps for each optimisation (in the first round of a search)
    #[arg(long, default_value = "10000")]
    steps: u64,
    /// The number of independent optimisations of each candidate
    #[arg(long, default_value = "4")]
    replicas: usize,
    /// The reduced pressure at the start of the compression
    #[arg(long, default_value = "0.1")]
    pressure_start: f64,
    /// The reduced pressure at the end of the compression
    #[arg(long, default_value = "1e6")]
    pressure_end: f64,
    /// The seed for the random number generator, chosen at random when not given
    #[arg(long)]
    seed: Option<u64>,
}

impl MonteCarloArgs {
    fn config(&self) -> MonteCarloConfig {
        MonteCarloConfig {
            steps: self.steps,
            pressure_start: self.pressure_start,
            pressure_end: self.pressure_end,
            ..Default::default()
        }
    }

    fn seed(&self) -> u64 {
        self.seed.unwrap_or_else(|| {
            let nanos = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos() as u64)
                .unwrap_or_default();
            Rng::seed_from_u64(nanos ^ std::process::id() as u64).next_u64()
        })
    }
}

/// Options for how long to search and which results to report
#[derive(Debug, Args)]
struct ReportArgs {
    /// Report every structure within this percentage of the best
    #[arg(long, default_value = "2")]
    within: f64,
    /// Report at most this many structures
    #[arg(long)]
    max_results: Option<usize>,
    /// List every candidate near the best, rather than grouping those which found the same
    /// crystal into distinct structures
    #[arg(long)]
    show_equivalent: bool,
    /// Finish after this many rounds without the best result improving
    #[arg(long, default_value = "2")]
    patience: usize,
    /// The maximum number of rounds
    #[arg(long, default_value = "8")]
    max_rounds: usize,
    /// The directory to write the results to
    #[arg(long, short, default_value = "packing-results")]
    output: PathBuf,
}

/// Options shared by the commands searching over plane groups
#[derive(Debug, Args)]
struct CandidateArgs {
    /// Comma separated plane groups to search, defaulting to all 17
    #[arg(long, value_delimiter = ',')]
    groups: Vec<String>,
    /// The maximum number of occupied Wyckoff sites
    #[arg(long, default_value = "1")]
    max_sites: usize,
    /// The maximum number of shapes in the unit cell
    #[arg(long, default_value = "12")]
    max_shapes: usize,
    #[command(flatten)]
    report: ReportArgs,
}

#[derive(Debug, Args)]
struct UnconstrainedArgs {
    #[command(flatten)]
    shape: ShapeArgs,
    #[command(flatten)]
    monte_carlo: MonteCarloArgs,
    #[command(flatten)]
    potential: PotentialArgs,
    /// Comma separated numbers of shapes in the periodic cell
    #[arg(long, value_delimiter = ',', default_value = "1,2,3,4")]
    shapes: Vec<usize>,
    #[command(flatten)]
    report: ReportArgs,
}

#[derive(Debug, Args)]
struct SearchArgs {
    #[command(flatten)]
    shape: ShapeArgs,
    #[command(flatten)]
    monte_carlo: MonteCarloArgs,
    #[command(flatten)]
    potential: PotentialArgs,
    #[command(flatten)]
    candidates: CandidateArgs,
    /// The factor by which the candidates are reduced, and steps increased, each round
    #[arg(long, default_value = "3")]
    eta: usize,
}

#[derive(Debug, Args)]
struct SurveyArgs {
    #[command(flatten)]
    shape: ShapeArgs,
    #[command(flatten)]
    monte_carlo: MonteCarloArgs,
    #[command(flatten)]
    potential: PotentialArgs,
    #[command(flatten)]
    candidates: CandidateArgs,
}

#[derive(Debug, Args)]
struct OptimiseArgs {
    #[command(flatten)]
    shape: ShapeArgs,
    #[command(flatten)]
    monte_carlo: MonteCarloArgs,
    #[command(flatten)]
    potential: PotentialArgs,
    /// The plane group, e.g. p2gg
    #[arg(long)]
    group: String,
    /// The occupied Wyckoff sites, e.g. "e" or "ab"
    #[arg(long)]
    sites: String,
    /// A structure to start the optimisation from
    #[arg(long)]
    start: Option<PathBuf>,
    /// The directory to write the results to
    #[arg(long, short, default_value = "packing-results")]
    output: PathBuf,
}

#[derive(Debug, Args)]
struct RenderArgs {
    #[command(flatten)]
    shape: ShapeArgs,
    /// The structure to render, as saved by the search or optimise commands
    input: PathBuf,
    /// The SVG file to write, defaulting to the input with an svg extension
    #[arg(long, short)]
    output: Option<PathBuf>,
    /// The number of unit cells to show either side of the central cell
    #[arg(long, default_value = "1")]
    cells: i32,
}

fn progress_bar(len: usize, message: String) -> ProgressBar {
    let bar = ProgressBar::new(len as u64);
    bar.set_style(
        ProgressStyle::with_template("{msg:30} [{bar:40}] {pos}/{len} {elapsed_precise}")
            .expect("Valid template")
            .progress_chars("=> "),
    );
    bar.set_message(message);
    bar
}

fn run_jobs(jobs: &[Job], message: String) -> Vec<crystal_packing::search::Outcome> {
    let bar = progress_bar(jobs.len(), message);
    let outcomes = jobs
        .par_iter()
        .map(|job| {
            let outcome = run_job(job);
            bar.inc(1);
            outcome
        })
        .collect();
    bar.finish();
    outcomes
}

fn write_structure(
    directory: &Path,
    name: &str,
    shape: &ShapeSpec,
    structure: &Structure,
) -> Result<(), Error> {
    let json = serde_json::to_string_pretty(&serde_json::json!({
        "shape": shape,
        "structure": structure,
    }))?;
    fs::write(directory.join(format!("{}.json", name)), json)?;
    svg::save(
        directory.join(format!("{}.svg", name)),
        &render_structure(shape, structure, 1)?,
    )?;
    Ok(())
}

fn search_config(
    monte_carlo: &MonteCarloArgs,
    potential: &PotentialArgs,
    report: &ReportArgs,
    strategy: Strategy,
    eta: usize,
) -> SearchConfig {
    SearchConfig {
        replicas: monte_carlo.replicas,
        steps: monte_carlo.steps,
        eta,
        strategy,
        patience: report.patience,
        max_rounds: report.max_rounds,
        within: report.within / 100.,
        seed: monte_carlo.seed(),
        monte_carlo: monte_carlo.config(),
        potential: potential.config(),
        ..Default::default()
    }
}

fn group_config(
    monte_carlo: &MonteCarloArgs,
    potential: &PotentialArgs,
    candidates: &CandidateArgs,
    strategy: Strategy,
    eta: usize,
) -> SearchConfig {
    SearchConfig {
        groups: candidates.groups.clone(),
        max_sites: candidates.max_sites,
        max_shapes: candidates.max_shapes,
        ..search_config(monte_carlo, potential, &candidates.report, strategy, eta)
    }
}

/// Run a search to completion, using every core
fn run_search(shape: &ShapeSpec, config: SearchConfig) -> Result<(Search, Duration), Error> {
    let seed = config.seed;
    let mut search = Search::new(shape.clone(), config)?;
    info!(
        "Searching {} candidates with seed {}",
        search.results.len(),
        seed
    );
    let start = Instant::now();
    while !search.is_finished() {
        let jobs = search.jobs();
        let message = format!(
            "Round {} ({:?}): {} candidates",
            search.round() + 1,
            search.phase(),
            search.num_active()
        );
        let outcomes = run_jobs(&jobs, message);
        search.submit(outcomes);
    }
    Ok((search, start.elapsed()))
}

fn energy_text(structure: Option<&Structure>) -> String {
    structure
        .and_then(|s| s.energy)
        .map_or(String::from("-"), |e| format!("{:.5}", e))
}

/// The difference of a score from the best, as a percentage
fn percent_below(score: f64, best: f64) -> f64 {
    100. * (best - score) / best.abs()
}

/// Print and save the structures within the tolerance of the best, best first.
///
/// By default candidates which found the same crystal are grouped into a single distinct
/// structure, listing the equivalent candidates.
fn report_near_best(search: &Search, shape: &ShapeSpec, report: &ReportArgs) -> Result<(), Error> {
    let within = report.within / 100.;
    let all = search.near_best(within);
    let distinct = search.distinct_near_best(within, None)?;
    // Each row is a structure with the candidates equivalent to it
    let rows: Vec<(&CandidateResult, Vec<&CandidateResult>)> = if report.show_equivalent {
        all.iter().map(|r| (*r, vec![])).collect()
    } else {
        distinct
            .iter()
            .map(|d| (d.result, d.equivalent.clone()))
            .collect()
    };
    let total = rows.len();
    let rows: Vec<_> = rows
        .into_iter()
        .take(report.max_results.unwrap_or(usize::MAX))
        .collect();
    let best = search.best().map_or(0., |s| s.score());
    let limited = if rows.len() < total {
        format!(", showing the best {}", rows.len())
    } else {
        String::new()
    };
    println!(
        "\n{} results within {}% of the best, {} distinct structure{}{}:",
        all.len(),
        report.within,
        distinct.len(),
        if distinct.len() == 1 { "" } else { "s" },
        limited
    );
    println!(
        "{:>4}  {:<16} {:>9} {:>10} {:>8} {:>7}  {}",
        "rank", "candidate", "phi", "energy", "below %", "shapes", "also found in"
    );
    let directory = report.output.join("near-best");
    if directory.exists() {
        fs::remove_dir_all(&directory)?;
    }
    fs::create_dir_all(&directory)?;
    let mut records = vec![];
    for (rank, (result, equivalent)) in rows.iter().enumerate() {
        let structure = result
            .best
            .as_ref()
            .expect("Near best results have a structure");
        let equivalent: Vec<String> = equivalent.iter().map(|r| r.candidate.label()).collect();
        println!(
            "{:>4}  {:<16} {:>9.5} {:>10} {:>8.3} {:>7}  {}",
            rank + 1,
            result.candidate.label(),
            structure.packing_fraction,
            energy_text(Some(structure)),
            percent_below(structure.score(), best),
            structure.num_shapes,
            equivalent.join(" ")
        );
        let name = format!(
            "{:02}-{}",
            rank + 1,
            result.candidate.label().replace(':', "-")
        );
        write_structure(&directory, &name, shape, structure)?;
        records.push(serde_json::json!({
            "candidate": result.candidate.label(),
            "structure": structure,
            "equivalent": equivalent,
        }));
    }
    println!();
    fs::write(
        report.output.join("near-best.json"),
        serde_json::to_string_pretty(&records)?,
    )?;
    Ok(())
}

fn print_best(best: &Structure, elapsed: Duration) {
    match best.energy {
        Some(energy) => println!(
            "\nLowest energy {:.6} (packing fraction {:.6}) in {} ({}) found in {:.1?}",
            energy,
            best.packing_fraction,
            best.group,
            best.sites_label(),
            elapsed
        ),
        None => println!(
            "\nBest packing fraction {:.6} in {} ({}) found in {:.1?}",
            best.packing_fraction,
            best.group,
            best.sites_label(),
            elapsed
        ),
    }
}

/// Write every candidate's best result as CSV
fn results_csv(search: &Search) -> String {
    let mut csv = String::from(
        "rank,group,sites,packing_fraction,energy,num_shapes,rounds,steps,a,b,angle\n",
    );
    for (rank, result) in search.ranking().iter().enumerate() {
        let best = result.best.as_ref();
        let (phi, shapes) = best.map_or((0., 0), |s| (s.packing_fraction, s.num_shapes));
        let (a, b, angle) = best.map_or((0., 0., 0.), |s| {
            (
                s.lattice.a().norm(),
                s.lattice.b().norm(),
                s.lattice.angle().to_degrees(),
            )
        });
        csv.push_str(&format!(
            "{},{},{},{:.8},{},{},{},{},{:.8},{:.8},{:.6}\n",
            rank + 1,
            result.candidate.group,
            result.candidate.sites.iter().collect::<String>(),
            phi,
            best.and_then(|s| s.energy)
                .map_or(String::new(), |e| format!("{:.8}", e)),
            shapes,
            result.rounds,
            result.steps,
            a,
            b,
            angle
        ));
    }
    csv
}

/// Write the files common to every search, and report the best structures
fn finish_report(
    search: &Search,
    shape: &ShapeSpec,
    report: &ReportArgs,
    elapsed: Duration,
) -> Result<(), Error> {
    let output = &report.output;
    fs::write(output.join("results.csv"), results_csv(search))?;
    fs::write(
        output.join("results.json"),
        serde_json::to_string_pretty(search)?,
    )?;
    report_near_best(search, shape, report)?;
    if let Some(best) = search.best() {
        write_structure(output, "best", shape, best)?;
        print_best(best, elapsed);
    }
    println!("Results written to {}", output.display());
    Ok(())
}

fn search(args: SearchArgs) -> Result<(), Error> {
    let shape = args.shape.spec()?;
    let config = group_config(
        &args.monte_carlo,
        &args.potential,
        &args.candidates,
        Strategy::Halving,
        args.eta,
    );
    let (search, elapsed) = run_search(&shape, config)?;
    fs::create_dir_all(&args.candidates.report.output)?;
    finish_report(&search, &shape, &args.candidates.report, elapsed)
}

fn unconstrained(args: UnconstrainedArgs) -> Result<(), Error> {
    let shape = args.shape.spec()?;
    let config = SearchConfig {
        unconstrained: args.shapes.clone(),
        ..search_config(
            &args.monte_carlo,
            &args.potential,
            &args.report,
            Strategy::Exhaustive,
            3,
        )
    };
    let (search, elapsed) = run_search(&shape, config)?;
    let output = &args.report.output;
    fs::create_dir_all(output)?;
    let best = search.best().map_or(0., |s| s.score());
    println!(
        "{:>6} {:>9} {:>10} {:>8} {:>9} {:>9} {:>8}",
        "shapes", "phi", "energy", "below %", "a", "b", "angle"
    );
    let mut csv = String::from("shapes,packing_fraction,energy,percent_below_best,a,b,angle\n");
    for result in search.results.iter() {
        let n = result.candidate.sites.len();
        let Some(structure) = result.best.as_ref() else {
            println!("{:>6} {:>9}", n, "failed");
            continue;
        };
        let (a, b, angle) = (
            structure.lattice.a().norm(),
            structure.lattice.b().norm(),
            structure.lattice.reduced_angle_degrees(),
        );
        let below = percent_below(structure.score(), best);
        println!(
            "{:>6} {:>9.5} {:>10} {:>8.3} {:>9.4} {:>9.4} {:>8.2}",
            n,
            structure.packing_fraction,
            energy_text(Some(structure)),
            below,
            a,
            b,
            angle
        );
        csv.push_str(&format!(
            "{},{:.8},{},{:.6},{:.8},{:.8},{:.6}\n",
            n,
            structure.packing_fraction,
            structure
                .energy
                .map_or(String::new(), |e| format!("{:.8}", e)),
            below,
            a,
            b,
            angle
        ));
        write_structure(output, &format!("shapes-{}", n), &shape, structure)?;
    }
    fs::write(output.join("unconstrained.csv"), csv)?;
    finish_report(&search, &shape, &args.report, elapsed)
}

fn survey(args: SurveyArgs) -> Result<(), Error> {
    let shape = args.shape.spec()?;
    let config = group_config(
        &args.monte_carlo,
        &args.potential,
        &args.candidates,
        Strategy::Exhaustive,
        3,
    );
    let (search, elapsed) = run_search(&shape, config)?;
    let output = &args.candidates.report.output;
    let directory = output.join("groups");
    fs::create_dir_all(&directory)?;
    let best = search.best().map_or(0., |s| s.score());

    println!(
        "{:<6} {:<8} {:>9} {:>10} {:>8} {:>7} {:>11}",
        "group", "sites", "phi", "energy", "below %", "shapes", "candidates"
    );
    let mut csv = String::from(
        "group,sites,packing_fraction,energy,percent_below_best,num_shapes,candidates\n",
    );
    for (group, result) in search.best_per_group() {
        let candidates = search
            .results
            .iter()
            .filter(|r| r.candidate.group == group)
            .count();
        match result.and_then(|r| r.best.as_ref().map(|s| (r, s))) {
            Some((result, structure)) => {
                let sites: String = result.candidate.sites.iter().collect();
                let below = percent_below(structure.score(), best);
                println!(
                    "{:<6} {:<8} {:>9.5} {:>10} {:>8.3} {:>7} {:>11}",
                    group,
                    sites,
                    structure.packing_fraction,
                    energy_text(Some(structure)),
                    below,
                    structure.num_shapes,
                    candidates
                );
                csv.push_str(&format!(
                    "{},{},{:.8},{},{:.6},{},{}\n",
                    group,
                    sites,
                    structure.packing_fraction,
                    structure
                        .energy
                        .map_or(String::new(), |e| format!("{:.8}", e)),
                    below,
                    structure.num_shapes,
                    candidates
                ));
                write_structure(&directory, &group, &shape, structure)?;
            }
            None => {
                println!(
                    "{:<6} {:<8} {:>9} {:>10} {:>8} {:>7} {:>11}",
                    group, "-", "-", "-", "-", "-", candidates
                );
                csv.push_str(&format!("{},,,,,,{}\n", group, candidates));
            }
        }
    }
    fs::write(output.join("survey.csv"), csv)?;
    fs::write(
        output.join("survey.json"),
        serde_json::to_string_pretty(&search)?,
    )?;
    finish_report(&search, &shape, &args.candidates.report, elapsed)
}

fn read_structure(path: &Path) -> Result<Structure, Error> {
    let value: serde_json::Value = serde_json::from_str(&fs::read_to_string(path)?)?;
    // Accept either a bare structure, or the files written by this program
    let structure = value.get("structure").cloned().unwrap_or(value);
    Ok(serde_json::from_value(structure)?)
}

fn optimise(args: OptimiseArgs) -> Result<(), Error> {
    let shape = args.shape.spec()?;
    let candidate = Candidate {
        group: args.group.clone(),
        sites: args.sites.chars().collect(),
    };
    let start = args.start.as_deref().map(read_structure).transpose()?;
    let seed = args.monte_carlo.seed();
    let jobs: Vec<Job> = (0..args.monte_carlo.replicas)
        .map(|replica| Job {
            shape: shape.clone(),
            candidate: candidate.clone(),
            seed: Rng::derive_seed(seed, &[replica as u64]),
            config: args.monte_carlo.config(),
            start: start.clone(),
            potential: args.potential.config(),
        })
        .collect();
    let outcomes = run_jobs(&jobs, candidate.label());
    for outcome in outcomes.iter() {
        if let Some(error) = &outcome.error {
            bail!("Optimisation failed: {}", error);
        }
    }
    let best = outcomes
        .into_iter()
        .filter_map(|o| o.structure)
        .max_by(|a, b| a.score().total_cmp(&b.score()))
        .context("No structure was found")?;
    fs::create_dir_all(&args.output)?;
    write_structure(&args.output, "best", &shape, &best)?;
    match best.energy {
        Some(energy) => println!(
            "Lowest energy {:.6}, packing fraction {:.6}, written to {}",
            energy,
            best.packing_fraction,
            args.output.display()
        ),
        None => println!(
            "Best packing fraction {:.6}, written to {}",
            best.packing_fraction,
            args.output.display()
        ),
    }
    Ok(())
}

fn render(args: RenderArgs) -> Result<(), Error> {
    let value: serde_json::Value = serde_json::from_str(&fs::read_to_string(&args.input)?)?;
    // Use the shape stored alongside the structure when available
    let shape = match value.get("shape") {
        Some(shape) => serde_json::from_value(shape.clone())?,
        None => args.shape.spec()?,
    };
    let structure = read_structure(&args.input)?;
    let output = args
        .output
        .unwrap_or_else(|| args.input.with_extension("svg"));
    svg::save(&output, &render_structure(&shape, &structure, args.cells)?)?;
    println!("Written {}", output.display());
    Ok(())
}

fn main() -> Result<(), Error> {
    let cli = Cli::parse();
    let level = match cli.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    env_logger::Builder::new().filter_level(level).init();
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }
    match cli.command {
        Command::Search(args) => search(args),
        Command::Survey(args) => survey(args),
        Command::Unconstrained(args) => unconstrained(args),
        Command::Optimise(args) => optimise(args),
        Command::Render(args) => render(args),
    }
}
