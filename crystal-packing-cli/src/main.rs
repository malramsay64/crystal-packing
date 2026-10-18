//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Error};
use clap::{Args, Parser, Subcommand, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, LevelFilter};
use rayon::prelude::*;

use crystal_packing::engine::Structure;
use crystal_packing::monte_carlo::MonteCarloConfig;
use crystal_packing::render::render_structure;
use crystal_packing::search::{run_job, Candidate, Job, Search, SearchConfig, ShapeSpec};
use crystal_packing::Rng;

mod legacy;

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
    /// Optimise a single plane group and choice of Wyckoff sites
    Optimise(OptimiseArgs),
    /// Render a saved structure as an SVG image
    Render(RenderArgs),
    /// The original optimisation, which also supports the Lennard-Jones potential
    Legacy(legacy::LegacyArgs),
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

#[derive(Debug, Args)]
struct SearchArgs {
    #[command(flatten)]
    shape: ShapeArgs,
    #[command(flatten)]
    monte_carlo: MonteCarloArgs,
    /// Comma separated plane groups to search, defaulting to all 17
    #[arg(long, value_delimiter = ',')]
    groups: Vec<String>,
    /// The maximum number of occupied Wyckoff sites
    #[arg(long, default_value = "1")]
    max_sites: usize,
    /// The maximum number of shapes in the unit cell
    #[arg(long, default_value = "12")]
    max_shapes: usize,
    /// The factor by which the candidates are reduced, and steps increased, each round
    #[arg(long, default_value = "3")]
    eta: usize,
    /// The directory to write the results to
    #[arg(long, short, default_value = "packing-results")]
    output: PathBuf,
    /// The number of the best structures to render
    #[arg(long, default_value = "5")]
    top: usize,
}

#[derive(Debug, Args)]
struct OptimiseArgs {
    #[command(flatten)]
    shape: ShapeArgs,
    #[command(flatten)]
    monte_carlo: MonteCarloArgs,
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

fn search(args: SearchArgs) -> Result<(), Error> {
    let shape = args.shape.spec()?;
    let seed = args.monte_carlo.seed();
    let config = SearchConfig {
        groups: args.groups.clone(),
        max_sites: args.max_sites,
        max_shapes: args.max_shapes,
        replicas: args.monte_carlo.replicas,
        steps: args.monte_carlo.steps,
        eta: args.eta,
        seed,
        monte_carlo: args.monte_carlo.config(),
    };
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
            "Round {}: {} candidates",
            search.round() + 1,
            search.num_active()
        );
        let outcomes = run_jobs(&jobs, message);
        search.submit(outcomes);
    }
    let elapsed = start.elapsed();

    fs::create_dir_all(&args.output)?;
    let ranking = search.ranking();
    println!(
        "{:>4}  {:<16} {:>8} {:>8} {:>7}",
        "rank", "candidate", "phi", "shapes", "rounds"
    );
    let mut csv = String::from("rank,group,sites,packing_fraction,num_shapes,rounds,a,b,angle\n");
    for (rank, result) in ranking.iter().enumerate() {
        let (phi, shapes) = result
            .best
            .as_ref()
            .map_or((0., 0), |s| (s.packing_fraction, s.num_shapes));
        if rank < 20 {
            println!(
                "{:>4}  {:<16} {:>8.5} {:>8} {:>7}",
                rank + 1,
                result.candidate.label(),
                phi,
                shapes,
                result.rounds
            );
        }
        let (a, b, angle) = result.best.as_ref().map_or((0., 0., 0.), |s| {
            (
                s.lattice.a().norm(),
                s.lattice.b().norm(),
                s.lattice.angle().to_degrees(),
            )
        });
        csv.push_str(&format!(
            "{},{},{},{:.8},{},{},{:.8},{:.8},{:.6}\n",
            rank + 1,
            result.candidate.group,
            result.candidate.sites.iter().collect::<String>(),
            phi,
            shapes,
            result.rounds,
            a,
            b,
            angle
        ));
        if rank < args.top {
            if let Some(structure) = &result.best {
                let name = format!(
                    "rank{:02}-{}",
                    rank + 1,
                    result.candidate.label().replace(':', "-")
                );
                write_structure(&args.output, &name, &shape, structure)?;
            }
        }
    }
    fs::write(args.output.join("results.csv"), csv)?;
    fs::write(
        args.output.join("results.json"),
        serde_json::to_string_pretty(&search)?,
    )?;
    if let Some(best) = search.best() {
        write_structure(&args.output, "best", &shape, best)?;
        println!(
            "\nBest packing fraction {:.6} in {} ({}) found in {:.1?}",
            best.packing_fraction,
            best.group,
            best.sites_label(),
            elapsed
        );
    }
    println!("Results written to {}", args.output.display());
    Ok(())
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
        .max_by(|a, b| a.packing_fraction.total_cmp(&b.packing_fraction))
        .context("No structure was found")?;
    fs::create_dir_all(&args.output)?;
    write_structure(&args.output, "best", &shape, &best)?;
    println!(
        "Best packing fraction {:.6}, written to {}",
        best.packing_fraction,
        args.output.display()
    );
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
        Command::Optimise(args) => optimise(args),
        Command::Render(args) => render(args),
        Command::Legacy(args) => legacy::run(args),
    }
}
