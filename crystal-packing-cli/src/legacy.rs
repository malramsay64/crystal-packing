//
// legacy.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use std::convert::TryInto;
use std::fs::File;
use std::io::prelude::*;
use std::path;
use std::path::PathBuf;

use anyhow::{anyhow, bail, Error};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use log::{debug, info};
use rayon::prelude::*;

use crystal_packing::traits::*;
use crystal_packing::wallpaper::{WallpaperGroup, WallpaperGroups};
use crystal_packing::{
    LJShape2, LineShape, MCOptimiser, MolecularShape2, PackedState2, PotentialState2, Rng,
};

#[derive(ClapArgs, Debug, Clone, Copy)]
pub struct BuildOptimiser {
    /// The number of steps to run the Monte-Carlo Optimisation.
    #[arg(short, long, default_value = "100")]
    steps: u64,

    /// The initial value of the "temperature" for the simulation. This is an indicator of how bad
    /// a move can be yet still be accepted.
    #[arg(long, default_value = "0.1")]
    kt_start: f64,

    /// The initial value of the "temperature" for the simulation. This is an indicator of how bad
    /// a move can be yet still be accepted.
    #[arg(long)]
    kt_finish: Option<f64>,

    /// The ratio the temperature is reduced every inner_steps
    #[arg(long)]
    kt_ratio: Option<f64>,

    /// The maximum size of a Monte-Carlo move
    #[arg(long, default_value = "0.01")]
    max_step_size: f64,

    /// The number of steps to run before reducing the temperature. This also defines the number of
    /// steps before the step size is updated.
    #[arg(long, default_value = "1000")]
    inner_steps: u64,

    /// This option is skipped on the command line and filled in when setting up the iterations.
    #[arg(skip)]
    seed: Option<u64>,

    /// The minimum change of the score within an inner loop. This is the indicator of convergence
    /// which allows for an early exit.
    #[arg(long)]
    convergence: Option<f64>,
}

impl Default for BuildOptimiser {
    fn default() -> Self {
        Self {
            kt_start: 0.1,
            kt_finish: Some(0.001),
            kt_ratio: None,
            max_step_size: 0.01,
            steps: 1000,
            inner_steps: 1000,
            seed: None,
            convergence: None,
        }
    }
}

#[allow(dead_code)]
impl BuildOptimiser {
    pub fn kt_start(&mut self, kt_start: f64) -> &mut Self {
        self.kt_start = kt_start;
        self
    }

    pub fn kt_finish(&mut self, kt_finish: f64) -> &mut Self {
        self.kt_finish = Some(kt_finish);
        self
    }

    pub fn kt_ratio(&mut self, kt_ratio: Option<f64>) -> &mut Self {
        self.kt_ratio = kt_ratio;
        self
    }

    pub fn max_step_size(&mut self, max_step_size: f64) -> &mut Self {
        self.max_step_size = max_step_size;
        self
    }

    pub fn steps(&mut self, steps: u64) -> &mut Self {
        self.steps = steps;
        self
    }

    pub fn inner_steps(&mut self, inner_steps: u64) -> &mut Self {
        self.inner_steps = inner_steps;
        self
    }

    pub fn seed(&mut self, seed: u64) -> &mut Self {
        self.seed = Some(seed);
        self
    }

    pub fn convergence(&mut self, converged: Option<f64>) -> &mut Self {
        self.convergence = converged;
        self
    }

    pub fn build(&self) -> MCOptimiser {
        let kt_ratio = match (self.kt_ratio, self.kt_finish) {
            (Some(ratio), _) => 1. - ratio,
            (None, Some(finish)) => f64::powf(finish / self.kt_start, 1. / self.steps as f64),
            (None, None) => 0.1,
        };
        debug!("Setting kt_ratio to: {}", kt_ratio);
        let seed = match self.seed {
            None => Rng::seed_from_u64(
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map(|d| d.as_nanos() as u64)
                    .unwrap_or_default(),
            )
            .next_u64(),
            Some(x) => x,
        };

        MCOptimiser::new(
            self.kt_start,
            kt_ratio,
            self.max_step_size,
            self.steps,
            u64::min(self.inner_steps, self.steps),
            seed,
            self.convergence,
        )
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Force {
    #[value(name = "LJ", alias = "lj")]
    LJ,
    #[value(name = "Hard", alias = "hard")]
    Hard,
}

/// The original optimisation, which also supports the Lennard-Jones potential
#[derive(Debug, ClapArgs)]
pub struct LegacyArgs {
    #[command(subcommand)]
    shape: Shapes,

    /// The defining symmetry of the unit cell
    #[arg(value_parser = clap::builder::PossibleValuesParser::new(WallpaperGroups::variants()))]
    wallpaper: String,

    /// The potential which is being optimised
    #[arg(short, long, value_enum, default_value = "Hard")]
    potential: Force,

    /// Where to save the best packed structure
    #[arg(long)]
    outfile: PathBuf,

    /// An initial configuration which is the starting point for optimisation
    #[arg(long)]
    start_config: Option<PathBuf>,

    /// The number of independent starting configurations to optimise
    #[arg(long, default_value = "100")]
    replications: u64,

    #[command(flatten)]
    optimisation: BuildOptimiser,
}

#[derive(Debug, Subcommand)]
enum Shapes {
    #[command(name = "polygon")]
    Polygon {
        /// The number of equally spaced sides
        #[arg(long, default_value = "4")]
        sides: usize,
    },
    #[command(name = "trimer")]
    Trimer {
        /// The distance from the central particle to the outer particles
        #[arg(short, long, default_value = "1.")]
        distance: f64,
        /// The angle between the two outer particles in degrees
        #[arg(short, long, default_value = "120")]
        angle: f64,
        /// The radius of the outer small particles
        #[arg(short, long, default_value = "0.637556")]
        radius: f64,
    },
    #[command(name = "circle")]
    Circle {},
}

fn analyse_state(
    outfile: path::PathBuf,
    start_configs: u64,
    state: impl State,
    optimiser: &BuildOptimiser,
) -> Result<(), Error> {
    let final_state = (0..start_configs)
        .into_par_iter()
        // Create collection of quickly optimised initial states
        .map(|index| {
            let result = optimiser
                .clone()
                .steps(100)
                .kt_start(0.)
                .seed(index)
                .convergence(None)
                .build()
                .optimise_state(state.clone());
            (index, result)
        })
        // Perform Monte carlo optimisation
        .map(|(index, opt_state)| {
            let result = optimiser
                .clone()
                .seed(index)
                .build()
                .optimise_state(opt_state);
            (index, result)
        })
        // Final optimsation to help find the minimum
        .map(|(index, opt_state)| {
            optimiser
                .clone()
                .kt_start(0.)
                .seed(index)
                .build()
                .optimise_state(opt_state)
        })
        .max()
        .ok_or_else(|| anyhow!("Error in running optimisation."))?;

    info!(
        "Final score: {}",
        final_state
            .score()
            .ok_or_else(|| anyhow!("State has become corrupted"))?
    );

    let serialised = serde_json::to_string(&final_state)?;

    File::create(outfile.with_extension("json"))?.write_all(serialised.as_bytes())?;
    svg::save(outfile.with_extension("svg"), &final_state.as_svg())?;

    Ok(())
}

pub fn run(args: LegacyArgs) -> Result<(), Error> {
    let wallpaper: WallpaperGroups = args.wallpaper.parse()?;
    let wg: WallpaperGroup = wallpaper.try_into()?;

    match (args.shape, args.potential) {
        (
            Shapes::Trimer {
                distance,
                angle,
                radius,
            },
            Force::LJ,
        ) => analyse_state(
            args.outfile,
            args.replications,
            PotentialState2::from_group(LJShape2::from_trimer(radius, angle, distance), &wg)?,
            &args.optimisation,
        ),
        (
            Shapes::Trimer {
                distance,
                angle,
                radius,
            },
            Force::Hard,
        ) => analyse_state(
            args.outfile,
            args.replications,
            PackedState2::from_group(MolecularShape2::from_trimer(radius, angle, distance), &wg)?,
            &args.optimisation,
        ),
        (Shapes::Circle {}, Force::LJ) => analyse_state(
            args.outfile,
            args.replications,
            PotentialState2::from_group(LJShape2::circle(), &wg)?,
            &args.optimisation,
        ),
        (Shapes::Circle {}, Force::Hard) => analyse_state(
            args.outfile,
            args.replications,
            PackedState2::from_group(MolecularShape2::circle(), &wg)?,
            &args.optimisation,
        ),
        (Shapes::Polygon { sides }, Force::Hard) => analyse_state(
            args.outfile,
            args.replications,
            PackedState2::from_group(LineShape::polygon(sides)?, &wg)?,
            &args.optimisation,
        ),
        (Shapes::Polygon { .. }, Force::LJ) => {
            bail!("Polygon with a LJ potential is not yet implemented")
        }
    }
}
