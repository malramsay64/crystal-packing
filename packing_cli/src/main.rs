//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use std::fs::File;
use std::io::prelude::*;
use std::path;
use std::path::PathBuf;

use anyhow::{anyhow, bail, Error};
use log::{debug, info, LevelFilter};
use rayon::prelude::*;
use clap::arg_enum;
use structopt::StructOpt;
use rand::prelude::*;

use packing::traits::*;
use packing::wallpaper::{get_wallpaper_group, WallpaperGroups};
use packing::{
    MCOptimiser, LJShape2, LineShape, MolecularShape2, PackedState2, PotentialState2,
};

#[derive(StructOpt, Debug, Clone, Copy)]
pub struct BuildOptimiser {
    /// The number of steps to run the Monte-Carlo Optimisation.
    #[structopt(short, long, default_value = "100")]
    steps: u64,

    /// The initial value of the "temperature" for the simulation. This is an indicator of how bad
    /// a move can be yet still be accepted.
    #[structopt(long, default_value = "0.1")]
    kt_start: f64,

    /// The initial value of the "temperature" for the simulation. This is an indicator of how bad
    /// a move can be yet still be accepted.
    #[structopt(long)]
    kt_finish: Option<f64>,

    /// The ratio the temperature is reduced every inner_steps
    #[structopt(long)]
    kt_ratio: Option<f64>,

    /// The maximum size of a Monte-Carlo move
    #[structopt(long, default_value = "0.01")]
    max_step_size: f64,

    /// The number of steps to run before reducing the temperature. This also defines the number of
    /// steps before the step size is updated.
    #[structopt(long, default_value = "1000")]
    inner_steps: u64,

    /// This option is skipped on the command line and filled in when setting up the iterations.
    #[structopt(skip)]
    seed: Option<u64>,

    /// The minimum change of the score within an inner loop. This is the indicator of convergence
    /// which allows for an early exit.
    #[structopt(long)]
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
            None => rand_pcg::Pcg64Mcg::from_entropy().gen(),
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

arg_enum! {
    #[derive(Debug)]
    enum Force {
        LJ,
        Hard,
    }
}

#[derive(Debug, StructOpt)]
#[structopt(name = "packing")]
struct Args {
    /// Pass many times for more log output
    ///
    /// By default, it'll only report errors. Passing `-v` one time also prints
    /// warnings, `-vv` enables info logging, `-vvv` debug, and `-vvvv` trace.
    #[structopt(long, short, parse(from_occurrences))]
    verbosity: u8,

    #[structopt(subcommand)]
    shape: Shapes,

    /// The defining symmetry of the unit cell
    #[structopt(possible_values = &WallpaperGroups::variants())]
    wallpaper: WallpaperGroups,

    /// The potential which is being optimised
    #[structopt(short, long, possible_values = &Force::variants(), default_value = "Hard")]
    potential: Force,

    /// Where to save the best packed structure
    #[structopt(long, parse(from_os_str))]
    outfile: PathBuf,

    /// An initial configuration which is the starting point for optimisation
    #[structopt(long, parse(from_os_str))]
    start_config: Option<PathBuf>,

    /// The number of independent starting configurations to optimise
    #[structopt(long, default_value = "100")]
    replications: u64,

    #[structopt(flatten)]
    optimisation: BuildOptimiser,
}

#[derive(Debug, StructOpt)]
enum Shapes {
    #[structopt(name = "polygon")]
    Polygon {
        /// The number of equally spaced sides
        #[structopt(long, default_value = "4")]
        sides: usize,
    },
    #[structopt(name = "trimer")]
    Trimer {
        /// The distance from the central particle to the outer particles
        #[structopt(short, long, default_value = "1.")]
        distance: f64,
        /// The angle between the two outer particles in degrees
        #[structopt(short, long, default_value = "120")]
        angle: f64,
        /// The radius of the outer small particles
        #[structopt(short, long, default_value = "0.637556")]
        radius: f64,
    },
    #[structopt(name = "circle")]
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

    File::create(outfile.with_extension("json"))?.write_all(&serialised.as_bytes())?;
    svg::save(outfile.with_extension("svg"), &final_state.as_svg())?;

    Ok(())
}

#[paw::main]
fn main(args: Args) -> Result<(), Error> {
    let log_level = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        2 => LevelFilter::Trace,
        _ => LevelFilter::Trace,
    };
    env_logger::Builder::new().filter_level(log_level).init();

    debug!("Logging Level: {}", log_level);

    let wg = get_wallpaper_group(args.wallpaper)?;

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
