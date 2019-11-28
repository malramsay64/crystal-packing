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
use serde_json;
use structopt::clap::arg_enum;
use structopt::StructOpt;

use packing::traits::*;
use packing::wallpaper::{get_wallpaper_group, WallpaperGroups};
use packing::{
    BuildOptimiser, LJShape2, LineShape, MolecularShape2, PackedState2, PotentialState2,
};

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
                .steps(1000)
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

    let final_score = final_state
        .score()
        .ok_or_else(|| anyhow!("State has become corrupted"))?;
    info!("Final score: {}", final_score);

    let serialised = serde_json::to_string(&final_state)?;

    File::create(outfile.clone().with_extension("json"))?.write_all(&serialised.as_bytes())?;
    svg::save(outfile.clone().with_extension("svg"), &final_state.as_svg())?;

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
