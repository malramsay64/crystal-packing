//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

// This is for the clap::arg_enum macro
#![allow(deprecated)]

use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;

use clap::{arg_enum, App, Arg, _clap_count_exprs, value_t};
use log::{debug, info};
use rayon::prelude::*;
use serde_json;
use simplelog::{Config, LevelFilter, TermLogger};

use packing::traits::*;
use packing::wallpaper::{WallpaperGroup, WallpaperGroups};
use packing::{
    BuildOptimiser, LJShape2, LineShape, MolecularShape2, PackedState2, PotentialState2,
};

struct CLIOptions {
    shape: ShapeTypes,
    group: WallpaperGroup,
    steps: u64,
    num_sides: usize,
    num_start_configs: u64,
    log_level: LevelFilter,
}

arg_enum! {
    #[derive(Debug)]
    pub enum ShapeTypes {
        Polygon,
        Trimer,
        Circle,
        LJTrimer,
    }
}

fn parallel_opt(
    steps: u64,
    start_configs: u64,
    state: impl State,
) -> Result<impl State, &'static str> {
    let opt = *BuildOptimiser::default().steps(steps);
    let final_state = (0..start_configs)
        .into_par_iter()
        .map(|_| opt.build().optimise_state(state.clone()))
        .max()
        .unwrap()?;
    info!("Final score: {}", final_state.score().unwrap());
    Ok(final_state)
}

fn cli() -> CLIOptions {
    let matches = App::new("packing")
        .version("0.1.0")
        .author("Malcolm Ramsay <malramsay64@gmail.com")
        .about("Find best tilings of 2d shapes")
        .arg(
            Arg::with_name("verbosity")
                .short("v")
                .multiple(true)
                .long("verbose"),
        )
        .arg(
            Arg::with_name("quiet")
                .short("q")
                .multiple(true)
                .long("quiet"),
        )
        .arg(
            Arg::with_name("wallpaper_group")
                .possible_values(&WallpaperGroups::variants())
                .required(true),
        )
        .arg(
            Arg::with_name("shape")
                .possible_values(&ShapeTypes::variants())
                .required(true),
        )
        .arg(
            Arg::with_name("sides")
                .long("--num-sides")
                .takes_value(true)
                .default_value("4"),
        )
        .arg(
            Arg::with_name("steps")
                .short("-s")
                .long("--steps")
                .takes_value(true)
                .default_value("100"),
        )
        .arg(
            Arg::with_name("replications")
                .long("--replications")
                .takes_value(true)
                .default_value("32"),
        )
        .get_matches();

    let wg = value_t!(matches.value_of("wallpaper_group"), WallpaperGroups).unwrap();
    info!("Using Wallpaper Group: {}", wg);

    let shape = value_t!(matches.value_of("shape"), ShapeTypes).unwrap();
    let num_sides = matches.value_of("sides").unwrap().parse().unwrap();

    let group = packing::wallpaper::get_wallpaper_group(wg).unwrap();

    let steps: u64 = matches.value_of("steps").unwrap().parse().unwrap();
    let num_start_configs: u64 = matches.value_of("replications").unwrap().parse().unwrap();

    let log_level =
        match matches.occurrences_of("verbosity") as i64 - matches.occurrences_of("quiet") as i64 {
            x if x <= -2 => LevelFilter::Error,
            -1 => LevelFilter::Warn,
            0 => LevelFilter::Info,
            1 => LevelFilter::Debug,
            x if 2 <= x => LevelFilter::Trace,
            _ => unreachable!(),
        };

    CLIOptions {
        shape,
        group,
        steps,
        num_sides,
        num_start_configs,
        log_level,
    }
}

fn main() -> Result<(), &'static str> {
    let options = cli();

    TermLogger::init(options.log_level, Config::default()).unwrap();
    debug!("Logging Level: {}", options.log_level);

    let mut file = File::create("structure.json").unwrap();

    match options.shape {
        ShapeTypes::Polygon => {
            let shape = LineShape::from_radial("Polygon", vec![1.; options.num_sides])?;
            let state = parallel_opt(
                options.steps,
                options.num_start_configs,
                PackedState2::from_group(shape, &options.group),
            );
            let serialised = serde_json::to_string(&state).unwrap();
            file.write_all(&serialised.as_bytes()).unwrap();
        }
        ShapeTypes::Trimer => {
            let shape = MolecularShape2::from_trimer(0.637_556, 2. * PI / 3., 1.);
            let state = parallel_opt(
                options.steps,
                options.num_start_configs,
                PackedState2::from_group(shape, &options.group),
            );
            let serialised = serde_json::to_string(&state).unwrap();
            file.write_all(&serialised.as_bytes()).unwrap();
        }
        ShapeTypes::Circle => {
            let shape = MolecularShape2::circle();
            let state = parallel_opt(
                options.steps,
                options.num_start_configs,
                PackedState2::from_group(shape, &options.group),
            );
            let serialised = serde_json::to_string(&state).unwrap();
            file.write_all(&serialised.as_bytes()).unwrap();
        }
        ShapeTypes::LJTrimer => {
            let shape = LJShape2::from_trimer(0.637_556, 2. * PI / 3., 1.);
            let state = parallel_opt(
                options.steps,
                options.num_start_configs,
                PotentialState2::from_group(shape, &options.group),
            );
            let serialised = serde_json::to_string(&state).unwrap();
            file.write_all(&serialised.as_bytes()).unwrap();
        }
    };
    Ok(())
}
