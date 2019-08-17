//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

// This is for the clap::arg_enum macro
// #![allow(deprecated)]

use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;
use std::{fs, path};

use clap::{value_t, value_t_or_exit, App, Arg};
use log::{debug, info};
use rayon::prelude::*;
use serde_json;
use simplelog::{Config, LevelFilter, TermLogger, TerminalMode};

use packing::traits::*;
use packing::wallpaper::{WallpaperGroup, WallpaperGroups};
use packing::{
    BuildOptimiser, LJShape2, LineShape, MolecularShape2, PackedState2, PotentialState2,
};

struct CLIOptions {
    state: StateTypes,
    steps: u64,
    num_start_configs: u64,
    log_level: LevelFilter,
    outfile: path::PathBuf,
    max_step_size: f64,
    kt: Option<f64>,
    kt_ratio: Option<f64>,
}

enum StateTypes {
    Polygon(PackedState2<LineShape>),
    Trimer(PackedState2<MolecularShape2>),
    Circle(PackedState2<MolecularShape2>),
    LJTrimer(PotentialState2<LJShape2>),
}

fn cli() -> CLIOptions {
    let matches = App::new("packing")
        .version("0.3.0")
        .author("Malcolm Ramsay <malramsay64@gmail.com")
        .about("Find best tilings of 2d shapes")
        .arg(
            Arg::with_name("verbosity")
                .short("v")
                .multiple(true)
                .long("verbose")
                .help("Output more debugging information to console.")
        )
        .arg(
            Arg::with_name("quiet")
                .short("q")
                .multiple(true)
                .long("quiet")
                .help("Suppress output to terminal.")
        )
        .arg(Arg::with_name("shape")
             .possible_values(&["trimer", "polygon", "circle"])
             .required(true)
             .help("Shape used for the optimisation.")
             )
        .arg(
            Arg::with_name("wallpaper")
                .possible_values(&WallpaperGroups::variants())
                .case_insensitive(true)
                .required(true)
                .help("The wallpaper group defining the symmetry of the unit cell.")
        )
        .arg(
            Arg::with_name("lj")
            .long("lj")
            .help("Pass to use the lennard jones potential instead of hard discs.")
        )
        .arg(
            Arg::with_name("radius")
            .long("radius")
            .help("The radius of the smaller particles in the trimer molecule.")
            .default_value("0.637556")
            .takes_value(true),
            )
        .arg(
            Arg::with_name("distance")
            .long("distance")
            .default_value("1")
            .help("The distance from the center of the large particle to the center of the small particles in the trimer molecule.")
            .takes_value(true),
            )
        .arg(
            Arg::with_name("angle")
            .long("angle")
            .default_value("120")
            .help("The angle between the two small particles in the trimer molecule.")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("sides")
            .long("--num-sides")
            .takes_value(true)
            .default_value("4")
            .help("The number of sides which should be used with the polygon shape."),
            )
        .arg(Arg::with_name("read_config")
             .long("read-config")
             .takes_value(true)
             .help("Read the configuration from a file. The shape has to be the same as in the configuration.")
        )
        .arg(Arg::with_name("outfile")
             .long("outfile")
             .takes_value(true)
             .default_value("structure.json")
             .help("Filename to output best packes structure to.")
             )
        .arg(
            Arg::with_name("steps")
                .short("-s")
                .long("--steps")
                .takes_value(true)
                .default_value("100")
                .help("The number of steps to use for the optimisation.")
        )
        .arg(
            Arg::with_name("kt")
                .long("--kt")
                .takes_value(true)
                .help("The initial temperature which guides the acceptance criterita.")
        )
        .arg(
            Arg::with_name("kt_ratio")
                .long("--ratio")
                .default_value("0.1")
                .takes_value(true)
                .help("The multiplicative factor which reduces the temperature.")
        )
        .arg(
            Arg::with_name("max_step_size")
                .long("--max-step")
                .takes_value(true)
                .default_value("1.0")
                .help("The number of steps to use for the optimisation.")
        )
        .arg(
            Arg::with_name("replications")
                .long("--replications")
                .takes_value(true)
                .default_value("32")
                .help("The number of replications of the optimisation simulation to run."),
        )
        .get_matches();

    println!("Starting Analysis");
    let state = if let Some(fname) = matches.value_of("read_config") {
        println!("Reading config from: {}", fname);
        let serialised = fs::read_to_string(fname).unwrap();
        match (matches.value_of("shape"), matches.is_present("lj")) {
            (Some("trimer"), false) => {
                StateTypes::Trimer(serde_json::from_str(&serialised).unwrap())
            }
            (Some("trimer"), true) => {
                StateTypes::LJTrimer(serde_json::from_str(&serialised).unwrap())
            }
            (Some("polygon"), _) => StateTypes::Polygon(serde_json::from_str(&serialised).unwrap()),
            _ => panic!(),
        }
    } else {
        let wg = value_t_or_exit!(matches.value_of("wallpaper"), WallpaperGroups);
        println!("Wallpaper Group: {}", wg);
        let group: WallpaperGroup = packing::wallpaper::get_wallpaper_group(wg).unwrap();
        match matches.value_of("shape") {
            Some("circle") => {
                StateTypes::Circle(PackedState2::from_group(MolecularShape2::circle(), &group))
            }
            Some("trimer") => {
                let radius: f64 = value_t_or_exit!(matches, "radius", f64);
                let distance: f64 = value_t_or_exit!(matches, "distance", f64);
                let angle: f64 = value_t_or_exit!(matches, "angle", f64);
                let angle = angle * PI / 180.;

                if matches.is_present("lj") {
                    StateTypes::LJTrimer(PotentialState2::from_group(
                        LJShape2::from_trimer(radius, angle, distance),
                        &group,
                    ))
                } else {
                    StateTypes::Trimer(PackedState2::from_group(
                        MolecularShape2::from_trimer(radius, angle, distance),
                        &group,
                    ))
                }
            }
            Some("polygon") => {
                let sides: usize = value_t_or_exit!(matches, "sides", usize);
                StateTypes::Polygon(PackedState2::from_group(
                    LineShape::from_radial("Polygon", vec![1.; sides]).unwrap(),
                    &group,
                ))
            }
            _ => panic!(),
        }
    };

    let steps: u64 = value_t_or_exit!(matches, "steps", u64);
    let max_step_size: f64 = value_t_or_exit!(matches, "max_step_size", f64);
    let kt_ratio: Option<f64> = value_t!(matches, "kt_ratio", f64).ok();
    let kt: Option<f64> = value_t!(matches, "kt", f64).ok();
    let num_start_configs: u64 = value_t_or_exit!(matches, "replications", u64);
    let outfile = value_t_or_exit!(matches, "outfile", path::PathBuf);

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
        state,
        steps,
        num_start_configs,
        log_level,
        outfile,
        max_step_size,
        kt,
        kt_ratio,
    }
}

fn analyse_state(
    outfile: path::PathBuf,
    start_configs: u64,
    state: impl State,
    optimiser: &BuildOptimiser,
) -> Result<(), &'static str> {
    let final_state = (0..start_configs)
        .into_par_iter()
        // Create collection of quickly optimised initial states
        .map(|index| {
            let result = optimiser
                .clone()
                .steps(1000)
                .kt_start(0.)
                .seed(index)
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
            let result = optimiser
                .clone()
                .kt_start(0.)
                .seed(index)
                .build()
                .optimise_state(opt_state);
            result
        })
        .max();

    if let Some(s) = final_state {
        if let Ok(score) = s.score() {
            info!("Final score: {}", score);

            let serialised = serde_json::to_string(&s).unwrap();

            let mut state_path = outfile.clone();
            state_path.set_extension("json");
            let mut state_file = File::create(state_path).unwrap();
            state_file.write_all(&serialised.as_bytes()).unwrap();
            svg::save("test.svg", &s.as_svg()).unwrap();
            return Ok(());
        }
    }
    Err("Final state contains invalid configuration.")
}

fn main() -> Result<(), &'static str> {
    let options = cli();
    TermLogger::init(options.log_level, Config::default(), TerminalMode::Mixed).unwrap();

    debug!("Logging Level: {}", options.log_level);

    let mut opt = *BuildOptimiser::default()
        .steps(options.steps)
        .max_step_size(options.max_step_size)
        .kt_ratio(options.kt_ratio);
    if let Some(kt) = options.kt {
        opt.kt_start(kt);
    }

    match options.state {
        StateTypes::Circle(state) | StateTypes::Trimer(state) => {
            analyse_state(options.outfile, options.num_start_configs, state, &opt)?
        }
        StateTypes::Polygon(state) => {
            analyse_state(options.outfile, options.num_start_configs, state, &opt)?
        }
        StateTypes::LJTrimer(state) => {
            analyse_state(options.outfile, options.num_start_configs, state, &opt)?
        }
    }

    Ok(())
}
