//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

// This is for the clap::arg_enum macro
#![allow(deprecated)]

use std::f64::consts::PI;
use std::fs;
use std::fs::File;
use std::io::prelude::*;

use clap::{value_t, App, Arg, SubCommand};
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
    state: StateTypes,
    steps: u64,
    num_start_configs: u64,
    log_level: LevelFilter,
}

enum StateTypes {
    Polygon(PackedState2<LineShape>),
    Trimer(PackedState2<MolecularShape2>),
    Circle(PackedState2<MolecularShape2>),
    LJTrimer(PotentialState2<LJShape2>),
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
                .case_insensitive(true)
                .required_unless("start_config"),
        )
        .subcommand(
            SubCommand::with_name("circle")
            )
        .subcommand(
            SubCommand::with_name("trimer")
                .arg(
                    Arg::with_name("interaction")
                        .long("interaction")
                        .takes_value(true)
                        .possible_values(&["lj", "packing"])
                        .required(true),
                )
                .arg(
                    Arg::with_name("radius")
                    .long("radius")
                    .help("The radius of the smaller particles.")
                    .default_value("0.637556")
                    .takes_value(true),
                    )
                .arg(
                    Arg::with_name("distance")
                    .long("distance")
                    .default_value("1")
                    .help("The distance from the center of the large particle to the center of the small particles.")
                    .takes_value(true),
                    )
                .arg(
                    Arg::with_name("angle")
                    .long("angle")
                    .default_value("120")
                    .help("The angle between the two small particles.")
                    .takes_value(true)
                    )
        )
        .subcommand(
            SubCommand::with_name("polygon").arg(
                Arg::with_name("sides")
                    .long("--num-sides")
                    .takes_value(true)
                    .default_value("4")
                    .help("The number of sides which should be used with the polygon shape."),
            ),
        )
        .subcommand(
            SubCommand::with_name("config")
            .arg(Arg::with_name("shape").possible_values(&["trimer", "ljtrimer", "polygon"]))
            .arg(Arg::with_name("filename"))
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
                .default_value("32")
                .help("The number of optimisation simulations to run."),
        )
        .get_matches();

    let wg = value_t!(matches.value_of("wallpaper_group"), WallpaperGroups).ok();
    let group: Option<WallpaperGroup> = match wg {
        Some(g) => Some(packing::wallpaper::get_wallpaper_group(g).unwrap()),
        None => None,
    };

    let state = match matches.subcommand() {
        ("circle", _) => StateTypes::Circle(PackedState2::from_group(
            MolecularShape2::circle(),
            &group.unwrap(),
        )),
        ("trimer", Some(matches)) => {
            let radius: f64 = matches
                .value_of("radius")
                .map(|x| x.parse().unwrap())
                .unwrap();
            let distance: f64 = matches
                .value_of("distance")
                .map(|x| x.parse().unwrap())
                .unwrap();
            let angle: f64 = matches
                .value_of("angle")
                .map(|x| x.parse().unwrap())
                .unwrap();
            let angle = angle * PI / 180.;

            match matches.value_of("interaction") {
                Some("packing") => StateTypes::Trimer(PackedState2::from_group(
                    MolecularShape2::from_trimer(radius, angle, distance),
                    &group.unwrap(),
                )),
                Some("lj") => StateTypes::LJTrimer(PotentialState2::from_group(
                    LJShape2::from_trimer(radius, angle, distance),
                    &group.unwrap(),
                )),
                _ => panic!(),
            }
        }
        ("polygon", Some(matches)) => {
            let sides: u64 = matches
                .value_of("num_sides")
                .map(|x| x.parse().unwrap())
                .unwrap();
            StateTypes::Polygon(PackedState2::from_group(
                LineShape::from_radial("Polygon", vec![1.; sides as usize]).unwrap(),
                &group.unwrap(),
            ))
        }
        ("config", Some(matches)) => {
            let serialised = fs::read_to_string("").unwrap();
            match matches.value_of("shape") {
                Some("trimer") => StateTypes::Trimer(serde_json::from_str(&serialised).unwrap()),
                Some("polygon") => StateTypes::Polygon(serde_json::from_str(&serialised).unwrap()),
                Some("ljtrimer") => {
                    StateTypes::LJTrimer(serde_json::from_str(&serialised).unwrap())
                }
                _ => panic!(),
            }
        }
        _ => panic!(),
    };

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
        state,
        steps,
        num_start_configs,
        log_level,
    }
}

fn main() -> Result<(), &'static str> {
    let options = cli();

    TermLogger::init(options.log_level, Config::default()).unwrap();
    debug!("Logging Level: {}", options.log_level);

    let mut file = File::create("structure.json").unwrap();

    match options.state {
        StateTypes::Circle(state) | StateTypes::Trimer(state) => {
            let s = parallel_opt(options.steps, options.num_start_configs, state.clone());
            let serialised = serde_json::to_string(&s).unwrap();
            file.write_all(&serialised.as_bytes()).unwrap();
        }
        StateTypes::Polygon(state) => {
            let s = parallel_opt(options.steps, options.num_start_configs, state.clone());
            let serialised = serde_json::to_string(&s).unwrap();
            file.write_all(&serialised.as_bytes()).unwrap();
        }
        StateTypes::LJTrimer(state) => {
            let s = parallel_opt(options.steps, options.num_start_configs, state.clone());
            let serialised = serde_json::to_string(&s).unwrap();
            file.write_all(&serialised.as_bytes()).unwrap();
        }
    }

    Ok(())
}
