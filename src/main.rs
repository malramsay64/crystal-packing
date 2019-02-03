//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

// This is for the clap::arg_enum macro
#![allow(deprecated)]

#[macro_use]
extern crate clap;
extern crate simplelog;
#[macro_use]
extern crate log;
extern crate rayon;

use clap::{App, Arg};
use packing;
#[allow(unused_imports)]
use packing::shape::{Atom, LineShape, MolecularShape, Shape};
use packing::wallpaper::{WallpaperGroup, WallpaperGroups};
use packing::PackedState;
use rayon::prelude::*;
use simplelog::{Config, LevelFilter, TermLogger};
use std::f64::consts::PI;

struct CLIOptions {
    shape: ShapeTypes,
    group: WallpaperGroup,
    steps: u64,
    num_sides: usize,
    log_level: LevelFilter,
}

arg_enum! {
    #[derive(Debug)]
    pub enum ShapeTypes {
        polygon,
        trimer,
        circle,
    }
}

fn get_packed_state<T>(options: CLIOptions, shape: T) -> Result<PackedState<T>, &'static str>
where
    T: Shape + Send + Sync,
{
    let wallpaper = packing::Wallpaper::new(&options.group);
    let isopointal = &[packing::WyckoffSite::new(options.group)];

    let state = PackedState::initialise(shape.clone(), wallpaper.clone(), isopointal);
    if state.check_intersection() {
        panic!("Initial state has intersetions...exiting.");
    }

    info!(
        "Init packing fraction: {}",
        state.packing_fraction().unwrap()
    );

    let mut vars = packing::MCVars::default();
    vars.steps = options.steps;
    vars.num_start_configs = 32;
    // Remove mutability
    let vars = vars;

    let final_state = (0..vars.num_start_configs)
        .into_par_iter()
        .map(|_| {
            let mut state =
                packing::PackedState::initialise(shape.clone(), wallpaper.clone(), isopointal);
            packing::monte_carlo_best_packing(&vars, &mut state)
        })
        .max()
        .unwrap()?;

    info!(
        "Cell Area: {}, Shape Area: {}, Num Shapes: {}",
        final_state.cell.area(),
        shape.area(),
        final_state.total_shapes(),
    );

    final_state.to_figure("test.txt");

    info!(
        "Final packing fraction: {}",
        final_state.packing_fraction().unwrap()
    );
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
        .get_matches();

    let wg = value_t!(matches.value_of("wallpaper_group"), WallpaperGroups).unwrap();
    println!("Using Wallpaper Group: {}", wg);

    let shape = value_t!(matches.value_of("shape"), ShapeTypes).unwrap();
    let num_sides = matches.value_of("sides").unwrap().parse().unwrap();

    let group = packing::wallpaper::get_wallpaper_group(wg).unwrap();

    let steps: u64 = matches.value_of("steps").unwrap().parse().unwrap();

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
        log_level,
    }
}

fn main() -> Result<(), &'static str> {
    let options = cli();

    TermLogger::init(options.log_level, Config::default()).unwrap();
    debug!("Logging Level: {}", options.log_level);

    match options.shape {
        ShapeTypes::polygon => {
            let shape = LineShape::from_radial("Polygon", vec![1.; options.num_sides])?;
            get_packed_state(options, shape).unwrap();
        }
        ShapeTypes::trimer => {
            let shape = MolecularShape::from_trimer(0.637_556, 2. * PI / 3., 1.);
            get_packed_state(options, shape).unwrap();
        }
        ShapeTypes::circle => {
            let shape = MolecularShape::circle();
            get_packed_state(options, shape).unwrap();
        }
    };
    Ok(())
}
