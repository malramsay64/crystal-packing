//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

// This is for the clap::arg_enum macro
#![allow(deprecated)]

use std::f64::consts::PI;

use clap::{arg_enum, App, Arg, _clap_count_exprs, value_t};
use log::{debug, info};
use rayon::prelude::*;
use simplelog::{Config, LevelFilter, TermLogger};

use packing::traits::*;
use packing::wallpaper::{Wallpaper, WallpaperGroup, WallpaperGroups, WyckoffSite};
use packing::{
    monte_carlo_best_packing, Cell2, LineShape, MCVars, MolecularShape2, OccupiedSite, PackedState,
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
    }
}

fn get_packed_state<S, C, T>(
    options: CLIOptions,
    shape: S,
) -> Result<PackedState<S, C, T>, &'static str>
where
    S: Shape + Send + Sync,
    C: Cell<Transform = S::Transform> + Send + Sync,
    T: Site<Transform = S::Transform> + Send + Sync,
{
    let wallpaper = Wallpaper::new(&options.group);
    let isopointal = &[WyckoffSite::new(options.group)];

    let state = PackedState::<S, C, T>::initialise(shape.clone(), wallpaper.clone(), isopointal);
    match state.score() {
        Err(_) => panic!("Initial state has intersections...exiting."),
        Ok(x) => info!("Init packing fraction: {}", x),
    };

    let mut vars = MCVars::default();
    vars.steps = options.steps;
    vars.num_start_configs = options.num_start_configs;
    // Remove mutability
    let vars = vars;

    let final_state = (0..vars.num_start_configs)
        .into_par_iter()
        .map(|_| {
            let state =
                PackedState::<S, C, T>::initialise(shape.clone(), wallpaper.clone(), isopointal);
            monte_carlo_best_packing(&vars, state)
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

    info!("Final packing fraction: {}", final_state.score().unwrap());
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
    println!("Using Wallpaper Group: {}", wg);

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

    match options.shape {
        ShapeTypes::Polygon => {
            let shape = LineShape::from_radial("Polygon", vec![1.; options.num_sides])?;
            get_packed_state::<LineShape, Cell2, OccupiedSite>(options, shape).unwrap();
        }
        ShapeTypes::Trimer => {
            let shape = MolecularShape2::from_trimer(0.637_556, 2. * PI / 3., 1.);
            get_packed_state::<MolecularShape2, Cell2, OccupiedSite>(options, shape).unwrap();
        }
        ShapeTypes::Circle => {
            let shape = MolecularShape2::circle();
            get_packed_state::<MolecularShape2, Cell2, OccupiedSite>(options, shape).unwrap();
        }
    };
    Ok(())
}
