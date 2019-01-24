//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//
#[macro_use]
extern crate log;
extern crate env_logger;

use packing;

fn main() {
    env_logger::init();

    debug!("test message");
    let square = packing::RadialShape {
        name: String::from("Square"),
        radial_points: vec![1., 1., 1., 1.],
        rotational_symmetries: 4,
        mirrors: 4,
    };

    let wallpaper = packing::Wallpaper {
        name: String::from("p2"),
        family: packing::CrystalFamily::Monoclinic,
    };

    let isopointal = &[packing::WyckoffSite {
        letter: 'd',
        symmetries: vec![
            packing::SymmetryTransform::new("x,y"),
            packing::SymmetryTransform::new("-x,-y"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let x = isopointal[0].clone();
    for sym in x.symmetries {
        println!("{:?}", sym);
    }

    let mut state = packing::PackedState::initialise(square, wallpaper, isopointal);

    for position in state.cartesian_positions() {
        println!("{:?}", position);
    }

    println!("Initial state intersects: {}", state.check_intersection());

    println!(
        "Cell Area: {}, Shape Area: {}",
        state.cell.area(),
        state.shape.area()
    );
    println!("Init packing fraction: {}", state.packing_fraction());

    let vars = packing::MCVars {
        kt_start: 0.1,
        kt_finish: 0.001,
        max_step_size: 0.1,
        num_start_configs: 1,
        steps: 1000,
    };

    let final_state = packing::monte_carlo_best_packing(&vars, &mut state);

    println!(
        "Cell Area: {}, Shape Area: {}",
        state.cell.area(),
        state.shape.area()
    );
    println!("{:?}", state.cell);

    println!("Final packing fraction: {}", final_state.packing_fraction());
}
