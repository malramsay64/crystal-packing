//
// packing.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use packing;
use packing::wallpaper::Wallpaper;
use packing::wallpaper::WyckoffSite;

use packing::U2::{Cell2, CrystalFamily, LineShape, OccupiedSite, Transform2};
#[allow(unused_imports)]
use packing::{monte_carlo_best_packing, FromSymmetry, MCVars, PackedState, Shape};

#[test]
fn test_packing_improves() -> Result<(), &'static str> {
    let square = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();

    let wallpaper = Wallpaper {
        name: String::from("p2"),
        family: CrystalFamily::Monoclinic,
    };

    let isopointal = &[WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform2::from_operations("x,y"),
            Transform2::from_operations("-x,-y"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let state =
        PackedState::<LineShape, Cell2, OccupiedSite>::initialise(square, wallpaper, isopointal);

    let init_packing = state.packing_fraction();

    let vars = MCVars {
        kt_start: 0.1,
        kt_finish: 0.001,
        max_step_size: 0.1,
        num_start_configs: 1,
        steps: 100,
        seed: Some(0),
    };

    let final_state = monte_carlo_best_packing(&vars, state);

    let final_packing = final_state?.packing_fraction();

    assert!(init_packing < final_packing);

    Ok(())
}

#[test]
fn test_periodic_overlap() -> Result<(), &'static str> {
    let square = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();

    let wallpaper = Wallpaper {
        name: String::from("p1"),
        family: CrystalFamily::Orthorhombic,
    };

    let isopointal = &[WyckoffSite {
        letter: 'a',
        symmetries: vec![Transform2::from_operations("x,y")],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let state =
        PackedState::<LineShape, Cell2, OccupiedSite>::initialise(square, wallpaper, isopointal);

    let init_packing = state.packing_fraction();

    let vars = MCVars {
        kt_start: 0.1,
        kt_finish: 0.001,
        max_step_size: 0.1,
        num_start_configs: 1,
        steps: 100,
        seed: Some(0),
    };

    let final_state = monte_carlo_best_packing(&vars, state);

    let final_packing = final_state?.packing_fraction();

    assert!(init_packing < final_packing);

    Ok(())
}
