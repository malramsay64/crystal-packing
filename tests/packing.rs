//
// packing.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use packing;
#[allow(unused_imports)]
use packing::shape::{LineShape, Shape};
use packing::symmetry::Transform;

#[test]
fn test_packing_improves() -> Result<(), &'static str> {
    let square = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();

    let wallpaper = packing::Wallpaper {
        name: String::from("p2"),
        family: packing::CrystalFamily::Monoclinic,
    };

    let isopointal = &[packing::WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform::from_operations("x,y"),
            Transform::from_operations("-x,-y"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let mut state = packing::PackedState::initialise(square, wallpaper, isopointal);

    let init_packing = state.packing_fraction();

    let vars = packing::MCVars {
        kt_start: 0.1,
        kt_finish: 0.001,
        max_step_size: 0.1,
        num_start_configs: 1,
        steps: 100,
        seed: Some(0),
    };

    let final_state = packing::monte_carlo_best_packing(&vars, &mut state);

    let final_packing = final_state?.packing_fraction();

    assert!(init_packing < final_packing);

    Ok(())
}
