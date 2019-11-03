//
// packing.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use anyhow::{anyhow, Error};
use packing;
use packing::traits::*;
use packing::wallpaper::Wallpaper;
use packing::wallpaper::WyckoffSite;
use packing::{
    BuildOptimiser, Cell2, CrystalFamily, LineShape, OccupiedSite, PackedState, Transform2,
};

#[test]
fn test_packing_improves() -> Result<(), Error> {
    let square = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();

    let wallpaper = Wallpaper {
        name: String::from("p2"),
        family: CrystalFamily::Monoclinic,
    };

    let isopointal = &[WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform2::from_operations("x,y").unwrap(),
            Transform2::from_operations("-x,-y").unwrap(),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let state =
        PackedState::<LineShape, Cell2, OccupiedSite>::initialise(square, wallpaper, isopointal);

    let init_packing = state
        .score()
        .ok_or_else(|| anyhow!("Invalid initial state"))?;

    let opt = BuildOptimiser::default()
        .seed(0)
        .steps(1000)
        .kt_start(0.)
        .kt_ratio(Some(0.))
        .max_step_size(0.001)
        .build();

    let final_state = opt.optimise_state(state);

    let final_packing = final_state
        .score()
        .ok_or_else(|| anyhow!("Invalid final state"))?;

    println!("Init Score: {} Final score {}", init_packing, final_packing);
    assert!(init_packing < final_packing);

    Ok(())
}
