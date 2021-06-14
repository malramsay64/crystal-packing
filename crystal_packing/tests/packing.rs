//
// packing.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use anyhow::{anyhow, Error};

use crystal_packing::traits::*;
use crystal_packing::wallpaper::Wallpaper;
use crystal_packing::wallpaper::WyckoffSite;
use crystal_packing::{CrystalFamily, LineShape, MCOptimiser, PackedState, Transform2};

#[test]
fn test_packing_improves() -> Result<(), Error> {
    let square = LineShape::from_radial("Square", vec![1., 1., 1., 1.])?;

    let wallpaper = Wallpaper {
        name: String::from("p2"),
        family: CrystalFamily::Monoclinic,
    };

    let isopointal = &[WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform2::from_operations("x,y")?,
            Transform2::from_operations("-x,-y")?,
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let state = PackedState::<LineShape>::initialise(square, wallpaper, isopointal);

    let init_packing = state
        .score()
        .ok_or_else(|| anyhow!("Invalid initial state"))?;

    let opt = MCOptimiser::new(0., 0., 0.001, 1000, 100, 0, None);

    let final_state = opt.optimise_state(state);

    let final_packing = final_state
        .score()
        .ok_or_else(|| anyhow!("Invalid final state"))?;

    println!("Init Score: {} Final score {}", init_packing, final_packing);
    assert!(init_packing < final_packing);

    Ok(())
}
