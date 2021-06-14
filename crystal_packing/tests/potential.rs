//
// potential.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use anyhow::{anyhow, Error};

use crystal_packing::traits::*;
use crystal_packing::wallpaper::Wallpaper;
use crystal_packing::wallpaper::WyckoffSite;
use crystal_packing::{CrystalFamily, LJShape2, MCOptimiser, PotentialState, Transform2};

#[test]
fn test_score_improves() -> Result<(), Error> {
    let square = LJShape2::from_trimer(0.63, 120., 1.);

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

    let state = PotentialState::<LJShape2>::initialise(square, wallpaper, isopointal);

    let init_score = state
        .score()
        .ok_or_else(|| anyhow!("Initial score is invalid"))?;

    let opt = MCOptimiser::new(0., 0., 0.001, 1000, 100, 0, None);

    let final_state = opt.optimise_state(state);

    let final_score = final_state
        .score()
        .ok_or_else(|| anyhow!("Final score is invalid"))?;

    println!("Init Score: {}, Final Score: {}", init_score, final_score);
    assert!(init_score < final_score);

    Ok(())
}
