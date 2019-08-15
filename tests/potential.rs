//
// potential.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#[allow(unused_imports)]
use itertools::Itertools;
#[allow(unused_imports)]
use serde::Deserialize;

use packing;
#[allow(unused_imports)]
use packing::traits::*;
use packing::wallpaper::Wallpaper;
use packing::wallpaper::WyckoffSite;
use packing::{
    BuildOptimiser, Cell2, CrystalFamily, FromSymmetry, LJShape2, OccupiedSite, PotentialState,
    Transform2,
};

#[test]
fn test_score_improves() -> Result<(), &'static str> {
    let square = LJShape2::from_trimer(0.63, 120., 1.);

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
        PotentialState::<LJShape2, Cell2, OccupiedSite>::initialise(square, wallpaper, isopointal);

    let init_score = state.score()?;

    let opt = BuildOptimiser::default().seed(0).build();

    let final_state = opt.optimise_state(state);

    let final_score = final_state.score()?;

    println!("Init Score: {}, Final Score: {}", init_score, final_score);
    assert!(init_score < final_score);

    Ok(())
}
