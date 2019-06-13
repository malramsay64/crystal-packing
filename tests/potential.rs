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
    monte_carlo_best_packing, Cell2, CrystalFamily, FromSymmetry, LJShape2, MCVars, OccupiedSite,
    PotentialState, Transform2,
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

    let vars = MCVars {
        kt_start: 0.1,
        kt_finish: 0.001,
        max_step_size: 0.1,
        num_start_configs: 1,
        steps: 100,
        seed: Some(0),
    };

    let final_state = monte_carlo_best_packing(&vars, state);

    let final_score = final_state?.score()?;

    println!("Init Score: {}, Final Score: {}", init_score, final_score);
    assert!(init_score < final_score);

    Ok(())
}
