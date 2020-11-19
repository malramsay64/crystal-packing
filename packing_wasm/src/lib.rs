use cfg_if::cfg_if;
use wasm_bindgen::prelude::*;

use packing::wallpaper::{Wallpaper, WyckoffSite};
use packing::*;

pub mod optimiser;
pub mod state;

cfg_if! {
    // When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
    // allocator.
    if #[cfg(feature = "wee_alloc")] {
        extern crate wee_alloc;
        #[global_allocator]
        static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;
    }
}

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}

#[wasm_bindgen]
pub fn greet() {
    alert("Hello, packing-wasm!");
}

#[wasm_bindgen]
pub fn setup_opt(kt: f64, step_size: f64, steps: u64) -> optimiser::Optimiser {
    optimiser::Optimiser::new(kt, step_size, steps, 0)
}

#[wasm_bindgen]
pub fn setup_state() -> state::JSState {
    let trimer = MolecularShape2::from_trimer(0.7, 120., 1.);

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

    state::JSState::new(PackedState::<MolecularShape2>::initialise(
        trimer, wallpaper, isopointal,
    ))
}
