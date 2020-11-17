extern crate cfg_if;
extern crate wasm_bindgen;

use cfg_if::cfg_if;
use wasm_bindgen::prelude::*;

use packing::*;
use packing::wallpaper::{Wallpaper, WyckoffSite};

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
extern {
    fn alert(s: &str);
}

#[wasm_bindgen]
pub fn greet() {
    alert("Hello, packing-wasm!");
}

#[wasm_bindgen]
pub struct State (PackedState<MolecularShape2>);

#[wasm_bindgen]
pub fn setup_state() -> State {
    let trimer = MolecularShape2::from_trimer(1., 0.7, 120.);

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

    State(PackedState::<MolecularShape2>::initialise(trimer, wallpaper, isopointal))
}

pub fn setup_opt() -> MCOptimiser {
    MCOptimiser::new(0., 0., 0.001, 1000, 100, 0, None)
}
