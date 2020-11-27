use cfg_if::cfg_if;
use std::convert::TryInto;
use std::str::FromStr;
use wasm_bindgen::prelude::*;

use packing::wallpaper::WallpaperGroups;
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
pub fn setup_state(radius: f64, angle: f64, distance: f64, wallpaper: String) -> state::JSState {
    let trimer = MolecularShape2::from_trimer(radius, angle, distance);

    let wallpaper: WallpaperGroup<'_> = WallpaperGroups::from_str(&wallpaper)
        .unwrap()
        .try_into()
        .expect("Shouldn't fail");

    state::JSState::new(PackedState::<MolecularShape2>::from_group(trimer, &wallpaper).unwrap())
}
