use packing::traits::{State, ToSVG};
use packing::*;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct JSState(PackedState<MolecularShape2>);

#[wasm_bindgen]
impl JSState {
    pub fn as_svg(&self) -> String {
        self.0.as_svg().to_string()
    }

    pub fn score(&self) -> Option<f64> {
        self.0.score()
    }
}

impl JSState {
    pub fn new(state: PackedState<MolecularShape2>) -> JSState {
        JSState(state)
    }

    pub fn generate_basis(&self) -> Vec<StandardBasis> {
        self.0.generate_basis()
    }
}
