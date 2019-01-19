//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use nalgebra::{Matrix2, Vector2};
use packing;

fn main() {
    let square = packing::Shape {
        name: String::from("Square"),
        radial_points: vec![1., 1., 1., 1.],
        rotational_symmetries: 4,
        mirrors: 4,
    };

    let wallpaper = packing::Wallpaper {
        name: String::from("p1"),
        family: packing::CrystalFamily::Monoclinic,
    };

    let isopointal = vec![packing::WyckoffSite {
        letter: 'a',
        symmetries: vec![packing::SymmetryTransform {
            rotation: Matrix2::new(1., 0., 1., 0.),
            translation: Vector2::new(0., 0.),
        }],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let state = packing::PackedState::initialise(square, wallpaper, isopointal, 0.1);

    assert_eq!(state.total_shapes(), 1);

    println!("Hello, world!");
}
