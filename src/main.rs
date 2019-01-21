//
// main.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use packing;

fn main() {
    let square = packing::Shape {
        name: String::from("Square"),
        radial_points: vec![1., 1., 1., 1.],
        rotational_symmetries: 4,
        mirrors: 4,
    };

    let wallpaper = packing::Wallpaper {
        name: String::from("p2mg"),
        family: packing::CrystalFamily::Monoclinic,
    };

    let isopointal = vec![packing::WyckoffSite {
        letter: 'd',
        symmetries: vec![
            packing::SymmetryTransform::new("x,y"),
            packing::SymmetryTransform::new("-x,-y"),
            packing::SymmetryTransform::new("-x+1/2,y"),
            packing::SymmetryTransform::new("x+1/2,-y"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    let state = packing::PackedState::initialise(square, wallpaper, isopointal);

    assert_eq!(state.total_shapes(), 4);
    println!(
        "Cell Area: {}, Shape Area: {}",
        state.cell.area(),
        state.shape.area()
    );
    println!("{}", state.packing_fraction());
}
