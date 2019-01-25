//
// intersection.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#[macro_use]
extern crate bencher;

use bencher::Bencher;
use nalgebra::{IsometryMatrix2, Vector2};
use packing;

fn create_polygon(sides: usize) -> packing::RadialShape {
    packing::RadialShape {
        name: String::from("Polygon"),
        radial_points: vec![1.; sides],
        rotational_symmetries: sides as u64,
        mirrors: sides as u64,
    }
}

fn setup_state(points: usize) -> packing::PackedState {
    let shape = create_polygon(points);

    let wallpaper = packing::Wallpaper {
        name: String::from("p2"),
        family: packing::CrystalFamily::Monoclinic,
    };

    let isopointal = &[packing::WyckoffSite {
        letter: 'd',
        symmetries: vec![
            packing::SymmetryTransform::new("x,y"),
            packing::SymmetryTransform::new("-x,-y"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    packing::PackedState::initialise(shape, wallpaper, isopointal)
}

fn setup_shapes(
    shape: &packing::RadialShape,
    points: usize,
) -> (packing::ShapeInstance, packing::ShapeInstance) {
    let si1 = packing::shape::ShapeInstance {
        shape: &shape,
        isometry: IsometryMatrix2::new(Vector2::new(-2., 0.), 0.),
    };
    let si2 = packing::ShapeInstance {
        shape: &shape,
        isometry: IsometryMatrix2::new(Vector2::new(2., 0.), 0.),
    };

    (si1, si2)
}

fn state_check_intersection_impl(bench: &mut Bencher, points: usize) {
    let state = setup_state(points);
    bench.iter(|| state.check_intersection());
}

fn shape_check_intersection_impl(bench: &mut Bencher, points: usize) {
    let shape = create_polygon(points);
    let (shape_i1, shape_i2) = setup_shapes(&shape, points);
    bench.iter(|| shape_i1.intersects(&shape_i2));
}

fn state_check_intersection_1(bench: &mut Bencher) {
    state_check_intersection_impl(bench, 1);
}

fn state_check_intersection_2(bench: &mut Bencher) {
    state_check_intersection_impl(bench, 2);
}

fn state_check_intersection_4(bench: &mut Bencher) {
    state_check_intersection_impl(bench, 4);
}

fn state_check_intersection_8(bench: &mut Bencher) {
    state_check_intersection_impl(bench, 8);
}

fn state_check_intersection_16(bench: &mut Bencher) {
    state_check_intersection_impl(bench, 16);
}

fn state_check_intersection_32(bench: &mut Bencher) {
    state_check_intersection_impl(bench, 32);
}

fn shape_check_intersection_1(bench: &mut Bencher) {
    shape_check_intersection_impl(bench, 1);
}

fn shape_check_intersection_2(bench: &mut Bencher) {
    shape_check_intersection_impl(bench, 2);
}

fn shape_check_intersection_4(bench: &mut Bencher) {
    shape_check_intersection_impl(bench, 4);
}

fn shape_check_intersection_8(bench: &mut Bencher) {
    shape_check_intersection_impl(bench, 8);
}

fn shape_check_intersection_16(bench: &mut Bencher) {
    shape_check_intersection_impl(bench, 16);
}

fn shape_check_intersection_32(bench: &mut Bencher) {
    shape_check_intersection_impl(bench, 32);
}

benchmark_group!(
    intersections,
    state_check_intersection_1,
    state_check_intersection_2,
    state_check_intersection_4,
    state_check_intersection_8,
    state_check_intersection_16,
    state_check_intersection_32,
    shape_check_intersection_1,
    shape_check_intersection_2,
    shape_check_intersection_4,
    shape_check_intersection_8,
    shape_check_intersection_16,
    shape_check_intersection_32,
);
benchmark_main!(intersections);
