//
// intersection.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;
use std::ops::Mul;

#[macro_use]
extern crate criterion;

use criterion::{Criterion, ParameterizedBenchmark};
use nalgebra::Vector2;

use packing::wallpaper::{Wallpaper, WyckoffSite};
use packing::U2::{Cell2, CrystalFamily, LineShape, OccupiedSite, Transform2};
use packing::{FromSymmetry, PackedState, Shape};

fn create_polygon(sides: usize) -> Result<LineShape, &'static str> {
    LineShape::from_radial("Polygon", vec![1.; sides])
}

fn setup_state(points: usize) -> PackedState<LineShape, Cell2, OccupiedSite> {
    let shape = create_polygon(points).unwrap();

    let wallpaper = Wallpaper {
        name: String::from("p2"),
        family: CrystalFamily::Monoclinic,
    };

    let isopointal = &[WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform2::from_operations("x,y"),
            Transform2::from_operations("-x,-y"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    PackedState::initialise(shape, wallpaper, isopointal)
}

fn state_check_intersection(c: &mut Criterion) {
    let parameters: Vec<usize> = (2..6).map(|x| 2_u64.pow(x) as usize).collect();
    let benchmark = ParameterizedBenchmark::new(
        "State Intersection Scaling",
        |b, &param| {
            let state = setup_state(param);
            b.iter(|| state.check_intersection())
        },
        parameters,
    );
    c.bench("test_bench_param", benchmark);
}

fn shape_check_intersection(c: &mut Criterion) {
    let parameters: Vec<usize> = (2..6).map(|x| 2_u64.pow(x) as usize).collect();

    let benchmark = ParameterizedBenchmark::new(
        "Shape Intersection Scaling",
        |b, &param| {
            let shape = create_polygon(param).unwrap();
            let si1 = shape.transform(&Transform2::new(Vector2::new(0.2, -5.3), PI / 3.));
            let si2 = shape.transform(&Transform2::new(Vector2::new(-0.2, 5.3), -PI / 3.));
            b.iter(|| si1.intersects(&si2))
        },
        parameters,
    );
    c.bench("test_bench_param", benchmark);
}

fn create_shape_instance(c: &mut Criterion) {
    let parameters: Vec<usize> = (2..6).map(|x| 2_u64.pow(x) as usize).collect();

    let benchmark = ParameterizedBenchmark::new(
        "Creation of ShapeInstance",
        |b, &param| {
            let shape = create_polygon(param).unwrap();
            let trans = Transform2::new(Vector2::new(0.2, -5.3), PI / 3.);
            b.iter(|| shape.transform(&trans))
        },
        parameters,
    );
    c.bench("create_shape_instance", benchmark);
}

criterion_group!(
    intersections,
    shape_check_intersection,
    create_shape_instance,
    state_check_intersection,
);
criterion_main!(intersections);
