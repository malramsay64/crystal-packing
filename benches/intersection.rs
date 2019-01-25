//
// intersection.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#[macro_use]
extern crate criterion;

use criterion::{Criterion, ParameterizedBenchmark};
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

fn setup_shapes(shape: &packing::RadialShape) -> (packing::ShapeInstance, packing::ShapeInstance) {
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

fn state_check_intersection(c: &mut Criterion) {
    let parameters: Vec<usize> = (0..4).map(|x| 2_u64.pow(x) as usize).collect();
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
    let parameters: Vec<usize> = (0..4).map(|x| 2_u64.pow(x) as usize).collect();

    let benchmark = ParameterizedBenchmark::new(
        "Shape Intersection Scaling",
        |b, &param| {
            let shape = create_polygon(param);
            let (si1, si2) = setup_shapes(&shape);
            b.iter(|| si1.intersects(&si2))
        },
        parameters,
    );
    c.bench("test_bench_param", benchmark);
}

criterion_group!(
    intersections,
    shape_check_intersection,
    state_check_intersection,
);
criterion_main!(intersections);
