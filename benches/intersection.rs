//
// intersection.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#[macro_use]
extern crate criterion;

use criterion::{Criterion, ParameterizedBenchmark};
use nalgebra::Vector2;
use packing::shape::{LineShape, Shape, ShapeInstance};
use packing::symmetry::Transform;

fn create_polygon(sides: usize) -> Result<LineShape, &'static str> {
    LineShape::from_radial("Polygon", vec![1.; sides])
}

fn setup_state(points: usize) -> packing::PackedState<LineShape> {
    let shape = create_polygon(points).unwrap();

    let wallpaper = packing::Wallpaper {
        name: String::from("p2"),
        family: packing::CrystalFamily::Monoclinic,
    };

    let isopointal = &[packing::WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform::from_operations("x,y"),
            Transform::from_operations("-x,-y"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    packing::PackedState::initialise(shape, wallpaper, isopointal)
}

fn setup_shapes<S>(shape: &S) -> (ShapeInstance<S::Component>, ShapeInstance<S::Component>)
where
    S: Shape,
{
    let si1 = ShapeInstance::from(shape, Transform::new(Vector2::new(-2., 0.), 0.));
    let si2 = ShapeInstance::from(shape, Transform::new(Vector2::new(2., 0.), 0.));

    (si1, si2)
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
