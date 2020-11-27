//
// intersection.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use anyhow::Error;
use criterion::BenchmarkId;
use criterion::{criterion_group, criterion_main, Criterion};

use packing::traits::*;
use packing::wallpaper::{Wallpaper, WyckoffSite};
use packing::{CrystalFamily, LineShape, MolecularShape2, OccupiedSite, PackedState, Transform2};

static BENCH_SIDES: &[usize] = &[4, 16, 64, 256];

/// Utility function to create a Packed State
///
/// This creates a packed state from the number of points used to create a shape.
///
fn create_packed_state(points: usize) -> Result<PackedState<LineShape>, Error> {
    let shape = LineShape::from_radial("Polygon", vec![1.; points])?;

    let wallpaper = Wallpaper {
        name: String::from("p2"),
        family: CrystalFamily::Monoclinic,
    };

    let isopointal = &[WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform2::from_operations("x,y")?,
            Transform2::from_operations("-x,-y")?,
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    }];

    Ok(PackedState::initialise(shape, wallpaper, isopointal))
}

fn state_check_intersection(c: &mut Criterion) {
    let mut group = c.benchmark_group("State Score");

    for &sides in BENCH_SIDES.iter() {
        group.bench_with_input(
            BenchmarkId::new("Polygon", sides),
            &create_packed_state(sides).expect("Creation of state failed"),
            |b, state| b.iter(|| state.score()),
        );
    }
    group.finish();
}

fn shape_check_intersection(c: &mut Criterion) {
    let mut group = c.benchmark_group("Shape Intersection");

    for &sides in BENCH_SIDES.iter() {
        group.bench_with_input(
            BenchmarkId::new("Polygon", sides),
            &LineShape::from_radial("Polygon", vec![1.; sides]).expect("Creation of shape failed"),
            |b, shape| {
                let si1 = shape.transform(&Transform2::new(PI / 3., (0.2, -5.3)));
                let si2 = shape.transform(&Transform2::new(-PI / 3., (-0.2, 5.3)));
                b.iter(|| si1.intersects(&si2))
            },
        );
    }
    group.bench_with_input(
        BenchmarkId::new("Molecule", 1),
        &MolecularShape2::circle(),
        |b, shape| {
            let si1 = shape.transform(&Transform2::new(PI / 3., (0.2, -5.3)));
            let si2 = shape.transform(&Transform2::new(-PI / 3., (-0.2, 5.3)));
            b.iter(|| si1.intersects(&si2))
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Molecule", 3),
        &MolecularShape2::from_trimer(0.637_556, 180., 1.0),
        |b, shape| {
            let si1 = shape.transform(&Transform2::new(PI / 3., (0.2, -5.3)));
            let si2 = shape.transform(&Transform2::new(-PI / 3., (-0.2, 5.3)));
            b.iter(|| si1.intersects(&si2))
        },
    );
    group.finish();
}

fn create_shape_instance(c: &mut Criterion) {
    let mut group = c.benchmark_group("Transform Shape");

    for &sides in BENCH_SIDES.iter() {
        group.bench_with_input(
            BenchmarkId::new("Polygon", sides),
            &LineShape::from_radial("Polygon", vec![1.; sides]).expect("Creation of shape failed"),
            |b, shape| {
                let trans = &Transform2::new(PI / 3., (0.2, -5.3));
                b.iter(|| shape.transform(trans))
            },
        );
    }
    group.bench_with_input(
        BenchmarkId::new("Molecule", 1),
        &MolecularShape2::circle(),
        |b, shape| {
            let trans = &Transform2::new(PI / 3., (0.2, -5.3));
            b.iter(|| shape.transform(trans))
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Molecule", 3),
        &MolecularShape2::from_trimer(0.637_556, 180., 1.0),
        |b, shape| {
            let trans = &Transform2::new(PI / 3., (0.2, -5.3));
            b.iter(|| shape.transform(trans))
        },
    );
    group.finish();
}

fn site_positions(c: &mut Criterion) {
    let site = OccupiedSite::from_wyckoff(&WyckoffSite {
        letter: 'd',
        symmetries: vec![
            Transform2::from_operations("x,y").expect("Transform is invalid"),
            Transform2::from_operations("-x,-y").expect("Transform is invalid"),
        ],
        num_rotations: 1,
        mirror_primary: false,
        mirror_secondary: false,
    });

    c.bench_function("Site Positions", |b| {
        b.iter(|| {
            for _ in site.positions() {
                criterion::black_box(0);
            }
        })
    });
}

fn state_modify_basis(c: &mut Criterion) {
    let state = create_packed_state(256).expect("Creation of state failed");
    let mut basis = state.generate_basis();

    c.bench_function("Modify Basis", |b| {
        b.iter(|| {
            for _ in 0..1000 {
                for value in basis.iter_mut() {
                    let val = value.get_value();
                    value.set_value(val + 0.1).unwrap();
                    value.set_value(val).unwrap();
                }
            }
        })
    });
}

criterion_group!(
    intersections,
    shape_check_intersection,
    create_shape_instance,
    state_check_intersection,
);

criterion_group!(general, site_positions, state_modify_basis);

criterion_main!(intersections, general);
