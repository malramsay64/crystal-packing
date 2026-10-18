//
// engine.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Benchmarks of the hard shape engine, compared with the original implementation.

use std::f64::consts::PI;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use crystal_packing::engine::Packing;
use crystal_packing::geometry::{Mat2, Vec2};
use crystal_packing::hard_shape::{HardShape, Molecule, Polygon};
use crystal_packing::monte_carlo::{MonteCarloConfig, Optimiser};
use crystal_packing::symmetry::PlaneGroup;
use crystal_packing::traits::{Intersect, Shape};
use crystal_packing::wallpaper::{Wallpaper, WyckoffSite};
use crystal_packing::{CrystalFamily, LineShape, MCOptimiser, PackedState, Rng, Transform2};

static SIDES: &[usize] = &[3, 5, 8, 16, 64];

/// Pairs of shapes close enough that the bounding circles overlap, the first pair overlapping and
/// the second separated.
/// The orientation and position of a shape
type Pose = (f64, Vec2);

const PAIRS: [(&str, Pose, Pose); 2] = [
    (
        "overlapping",
        (PI / 3., Vec2::new(0.2, -0.7)),
        (-PI / 3., Vec2::new(-0.2, 0.7)),
    ),
    (
        "separated",
        (0., Vec2::new(0., -0.95)),
        (PI, Vec2::new(0., 0.95)),
    ),
];

fn pair_overlap(c: &mut Criterion) {
    for (label, (angle_a, pos_a), (angle_b, pos_b)) in PAIRS {
        let mut group = c.benchmark_group(format!("Pair overlap {}", label));
        for &sides in SIDES {
            let shape = Polygon::regular(sides).unwrap();
            let mut a = Default::default();
            let mut b = Default::default();
            shape.place(&Mat2::rotation(angle_a), pos_a, &mut a);
            shape.place(&Mat2::rotation(angle_b), pos_b, &mut b);
            group.bench_with_input(BenchmarkId::new("SAT", sides), &sides, |bench, _| {
                bench.iter(|| Polygon::overlaps(black_box(&a), black_box(&b), Vec2::ZERO))
            });

            let legacy = LineShape::from_radial("Polygon", vec![1.; sides]).unwrap();
            let la = legacy.transform(&Transform2::new(angle_a, (pos_a.x, pos_a.y)));
            let lb = legacy.transform(&Transform2::new(angle_b, (pos_b.x, pos_b.y)));
            group.bench_with_input(BenchmarkId::new("Legacy", sides), &sides, |bench, _| {
                bench.iter(|| black_box(&la).intersects(black_box(&lb)))
            });
        }
        let disk = Molecule::disk();
        let mut a = Default::default();
        let mut b = Default::default();
        disk.place(&Mat2::IDENTITY, pos_a, &mut a);
        disk.place(&Mat2::IDENTITY, pos_b, &mut b);
        group.bench_function("Disk", |bench| {
            bench.iter(|| Molecule::overlaps(black_box(&a), black_box(&b), Vec2::ZERO))
        });
        group.finish();
    }
}

/// A dense packing of pentagons in p2, the starting point for the state benchmarks
fn dense_packing(sides: usize) -> Packing<Polygon> {
    let group = PlaneGroup::from_name("p2").unwrap();
    let mut rng = Rng::seed_from_u64(0);
    let packing = Packing::new(Polygon::regular(sides).unwrap(), group, &['e'], &mut rng).unwrap();
    let mut optimiser = Optimiser::new(packing, MonteCarloConfig::default(), 0);
    optimiser.run(20_000);
    let structure = optimiser.finish();
    Packing::from_structure(Polygon::regular(sides).unwrap(), &structure).unwrap()
}

fn legacy_state(sides: usize) -> PackedState<LineShape> {
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
    PackedState::initialise(
        LineShape::from_radial("Polygon", vec![1.; sides]).unwrap(),
        wallpaper,
        isopointal,
    )
}

fn state_overlap(c: &mut Criterion) {
    let mut group = c.benchmark_group("State overlap check");
    for &sides in &[5, 16] {
        let packing = dense_packing(sides);
        group.bench_with_input(BenchmarkId::new("Full", sides), &sides, |b, _| {
            b.iter(|| black_box(&packing).any_overlap())
        });
        group.bench_with_input(BenchmarkId::new("Site", sides), &sides, |b, _| {
            b.iter(|| black_box(&packing).site_overlaps(0))
        });
        let legacy = legacy_state(sides);
        group.bench_with_input(BenchmarkId::new("Legacy", sides), &sides, |b, _| {
            b.iter(|| crystal_packing::traits::State::score(black_box(&legacy)))
        });
    }
    group.finish();
}

fn monte_carlo(c: &mut Criterion) {
    let mut group = c.benchmark_group("Monte Carlo 10000 steps");
    group.sample_size(20);
    let config = MonteCarloConfig {
        steps: 10_000,
        ..Default::default()
    };
    group.bench_function("Pentagon p2", |b| {
        b.iter(|| {
            let group = PlaneGroup::from_name("p2").unwrap();
            let mut rng = Rng::seed_from_u64(0);
            let packing =
                Packing::new(Polygon::regular(5).unwrap(), group, &['e'], &mut rng).unwrap();
            let mut optimiser = Optimiser::new(packing, config, 0);
            optimiser.run(10_000)
        })
    });
    group.bench_function("Legacy pentagon p2", |b| {
        b.iter(|| {
            let optimiser = MCOptimiser::new(0.1, 0.9, 0.01, 10_000, 1000, 0, None);
            crystal_packing::traits::State::score(&optimiser.optimise_state(legacy_state(5)))
        })
    });
    group.finish();
}

criterion_group!(benches, pair_overlap, state_overlap, monte_carlo);
criterion_main!(benches);
