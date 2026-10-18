//
// regression.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Known densest packings which a full search over the plane groups must find.
//!
//! These take a few seconds each in release mode so are ignored by default. Run them with
//!
//! ```text
//! cargo test --release -p crystal_packing --test regression -- --ignored
//! ```

use std::f64::consts::PI;

use crystal_packing::search::{AnyShape, Search, SearchConfig, ShapeSpec};
use crystal_packing::validate::validate;

const TOLERANCE: f64 = 1e-3;

fn search(shape: ShapeSpec) -> f64 {
    let config = SearchConfig {
        seed: 1,
        ..Default::default()
    };
    let mut search = Search::new(shape.clone(), config).unwrap();
    let best = search.run().expect("No structure found").clone();
    // Every result is validated as part of the search, check the best again to be certain
    match shape.build().unwrap() {
        AnyShape::Polygon(s) => validate(&s, &best).unwrap(),
        AnyShape::Molecule(s) => validate(&s, &best).unwrap(),
    }
    println!(
        "{:?}: {:.6} in {} ({})",
        shape,
        best.packing_fraction,
        best.group,
        best.sites_label()
    );
    best.packing_fraction
}

fn assert_close(found: f64, expected: f64) {
    assert!(
        (found - expected).abs() < TOLERANCE,
        "Found {:.6}, expected {:.6}",
        found,
        expected
    );
    // The packing fraction can exceed the optimum only by the overlap tolerance
    assert!(
        found < expected + 1e-6,
        "Found {} exceeds {}",
        found,
        expected
    );
}

#[test]
#[ignore]
fn disk() {
    assert_close(search(ShapeSpec::Disk), PI / (2. * 3f64.sqrt()));
}

#[test]
#[ignore]
fn triangle() {
    assert_close(search(ShapeSpec::Polygon { sides: 3 }), 1.);
}

#[test]
#[ignore]
fn square() {
    assert_close(search(ShapeSpec::Polygon { sides: 4 }), 1.);
}

#[test]
#[ignore]
fn pentagon() {
    // The double lattice packing (Kuperberg & Kuperberg 1990)
    assert_close(
        search(ShapeSpec::Polygon { sides: 5 }),
        (5. - 5f64.sqrt()) / 3.,
    );
}

#[test]
#[ignore]
fn hexagon() {
    assert_close(search(ShapeSpec::Polygon { sides: 6 }), 1.);
}

#[test]
#[ignore]
fn heptagon() {
    // The best known packing is a double lattice (Kuperberg & Kuperberg 1990)
    assert_close(search(ShapeSpec::Polygon { sides: 7 }), 0.892690686);
}

#[test]
#[ignore]
fn trimer() {
    // No analytical result, this guards against regressions of the search
    let found = search(ShapeSpec::Trimer {
        radius: 0.637556,
        angle: 120.,
        distance: 1.,
    });
    assert!(found > 0.94, "{}", found);
}
