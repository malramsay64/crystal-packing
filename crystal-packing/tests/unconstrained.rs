//
// unconstrained.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Packing shapes with no symmetry constraints: N independent shapes in a periodic cell.
//!
//! Without symmetry the optimisation has to discover the relationships between shapes itself,
//! so these check it can find the known dense packings from independent shapes.
//!
//! ```text
//! cargo test --release -p crystal_packing --test unconstrained -- --ignored
//! ```

use crystal_packing::search::{AnyShape, Search, SearchConfig, ShapeSpec, Strategy};
use crystal_packing::validate::validate;

/// The best packing fraction for each number of shapes
fn unconstrained(shape: ShapeSpec, counts: Vec<usize>) -> Vec<(usize, f64)> {
    unconstrained_with(shape, counts, 0)
}

fn unconstrained_with(
    shape: ShapeSpec,
    counts: Vec<usize>,
    steps_per_shape: u64,
) -> Vec<(usize, f64)> {
    let config = SearchConfig {
        unconstrained: counts,
        strategy: Strategy::Exhaustive,
        steps_per_shape,
        seed: 2,
        ..Default::default()
    };
    let mut search = Search::new(shape.clone(), config).unwrap();
    search.run();
    let built = shape.build().unwrap();
    search
        .results
        .iter()
        .map(|result| {
            let best = result.best.as_ref().expect("Every count has a structure");
            match &built {
                AnyShape::Polygon(s) => validate(s, best).unwrap(),
                AnyShape::Molecule(s) => validate(s, best).unwrap(),
            }
            println!("{} shapes: {:.6}", best.num_shapes, best.packing_fraction);
            (best.num_shapes, best.packing_fraction)
        })
        .collect()
}

#[test]
#[ignore]
fn squares() {
    for (n, phi) in unconstrained(ShapeSpec::Polygon { sides: 4 }, vec![1, 2, 3]) {
        assert!(phi > 0.999, "{} squares: {}", n, phi);
    }
}

#[test]
#[ignore]
fn triangles() {
    // Triangles need two orientations, so a single triangle per cell can't fill space
    let results = unconstrained(ShapeSpec::Polygon { sides: 3 }, vec![1, 2]);
    assert!(results[0].1 < 0.7, "{}", results[0].1);
    assert!(results[1].1 > 0.999, "{}", results[1].1);
}

#[test]
#[ignore]
fn pentagons() {
    // Two pentagons find the double lattice packing without being told about the 2-fold rotation
    let results = unconstrained(ShapeSpec::Polygon { sides: 5 }, vec![1, 2]);
    let double_lattice = (5. - 5f64.sqrt()) / 3.;
    assert!(results[0].1 < double_lattice - 0.05, "{}", results[0].1);
    assert!(
        (results[1].1 - double_lattice).abs() < 1e-3,
        "{}",
        results[1].1
    );
}

#[test]
#[ignore]
fn disks() {
    for (n, phi) in unconstrained(ShapeSpec::Disk, vec![1, 2, 4]) {
        assert!((phi - 0.9069).abs() < 1e-3, "{} disks: {}", n, phi);
    }
}

#[test]
#[ignore]
fn many_pentagons() {
    // Four and six independent pentagons, which each need to find the pairing of the double lattice
    let double_lattice = (5. - 5f64.sqrt()) / 3.;
    for (n, phi) in unconstrained(ShapeSpec::Polygon { sides: 5 }, vec![4, 6]) {
        assert!(
            (phi - double_lattice).abs() < 1e-3,
            "{} pentagons: {}",
            n,
            phi
        );
    }
}

#[test]
#[ignore]
fn heptagons() {
    let results = unconstrained(ShapeSpec::Polygon { sides: 7 }, vec![2]);
    assert!(
        (results[0].1 - 0.892690686).abs() < 1e-3,
        "{}",
        results[0].1
    );
}

#[test]
#[ignore]
fn larger_systems() {
    // With a budget proportional to the number of shapes, larger systems still find the
    // densest packings
    let double_lattice = (5. - 5f64.sqrt()) / 3.;
    for (n, phi) in unconstrained_with(ShapeSpec::Polygon { sides: 5 }, vec![12], 2500) {
        assert!(
            (phi - double_lattice).abs() < 1e-3,
            "{} pentagons: {}",
            n,
            phi
        );
    }
    for (n, phi) in unconstrained_with(ShapeSpec::Disk, vec![24], 2500) {
        assert!((phi - 0.906900).abs() < 1e-3, "{} disks: {}", n, phi);
    }
}
