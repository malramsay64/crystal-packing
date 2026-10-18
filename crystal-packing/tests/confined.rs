//
// confined.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Packing a finite number of shapes in a container, compared with known optimal packings.
//!
//! ```text
//! cargo test --release -p crystal_packing --test confined -- --ignored
//! ```

use std::f64::consts::PI;

use crystal_packing::container::Container;
use crystal_packing::search::{AnyShape, Search, SearchConfig, ShapeSpec, Strategy};
use crystal_packing::validate::validate_confined;

/// The best packing fraction for each number of shapes in the container
fn confined(shape: ShapeSpec, container: Container, counts: Vec<usize>) -> Vec<f64> {
    let config = SearchConfig {
        unconstrained: counts,
        container: Some(container),
        strategy: Strategy::Exhaustive,
        steps_per_shape: 2500,
        seed: 5,
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
            assert_eq!(best.container.unwrap().kind, container);
            match &built {
                AnyShape::Polygon(s) => validate_confined(s, best).unwrap(),
                AnyShape::Molecule(s) => validate_confined(s, best).unwrap(),
            }
            println!(
                "{:?} {} in {:?}: {:.6}",
                best.shape, best.num_shapes, container, best.packing_fraction
            );
            best.packing_fraction
        })
        .collect()
}

fn assert_close(found: f64, expected: f64, label: &str) {
    assert!(
        (found - expected).abs() < 2e-3 && found < expected + 1e-6,
        "{}: found {:.6}, expected {:.6}",
        label,
        found,
        expected
    );
}

#[test]
#[ignore]
fn disks_in_a_circle() {
    // Known optima: 2 disks need a container of radius 2, 3 disks 1 + 2/√3, and 7 disks
    // (a hexagon around a central disk) radius 3
    let found = confined(ShapeSpec::Disk, Container::Circle, vec![2, 3, 7]);
    assert_close(found[0], 2. / 4., "2 disks");
    let r3 = 1. + 2. / 3f64.sqrt();
    assert_close(found[1], 3. / (r3 * r3), "3 disks");
    assert_close(found[2], 7. / 9., "7 disks");
}

#[test]
#[ignore]
fn disks_in_a_square() {
    // 2 disks along the diagonal of a square of side 2 + √2, and 4 disks in a square of side 4
    let found = confined(ShapeSpec::Disk, Container::Square, vec![2, 4]);
    let side = 2. + 2f64.sqrt();
    assert_close(found[0], 2. * PI / (side * side), "2 disks");
    assert_close(found[1], PI / 4., "4 disks");
}

#[test]
#[ignore]
fn squares() {
    let square = ShapeSpec::Polygon { sides: 4 };
    // Four squares fill a square, while two fill only half of the smallest square around them
    let found = confined(square.clone(), Container::Square, vec![2, 4]);
    assert_close(found[0], 0.5, "2 squares in a square");
    assert_close(found[1], 1., "4 squares in a square");
    // With a free aspect ratio two squares fill a rectangle
    let found = confined(square, Container::Rectangle, vec![2]);
    assert_close(found[0], 1., "2 squares in a rectangle");
}

#[test]
#[ignore]
fn triangle_in_a_circle() {
    // A unit triangle is inscribed in a unit circle
    let found = confined(ShapeSpec::Polygon { sides: 3 }, Container::Circle, vec![1]);
    assert_close(found[0], 3. * 3f64.sqrt() / 4. / PI, "1 triangle");
}
