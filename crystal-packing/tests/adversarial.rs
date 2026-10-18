//
// adversarial.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Shapes chosen to stress the optimisation: extreme aspect ratios, concave shapes, shapes with
//! no symmetry, shapes whose origin is far from their centre, and shapes with many vertices.
//!
//! Many have known answers because they tile the plane: every triangle and every quadrilateral
//! (convex or not) tiles the plane with two copies related by a 2-fold rotation (p2), and the
//! L-tromino and plus pentomino tile the plane by translation alone (p1).
//!
//! These take a few seconds each in release mode so are ignored by default. Run them with
//!
//! ```text
//! cargo test --release -p crystal_packing --test adversarial -- --ignored
//! ```

use crystal_packing::search::{AnyShape, Search, SearchConfig, ShapeSpec};
use crystal_packing::validate::validate;

fn vertices(name: &str, points: &[[f64; 2]]) -> ShapeSpec {
    ShapeSpec::Vertices {
        name: name.into(),
        vertices: points.to_vec(),
    }
}

/// Search the given groups (all groups when empty), returning the best packing fraction
fn search(shape: ShapeSpec, groups: &[&str]) -> f64 {
    let config = SearchConfig {
        groups: groups.iter().map(|g| g.to_string()).collect(),
        seed: 3,
        ..Default::default()
    };
    let mut search = Search::new(shape.clone(), config).unwrap();
    let best = search.run().expect("No structure found").clone();
    match shape.build().unwrap() {
        AnyShape::Polygon(s) => validate(&s, &best).unwrap(),
        AnyShape::Molecule(s) => validate(&s, &best).unwrap(),
    }
    println!(
        "{}: {:.6} in {} ({})",
        best.shape,
        best.packing_fraction,
        best.group,
        best.sites_label()
    );
    assert!(
        best.packing_fraction <= 1. + 1e-6,
        "Packing fraction above 1"
    );
    best.packing_fraction
}

fn assert_tiles(found: f64) {
    assert!(
        found > 0.995,
        "Expected the shape to tile the plane, found {}",
        found
    );
}

#[test]
#[ignore]
fn long_rectangle() {
    let shape = vertices(
        "rectangle",
        &[[-10., -0.5], [10., -0.5], [10., 0.5], [-10., 0.5]],
    );
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn irregular_triangle() {
    let shape = vertices("triangle", &[[0., 0.], [3., 0.4], [0.7, 1.9]]);
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn obtuse_sliver_triangle() {
    // A very flat triangle with an angle close to 180°
    let shape = vertices("sliver", &[[-2., 0.], [2., 0.], [0.3, 0.25]]);
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn irregular_quadrilateral() {
    let shape = vertices(
        "quadrilateral",
        &[[0., 0.], [2.2, -0.3], [2.6, 1.5], [-0.4, 1.1]],
    );
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn concave_dart() {
    // A concave quadrilateral, which still tiles the plane in p2
    let shape = vertices("dart", &[[0., 0.], [2., -1.], [0.6, 0.], [2., 1.]]);
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn l_tromino() {
    let shape = vertices(
        "L-tromino",
        &[[0., 0.], [2., 0.], [2., 1.], [1., 1.], [1., 2.], [0., 2.]],
    );
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn plus_pentomino() {
    let shape = vertices(
        "plus",
        &[
            [-0.5, -1.5],
            [0.5, -1.5],
            [0.5, -0.5],
            [1.5, -0.5],
            [1.5, 0.5],
            [0.5, 0.5],
            [0.5, 1.5],
            [-0.5, 1.5],
            [-0.5, 0.5],
            [-1.5, 0.5],
            [-1.5, -0.5],
            [-0.5, -0.5],
        ],
    );
    assert_tiles(search(shape, &["p1", "p2", "p4"]));
}

#[test]
#[ignore]
fn origin_at_vertex() {
    // The shape is rotated about, and placed on sites by, a vertex rather than its centre
    let shape = vertices("corner triangle", &[[0., 0.], [1., 0.], [0.5, 0.866]]);
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn origin_outside_shape() {
    // A square whose origin lies outside the shape entirely
    let shape = vertices("offset square", &[[2., 2.], [3., 2.], [3., 3.], [2., 3.]]);
    assert_tiles(search(shape, &["p1", "p2"]));
}

#[test]
#[ignore]
fn many_sided_polygon() {
    // Close to a disk, which packs at π/√12 = 0.9069
    let found = search(ShapeSpec::Polygon { sides: 64 }, &[]);
    assert!(found > 0.905, "{}", found);
}

#[test]
#[ignore]
fn pentagram() {
    // A concave star with deep notches, packed by interlocking points
    let points: Vec<[f64; 2]> = (0..10)
        .map(|i| {
            let angle = std::f64::consts::PI / 5. * i as f64;
            let r = if i % 2 == 0 { 1. } else { 0.382 };
            [r * angle.sin(), r * angle.cos()]
        })
        .collect();
    let found = search(vertices("pentagram", &points), &[]);
    assert!(found > 0.55, "{}", found);
}

#[test]
#[ignore]
fn overlapping_trimer() {
    // Atoms overlapping heavily, so the shape is nearly a single disk
    let found = search(
        ShapeSpec::Trimer {
            radius: 0.9,
            angle: 60.,
            distance: 0.3,
        },
        &[],
    );
    assert!(found > 0.85, "{}", found);
}
