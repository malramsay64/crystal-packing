//
// known_geometries.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! The optimisation finds the known densest arrangements, not just their packing fractions.
//!
//! Each known structure is constructed by hand, checked to be a valid packing with the expected
//! packing fraction, then compared with the structures found by a search. Crystals are compared
//! using the environment of every shape, so the search can find them in any plane group, cell or
//! orientation; arrangements in a container are compared using distances, so they can be rotated
//! or reflected.
//!
//! ```text
//! cargo test --release -p crystal_packing --test known_geometries -- --ignored
//! ```

use std::f64::consts::PI;

use crystal_packing::compare::{equivalent, equivalent_confined};
use crystal_packing::container::{Container, ContainerState};
use crystal_packing::engine::{Packing, SiteRecord, Structure};
use crystal_packing::hard_shape::{HardShape, Molecule, Polygon};
use crystal_packing::lattice::Lattice;
use crystal_packing::search::{Search, SearchConfig, ShapeSpec, Strategy};
use crystal_packing::symmetry::LatticeFamily;
use crystal_packing::validate::{validate, validate_confined};

fn site(x: f64, y: f64, angle: f64) -> SiteRecord {
    SiteRecord {
        letter: 'a',
        x,
        y,
        angle,
    }
}

/// A crystal with every site given, checked to be a valid packing
fn crystal<S: HardShape>(
    shape: &S,
    group: &str,
    sites: Vec<SiteRecord>,
    lattice: Lattice,
) -> Structure {
    let structure = Structure {
        shape: shape.name().into(),
        group: group.into(),
        sites,
        lattice,
        num_shapes: 0,
        packing_fraction: 0.,
        energy: None,
        container: None,
    };
    // Recalculate the number of shapes and packing fraction
    let structure = Packing::from_structure(shape.clone(), &structure)
        .unwrap()
        .structure();
    validate(shape, &structure).unwrap();
    structure
}

/// Shapes at the given positions in a container, checked to be a valid packing
fn confined<S: HardShape>(
    shape: &S,
    container: Container,
    width: f64,
    height: f64,
    sites: Vec<SiteRecord>,
) -> Structure {
    let state = ContainerState {
        kind: container,
        width,
        height,
    };
    let structure = Structure {
        shape: shape.name().into(),
        group: "finite".into(),
        num_shapes: sites.len(),
        packing_fraction: sites.len() as f64 * shape.area() / state.area(),
        sites,
        lattice: Lattice::from_parameters(LatticeFamily::Rectangular, &[width, height]),
        energy: None,
        container: Some(state),
    };
    validate_confined(shape, &structure).unwrap();
    structure
}

/// The structures found by a search of every plane group within `within` of the best
fn search_crystals(shape: ShapeSpec, within: f64) -> Vec<Structure> {
    let config = SearchConfig {
        seed: 11,
        ..Default::default()
    };
    let mut search = Search::new(shape, config).unwrap();
    search.run();
    search
        .near_best(within)
        .into_iter()
        .filter_map(|r| r.best.clone())
        .collect()
}

/// The best structure for each number of shapes in a container
fn search_container(shape: ShapeSpec, container: Container, counts: Vec<usize>) -> Vec<Structure> {
    let config = SearchConfig {
        unconstrained: counts,
        container: Some(container),
        strategy: Strategy::Exhaustive,
        steps_per_shape: 2500,
        seed: 11,
        ..Default::default()
    };
    let mut search = Search::new(shape, config).unwrap();
    search.run();
    search
        .results
        .iter()
        .filter_map(|r| r.best.clone())
        .collect()
}

fn assert_found<S: HardShape>(shape: &S, found: &[Structure], known: &Structure, label: &str) {
    let matches = found
        .iter()
        .filter(|s| equivalent(shape, s, known).unwrap())
        .count();
    let summary: Vec<String> = found
        .iter()
        .map(|s| format!("{}:{} {:.6}", s.group, s.sites_label(), s.packing_fraction))
        .collect();
    println!(
        "{}: {} of {} match ({})",
        label,
        matches,
        found.len(),
        summary.join(", ")
    );
    assert!(matches > 0, "{}: the known structure wasn't found", label);
}

#[test]
#[ignore]
fn hexagonal_disks() {
    // Disks on a hexagonal lattice with a spacing of a diameter
    let disk = Molecule::disk();
    let known = crystal(
        &disk,
        "p6mm",
        vec![site(0., 0., 0.)],
        Lattice::from_parameters(LatticeFamily::Hexagonal, &[2.]),
    );
    assert!((known.packing_fraction - PI / 12f64.sqrt()).abs() < 1e-12);
    // Every structure close to the best is the hexagonal packing
    let found = search_crystals(ShapeSpec::Disk, 1e-4);
    assert_found(&disk, &found, &known, "hexagonal disks");
    assert!(found.iter().all(|s| equivalent(&disk, s, &known).unwrap()));
}

#[test]
#[ignore]
fn square_lattice() {
    // Squares with a side of √2 (unit circumradius) on a square lattice, rotated to align with it
    let square = Polygon::regular(4).unwrap();
    let side = 2f64.sqrt();
    let known = crystal(
        &square,
        "p4mm",
        vec![site(0., 0., PI / 4.)],
        Lattice::from_parameters(LatticeFamily::Square, &[side]),
    );
    assert!((known.packing_fraction - 1.).abs() < 1e-12);
    // Rows of squares can also slide past each other, so other tilings are equally dense
    let found = search_crystals(ShapeSpec::Polygon { sides: 4 }, 1e-4);
    assert_found(&square, &found, &known, "square lattice");
}

#[test]
#[ignore]
fn triangle_tiling() {
    // Triangles with unit circumradius (side √3) alternately pointing up and down, an up and a
    // down triangle sharing an edge are related by a 2 fold rotation about its midpoint
    let triangle = Polygon::regular(3).unwrap();
    let s = 3f64.sqrt();
    let known = crystal(
        &triangle,
        "p2",
        // The up triangle has its centroid at (0, 1/2) in Cartesian coordinates
        vec![SiteRecord {
            letter: 'e',
            x: -1. / 6.,
            y: 1. / 3.,
            angle: 0.,
        }],
        Lattice::from_parameters(LatticeFamily::Oblique, &[s, s / 2., 1.5]),
    );
    assert!((known.packing_fraction - 1.).abs() < 1e-12);
    let found = search_crystals(ShapeSpec::Polygon { sides: 3 }, 1e-4);
    assert_found(&triangle, &found, &known, "triangle tiling");
}

#[test]
#[ignore]
fn honeycomb() {
    // Hexagons with unit circumradius on a hexagonal lattice with a spacing of √3, sharing edges
    let hexagon = Polygon::regular(6).unwrap();
    let known = crystal(
        &hexagon,
        "p6",
        vec![site(0., 0., 0.)],
        Lattice::from_parameters(LatticeFamily::Hexagonal, &[3f64.sqrt()]),
    );
    assert!((known.packing_fraction - 1.).abs() < 1e-12);
    // The honeycomb is the only way hexagons fill space
    let found = search_crystals(ShapeSpec::Polygon { sides: 6 }, 1e-4);
    assert_found(&hexagon, &found, &known, "honeycomb");
    assert!(found
        .iter()
        .all(|s| equivalent(&hexagon, s, &known).unwrap()));
}

fn assert_container_found<S: HardShape>(
    shape: &S,
    found: &Structure,
    known: &Structure,
    label: &str,
) {
    println!(
        "{}: found {:.6}, known {:.6}",
        label, found.packing_fraction, known.packing_fraction
    );
    assert!(
        equivalent_confined(shape, found, known).unwrap(),
        "{}: found a different arrangement\n{:?}",
        label,
        found
    );
}

#[test]
#[ignore]
fn disks_in_a_circle() {
    let disk = Molecule::disk();
    // Two disks side by side across the diameter of a circle of radius 2
    let two = confined(
        &disk,
        Container::Circle,
        4.,
        4.,
        vec![site(-1., 0., 0.), site(1., 0., 0.)],
    );
    // Three disks in a triangle, their centres 2/√3 from the centre
    let r = 2. / 3f64.sqrt();
    let three = confined(
        &disk,
        Container::Circle,
        2. * (1. + r),
        2. * (1. + r),
        (0..3)
            .map(|k| {
                let angle = PI / 2. + 2. * PI / 3. * k as f64;
                site(r * angle.cos(), r * angle.sin(), 0.)
            })
            .collect(),
    );
    // Seven disks, a hexagon around a central disk, in a circle of radius 3
    let seven = confined(
        &disk,
        Container::Circle,
        6.,
        6.,
        std::iter::once(site(0., 0., 0.))
            .chain((0..6).map(|k| {
                let angle = PI / 3. * k as f64;
                site(2. * angle.cos(), 2. * angle.sin(), 0.)
            }))
            .collect(),
    );
    let found = search_container(ShapeSpec::Disk, Container::Circle, vec![2, 3, 7]);
    assert_container_found(&disk, &found[0], &two, "2 disks in a circle");
    assert_container_found(&disk, &found[1], &three, "3 disks in a circle");
    assert_container_found(&disk, &found[2], &seven, "7 disks in a circle");
}

#[test]
#[ignore]
fn disks_in_a_square() {
    let disk = Molecule::disk();
    // Two disks in opposite corners of a square of side 2 + √2
    let offset = 1. / 2f64.sqrt();
    let side = 2. + 2f64.sqrt();
    let two = confined(
        &disk,
        Container::Square,
        side,
        side,
        vec![site(-offset, -offset, 0.), site(offset, offset, 0.)],
    );
    // Four disks in a square of side 4
    let four = confined(
        &disk,
        Container::Square,
        4.,
        4.,
        vec![
            site(-1., -1., 0.),
            site(1., -1., 0.),
            site(-1., 1., 0.),
            site(1., 1., 0.),
        ],
    );
    let found = search_container(ShapeSpec::Disk, Container::Square, vec![2, 4]);
    assert_container_found(&disk, &found[0], &two, "2 disks in a square");
    assert_container_found(&disk, &found[1], &four, "4 disks in a square");
}

#[test]
#[ignore]
fn squares_in_a_rectangle() {
    // Two squares of side √2 side by side, filling a rectangle
    let square = Polygon::regular(4).unwrap();
    let side = 2f64.sqrt();
    let known = confined(
        &square,
        Container::Rectangle,
        2. * side,
        side,
        vec![site(-side / 2., 0., PI / 4.), site(side / 2., 0., PI / 4.)],
    );
    assert!((known.packing_fraction - 1.).abs() < 1e-12);
    let found = search_container(
        ShapeSpec::Polygon { sides: 4 },
        Container::Rectangle,
        vec![2],
    );
    assert_container_found(&square, &found[0], &known, "2 squares in a rectangle");
}
