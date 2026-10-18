//
// validate.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Independent validation of packed structures.
//!
//! The validation shares as little code as possible with the optimisation. It rebuilds every shape
//! of a structure from its outline, and checks every pair of shapes in a supercell large enough to
//! contain all contacts using a brute force test (edges crossing or one shape containing the other)
//! rather than the separating axis theorem used in the engine. Shapes are shrunk by a small amount
//! so shapes in contact are not reported as overlapping.

use anyhow::{bail, Error};

use crate::engine::Structure;
use crate::geometry::{Mat2, Vec2};
use crate::hard_shape::{HardShape, Outline};
use crate::symmetry::{OrientationFreedom, PlaneGroup};

/// The relative amount each shape is shrunk by before checking for overlaps
const SHRINK: f64 = 1e-7;

fn segments_cross(p1: Vec2, p2: Vec2, q1: Vec2, q2: Vec2) -> bool {
    let d1 = (p2 - p1).cross(q1 - p1);
    let d2 = (p2 - p1).cross(q2 - p1);
    let d3 = (q2 - q1).cross(p1 - q1);
    let d4 = (q2 - q1).cross(p2 - q1);
    d1 * d2 < 0. && d3 * d4 < 0.
}

fn polygon_contains(polygon: &[Vec2], p: Vec2) -> bool {
    let n = polygon.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (polygon[i], polygon[(i + 1) % n]);
        if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x) {
            inside = !inside;
        }
    }
    inside
}

fn segment_distance(p: Vec2, a: Vec2, b: Vec2) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.norm_sq()).clamp(0., 1.);
    (a + ab * t - p).norm()
}

fn outlines_overlap(a: &Outline, b: &Outline) -> bool {
    match (a, b) {
        (Outline::Circle(ca, ra), Outline::Circle(cb, rb)) => (*ca - *cb).norm() < ra + rb,
        (Outline::Polygon(pa), Outline::Polygon(pb)) => {
            let (n, m) = (pa.len(), pb.len());
            for i in 0..n {
                for j in 0..m {
                    if segments_cross(pa[i], pa[(i + 1) % n], pb[j], pb[(j + 1) % m]) {
                        return true;
                    }
                }
            }
            pa.iter().any(|&p| polygon_contains(pb, p))
                || pb.iter().any(|&p| polygon_contains(pa, p))
        }
        (Outline::Polygon(p), Outline::Circle(c, r))
        | (Outline::Circle(c, r), Outline::Polygon(p)) => {
            let n = p.len();
            polygon_contains(p, *c)
                || (0..n).any(|i| segment_distance(*c, p[i], p[(i + 1) % n]) < *r)
        }
    }
}

fn shrink(outline: Outline, centre: Vec2) -> Outline {
    let factor = 1. - SHRINK;
    match outline {
        Outline::Circle(c, r) => Outline::Circle(centre + (c - centre) * factor, r * factor),
        Outline::Polygon(p) => Outline::Polygon(
            p.into_iter()
                .map(|v| centre + (v - centre) * factor)
                .collect(),
        ),
    }
}

fn translate(outline: &Outline, shift: Vec2) -> Outline {
    match outline {
        Outline::Circle(c, r) => Outline::Circle(*c + shift, *r),
        Outline::Polygon(p) => Outline::Polygon(p.iter().map(|&v| v + shift).collect()),
    }
}

fn bounding_radius(outlines: &[Outline], centre: Vec2) -> f64 {
    outlines
        .iter()
        .map(|o| match o {
            Outline::Circle(c, r) => (*c - centre).norm() + r,
            Outline::Polygon(p) => p.iter().map(|&v| (v - centre).norm()).fold(0., f64::max),
        })
        .fold(0., f64::max)
}

/// Check a structure is a valid packing of the shape, returning an error describing any problem.
pub fn validate<S: HardShape>(shape: &S, structure: &Structure) -> Result<(), Error> {
    let group = PlaneGroup::from_name(&structure.group)?;
    let lattice = &structure.lattice;
    if lattice.family != group.family || !lattice.is_valid() {
        bail!("The lattice is invalid for group {}", group.name);
    }
    let (a, b) = (lattice.a(), lattice.b());

    // Build every shape in the unit cell
    let mut shapes: Vec<(Vec2, Vec<Outline>)> = vec![];
    let mut num_shapes = 0;
    for record in structure.sites.iter() {
        let site = group.wyckoff_position(record.letter)?;
        if let Some(OrientationFreedom::Discrete(angles)) =
            site.orientation_freedom(shape.symmetry())
        {
            let period = if shape.symmetry().rotations > 0 {
                2. * std::f64::consts::PI / shape.symmetry().rotations as f64
            } else {
                2. * std::f64::consts::PI
            };
            let allowed = shape.symmetry().continuous
                || angles
                    .iter()
                    .any(|&t| crate::geometry::angle_distance(t, record.angle, period) < 1e-8);
            if !allowed {
                bail!(
                    "Site {} has an orientation which breaks its mirror symmetry",
                    record.letter
                );
            }
        } else if site.orientation_freedom(shape.symmetry()).is_none() {
            bail!("The shape can't occupy site {}", record.letter);
        }
        let params = Vec2::new(record.x, record.y);
        for (op, orientation) in site.orbit.iter().zip(site.orientations.iter()) {
            let f = op.apply(params);
            let f = Vec2::new(f.x - f.x.floor(), f.y - f.y.floor());
            let centre = a * f.x + b * f.y;
            let orientation = *orientation * Mat2::rotation(record.angle);
            let outlines = shape
                .outlines(&orientation, centre)
                .into_iter()
                .map(|o| shrink(o, centre))
                .collect();
            shapes.push((centre, outlines));
            num_shapes += 1;
        }
    }
    if num_shapes != structure.num_shapes {
        bail!(
            "Expected {} shapes, found {}",
            structure.num_shapes,
            num_shapes
        );
    }
    let phi = num_shapes as f64 * shape.area() / lattice.area();
    if (phi - structure.packing_fraction).abs() > 1e-9 {
        bail!(
            "Packing fraction {} doesn't match the structure {}",
            phi,
            structure.packing_fraction
        );
    }

    // Every pair of shapes, including periodic images, closer than the contact distance
    let radius = shapes
        .iter()
        .map(|(c, o)| bounding_radius(o, *c))
        .fold(0., f64::max);
    let contact = 2. * radius;
    // For each pair, enumerate the translations t = n u + m v of the reduced basis which bring
    // the shapes within the contact distance, using |m| ≤ |d + t| |u| / |u × v| and similarly
    // for n, where d is the separation of the shapes.
    let reduced = lattice.reduced_basis();
    let (u, v) = (reduced.column(0), reduced.column(1));
    let inverse = reduced.inverse();
    let range_n = (contact * v.norm() / lattice.area()).ceil() as i64 + 1;
    let range_m = (contact * u.norm() / lattice.area()).ceil() as i64 + 1;
    for (i, (ca, oa)) in shapes.iter().enumerate() {
        for (j, (cb, ob)) in shapes.iter().enumerate() {
            let f = inverse * (*cb - *ca);
            let (n0, m0) = (-f.x.round() as i64, -f.y.round() as i64);
            for m in (m0 - range_m)..=(m0 + range_m) {
                for n in (n0 - range_n)..=(n0 + range_n) {
                    let shift = u * n as f64 + v * m as f64;
                    if i == j && n == 0 && m == 0 {
                        continue;
                    }
                    if (*cb + shift - *ca).norm() >= contact {
                        continue;
                    }
                    for x in oa.iter() {
                        for y in ob.iter() {
                            if outlines_overlap(x, &translate(y, shift)) {
                                bail!(
                                    "Shapes {} and {} overlap with translation ({}, {})",
                                    i,
                                    j,
                                    n,
                                    m
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// Check a structure of shapes in a container is valid: every shape is inside the container and
/// no two shapes overlap. This shares no code with the optimisation beyond the shape outlines.
pub fn validate_confined<S: HardShape>(shape: &S, structure: &Structure) -> Result<(), Error> {
    let Some(container) = structure.container else {
        bail!("The structure has no container");
    };
    if structure.sites.len() != structure.num_shapes {
        bail!(
            "Expected {} shapes, found {}",
            structure.num_shapes,
            structure.sites.len()
        );
    }
    let phi = structure.num_shapes as f64 * shape.area() / container.area();
    if (phi - structure.packing_fraction).abs() > 1e-9 {
        bail!(
            "Packing fraction {} doesn't match the structure {}",
            phi,
            structure.packing_fraction
        );
    }
    let shapes: Vec<(Vec2, Vec<Outline>)> = structure
        .sites
        .iter()
        .map(|site| {
            let centre = Vec2::new(site.x, site.y);
            let outlines = shape
                .outlines(&Mat2::rotation(site.angle), centre)
                .into_iter()
                .map(|o| shrink(o, centre))
                .collect();
            (centre, outlines)
        })
        .collect();
    let (half_width, half_height) = (0.5 * container.width, 0.5 * container.height);
    let inside = |p: Vec2, r: f64| match container.kind {
        crate::container::Container::Circle => p.norm() + r <= half_width,
        _ => p.x.abs() + r <= half_width && p.y.abs() + r <= half_height,
    };
    for (i, (_, outlines)) in shapes.iter().enumerate() {
        for outline in outlines {
            let contained = match outline {
                Outline::Polygon(p) => p.iter().all(|&v| inside(v, 0.)),
                Outline::Circle(c, r) => inside(*c, *r),
            };
            if !contained {
                bail!("Shape {} is outside the container", i);
            }
        }
    }
    for (i, (_, a)) in shapes.iter().enumerate() {
        for (j, (_, b)) in shapes.iter().enumerate().skip(i + 1) {
            if a.iter().any(|x| b.iter().any(|y| outlines_overlap(x, y))) {
                bail!("Shapes {} and {} overlap", i, j);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::engine::Packing;
    use crate::hard_shape::{Molecule, Polygon};
    use crate::lattice::Lattice;
    use crate::monte_carlo::{MonteCarloConfig, Optimiser};
    use crate::random::Rng;
    use crate::symmetry::LatticeFamily;

    #[test]
    fn detects_overlap() {
        let square = Polygon::regular(4).unwrap();
        let group = PlaneGroup::from_name("p1").unwrap();
        let mut rng = Rng::seed_from_u64(0);
        let mut packing = Packing::new(square.clone(), group, &['a'], &mut rng).unwrap();
        assert!(validate(&square, &packing.structure()).is_ok());
        packing.set_lattice(Lattice::from_parameters(
            LatticeFamily::Oblique,
            &[1., 0., 1.],
        ));
        assert!(validate(&square, &packing.structure()).is_err());
    }

    #[test]
    fn optimised_structures_are_valid() {
        let cases: Vec<(&str, Vec<char>)> = vec![
            ("p2", vec!['e']),
            ("p2gg", vec!['c']),
            ("p1", vec!['a', 'a']),
            ("p2mg", vec!['d']),
            ("p6", vec!['d']),
        ];
        for (group, sites) in cases {
            let shape = Polygon::regular(5).unwrap();
            let g = PlaneGroup::from_name(group).unwrap();
            let mut rng = Rng::seed_from_u64(1);
            let packing = Packing::new(shape.clone(), g, &sites, &mut rng).unwrap();
            let config = MonteCarloConfig {
                steps: 5000,
                ..Default::default()
            };
            let mut optimiser = Optimiser::new(packing, config, 1);
            optimiser.run(5000);
            let structure = optimiser.finish();
            validate(&shape, &structure).unwrap();
        }
        let disk = Molecule::disk();
        let g = PlaneGroup::from_name("p2").unwrap();
        let mut rng = Rng::seed_from_u64(1);
        let packing = Packing::new(disk.clone(), g, &['e'], &mut rng).unwrap();
        let mut optimiser = Optimiser::new(packing, MonteCarloConfig::default(), 1);
        optimiser.run(20_000);
        validate(&disk, &optimiser.finish()).unwrap();
    }
}
