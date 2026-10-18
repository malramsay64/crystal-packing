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

use crate::container::Container;
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

/// Whether two outlines overlap: edges properly crossing, one containing a vertex of the other,
/// or circles closer than the sum of their radii. Polygons with a vertex on the edge of the other
/// can count as overlapping, which is why the validation shrinks the shapes first.
pub(crate) fn outlines_overlap(a: &Outline, b: &Outline) -> bool {
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

/// A brute force overlap check of two shapes from their outlines, the reference for the tests of
/// the overlap checks used in the optimisation
#[cfg(test)]
pub(crate) fn brute_force_overlap<S: HardShape>(
    shape: &S,
    a: (&Mat2, Vec2),
    b: (&Mat2, Vec2),
) -> bool {
    let a = shape.outlines(a.0, a.1);
    let b = shape.outlines(b.0, b.1);
    a.iter().any(|x| b.iter().any(|y| outlines_overlap(x, y)))
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

/// The area of the convex hull of a set of points, found by gift wrapping
fn wrapped_hull_area(points: &[Vec2]) -> f64 {
    let Some(start) = points
        .iter()
        .copied()
        .min_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)))
    else {
        return 0.;
    };
    let mut hull = vec![start];
    let mut current = start;
    loop {
        // The next point has every other point to its left
        let mut next = current;
        for &p in points {
            let turn = (next - current).cross(p - current);
            if next == current
                || turn < 0.
                || (turn == 0. && (p - current).norm_sq() > (next - current).norm_sq())
            {
                next = p;
            }
        }
        if next == start || next == current || hull.len() > points.len() {
            break;
        }
        hull.push(next);
        current = next;
    }
    let n = hull.len();
    0.5 * (0..n)
        .map(|i| hull[i].cross(hull[(i + 1) % n]))
        .sum::<f64>()
}

/// The area of the convex hull of shapes, with each circle replaced by many points around it
fn outline_hull_area(outlines: &[Outline]) -> f64 {
    let points: Vec<Vec2> = outlines
        .iter()
        .flat_map(|o| match o {
            Outline::Polygon(p) => p.clone(),
            Outline::Circle(c, r) => {
                let n = 2048;
                let radius = r / (std::f64::consts::PI / n as f64).cos();
                (0..n)
                    .map(|k| {
                        let t = 2. * std::f64::consts::PI * k as f64 / n as f64;
                        *c + Vec2::new(t.cos(), t.sin()) * radius
                    })
                    .collect()
            }
        })
        .collect();
    wrapped_hull_area(&points)
}

/// Check a structure of shapes in a container is valid: every shape is inside the container and
/// no two shapes overlap, including the periodic images of a strip. For a cluster, the area of the
/// convex hull is checked. This shares no code with the optimisation beyond the shape outlines.
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
    if container.kind == Container::Cluster {
        let outlines: Vec<Outline> = structure
            .sites
            .iter()
            .flat_map(|site| shape.outlines(&Mat2::rotation(site.angle), Vec2::new(site.x, site.y)))
            .collect();
        let hull = outline_hull_area(&outlines);
        if (hull - container.area()).abs() > 1e-4 * hull {
            bail!(
                "The area of the convex hull {} doesn't match the structure {}",
                hull,
                container.area()
            );
        }
    }
    let (half_width, half_height) = (0.5 * container.width, 0.5 * container.height);
    let inside = |p: Vec2, r: f64| match container.kind {
        Container::Circle => p.norm() + r <= half_width,
        Container::Square | Container::Rectangle => {
            p.x.abs() + r <= half_width && p.y.abs() + r <= half_height
        }
        Container::Strip => p.y.abs() + r <= half_height,
        Container::Cluster => true,
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
    // The periodic images along x in a strip
    let images = if container.kind == Container::Strip {
        let reach = shapes
            .iter()
            .map(|(c, o)| bounding_radius(o, *c) + c.x.abs())
            .fold(0., f64::max);
        (2. * reach / container.width).ceil() as i64 + 1
    } else {
        0
    };
    for (i, (_, a)) in shapes.iter().enumerate() {
        for (j, (_, b)) in shapes.iter().enumerate().skip(i) {
            for k in -images..=images {
                if j == i && k <= 0 {
                    continue;
                }
                let shift = Vec2::new(k as f64 * container.width, 0.);
                if a.iter()
                    .any(|x| b.iter().any(|y| outlines_overlap(x, &translate(y, shift))))
                {
                    bail!("Shapes {} and {} overlap", i, j);
                }
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
    use proptest::prelude::*;

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

    /// A packing at a random state, compressed by `scale` so overlaps are likely
    fn compressed<S: HardShape>(
        shape: S,
        group: &str,
        sites: &[char],
        seed: u64,
        scale: f64,
    ) -> Packing<S> {
        let group = PlaneGroup::from_name(group).unwrap();
        let mut rng = Rng::seed_from_u64(seed);
        let mut packing = Packing::new(shape, group, sites, &mut rng).unwrap();
        let lattice = packing.lattice().scaled(scale);
        packing.set_lattice(lattice);
        packing
    }

    fn cases() -> impl Strategy<Value = (&'static str, Vec<char>)> {
        prop::sample::select(vec![
            ("p1", vec!['a']),
            ("p1", vec!['a', 'a']),
            ("p2", vec!['e']),
            ("p2", vec!['e', 'e']),
            ("p1g1", vec!['a']),
            ("p2gg", vec!['c']),
            ("p2mg", vec!['d']),
            ("p3", vec!['d']),
            ("p4", vec!['d']),
            ("p6", vec!['d']),
        ])
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(200))]

        // The symmetry reduced overlap checks of the engine agree with the brute force checks of
        // every shape in a supercell, for convex, concave and molecular shapes
        #[test]
        fn engine_matches_validation(
            (group, sites) in cases(),
            seed in 0u64..10_000,
            scale in 0.25f64..1.,
        ) {
            let pentagon = compressed(Polygon::regular(5).unwrap(), group, &sites, seed, scale);
            let valid = validate(pentagon.shape(), &pentagon.structure());
            prop_assert_eq!(pentagon.any_overlap(), valid.is_err(), "{:?}", valid);

            let star = Polygon::from_radial("star", &[1., 0.4, 0.9, 0.5, 1., 0.3]).unwrap();
            let star = compressed(star, group, &sites, seed, scale);
            let valid = validate(star.shape(), &star.structure());
            prop_assert_eq!(star.any_overlap(), valid.is_err(), "{:?}", valid);

            let trimer = Molecule::trimer(0.7, 100., 1.2).unwrap();
            let trimer = compressed(trimer, group, &sites, seed, scale);
            let valid = validate(trimer.shape(), &trimer.structure());
            prop_assert_eq!(trimer.any_overlap(), valid.is_err(), "{:?}", valid);
        }
    }

    #[test]
    fn outline_overlaps() {
        let square = |x: f64, y: f64, side: f64| {
            Outline::Polygon(vec![
                Vec2::new(x, y),
                Vec2::new(x + side, y),
                Vec2::new(x + side, y + side),
                Vec2::new(x, y + side),
            ])
        };
        // Crossing edges, containment either way, and separate shapes
        assert!(outlines_overlap(&square(0., 0., 1.), &square(0.5, 0.5, 1.)));
        assert!(outlines_overlap(&square(0., 0., 3.), &square(1., 1., 1.)));
        assert!(outlines_overlap(&square(1., 1., 1.), &square(0., 0., 3.)));
        assert!(!outlines_overlap(&square(0., 0., 1.), &square(2., 0., 1.)));
        // Touching shapes don't overlap once shrunk about their centres, but do when pushed
        // together by much more than the shrinking
        let shrunk = |x: f64, y: f64| shrink(square(x, y, 1.), Vec2::new(x + 0.5, y + 0.5));
        assert!(!outlines_overlap(&shrunk(0., 0.), &shrunk(1., 0.)));
        assert!(!outlines_overlap(&shrunk(0., 0.), &shrunk(1., 1.)));
        assert!(!outlines_overlap(&shrunk(0., 0.), &shrunk(1., 0.5)));
        assert!(outlines_overlap(&shrunk(0., 0.), &shrunk(1. - 1e-5, 0.5)));
        // Circles closer than the sum of their radii, touching and apart
        let circle = |x: f64, r: f64| Outline::Circle(Vec2::new(x, 0.), r);
        assert!(outlines_overlap(&circle(0., 1.), &circle(1.5, 1.)));
        assert!(!outlines_overlap(&circle(0., 1.), &circle(2., 1.)));
        // A circle against a polygon: inside, crossing an edge, and outside
        assert!(outlines_overlap(&square(-1., -1., 2.), &circle(0., 0.1)));
        assert!(outlines_overlap(&square(0., 0., 1.), &circle(-0.5, 0.6)));
        assert!(!outlines_overlap(&square(0., 0., 1.), &circle(-0.5, 0.4)));
    }

    #[test]
    fn hull_by_gift_wrapping() {
        let points: Vec<Vec2> = [(0., 0.), (2., 0.), (2., 2.), (0., 2.), (1., 1.), (1., 0.)]
            .iter()
            .map(|&(x, y)| Vec2::new(x, y))
            .collect();
        assert!((wrapped_hull_area(&points) - 4.).abs() < 1e-12);
        assert_eq!(wrapped_hull_area(&[]), 0.);
        assert_eq!(wrapped_hull_area(&points[..1]), 0.);
        assert_eq!(wrapped_hull_area(&points[..2]), 0.);
    }
}
