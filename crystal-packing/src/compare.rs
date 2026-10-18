//
// compare.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Recognising when two structures are the same crystal.
//!
//! The same crystal can be found in many plane groups (a structure with p2gg symmetry also has p2
//! and p1 symmetry) and described by different unit cells, origins and orientations. Structures are
//! compared using the local environment of every shape, which doesn't depend on any of these
//! choices.

use anyhow::{bail, Error};

use crate::container::Container;
use crate::engine::Structure;
use crate::geometry::{Mat2, Vec2};
use crate::hard_shape::{HardShape, Outline};
use crate::symmetry::PlaneGroup;

/// Relative tolerance on the packing fraction for structures to be the same
const FRACTION_TOLERANCE: f64 = 1e-4;
/// Tolerance on distances, relative to the size of the shape. Optimised structures which haven't
/// fully converged can differ by around 1% of the shape size, while distinct structures with
/// packing fractions within the tolerance above differ by much more.
const DISTANCE_TOLERANCE: f64 = 1e-2;

/// A description of a structure which doesn't depend on its symmetry group, unit cell, origin or
/// orientation.
///
/// The environment of a shape is the sorted list of distances from each of its points (the
/// vertices of a polygon or the centres of the disks in a molecule) to the points of every other
/// shape nearby. A structure is described by its packing fraction and the distinct environments
/// of its shapes, along with the fraction of shapes having each environment. Distances don't
/// change under reflection, so a structure and its mirror image have the same fingerprint.
#[derive(Debug, Clone, PartialEq)]
pub struct Fingerprint {
    packing_fraction: f64,
    /// The length scale for comparing distances
    scale: f64,
    environments: Vec<(Vec<f64>, f64)>,
    /// Distances are compared up to this limit
    limit: f64,
}

fn points(outlines: &[Outline]) -> Vec<Vec2> {
    outlines
        .iter()
        .flat_map(|o| match o {
            Outline::Polygon(p) => p.clone(),
            Outline::Circle(c, _) => vec![*c],
        })
        .collect()
}

/// Compare two sorted lists of distances up to `limit`.
///
/// Distances within the tolerance of the limit are ignored, since a tiny difference between two
/// structures can move them either side of it.
fn same_environment(a: &[f64], b: &[f64], tolerance: f64, limit: f64) -> bool {
    let a: Vec<f64> = a
        .iter()
        .copied()
        .filter(|&d| d < limit - tolerance)
        .collect();
    let b: Vec<f64> = b
        .iter()
        .copied()
        .filter(|&d| d < limit - tolerance)
        .collect();
    let extra = |longer: &[f64], shorter: &[f64]| {
        // Entries beyond the shorter list must be near the limit
        longer[shorter.len()..]
            .iter()
            .all(|&d| d > limit - 3. * tolerance)
    };
    let n = a.len().min(b.len());
    let lengths_match = if a.len() > b.len() {
        extra(&a, &b)
    } else {
        extra(&b, &a)
    };
    lengths_match
        && a[..n]
            .iter()
            .zip(b[..n].iter())
            .all(|(x, y)| (x - y).abs() < tolerance)
}

impl Fingerprint {
    pub fn new<S: HardShape>(shape: &S, structure: &Structure) -> Result<Fingerprint, Error> {
        let group = PlaneGroup::from_name(&structure.group)?;
        let lattice = &structure.lattice;
        let scale = shape.circumradius();

        // The centre and points of every shape in the unit cell
        let mut shapes: Vec<(Vec2, Vec<Vec2>)> = vec![];
        for record in structure.sites.iter() {
            let site = group.wyckoff_position(record.letter)?;
            for (op, orientation) in site.orbit.iter().zip(site.orientations.iter()) {
                let f = op.apply(Vec2::new(record.x, record.y));
                let centre = lattice.to_cartesian(Vec2::new(f.x - f.x.floor(), f.y - f.y.floor()));
                let orientation = *orientation * Mat2::rotation(record.angle);
                shapes.push((centre, points(&shape.outlines(&orientation, centre))));
            }
        }

        // Distances between points up to the cutoff; those close to the cutoff are dropped since
        // a small change could move them either side of it.
        let cutoff = 2.4 * scale;
        let keep = cutoff - 0.1 * scale;
        let reduced = lattice.reduced_basis();
        let (u, v) = (reduced.column(0), reduced.column(1));
        let inverse = reduced.inverse();
        let reach = cutoff + 2. * scale;
        let range_n = (reach * v.norm() / lattice.area()).ceil() as i64 + 1;
        let range_m = (reach * u.norm() / lattice.area()).ceil() as i64 + 1;

        let mut environments: Vec<(Vec<f64>, f64)> = vec![];
        let weight = 1. / shapes.len() as f64;
        for (i, (ca, pa)) in shapes.iter().enumerate() {
            let mut distances = vec![];
            for (j, (cb, pb)) in shapes.iter().enumerate() {
                let f = inverse * (*cb - *ca);
                let (n0, m0) = (-f.x.round() as i64, -f.y.round() as i64);
                for m in (m0 - range_m)..=(m0 + range_m) {
                    for n in (n0 - range_n)..=(n0 + range_n) {
                        if i == j && n == 0 && m == 0 {
                            continue;
                        }
                        let t = u * n as f64 + v * m as f64;
                        if (*cb + t - *ca).norm() > reach {
                            continue;
                        }
                        for &x in pa.iter() {
                            for &y in pb.iter() {
                                let d = (y + t - x).norm();
                                if d < cutoff {
                                    distances.push(d);
                                }
                            }
                        }
                    }
                }
            }
            distances.sort_by(f64::total_cmp);
            let tolerance = DISTANCE_TOLERANCE * scale;
            match environments
                .iter_mut()
                .find(|(e, _)| same_environment(e, &distances, tolerance, keep))
            {
                Some((_, fraction)) => *fraction += weight,
                None => environments.push((distances, weight)),
            }
        }
        Ok(Fingerprint {
            packing_fraction: structure.packing_fraction,
            scale,
            environments,
            limit: keep,
        })
    }

    /// A fingerprint which doesn't match any other
    pub fn unique() -> Fingerprint {
        Fingerprint {
            packing_fraction: f64::NAN,
            scale: 1.,
            environments: vec![],
            limit: 0.,
        }
    }

    /// Whether two fingerprints describe the same crystal
    pub fn matches(&self, other: &Fingerprint) -> bool {
        let phi = self.packing_fraction.abs().max(1e-12);
        if (self.packing_fraction - other.packing_fraction).abs() > FRACTION_TOLERANCE * phi {
            return false;
        }
        if self.environments.len() != other.environments.len() {
            return false;
        }
        // Distances near the cutoff can fall either side of it in nearly identical structures,
        // so the comparison allows a slightly looser tolerance than the grouping.
        let tolerance = 2. * DISTANCE_TOLERANCE * self.scale;
        self.environments.iter().all(|(environment, fraction)| {
            other.environments.iter().any(|(e, f)| {
                (fraction - f).abs() < 1e-9
                    && same_environment(environment, e, tolerance, self.limit)
            })
        })
    }
}

/// Whether two structures of a shape are the same crystal
pub fn equivalent<S: HardShape>(shape: &S, a: &Structure, b: &Structure) -> Result<bool, Error> {
    Ok(Fingerprint::new(shape, a)?.matches(&Fingerprint::new(shape, b)?))
}

/// Whether two structures of shapes in a container are the same arrangement.
///
/// The containers need to be the same, and the arrangements are compared using the distances of
/// the points of the shapes from the centre of the container and from each other, which don't
/// change when the arrangement is rotated or reflected, or the shapes are listed in another order.
/// In a strip, the distance from the centre line is used, with the distances between points
/// measured to the nearest periodic image, so moving the shapes along the strip doesn't matter.
/// A cluster uses the distances from the centre of its points.
pub fn equivalent_confined<S: HardShape>(
    shape: &S,
    a: &Structure,
    b: &Structure,
) -> Result<bool, Error> {
    let (Some(ca), Some(cb)) = (a.container, b.container) else {
        bail!("Both structures need to be in a container");
    };
    let scale = shape.circumradius();
    let tolerance = DISTANCE_TOLERANCE * scale;
    let mut size_a = [ca.width, ca.height];
    let mut size_b = [cb.width, cb.height];
    // The period and width of a strip can't be swapped
    if ca.kind != Container::Strip {
        size_a.sort_by(f64::total_cmp);
        size_b.sort_by(f64::total_cmp);
    }
    if ca.kind != cb.kind
        || a.num_shapes != b.num_shapes
        || size_a
            .iter()
            .zip(size_b.iter())
            .any(|(x, y)| (x - y).abs() > tolerance)
    {
        return Ok(false);
    }
    let features = |structure: &Structure| {
        let points: Vec<Vec2> = structure
            .sites
            .iter()
            .flat_map(|site| {
                let centre = Vec2::new(site.x, site.y);
                points(&shape.outlines(&Mat2::rotation(site.angle), centre))
            })
            .collect();
        let kind = ca.kind;
        let period = structure.container.map_or(0., |c| c.width);
        let centre =
            points.iter().fold(Vec2::ZERO, |a, &p| a + p) * (1. / points.len().max(1) as f64);
        let mut radial: Vec<f64> = points
            .iter()
            .map(|p| match kind {
                Container::Strip => p.y.abs(),
                Container::Cluster => (*p - centre).norm(),
                _ => p.norm(),
            })
            .collect();
        let separation = |p: Vec2, q: Vec2| {
            let mut d = q - p;
            if kind == Container::Strip {
                d.x -= period * (d.x / period).round();
            }
            d.norm()
        };
        let mut pairs: Vec<f64> = points
            .iter()
            .enumerate()
            .flat_map(|(i, p)| points.iter().skip(i + 1).map(move |q| separation(*p, *q)))
            .collect();
        radial.sort_by(f64::total_cmp);
        pairs.sort_by(f64::total_cmp);
        (radial, pairs)
    };
    let close = |x: &[f64], y: &[f64]| {
        x.len() == y.len()
            && x.iter()
                .zip(y.iter())
                .all(|(p, q)| (p - q).abs() < tolerance)
    };
    let (radial_a, pairs_a) = features(a);
    let (radial_b, pairs_b) = features(b);
    Ok(close(&radial_a, &radial_b) && close(&pairs_a, &pairs_b))
}

#[cfg(test)]
mod confined_test {
    use super::*;
    use crate::container::{Container, ContainerState};
    use crate::engine::SiteRecord;
    use crate::hard_shape::Molecule;
    use crate::lattice::Lattice;
    use crate::symmetry::LatticeFamily;

    /// Seven unit disks in a circle of radius 3: one at the centre and six around it
    fn seven_disks(rotation: f64, reflect: bool) -> Structure {
        let mut sites = vec![SiteRecord {
            letter: 'a',
            x: 0.,
            y: 0.,
            angle: 0.,
        }];
        for k in 0..6 {
            let angle = rotation + std::f64::consts::PI / 3. * k as f64;
            let y = 2. * angle.sin();
            sites.push(SiteRecord {
                letter: 'a',
                x: 2. * angle.cos(),
                y: if reflect { -y } else { y },
                angle: 0.,
            });
        }
        Structure {
            shape: "Disk".into(),
            group: "finite".into(),
            sites,
            lattice: Lattice::from_parameters(LatticeFamily::Rectangular, &[6., 6.]),
            num_shapes: 7,
            packing_fraction: 7. / 9.,
            energy: None,
            container: Some(ContainerState {
                kind: Container::Circle,
                width: 6.,
                height: 6.,
            }),
        }
    }

    #[test]
    fn rotated_and_reflected_are_equivalent() {
        let disk = Molecule::disk();
        let a = seven_disks(0., false);
        let b = seven_disks(0.37, true);
        // The shapes can be listed in any order
        let mut c = seven_disks(1.1, false);
        c.sites.reverse();
        assert!(equivalent_confined(&disk, &a, &b).unwrap());
        assert!(equivalent_confined(&disk, &a, &c).unwrap());
    }

    /// Two disks in a zig-zag in a strip, shifted along the strip by `shift`
    fn zig_zag(shift: f64, flip: bool) -> Structure {
        let s3 = 3f64.sqrt();
        let site = |x: f64, y: f64| SiteRecord {
            letter: 'a',
            x: x + shift,
            y: if flip { -y } else { y },
            angle: 0.,
        };
        Structure {
            shape: "Disk".into(),
            group: "finite".into(),
            sites: vec![site(0., -0.5 * s3), site(1., 0.5 * s3)],
            lattice: Lattice::from_parameters(LatticeFamily::Rectangular, &[2., 2. + s3]),
            num_shapes: 2,
            packing_fraction: std::f64::consts::PI / (2. + s3),
            energy: None,
            container: Some(ContainerState {
                kind: Container::Strip,
                width: 2.,
                height: 2. + s3,
            }),
        }
    }

    #[test]
    fn strip_shifted_and_reflected_are_equivalent() {
        let disk = Molecule::disk();
        let a = zig_zag(0., false);
        assert!(equivalent_confined(&disk, &a, &zig_zag(0.7, false)).unwrap());
        assert!(equivalent_confined(&disk, &a, &zig_zag(-0.3, true)).unwrap());
        // Disks in a single row are different, even with the same container
        let mut b = a.clone();
        b.sites[1].y = 0.5;
        assert!(!equivalent_confined(&disk, &a, &b).unwrap());
    }

    #[test]
    fn cluster_moved_and_rotated_is_equivalent() {
        let disk = Molecule::disk();
        let s3 = 3f64.sqrt();
        let cluster = |centres: &[(f64, f64)], rotation: f64, offset: (f64, f64)| {
            let m = crate::geometry::Mat2::rotation(rotation);
            let sites = centres
                .iter()
                .map(|&(x, y)| {
                    let p = m * Vec2::new(x, y);
                    SiteRecord {
                        letter: 'a',
                        x: p.x + offset.0,
                        y: p.y + offset.1,
                        angle: 0.,
                    }
                })
                .collect();
            let area = 2. * s3 + 8. + std::f64::consts::PI;
            Structure {
                shape: "Disk".into(),
                group: "finite".into(),
                sites,
                lattice: Lattice::from_parameters(
                    LatticeFamily::Rectangular,
                    &[area.sqrt(), area.sqrt()],
                ),
                num_shapes: 4,
                packing_fraction: 4. * std::f64::consts::PI / area,
                energy: None,
                container: Some(ContainerState {
                    kind: Container::Cluster,
                    width: area.sqrt(),
                    height: area.sqrt(),
                }),
            }
        };
        let rhombus = [(0., 0.), (2., 0.), (1., s3), (3., s3)];
        let a = cluster(&rhombus, 0., (0., 0.));
        let b = cluster(&rhombus, 1.3, (5., -2.));
        assert!(equivalent_confined(&disk, &a, &b).unwrap());
        // The mirror image
        let other = [(0., 0.), (2., 0.), (-1., s3), (1., s3)];
        let c = cluster(&other, 0., (0., 0.));
        assert!(equivalent_confined(&disk, &a, &c).unwrap());
        let bent = [(0., 0.), (2., 0.), (1., s3), (3., s3 - 0.1)];
        let d = cluster(&bent, 0., (0., 0.));
        assert!(!equivalent_confined(&disk, &a, &d).unwrap());
    }

    #[test]
    fn different_arrangements_differ() {
        let disk = Molecule::disk();
        let a = seven_disks(0., false);
        let mut b = a.clone();
        // Move one outer disk inwards slightly
        b.sites[3].x *= 0.97;
        b.sites[3].y *= 0.97;
        assert!(!equivalent_confined(&disk, &a, &b).unwrap());
        // A different container
        let mut c = a.clone();
        c.container.as_mut().unwrap().kind = Container::Square;
        assert!(!equivalent_confined(&disk, &a, &c).unwrap());
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::engine::{Packing, Structure};
    use crate::hard_shape::{HardShape, Molecule, Polygon};
    use crate::lattice::Lattice;
    use crate::monte_carlo::{MonteCarloConfig, Optimiser};
    use crate::random::Rng;
    use crate::symmetry::{LatticeFamily, PlaneGroup};

    fn optimise<S: HardShape>(shape: S, group: &str, sites: &[char], seed: u64) -> Structure {
        let group = PlaneGroup::from_name(group).unwrap();
        let mut rng = Rng::seed_from_u64(seed);
        let packing = Packing::new(shape, group, sites, &mut rng).unwrap();
        let mut optimiser = Optimiser::new(packing, MonteCarloConfig::default(), seed);
        optimiser.run(20_000);
        optimiser.finish()
    }

    /// The best of several optimisations, so the densest structure is found
    fn best<S: HardShape>(shape: S, group: &str, sites: &[char]) -> Structure {
        (0..6)
            .map(|seed| optimise(shape.clone(), group, sites, seed))
            .max_by(|a, b| a.packing_fraction.total_cmp(&b.packing_fraction))
            .unwrap()
    }

    #[test]
    fn pentagon_double_lattice_in_different_groups() {
        let pentagon = Polygon::regular(5).unwrap();
        let p2 = best(pentagon.clone(), "p2", &['e']);
        let p2gg = best(pentagon.clone(), "p2gg", &['c']);
        let p1 = best(pentagon.clone(), "p1", &['a']);
        assert!((p2.packing_fraction - p2gg.packing_fraction).abs() < 1e-4);
        assert!(equivalent(&pentagon, &p2, &p2gg).unwrap());
        assert!(equivalent(&pentagon, &p2gg, &p2).unwrap());
        // The lattice packing is a different, less dense structure
        assert!(!equivalent(&pentagon, &p2, &p1).unwrap());
    }

    #[test]
    fn hexagonal_disks_in_different_groups() {
        let disk = Molecule::disk();
        let p6mm = best(disk.clone(), "p6mm", &['a']);
        let p1 = best(disk.clone(), "p1", &['a', 'a']);
        assert!(equivalent(&disk, &p6mm, &p1).unwrap());
    }

    #[test]
    fn invariant_to_cell_and_origin() {
        let pentagon = Polygon::regular(5).unwrap();
        let structure = best(pentagon.clone(), "p2", &['e']);
        let mut packing = Packing::from_structure(pentagon.clone(), &structure).unwrap();

        // A different basis for the same lattice, b' = b + 3a, where the fractional coordinates
        // become x' = x - 3y, followed by reduction which also rotates the structure
        let p = packing.lattice().parameters();
        packing.set_lattice(Lattice::from_parameters(
            LatticeFamily::Oblique,
            &[p[0], p[1] + 3. * p[0], p[2]],
        ));
        let mut state = packing.site_state(0);
        state.position.x -= 3. * state.position.y;
        packing.set_site(0, state);
        packing.reduce_lattice();
        let rebased = packing.structure();
        assert!(equivalent(&pentagon, &structure, &rebased).unwrap());

        // Shifting every shape by the same amount
        let mut shifted = structure.clone();
        for site in shifted.sites.iter_mut() {
            site.x += 0.123;
            site.y -= 0.31;
        }
        // In p2 the shift changes the structure unless applied to the whole crystal, so use p1
        let p1 = best(pentagon.clone(), "p1", &['a', 'a']);
        let mut moved = p1.clone();
        for site in moved.sites.iter_mut() {
            site.x += 0.123;
            site.y -= 0.31;
        }
        assert!(equivalent(&pentagon, &p1, &moved).unwrap());

        // Rotating every shape and the cell together
        let mut rotated = p1.clone();
        for site in rotated.sites.iter_mut() {
            site.angle += 0.7;
        }
        let a = rotated.lattice.a();
        let b = rotated.lattice.b();
        let rotation = crate::geometry::Mat2::rotation(0.7);
        let (ra, rb) = (rotation * a, rotation * b);
        // Express the rotated cell in upper triangular form by rotating back the basis, which is
        // the same as rotating the shapes the other way
        let angle = ra.angle();
        let back = crate::geometry::Mat2::rotation(-angle);
        let (ua, ub) = (back * ra, back * rb);
        for site in rotated.sites.iter_mut() {
            site.angle -= angle;
        }
        rotated.lattice = Lattice::from_parameters(LatticeFamily::Oblique, &[ua.x, ub.x, ub.y]);
        assert!(equivalent(&pentagon, &p1, &rotated).unwrap());
    }

    #[test]
    fn different_densities_differ() {
        let square = Polygon::regular(4).unwrap();
        let dense = best(square.clone(), "p1", &['a']);
        let mut loose = dense.clone();
        loose.lattice = loose.lattice.scaled(1.1);
        loose.packing_fraction /= 1.21;
        assert!(!equivalent(&square, &dense, &loose).unwrap());
    }

    #[test]
    fn mirror_image_of_achiral_shape() {
        // The square lattice in p4mm (with mirrors) is the same as in p4 (without). A single
        // square in p1 isn't used, since it can also tile in sheared rows.
        let square = Polygon::regular(4).unwrap();
        let p4mm = best(square.clone(), "p4mm", &['a']);
        let p4 = best(square.clone(), "p4", &['a']);
        assert!(p4mm.packing_fraction > 0.999 && p4.packing_fraction > 0.999);
        assert!(equivalent(&square, &p4mm, &p4).unwrap());
    }
}
