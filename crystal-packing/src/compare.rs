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

use anyhow::Error;

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
