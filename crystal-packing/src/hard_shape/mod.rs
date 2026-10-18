//
// mod.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Hard shapes for the packing engine.
//!
//! A [`HardShape`] is defined in its own body frame with the origin as the point which is placed on
//! a Wyckoff site. The engine positions copies of the shape with an orientation matrix (which may
//! include a reflection) and a translation, caching the result as the shape's `Placed` type. Overlap
//! checks then operate on two placed shapes, with an additional translation applied to the second
//! so periodic images don't need to be stored separately.

use std::f64::consts::PI;
use std::fmt;

use serde::{Deserialize, Serialize};

use crate::geometry::{angle_distance, wrap_angle, Mat2, Vec2};

mod molecule;
mod polygon;

pub use molecule::{Disk, Molecule, PlacedMolecule};
pub use polygon::Polygon;

/// The relative tolerance used when deciding whether two shapes are touching or overlapping.
///
/// Shapes which overlap by less than this fraction of their size are considered to be touching,
/// allowing the optimisation to converge onto contacts without numerical noise causing rejections.
pub const OVERLAP_TOLERANCE: f64 = 1e-10;

/// The point symmetries of a shape about its origin.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ShapeSymmetry {
    /// The order of the rotational symmetry, a value of 1 is no rotational symmetry.
    pub rotations: u32,
    /// The angles of the mirror lines in the body frame, in the range [0, π)
    pub mirrors: Vec<f64>,
    /// The shape is invariant to all rotations and reflections (a disk).
    pub continuous: bool,
}

impl ShapeSymmetry {
    pub fn none() -> Self {
        ShapeSymmetry {
            rotations: 1,
            mirrors: vec![],
            continuous: false,
        }
    }

    pub fn continuous() -> Self {
        ShapeSymmetry {
            rotations: 0,
            mirrors: vec![],
            continuous: true,
        }
    }

    /// Detect the symmetries of a set of points with an associated label (e.g. a radius).
    ///
    /// The shape is symmetric to an operation when every point maps onto a point with the same
    /// label. Candidate mirror lines are the lines through each point and through the midpoint of
    /// each pair of points with the same label, which covers all possible mirrors of a finite set.
    pub fn detect(points: &[(Vec2, f64)]) -> Self {
        let scale = points
            .iter()
            .map(|(p, _)| p.norm())
            .fold(0., f64::max)
            .max(1e-12);
        let tol = 1e-8 * scale;

        let maps_onto_self = |transform: &Mat2| {
            points.iter().all(|(p, label)| {
                let q = *transform * *p;
                points
                    .iter()
                    .any(|(r, l)| (*r - q).norm() < tol && (l - label).abs() < tol)
            })
        };

        let max_order = points.len().max(1) as u32;
        let rotations = (2..=max_order)
            .rev()
            .find(|&k| maps_onto_self(&Mat2::rotation(2. * PI / k as f64)))
            .unwrap_or(1);

        let mut candidates: Vec<f64> = vec![];
        for (i, (p, label)) in points.iter().enumerate() {
            if p.norm() > tol {
                candidates.push(p.angle());
            }
            for (q, other) in points.iter().skip(i + 1) {
                if (label - other).abs() < tol {
                    let mid = (*p + *q) * 0.5;
                    if mid.norm() > tol {
                        candidates.push(mid.angle());
                    } else {
                        // The points are opposite each other, the mirror is perpendicular
                        candidates.push((*q - *p).perp().angle());
                    }
                }
            }
        }

        let mut mirrors: Vec<f64> = vec![];
        for angle in candidates {
            let angle = wrap_angle(angle, PI);
            if mirrors.iter().any(|&m| angle_distance(m, angle, PI) < 1e-8) {
                continue;
            }
            if maps_onto_self(&Mat2::reflection(angle)) {
                mirrors.push(angle);
            }
        }
        mirrors.sort_by(f64::total_cmp);

        ShapeSymmetry {
            rotations,
            mirrors,
            continuous: false,
        }
    }
}

/// An outline used for rendering shapes
#[derive(Debug, Clone, PartialEq)]
pub enum Outline {
    Polygon(Vec<Vec2>),
    Circle(Vec2, f64),
}

/// A shape which can be packed by the hard particle engine.
pub trait HardShape: Clone + Send + Sync + fmt::Debug {
    /// The cached representation of the shape at a position and orientation
    type Placed: Clone + Default + Send + Sync + fmt::Debug;

    fn name(&self) -> &str;

    fn area(&self) -> f64;

    /// The radius of the smallest origin-centred circle containing the shape
    fn circumradius(&self) -> f64;

    /// The radius of the largest origin-centred circle contained within the shape.
    ///
    /// Two shapes whose centres are closer than the sum of their inradii must overlap.
    fn inradius(&self) -> f64;

    fn symmetry(&self) -> &ShapeSymmetry;

    /// Place the shape with the given orientation and position, storing the result in `out`.
    ///
    /// This reuses the allocations in `out` so it can be called in the hot loop.
    fn place(&self, orientation: &Mat2, position: Vec2, out: &mut Self::Placed);

    /// Check whether two placed shapes overlap, translating `b` by `shift`.
    fn overlaps(a: &Self::Placed, b: &Self::Placed, shift: Vec2) -> bool;

    /// The furthest extent of a placed shape in a direction, max(p · direction) over the shape,
    /// where the direction is a unit vector.
    fn support(placed: &Self::Placed, direction: Vec2) -> f64;

    /// The smallest box containing a placed shape, as its lower left and upper right corners
    fn bounds(placed: &Self::Placed) -> (Vec2, Vec2) {
        let (x, y) = (Vec2::new(1., 0.), Vec2::new(0., 1.));
        (
            Vec2::new(-Self::support(placed, -x), -Self::support(placed, -y)),
            Vec2::new(Self::support(placed, x), Self::support(placed, y)),
        )
    }

    /// The largest distance from a point to any point of a placed shape
    fn farthest(placed: &Self::Placed, point: Vec2) -> f64;

    /// Add circles to `out` whose convex hull is the convex hull of the placed shape, as
    /// (centre, radius), where the vertices of a polygon have a radius of zero
    fn hull_circles(placed: &Self::Placed, out: &mut Vec<(Vec2, f64)>);

    /// The outlines of the shape at a position and orientation, used for rendering and for the
    /// independent validation of structures.
    fn outlines(&self, orientation: &Mat2, position: Vec2) -> Vec<Outline>;
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn regular(n: usize, offset: f64) -> Vec<(Vec2, f64)> {
        (0..n)
            .map(|i| (Vec2::from_angle(offset + 2. * PI * i as f64 / n as f64), 0.))
            .collect()
    }

    #[test]
    fn detect_square() {
        let sym = ShapeSymmetry::detect(&regular(4, 0.));
        assert_eq!(sym.rotations, 4);
        assert_eq!(sym.mirrors.len(), 4);
        assert_abs_diff_eq!(sym.mirrors[0], 0.);
        assert_abs_diff_eq!(sym.mirrors[1], PI / 4., epsilon = 1e-12);
    }

    #[test]
    fn detect_triangle() {
        let sym = ShapeSymmetry::detect(&regular(3, PI / 2.));
        assert_eq!(sym.rotations, 3);
        assert_eq!(sym.mirrors.len(), 3);
    }

    #[test]
    fn detect_pentagon() {
        let sym = ShapeSymmetry::detect(&regular(5, 0.3));
        assert_eq!(sym.rotations, 5);
        assert_eq!(sym.mirrors.len(), 5);
    }

    #[test]
    fn detect_asymmetric() {
        let points = vec![
            (Vec2::new(1., 0.), 0.),
            (Vec2::new(0., 2.), 0.),
            (Vec2::new(-1., -0.5), 0.),
        ];
        let sym = ShapeSymmetry::detect(&points);
        assert_eq!(sym.rotations, 1);
        assert!(sym.mirrors.is_empty());
    }

    #[test]
    fn detect_rectangle() {
        let points = vec![
            (Vec2::new(2., 1.), 0.),
            (Vec2::new(-2., 1.), 0.),
            (Vec2::new(-2., -1.), 0.),
            (Vec2::new(2., -1.), 0.),
        ];
        let sym = ShapeSymmetry::detect(&points);
        assert_eq!(sym.rotations, 2);
        assert_eq!(sym.mirrors.len(), 2);
        assert_abs_diff_eq!(sym.mirrors[0], 0.);
        assert_abs_diff_eq!(sym.mirrors[1], PI / 2., epsilon = 1e-12);
    }
}
