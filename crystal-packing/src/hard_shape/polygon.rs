//
// polygon.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::f64::consts::PI;

use anyhow::{bail, Error};
use serde::{Deserialize, Serialize};

use super::{HardShape, Outline, ShapeSymmetry, OVERLAP_TOLERANCE};
use crate::geometry::{Mat2, Vec2};

/// A simple polygon, which can be either convex or concave.
///
/// The polygon is decomposed into convex pieces on construction, with overlaps detected using the
/// separating axis theorem on each pair of pieces. For a convex polygon there is a single piece.
///
/// ```
/// use crystal_packing::hard_shape::{HardShape, Polygon};
/// let square = Polygon::regular(4).unwrap();
/// assert!((square.area() - 2.).abs() < 1e-12);
/// assert_eq!(square.symmetry().rotations, 4);
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Polygon {
    name: String,
    /// The vertices in anticlockwise order
    vertices: Vec<Vec2>,
    /// The indices of the vertices making up each convex piece
    pieces: Vec<Vec<usize>>,
    area: f64,
    circumradius: f64,
    inradius: f64,
    symmetry: ShapeSymmetry,
}

/// A convex piece of a polygon in world coordinates
#[derive(Debug, Clone, Default)]
pub struct PlacedPiece {
    vertices: Vec<Vec2>,
    /// Outward unit normals of each edge, used as the separating axes
    axes: Vec<Vec2>,
    /// The projection of each edge onto its normal, which is the furthest extent of the piece
    /// along that axis.
    offsets: Vec<f64>,
    /// The centre and radius of a circle enclosing the piece
    centre: Vec2,
    radius: f64,
}

#[derive(Debug, Clone, Default)]
pub struct PlacedPolygon {
    pieces: Vec<PlacedPiece>,
    scale: f64,
}

fn signed_area(vertices: &[Vec2]) -> f64 {
    let n = vertices.len();
    (0..n)
        .map(|i| vertices[i].cross(vertices[(i + 1) % n]))
        .sum::<f64>()
        * 0.5
}

fn point_in_triangle(p: Vec2, a: Vec2, b: Vec2, c: Vec2) -> bool {
    let d1 = (b - a).cross(p - a);
    let d2 = (c - b).cross(p - b);
    let d3 = (a - c).cross(p - c);
    d1 >= 0. && d2 >= 0. && d3 >= 0.
}

fn is_convex(vertices: &[Vec2]) -> bool {
    let n = vertices.len();
    (0..n).all(|i| {
        let a = vertices[i];
        let b = vertices[(i + 1) % n];
        let c = vertices[(i + 2) % n];
        (b - a).cross(c - b) >= -1e-12
    })
}

/// Triangulate a simple anticlockwise polygon using ear clipping
fn triangulate(vertices: &[Vec2]) -> Result<Vec<Vec<usize>>, Error> {
    let mut remaining: Vec<usize> = (0..vertices.len()).collect();
    let mut triangles = vec![];
    while remaining.len() > 3 {
        let n = remaining.len();
        let ear = (0..n).find(|&i| {
            let (ia, ib, ic) = (
                remaining[(i + n - 1) % n],
                remaining[i],
                remaining[(i + 1) % n],
            );
            let (a, b, c) = (vertices[ia], vertices[ib], vertices[ic]);
            if (b - a).cross(c - b) <= 0. {
                return false;
            }
            remaining
                .iter()
                .filter(|&&j| j != ia && j != ib && j != ic)
                .all(|&j| !point_in_triangle(vertices[j], a, b, c))
        });
        match ear {
            Some(i) => {
                triangles.push(vec![
                    remaining[(i + n - 1) % n],
                    remaining[i],
                    remaining[(i + 1) % n],
                ]);
                remaining.remove(i);
            }
            None => bail!("Unable to triangulate polygon, is it self intersecting?"),
        }
    }
    triangles.push(remaining);
    Ok(triangles)
}

fn segment_distance(p: Vec2, a: Vec2, b: Vec2) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.norm_sq()).clamp(0., 1.);
    (a + ab * t - p).norm()
}

fn contains_point(vertices: &[Vec2], p: Vec2) -> bool {
    let n = vertices.len();
    let mut inside = false;
    for i in 0..n {
        let a = vertices[i];
        let b = vertices[(i + 1) % n];
        if (a.y > p.y) != (b.y > p.y) {
            let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if p.x < x {
                inside = !inside;
            }
        }
    }
    inside
}

impl Polygon {
    /// Create a polygon from a list of vertices, given in either winding order.
    ///
    /// The origin of the coordinate system is the point placed on the Wyckoff site and the point
    /// about which the shape is rotated.
    pub fn from_vertices(name: &str, vertices: Vec<Vec2>) -> Result<Polygon, Error> {
        if vertices.len() < 3 {
            bail!("A polygon requires at least 3 vertices")
        }
        let mut vertices = vertices;
        let mut area = signed_area(&vertices);
        if area.abs() < 1e-14 {
            bail!("The polygon has no area")
        }
        if area < 0. {
            vertices.reverse();
            area = -area;
        }
        let pieces = if is_convex(&vertices) {
            vec![(0..vertices.len()).collect()]
        } else {
            triangulate(&vertices)?
        };
        let circumradius = vertices.iter().map(|v| v.norm()).fold(0., f64::max);
        let n = vertices.len();
        let inradius = if contains_point(&vertices, Vec2::ZERO) {
            (0..n)
                .map(|i| segment_distance(Vec2::ZERO, vertices[i], vertices[(i + 1) % n]))
                .fold(f64::INFINITY, f64::min)
        } else {
            0.
        };
        let points: Vec<_> = vertices.iter().map(|&v| (v, 0.)).collect();
        let symmetry = ShapeSymmetry::detect(&points);
        Ok(Polygon {
            name: name.to_string(),
            vertices,
            pieces,
            area,
            circumradius,
            inradius,
            symmetry,
        })
    }

    /// Create a shape from points equally spaced in angle with the given distances from the origin.
    ///
    /// The first point lies on the positive y axis, with subsequent points moving clockwise.
    pub fn from_radial(name: &str, radii: &[f64]) -> Result<Polygon, Error> {
        if radii.len() < 3 {
            bail!("The number of points provided is too few to create a 2D shape.")
        }
        let dtheta = 2. * PI / radii.len() as f64;
        let vertices = radii
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let angle = i as f64 * dtheta;
                Vec2::new(r * libm::sin(angle), r * libm::cos(angle))
            })
            .collect();
        Polygon::from_vertices(name, vertices)
    }

    /// A regular polygon with vertices a unit distance from the centre
    pub fn regular(sides: usize) -> Result<Polygon, Error> {
        Polygon::from_radial(&format!("Polygon-{}", sides), &vec![1.; sides])
    }

    pub fn vertices(&self) -> &[Vec2] {
        &self.vertices
    }

    pub fn is_convex(&self) -> bool {
        self.pieces.len() == 1
    }
}

/// Check for the overlap of two convex pieces using the separating axis theorem
#[inline]
fn pieces_overlap(a: &PlacedPiece, b: &PlacedPiece, shift: Vec2, tol: f64) -> bool {
    // For convex pieces a separating axis, if one exists, is normal to an edge of one of the
    // pieces. The furthest extent of a piece along one of its own outward edge normals is the
    // edge itself, so the axis separates the pieces when every vertex of the other piece is beyond
    // the edge. Checking stops at the first vertex which isn't, which makes the common case of an
    // axis not separating the pieces cheap.
    for (axis, offset) in a.axes.iter().zip(a.offsets.iter()) {
        let limit = offset - tol - axis.dot(shift);
        if b.vertices.iter().all(|v| axis.dot(*v) >= limit) {
            return false;
        }
    }
    for (axis, offset) in b.axes.iter().zip(b.offsets.iter()) {
        let limit = offset + axis.dot(shift) - tol;
        if a.vertices.iter().all(|v| axis.dot(*v) >= limit) {
            return false;
        }
    }
    true
}

impl HardShape for Polygon {
    type Placed = PlacedPolygon;

    fn name(&self) -> &str {
        &self.name
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn circumradius(&self) -> f64 {
        self.circumradius
    }

    fn inradius(&self) -> f64 {
        self.inradius
    }

    fn symmetry(&self) -> &ShapeSymmetry {
        &self.symmetry
    }

    fn place(&self, orientation: &Mat2, position: Vec2, out: &mut PlacedPolygon) {
        out.pieces.resize_with(self.pieces.len(), Default::default);
        out.scale = self.circumradius;
        for (indices, piece) in self.pieces.iter().zip(out.pieces.iter_mut()) {
            piece.vertices.clear();
            piece.vertices.extend(
                indices
                    .iter()
                    .map(|&i| *orientation * self.vertices[i] + position),
            );
            let n = piece.vertices.len();
            // The vertices are anticlockwise, unless the orientation includes a reflection
            let sign = if orientation.det() > 0. { -1. } else { 1. };
            piece.axes.clear();
            piece.axes.extend((0..n).map(|i| {
                (piece.vertices[(i + 1) % n] - piece.vertices[i])
                    .perp()
                    .normalised()
                    * sign
            }));
            piece.offsets.clear();
            piece.offsets.extend(
                piece
                    .axes
                    .iter()
                    .zip(piece.vertices.iter())
                    .map(|(axis, v)| axis.dot(*v)),
            );
            let centre =
                piece.vertices.iter().fold(Vec2::ZERO, |acc, &v| acc + v) * (1. / n as f64);
            piece.radius = piece
                .vertices
                .iter()
                .map(|&v| (v - centre).norm())
                .fold(0., f64::max);
            piece.centre = centre;
        }
    }

    fn overlaps(a: &PlacedPolygon, b: &PlacedPolygon, shift: Vec2) -> bool {
        let tol = OVERLAP_TOLERANCE * (a.scale + b.scale);
        for pa in a.pieces.iter() {
            for pb in b.pieces.iter() {
                let r = pa.radius + pb.radius;
                if (pb.centre + shift - pa.centre).norm_sq() >= r * r {
                    continue;
                }
                if pieces_overlap(pa, pb, shift, tol) {
                    return true;
                }
            }
        }
        false
    }

    fn outlines(&self, orientation: &Mat2, position: Vec2) -> Vec<Outline> {
        vec![Outline::Polygon(
            self.vertices
                .iter()
                .map(|&v| *orientation * v + position)
                .collect(),
        )]
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    fn placed(shape: &Polygon, angle: f64, position: Vec2) -> PlacedPolygon {
        let mut out = PlacedPolygon::default();
        shape.place(&Mat2::rotation(angle), position, &mut out);
        out
    }

    /// An independent, brute force overlap check for two simple polygons.
    ///
    /// Two polygons overlap if any of their edges properly cross or one contains a vertex of
    /// the other.
    fn brute_force_overlap(a: &[Vec2], b: &[Vec2]) -> bool {
        let cross = |p1: Vec2, p2: Vec2, q1: Vec2, q2: Vec2| {
            let d1 = (p2 - p1).cross(q1 - p1);
            let d2 = (p2 - p1).cross(q2 - p1);
            let d3 = (q2 - q1).cross(p1 - q1);
            let d4 = (q2 - q1).cross(p2 - q1);
            d1 * d2 < 0. && d3 * d4 < 0.
        };
        let n = a.len();
        let m = b.len();
        for i in 0..n {
            for j in 0..m {
                if cross(a[i], a[(i + 1) % n], b[j], b[(j + 1) % m]) {
                    return true;
                }
            }
        }
        contains_point(a, b[0]) || contains_point(b, a[0])
    }

    fn world(shape: &Polygon, angle: f64, position: Vec2) -> Vec<Vec2> {
        let r = Mat2::rotation(angle);
        shape.vertices().iter().map(|&v| r * v + position).collect()
    }

    #[test]
    fn square_area_and_radii() {
        let square = Polygon::regular(4).unwrap();
        assert_abs_diff_eq!(square.area(), 2., epsilon = 1e-12);
        assert_abs_diff_eq!(square.circumradius(), 1., epsilon = 1e-12);
        assert_abs_diff_eq!(square.inradius(), 1. / 2f64.sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn regular_area() {
        for n in 3..12 {
            let shape = Polygon::regular(n).unwrap();
            let expected = 0.5 * n as f64 * (2. * PI / n as f64).sin();
            assert_abs_diff_eq!(shape.area(), expected, epsilon = 1e-12);
            assert_eq!(shape.symmetry().rotations, n as u32);
            assert_eq!(shape.symmetry().mirrors.len(), n);
        }
    }

    #[test]
    fn touching_squares_do_not_overlap() {
        // Squares of side sqrt(2) rotated to be axis aligned
        let square = Polygon::regular(4).unwrap();
        let a = placed(&square, PI / 4., Vec2::ZERO);
        let b = placed(&square, PI / 4., Vec2::ZERO);
        let side = 2f64.sqrt();
        assert!(!Polygon::overlaps(&a, &b, Vec2::new(side, 0.)));
        assert!(!Polygon::overlaps(&a, &b, Vec2::new(side, side)));
        assert!(Polygon::overlaps(&a, &b, Vec2::new(side * 0.999, 0.)));
    }

    #[test]
    fn identical_position_overlaps() {
        let square = Polygon::regular(4).unwrap();
        let a = placed(&square, 0., Vec2::ZERO);
        assert!(Polygon::overlaps(&a, &a, Vec2::ZERO));
    }

    #[test]
    fn containment_overlaps() {
        // The line based check missed a shape completely contained within another
        let big = Polygon::from_radial("big", &[2.; 4]).unwrap();
        let small = Polygon::from_radial("small", &[0.5; 4]).unwrap();
        let a = placed(&big, 0., Vec2::ZERO);
        let b = placed(&small, 0.3, Vec2::ZERO);
        assert!(Polygon::overlaps(&a, &b, Vec2::new(0.1, 0.1)));
        assert!(Polygon::overlaps(&b, &a, Vec2::new(0.1, 0.1)));
    }

    #[test]
    fn collinear_edge_overlap() {
        // Two axis aligned squares sharing part of an edge line but overlapping in area
        let square = Polygon::regular(4).unwrap();
        let a = placed(&square, PI / 4., Vec2::ZERO);
        assert!(Polygon::overlaps(&a, &a, Vec2::new(0.5, 0.)));
    }

    #[test]
    fn concave_shape() {
        // An L shape, which is decomposed into triangles
        let l = Polygon::from_vertices(
            "L",
            vec![
                Vec2::new(0., 0.),
                Vec2::new(2., 0.),
                Vec2::new(2., 1.),
                Vec2::new(1., 1.),
                Vec2::new(1., 2.),
                Vec2::new(0., 2.),
            ],
        )
        .unwrap();
        assert!(!l.is_convex());
        assert_abs_diff_eq!(l.area(), 3., epsilon = 1e-12);
        let a = placed(&l, 0., Vec2::ZERO);
        // A small square in the notch of the L doesn't overlap
        let small = Polygon::from_vertices(
            "s",
            vec![
                Vec2::new(0., 0.),
                Vec2::new(0.8, 0.),
                Vec2::new(0.8, 0.8),
                Vec2::new(0., 0.8),
            ],
        )
        .unwrap();
        let b = placed(&small, 0., Vec2::new(1.1, 1.1));
        assert!(!Polygon::overlaps(&a, &b, Vec2::ZERO));
        assert!(Polygon::overlaps(&a, &b, Vec2::new(-0.2, 0.)));
    }

    proptest! {
        #[test]
        fn sat_matches_brute_force(
            sides_a in 3usize..9,
            sides_b in 3usize..9,
            angle_a in 0. ..(2. * PI),
            angle_b in 0. ..(2. * PI),
            x in -2.5f64..2.5,
            y in -2.5f64..2.5,
        ) {
            let a = Polygon::regular(sides_a).unwrap();
            let b = Polygon::regular(sides_b).unwrap();
            let pa = placed(&a, angle_a, Vec2::ZERO);
            let pb = placed(&b, angle_b, Vec2::ZERO);
            let shift = Vec2::new(x, y);
            let expected = brute_force_overlap(
                &world(&a, angle_a, Vec2::ZERO),
                &world(&b, angle_b, shift),
            );
            prop_assert_eq!(Polygon::overlaps(&pa, &pb, shift), expected);
        }

        #[test]
        fn sat_matches_brute_force_reflected(
            angle_a in 0. ..(2. * PI),
            angle_b in 0. ..(2. * PI),
            mirror in 0. ..PI,
            x in -2.5f64..2.5,
            y in -2.5f64..2.5,
        ) {
            // An irregular shape, reflected so the winding order of the vertices is reversed
            let shape = Polygon::from_radial("irregular", &[1., 0.6, 0.9, 0.4, 0.8]).unwrap();
            let reflected = Mat2::reflection(mirror) * Mat2::rotation(angle_b);
            let pa = placed(&shape, angle_a, Vec2::ZERO);
            let mut pb = PlacedPolygon::default();
            shape.place(&reflected, Vec2::ZERO, &mut pb);
            let shift = Vec2::new(x, y);
            let world_b: Vec<Vec2> = shape.vertices().iter().map(|&v| reflected * v + shift).collect();
            let expected = brute_force_overlap(&world(&shape, angle_a, Vec2::ZERO), &world_b);
            prop_assert_eq!(Polygon::overlaps(&pa, &pb, shift), expected);
        }

        #[test]
        fn overlap_is_symmetric(
            angle_a in 0. ..(2. * PI),
            angle_b in 0. ..(2. * PI),
            x in -2.5f64..2.5,
            y in -2.5f64..2.5,
        ) {
            let shape = Polygon::regular(5).unwrap();
            let pa = placed(&shape, angle_a, Vec2::ZERO);
            let pb = placed(&shape, angle_b, Vec2::ZERO);
            let shift = Vec2::new(x, y);
            prop_assert_eq!(
                Polygon::overlaps(&pa, &pb, shift),
                Polygon::overlaps(&pb, &pa, -shift)
            );
        }

        #[test]
        fn overlap_is_translation_invariant(
            angle_a in 0. ..(2. * PI),
            angle_b in 0. ..(2. * PI),
            x in -2.5f64..2.5,
            y in -2.5f64..2.5,
            tx in -10f64..10.,
            ty in -10f64..10.,
        ) {
            let shape = Polygon::regular(3).unwrap();
            let t = Vec2::new(tx, ty);
            let pa = placed(&shape, angle_a, Vec2::ZERO);
            let pb = placed(&shape, angle_b, Vec2::new(x, y));
            let qa = placed(&shape, angle_a, t);
            let qb = placed(&shape, angle_b, Vec2::new(x, y) + t);
            prop_assert_eq!(
                Polygon::overlaps(&pa, &pb, Vec2::ZERO),
                Polygon::overlaps(&qa, &qb, Vec2::ZERO)
            );
        }
    }
}
