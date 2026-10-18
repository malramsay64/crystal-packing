//
// render.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Rendering packed structures as SVG images.

use svg::node::element::{Circle, Group, Line, Polygon as SvgPolygon, Rectangle};
use svg::Document;

use anyhow::Error;

use crate::container::{Confined, Container};
use crate::engine::{Packing, Structure};
use crate::geometry::{convex_hull, Vec2};
use crate::hard_shape::{HardShape, Outline};
use crate::search::{AnyShape, ShapeSpec};

/// Colours for the shapes on each occupied site
const PALETTE: [&str; 6] = [
    "#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2", "#b279a2",
];

fn points(vertices: &[Vec2]) -> String {
    vertices
        .iter()
        .map(|v| format!("{:.6},{:.6}", v.x, -v.y))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Render a packing, showing `cells` unit cells either side of the central cell.
pub fn render<S: HardShape>(packing: &Packing<S>, cells: i32) -> Document {
    let lattice = packing.lattice();
    let scale = packing.shape().circumradius();
    let stroke = format!("{:.6}", 0.02 * scale);

    let mut shapes = Group::new()
        .set("stroke", "#222")
        .set("stroke-width", stroke.clone());
    let mut min = Vec2::new(f64::INFINITY, f64::INFINITY);
    let mut max = Vec2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for (site, centre, orientation, _) in packing.placed_shapes(cells) {
        let colour = PALETTE[site % PALETTE.len()];
        for outline in packing.shape().outlines(&orientation, centre) {
            match outline {
                Outline::Polygon(vertices) => {
                    for v in vertices.iter() {
                        min = Vec2::new(min.x.min(v.x), min.y.min(v.y));
                        max = Vec2::new(max.x.max(v.x), max.y.max(v.y));
                    }
                    shapes = shapes.add(
                        SvgPolygon::new()
                            .set("points", points(&vertices))
                            .set("fill", colour),
                    );
                }
                Outline::Circle(c, r) => {
                    min = Vec2::new(min.x.min(c.x - r), min.y.min(c.y - r));
                    max = Vec2::new(max.x.max(c.x + r), max.y.max(c.y + r));
                    shapes = shapes.add(
                        Circle::new()
                            .set("cx", format!("{:.6}", c.x))
                            .set("cy", format!("{:.6}", -c.y))
                            .set("r", format!("{:.6}", r))
                            .set("fill", colour),
                    );
                }
            }
        }
    }

    let corners = [
        Vec2::ZERO,
        lattice.a(),
        lattice.a() + lattice.b(),
        lattice.b(),
    ];
    let cell = SvgPolygon::new()
        .set("points", points(&corners))
        .set("fill", "none")
        .set("stroke", "#000")
        .set("stroke-width", format!("{:.6}", 0.06 * scale));

    let margin = 0.2 * scale;
    Document::new()
        .set(
            "viewBox",
            (
                min.x - margin,
                -max.y - margin,
                max.x - min.x + 2. * margin,
                max.y - min.y + 2. * margin,
            ),
        )
        .set("width", 800)
        .add(shapes)
        .add(cell)
}

/// The points of outlines, with circles replaced by polygons, for finding their extent
fn outline_points(outlines: &[Outline]) -> Vec<Vec2> {
    outlines
        .iter()
        .flat_map(|o| match o {
            Outline::Polygon(p) => p.clone(),
            Outline::Circle(c, r) => (0..64)
                .map(|k| *c + Vec2::from_angle(std::f64::consts::PI * k as f64 / 32.) * *r)
                .collect(),
        })
        .collect()
}

/// Render shapes in a container.
///
/// A strip is drawn with the periodic images either side of the period, which is marked with a
/// dashed outline, and a cluster with the convex hull of its shapes.
pub fn render_confined<S: HardShape>(confined: &Confined<S>) -> Document {
    let scale = confined.shape().circumradius();
    let size = confined.size();
    let container = confined.container();
    let images: &[f64] = if container.is_periodic() {
        &[-1., 0., 1.]
    } else {
        &[0.]
    };
    let mut shapes = Group::new()
        .set("stroke", "#222")
        .set("stroke-width", format!("{:.6}", 0.02 * scale));
    let mut all_outlines = vec![];
    for &k in images {
        let shift = Vec2::new(k * size.width, 0.);
        for (index, centre, orientation, _) in confined.placed_shapes() {
            let colour = PALETTE[index % PALETTE.len()];
            let opacity = if k == 0. { "1" } else { "0.45" };
            for outline in confined.shape().outlines(&orientation, centre + shift) {
                shapes = match &outline {
                    Outline::Polygon(vertices) => shapes.add(
                        SvgPolygon::new()
                            .set("points", points(vertices))
                            .set("fill", colour)
                            .set("fill-opacity", opacity),
                    ),
                    Outline::Circle(c, r) => shapes.add(
                        Circle::new()
                            .set("cx", format!("{:.6}", c.x))
                            .set("cy", format!("{:.6}", -c.y))
                            .set("r", format!("{:.6}", r))
                            .set("fill", colour)
                            .set("fill-opacity", opacity),
                    ),
                };
                all_outlines.push(outline);
            }
        }
    }
    let (w, h) = (size.width, size.height);
    let stroke = format!("{:.6}", 0.06 * scale);
    let margin = 0.2 * scale;
    let extent = outline_points(&all_outlines);
    let (x0, x1, y0, y1) = match container {
        Container::Strip => (-1.5 * w, 1.5 * w, -0.5 * h, 0.5 * h),
        Container::Cluster => extent.iter().fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |(a, b, c, d), p| (a.min(p.x), b.max(p.x), c.min(p.y), d.max(p.y)),
        ),
        _ => (-0.5 * w, 0.5 * w, -0.5 * h, 0.5 * h),
    };
    let document = Document::new()
        .set(
            "viewBox",
            (
                x0 - margin,
                -y1 - margin,
                x1 - x0 + 2. * margin,
                y1 - y0 + 2. * margin,
            ),
        )
        .set("width", 800)
        .add(shapes);
    let outline = |x: f64, y: f64, width: f64, height: f64| {
        Rectangle::new()
            .set("x", format!("{:.6}", x))
            .set("y", format!("{:.6}", y))
            .set("width", format!("{:.6}", width))
            .set("height", format!("{:.6}", height))
            .set("fill", "none")
            .set("stroke", "#000")
            .set("stroke-width", stroke.clone())
    };
    match container {
        Container::Circle => document.add(
            Circle::new()
                .set("cx", 0)
                .set("cy", 0)
                .set("r", format!("{:.6}", 0.5 * w))
                .set("fill", "none")
                .set("stroke", "#000")
                .set("stroke-width", stroke.clone()),
        ),
        Container::Square | Container::Rectangle => document.add(outline(-0.5 * w, -0.5 * h, w, h)),
        Container::Strip => {
            let wall = |y: f64| {
                Line::new()
                    .set("x1", format!("{:.6}", x0 - margin))
                    .set("x2", format!("{:.6}", x1 + margin))
                    .set("y1", format!("{:.6}", y))
                    .set("y2", format!("{:.6}", y))
                    .set("stroke", "#000")
                    .set("stroke-width", stroke.clone())
            };
            document.add(wall(-0.5 * h)).add(wall(0.5 * h)).add(
                outline(-0.5 * w, -0.5 * h, w, h)
                    .set("stroke-dasharray", format!("{:.6}", 0.2 * scale))
                    .set("stroke-width", format!("{:.6}", 0.03 * scale)),
            )
        }
        Container::Cluster => document.add(
            SvgPolygon::new()
                .set("points", points(&convex_hull(&extent)))
                .set("fill", "none")
                .set("stroke", "#000")
                .set("stroke-dasharray", format!("{:.6}", 0.2 * scale))
                .set("stroke-width", format!("{:.6}", 0.03 * scale)),
        ),
    }
}

/// Render a saved structure of the described shape
pub fn render_structure(
    shape: &ShapeSpec,
    structure: &Structure,
    cells: i32,
) -> Result<Document, Error> {
    if structure.container.is_some() {
        return Ok(match shape.build()? {
            AnyShape::Polygon(s) => render_confined(&Confined::from_structure(s, structure)?),
            AnyShape::Molecule(s) => render_confined(&Confined::from_structure(s, structure)?),
        });
    }
    Ok(match shape.build()? {
        AnyShape::Polygon(s) => render(&Packing::from_structure(s, structure)?, cells),
        AnyShape::Molecule(s) => render(&Packing::from_structure(s, structure)?, cells),
    })
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::hard_shape::{Molecule, Polygon};
    use crate::random::Rng;
    use crate::symmetry::PlaneGroup;

    #[test]
    fn render_polygons() {
        let group = PlaneGroup::from_name("p2").unwrap();
        let mut rng = Rng::seed_from_u64(0);
        let packing = Packing::new(Polygon::regular(5).unwrap(), group, &['e'], &mut rng).unwrap();
        let svg = render(&packing, 1).to_string();
        // Two shapes in each of 9 cells, plus the cell outline
        assert_eq!(svg.matches("<polygon").count(), 19);
    }

    #[test]
    fn render_container() {
        let mut rng = Rng::seed_from_u64(0);
        let confined = Confined::new(Molecule::disk(), Container::Circle, 3, &mut rng).unwrap();
        let svg = render_confined(&confined).to_string();
        // Three disks and the container
        assert_eq!(svg.matches("<circle").count(), 4);
        let confined =
            Confined::new(Polygon::regular(4).unwrap(), Container::Square, 2, &mut rng).unwrap();
        let svg = render_confined(&confined).to_string();
        assert_eq!(svg.matches("<polygon").count(), 2);
        assert_eq!(svg.matches("<rect").count(), 1);
    }

    #[test]
    fn render_strip_and_cluster() {
        let mut rng = Rng::seed_from_u64(0);
        let confined = Confined::new(Molecule::disk(), Container::Strip, 2, &mut rng).unwrap();
        let svg = render_confined(&confined).to_string();
        // Two disks and their images either side, the walls and the period
        assert_eq!(svg.matches("<circle").count(), 6);
        assert_eq!(svg.matches("<line").count(), 2);
        assert_eq!(svg.matches("<rect").count(), 1);
        let confined = Confined::new(
            Polygon::regular(4).unwrap(),
            Container::Cluster,
            3,
            &mut rng,
        )
        .unwrap();
        let svg = render_confined(&confined).to_string();
        // Three squares and the hull
        assert_eq!(svg.matches("<polygon").count(), 4);
    }

    #[test]
    fn render_disks() {
        let group = PlaneGroup::from_name("p6mm").unwrap();
        let mut rng = Rng::seed_from_u64(0);
        let packing = Packing::new(Molecule::disk(), group, &['a'], &mut rng).unwrap();
        let svg = render(&packing, 1).to_string();
        assert_eq!(svg.matches("<circle").count(), 9);
    }
}
