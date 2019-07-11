//
// svg.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use svg::node::element;
use svg::Document;

use crate::traits::*;
use crate::*;

impl ToSVG for Line2 {
    type Value = element::Path;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        let data = element::path::Data::new()
            .move_to((self.start.x * scaling, self.start.y * scaling))
            .line_by((self.end.x * scaling, self.end.y * scaling));
        element::Path::new().set("d", data)
    }
}

impl ToSVG for LJ2 {
    type Value = element::Circle;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        element::Circle::new()
            .set("r", self.sigma * scaling / 2.)
            .set("cx", self.position.x * scaling)
            .set("cy", self.position.y * scaling)
            .set("fill", "green")
    }
}

impl ToSVG for LJShape2 {
    type Value = element::Group;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        let mut smol = element::Group::new();
        for item in self {
            smol = smol.add(item.as_svg(scaling))
        }
        smol
    }
}

impl ToSVG for LineShape {
    type Value = element::Group;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        let mut data = element::path::Data::new();
        for item in self {
            data = data.line_by((item.end.x * scaling, item.end.y * scaling));
        }
        element::Group::new().add(element::Path::new().set("d", data))
    }
}

impl ToSVG for Cell2 {
    type Value = element::Group;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        let corners = self.get_corners();
        let cell_data = element::path::Data::new()
            .move_to((corners[0].x * scaling, corners[0].y * scaling))
            .line_to((corners[1].x * scaling, corners[1].y * scaling))
            .line_to((corners[2].x * scaling, corners[2].y * scaling))
            .line_to((corners[3].x * scaling, corners[3].y * scaling))
            .close();

        element::Group::new().add(
            element::Path::new()
                .set("fill", "None")
                .set("stroke", "grey")
                .set("stroke-width", 0.1 * scaling)
                .set("d", cell_data),
        )
    }
}

impl ToSVG for Atom2 {
    type Value = element::Circle;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        element::Circle::new()
            .set("r", self.radius * scaling)
            .set("cx", self.position.x * scaling)
            .set("cy", self.position.y * scaling)
            .set("fill", "green")
    }
}

impl ToSVG for MolecularShape2 {
    type Value = element::Group;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        let mut smol = element::Group::new();
        for item in self {
            smol = smol.add(item.as_svg(scaling))
        }
        smol
    }
}

impl<S, C, T> ToSVG for PackedState<S, C, T>
where
    S: Shape + Intersect,
    C: Cell,
    T: Site,
{
    type Value = Document;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        let padding = self.shape.enclosing_radius();
        let viewbox = self
            .cell
            .get_corners()
            .iter()
            .fold((0., 0., 0., 0.), |acc, p| {
                (
                    f64::min((p.x - padding) * scaling, acc.0),
                    f64::min((p.y - padding) * scaling, acc.1),
                    f64::max(2. * (p.x + padding) * scaling, acc.2),
                    f64::max(2. * (p.y + padding) * scaling, acc.3),
                )
            });
        let mut doc = Document::new().set("viewBox", viewbox).add(
            element::Definitions::new()
                .add(self.cell.as_svg(scaling).set("id", "cell"))
                .add(self.shape.as_svg(scaling).set("id", "mol")),
        );
        doc = doc.add(
            element::Use::new()
                .set("href", "#cell")
                .set("transform", format!("translate({}, {})", 0, 0)),
        );
        for position in self.cartesian_positions() {
            doc = doc.add(element::Use::new().set("href", "#mol").set(
                "transform",
                format!(
                    "rotate({0}, {1}, {2}) translate({1}, {2})",
                    position.rotation.angle() * 180. / std::f64::consts::PI,
                    (position.translation.vector.x) * scaling,
                    (position.translation.vector.y) * scaling
                ),
            ));
        }
        doc
    }
}

impl<S, C, T> ToSVG for PotentialState<S, C, T>
where
    S: Shape + Potential,
    C: Cell,
    T: Site,
{
    type Value = Document;

    fn as_svg(&self, scaling: f64) -> Self::Value {
        let padding = self.shape.enclosing_radius();
        let viewbox = self
            .cell
            .get_corners()
            .iter()
            .fold((0., 0., 0., 0.), |acc, p| {
                (
                    f64::min((p.x - padding) * scaling, acc.0),
                    f64::min((p.y - padding) * scaling, acc.1),
                    f64::max(2. * (p.x + padding) * scaling, acc.2),
                    f64::max(2. * (p.y + padding) * scaling, acc.3),
                )
            });
        let mut doc = Document::new().set("viewBox", viewbox).add(
            element::Definitions::new()
                .add(self.cell.as_svg(scaling).set("id", "cell"))
                .add(self.shape.as_svg(scaling).set("id", "mol")),
        );
        doc = doc.add(
            element::Use::new()
                .set("href", "#cell")
                .set("transform", format!("translate({}, {})", 0, 0)),
        );
        for position in self.cartesian_positions() {
            doc = doc.add(element::Use::new().set("href", "#mol").set(
                "transform",
                format!(
                    "rotate({0}, {1}, {2}) translate({1}, {2})",
                    position.rotation.angle() * 180. / std::f64::consts::PI,
                    (position.translation.vector.x) * scaling,
                    (position.translation.vector.y) * scaling
                ),
            ));
        }
        doc
    }
}
