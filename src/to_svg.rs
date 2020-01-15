//
// svg.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use nalgebra::Matrix3;
use svg::node::element;
use svg::Document;

use crate::traits::*;
use crate::*;

impl ToSVG for Line2 {
    type Value = element::Path;

    fn as_svg(&self) -> Self::Value {
        let data = element::path::Data::new()
            .move_to((self.start.x, self.start.y))
            .line_by((self.end.x, self.end.y));
        element::Path::new().set("d", data)
    }
}

impl ToSVG for LJ2 {
    type Value = element::Circle;

    fn as_svg(&self) -> Self::Value {
        element::Circle::new()
            .set("r", self.sigma / 2.)
            .set("cx", self.position.x)
            .set("cy", self.position.y)
    }
}

impl ToSVG for LJShape2 {
    type Value = element::Group;

    fn as_svg(&self) -> Self::Value {
        let mut smol = element::Group::new();
        for item in self {
            smol = smol.add(item.as_svg())
        }
        smol
    }
}

impl ToSVG for LineShape {
    type Value = element::Group;

    fn as_svg(&self) -> Self::Value {
        let start = self.items[0].start;
        let mut data = element::path::Data::new().move_to((start.x, start.y));

        for item in self {
            data = data.line_to((item.end.x, item.end.y));
        }
        element::Group::new().add(element::Path::new().set("d", data.close()))
    }
}

impl ToSVG for Cell2 {
    type Value = element::Group;

    fn as_svg(&self) -> Self::Value {
        let corners = self.get_corners();
        let cell_data = element::path::Data::new()
            .move_to((corners[0].x, corners[0].y))
            .line_to((corners[1].x, corners[1].y))
            .line_to((corners[2].x, corners[2].y))
            .line_to((corners[3].x, corners[3].y))
            .close();

        element::Group::new().add(
            element::Path::new()
                .set("fill", "None")
                .set("stroke", "grey")
                .set("stroke-width", 0.1)
                .set("d", cell_data),
        )
    }
}

impl ToSVG for Atom2 {
    type Value = element::Circle;

    fn as_svg(&self) -> Self::Value {
        element::Circle::new()
            .set("r", self.radius)
            .set("cx", self.position.x)
            .set("cy", self.position.y)
    }
}

impl ToSVG for MolecularShape2 {
    type Value = element::Group;

    fn as_svg(&self) -> Self::Value {
        let mut smol = element::Group::new();
        for item in self {
            smol = smol.add(item.as_svg())
        }
        smol
    }
}

impl ToSVG for Transform2 {
    type Value = element::Use;

    fn as_svg(&self) -> Self::Value {
        let matrix: Matrix3<f64> = self.clone().into();
        element::Use::new().set(
            "transform",
            format!(
                "matrix({0} {1} {2} {3} {4} {5})",
                matrix[(0, 0)],
                matrix[(1, 0)],
                matrix[(0, 1)],
                matrix[(1, 1)],
                matrix[(0, 2)],
                matrix[(1, 2)],
            ),
        )
    }
}

impl<S> ToSVG for PotentialState<S>
where
    S: Shape + Potential,
{
    type Value = Document;

    fn as_svg(&self) -> Self::Value {
        let padding = self.shape.enclosing_radius();
        let viewbox =
            self.cell
                .get_corners()
                .iter()
                .map(|p| p * 3.)
                .fold((0., 0., 0., 0.), |acc, p| {
                    (
                        f64::min(p.x - padding, acc.0),
                        f64::min(p.y - padding, acc.1),
                        f64::max(2. * (p.x + padding), acc.2),
                        f64::max(2. * (p.y + padding), acc.3),
                    )
                });
        let mut doc = Document::new().set("viewBox", viewbox).add(
            element::Definitions::new()
                .add(self.cell.as_svg().set("id", "cell"))
                .add(self.shape.as_svg().set("id", "mol")),
        );
        for transform in self.cell.periodic_images(Transform2::identity(), 1, true) {
            doc = doc.add(transform.as_svg().set("href", "#cell"));
        }

        for position in self.relative_positions() {
            let transform = self.cell.to_cartesian_isometry(position);
            doc = doc.add(transform.as_svg().set("href", "#mol").set("fill", "blue"));
            for periodic in self.cell.periodic_images(position, 1, false) {
                doc = doc.add(periodic.as_svg().set("href", "#mol").set("fill", "green"));
            }
        }
        doc
    }
}

impl<S> ToSVG for PackedState<S>
where
    S: Shape + Intersect,
{
    type Value = Document;

    fn as_svg(&self) -> Self::Value {
        let padding = self.shape.enclosing_radius();
        let viewbox =
            self.cell
                .get_corners()
                .iter()
                .map(|p| p * 3.)
                .fold((0., 0., 0., 0.), |acc, p| {
                    (
                        f64::min(p.x - padding, acc.0),
                        f64::min(p.y - padding, acc.1),
                        f64::max(2. * (p.x + padding), acc.2),
                        f64::max(2. * (p.y + padding), acc.3),
                    )
                });
        let mut doc = Document::new().set("viewBox", viewbox).add(
            element::Definitions::new()
                .add(self.cell.as_svg().set("id", "cell"))
                .add(self.shape.as_svg().set("id", "mol")),
        );
        for transform in self.cell.periodic_images(Transform2::identity(), 1, true) {
            doc = doc.add(transform.as_svg().set("href", "#cell"));
        }
        for position in self.relative_positions() {
            let matrix = self.cell.to_cartesian_isometry(position);
            doc = doc.add(matrix.as_svg().set("href", "#mol").set("fill", "blue"));
            for transform in self.cell.periodic_images(position, 1, false) {
                doc = doc.add(transform.as_svg().set("href", "#mol").set("fill", "green"));
            }
        }
        doc
    }
}
