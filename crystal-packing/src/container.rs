//
// container.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Packing a finite number of shapes inside a container with hard walls.
//!
//! Rather than an infinite crystal, a fixed number of shapes are placed inside a square, a
//! rectangle (with a free aspect ratio) or a circle, which is compressed around them. The packing
//! fraction is the area of the shapes divided by the area of the container. The same moves and
//! compression schedule are used as for crystals, with the cell moves replaced by changes of the
//! container size.
//!
//! Two containers have fewer walls. A strip is periodic in one direction: the shapes are repeated
//! along x with a period equal to the width, between walls at the top and bottom, so the packing
//! fraction is that of an infinite channel. A cluster has no walls or periodic boundaries at all,
//! and the packing fraction is the area of the shapes divided by the area of their convex hull.

use std::cell::Cell;
use std::f64::consts::PI;

use anyhow::{bail, Error};
use serde::{Deserialize, Serialize};

use crate::cells::CellList;
use crate::engine::{SiteRecord, Structure};
use crate::geometry::{convex_hull, perimeter, polygon_area, wrap_angle, Mat2, Vec2};
use crate::hard_shape::{HardShape, Outline, OVERLAP_TOLERANCE};
use crate::lattice::Lattice;
use crate::monte_carlo::{MonteCarloConfig, PressureSchedule, Progress, StepSize};
use crate::random::Rng;
use crate::symmetry::LatticeFamily;

/// The shape of a container
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Container {
    Square,
    /// A rectangle with a free aspect ratio
    Rectangle,
    Circle,
    /// Periodic along x with the width as the period, between walls at the top and bottom
    Strip,
    /// No walls or periodic boundaries, with the area being the convex hull of the shapes
    Cluster,
}

/// The area of the convex hull of circles of different sizes.
///
/// The hull is found from its support function h(θ), the furthest extent of any circle in the
/// direction θ, which is made of sinusoids c·u(θ) + r, one for each circle on the hull. Sweeping
/// θ around, the next circle on the hull is the first whose sinusoid rises above the current one,
/// and the area is ½∮ h² - h'² dθ, which is integrated exactly over each part of the sweep.
fn circle_hull_area(circles: &[(Vec2, f64)]) -> f64 {
    let support = |(c, r): (Vec2, f64), t: f64| c.x * libm::cos(t) + c.y * libm::sin(t) + r;
    // The circle on the hull at θ = 0, breaking ties by the one which stays furthest after it
    let Some(mut current) = (0..circles.len()).max_by(|&a, &b| {
        support(circles[a], 0.)
            .total_cmp(&support(circles[b], 0.))
            .then(circles[a].0.y.total_cmp(&circles[b].0.y))
    }) else {
        return 0.;
    };
    // The integral of h² - h'² for a single circle between two angles
    let integral = |(c, r): (Vec2, f64), t0: f64, t1: f64| {
        let (a, b) = (c.x, c.y);
        0.5 * (a * a - b * b) * (libm::sin(2. * t1) - libm::sin(2. * t0))
            - a * b * (libm::cos(2. * t1) - libm::cos(2. * t0))
            + 2. * r * (a * (libm::sin(t1) - libm::sin(t0)) - b * (libm::cos(t1) - libm::cos(t0)))
            + r * r * (t1 - t0)
    };
    let mut t = 0.;
    let mut total = 0.;
    for _ in 0..4 * circles.len() + 4 {
        let (ci, ri) = circles[current];
        // The first angle after t at which another circle rises above the current one, where
        // (c_j - c_i)·u(θ) = r_i - r_j with the difference increasing
        let mut next = (2. * PI, current);
        for (j, &(cj, rj)) in circles.iter().enumerate() {
            let d = cj - ci;
            let distance = d.norm();
            let k = ri - rj;
            if j == current || distance <= k.abs() {
                continue;
            }
            let crossing = d.angle() - libm::acos(k / distance);
            let delta = wrap_angle(crossing - t, 2. * PI);
            if delta > 1e-12 && t + delta < next.0 {
                next = (t + delta, j);
            }
        }
        total += integral(circles[current], t, next.0);
        t = next.0;
        current = next.1;
        if t >= 2. * PI {
            break;
        }
    }
    0.5 * total
}

/// The area of the convex hull of a set of circles, as (centre, radius).
///
/// The hull of circles of equal radius r is the hull of their centres grown by r, with an area of
/// A + P r + π r², from the area A and perimeter P of the hull of the centres, which also covers
/// polygons with r = 0. Circles of different sizes use the support function of the hull. Both are
/// exact.
pub fn circles_hull_area(circles: &[(Vec2, f64)]) -> f64 {
    let r0 = circles.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
    let r1 = circles.iter().map(|c| c.1).fold(0., f64::max);
    if circles.is_empty() {
        0.
    } else if r1 - r0 <= 1e-12 * r1 {
        let centres: Vec<Vec2> = circles.iter().map(|c| c.0).collect();
        let hull = convex_hull(&centres);
        polygon_area(&hull) + perimeter(&hull) * r0 + PI * r0 * r0
    } else {
        circle_hull_area(circles)
    }
}

/// The area of the convex hull of a set of outlines
pub fn hull_area(outlines: &[Outline]) -> f64 {
    let circles: Vec<(Vec2, f64)> = outlines
        .iter()
        .flat_map(|o| match o {
            Outline::Polygon(p) => p.iter().map(|&v| (v, 0.)).collect(),
            Outline::Circle(c, r) => vec![(*c, *r)],
        })
        .collect();
    circles_hull_area(&circles)
}

/// The size of a container; a circle has a diameter equal to the width and height
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ContainerSize {
    pub width: f64,
    pub height: f64,
}

/// A container and its size, as saved with a structure
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ContainerState {
    pub kind: Container,
    pub width: f64,
    pub height: f64,
}

impl ContainerState {
    pub fn size(&self) -> ContainerSize {
        ContainerSize {
            width: self.width,
            height: self.height,
        }
    }

    pub fn area(&self) -> f64 {
        self.kind.area(&self.size())
    }
}

impl Container {
    pub fn area(&self, size: &ContainerSize) -> f64 {
        match self {
            Container::Square | Container::Rectangle | Container::Strip => size.width * size.height,
            Container::Circle => PI * size.width * size.width / 4.,
            // A cluster is saved with the side of a square of the same area as its hull
            Container::Cluster => size.width * size.height,
        }
    }

    /// Keep the size consistent with the container's shape
    fn constrain(&self, size: ContainerSize) -> ContainerSize {
        match self {
            Container::Rectangle | Container::Strip => size,
            Container::Square | Container::Circle | Container::Cluster => ContainerSize {
                width: size.width,
                height: size.width,
            },
        }
    }

    /// Whether a placed shape is inside the container, allowing for the overlap tolerance
    fn contains<S: HardShape>(&self, size: &ContainerSize, placed: &S::Placed, scale: f64) -> bool {
        let tolerance = OVERLAP_TOLERANCE * scale;
        match self {
            Container::Square | Container::Rectangle => {
                let (w, h) = (0.5 * size.width + tolerance, 0.5 * size.height + tolerance);
                let (lo, hi) = S::bounds(placed);
                -lo.x <= w && hi.x <= w && -lo.y <= h && hi.y <= h
            }
            Container::Circle => S::farthest(placed, Vec2::ZERO) <= 0.5 * size.width + tolerance,
            Container::Strip => {
                let h = 0.5 * size.height + tolerance;
                let (lo, hi) = S::bounds(placed);
                -lo.y <= h && hi.y <= h
            }
            Container::Cluster => true,
        }
    }

    /// Whether the area of the container depends on the positions of the shapes
    pub fn is_cluster(&self) -> bool {
        *self == Container::Cluster
    }

    /// Whether the container is periodic along x
    pub fn is_periodic(&self) -> bool {
        *self == Container::Strip
    }
}

/// The state before the last change of the shapes in a container, so it can be undone, with the
/// area of a cluster before it
#[derive(Debug, Clone)]
enum Undo {
    None,
    Shape(usize, Vec2, f64, Option<f64>),
    /// The size, the area of a cluster, and whether the cell list was rebuilt
    Size(ContainerSize, Option<f64>, bool),
}

/// A finite number of shapes inside a container.
///
/// The positions of the shapes are stored relative to the size of the container, as fractions of
/// its half width and half height, so changing the size of the container scales the positions.
/// In a strip the x coordinate is wrapped into the period, and the size of a cluster only sets
/// the scale of the positions.
#[derive(Debug, Clone)]
pub struct Confined<S: HardShape> {
    shape: S,
    container: Container,
    size: ContainerSize,
    positions: Vec<Vec2>,
    angles: Vec<f64>,
    placed: Vec<S::Placed>,
    contact: f64,
    core: f64,
    /// The area of the convex hull of a cluster, when it is up to date
    hull: Cell<Option<f64>>,
    /// The placed shapes before the last change, swapped back to undo it
    saved: Vec<S::Placed>,
    undo: Undo,
    /// The shapes by position, in a grid covering `extent` of their positions relative to the
    /// container, so only nearby shapes are checked
    cells: CellList,
    extent: (f64, f64),
    /// The cell list before it was last rebuilt, swapped back to undo a change of size
    saved_cells: (CellList, (f64, f64)),
}

impl<S: HardShape> Confined<S> {
    /// Place `n` shapes at random positions in a container large enough that they don't overlap
    pub fn new(
        shape: S,
        container: Container,
        n: usize,
        rng: &mut Rng,
    ) -> Result<Confined<S>, Error> {
        if n == 0 {
            bail!("At least one shape is needed");
        }
        let mut positions = vec![];
        while positions.len() < n {
            let p = Vec2::new(2. * rng.centered(), 2. * rng.centered());
            // Inside the unit circle so the position is valid for every container
            if p.norm() < 0.9 {
                positions.push(p * 0.8);
            }
        }
        let angles = (0..n).map(|_| 2. * PI * rng.uniform()).collect();
        let area = n as f64 * shape.area() / 0.05;
        let width = (area
            / container.area(&ContainerSize {
                width: 1.,
                height: 1.,
            }))
        .sqrt();
        let contact = 2. * shape.circumradius();
        let core = 2. * shape.inradius() * (1. - 1e-8);
        let mut confined = Confined {
            shape,
            container,
            size: ContainerSize {
                width,
                height: width,
            },
            positions,
            angles,
            placed: vec![S::Placed::default(); n],
            contact,
            core,
            hull: Cell::new(None),
            saved: vec![S::Placed::default(); n],
            undo: Undo::None,
            cells: CellList::default(),
            extent: (2., 2.),
            saved_cells: (CellList::default(), (2., 2.)),
        };
        confined.set_size(confined.size);
        for _ in 0..200 {
            if !confined.any_overlap() {
                return Ok(confined);
            }
            let size = confined.size;
            confined.set_size(ContainerSize {
                width: size.width * 1.5,
                height: size.height * 1.5,
            });
        }
        bail!("Unable to find a starting configuration without overlaps")
    }

    /// Recreate a confined structure from a saved structure
    pub fn from_structure(shape: S, structure: &Structure) -> Result<Confined<S>, Error> {
        let Some(state) = structure.container else {
            bail!("The structure has no container");
        };
        let n = structure.sites.len();
        let mut confined = Confined {
            contact: 2. * shape.circumradius(),
            core: 2. * shape.inradius() * (1. - 1e-8),
            shape,
            container: state.kind,
            size: state.size(),
            positions: vec![Vec2::ZERO; n],
            angles: vec![0.; n],
            placed: vec![S::Placed::default(); n],
            hull: Cell::new(None),
            saved: vec![S::Placed::default(); n],
            undo: Undo::None,
            cells: CellList::default(),
            extent: (2., 2.),
            saved_cells: (CellList::default(), (2., 2.)),
        };
        for (i, site) in structure.sites.iter().enumerate() {
            confined.set_shape(i, Vec2::new(site.x, site.y), site.angle);
        }
        confined.rebuild_cells();
        Ok(confined)
    }

    pub fn shape(&self) -> &S {
        &self.shape
    }

    pub fn container(&self) -> Container {
        self.container
    }

    pub fn size(&self) -> ContainerSize {
        self.size
    }

    pub fn num_shapes(&self) -> usize {
        self.positions.len()
    }

    fn to_cartesian(&self, fractional: Vec2) -> Vec2 {
        Vec2::new(
            0.5 * fractional.x * self.size.width,
            0.5 * fractional.y * self.size.height,
        )
    }

    pub fn position(&self, index: usize) -> Vec2 {
        self.to_cartesian(self.positions[index])
    }

    pub fn angle(&self, index: usize) -> f64 {
        self.angles[index]
    }

    fn place(&mut self, index: usize) {
        let position = self.position(index);
        let orientation = Mat2::rotation(self.angles[index]);
        self.shape
            .place(&orientation, position, &mut self.placed[index]);
        self.hull.set(None);
    }

    /// Move a shape to a position (in Cartesian coordinates) and orientation. This can be undone
    /// with [`Confined::undo`].
    pub fn set_shape(&mut self, index: usize, position: Vec2, angle: f64) {
        self.undo = Undo::Shape(
            index,
            self.positions[index],
            self.angles[index],
            self.hull.get(),
        );
        std::mem::swap(&mut self.placed[index], &mut self.saved[index]);
        let mut x = 2. * position.x / self.size.width;
        if self.container.is_periodic() {
            x -= 2. * libm::floor(0.5 * (x + 1.));
        }
        self.positions[index] = Vec2::new(x, 2. * position.y / self.size.height);
        self.angles[index] = wrap_angle(angle, 2. * PI);
        self.place(index);
        self.update_cell(index);
    }

    /// Change the size of the container, scaling the positions of the shapes with it. This can be
    /// undone with [`Confined::undo`].
    pub fn set_size(&mut self, size: ContainerSize) {
        self.undo = Undo::Size(self.size, self.hull.get(), false);
        std::mem::swap(&mut self.placed, &mut self.saved);
        self.size = self.container.constrain(size);
        for index in 0..self.positions.len() {
            self.place(index);
        }
        if self.refresh_cells() {
            if let Undo::Size(.., rebuilt) = &mut self.undo {
                *rebuilt = true;
            }
        }
    }

    /// Keep the cell list up to date after a change of size, returning whether it was rebuilt,
    /// with the previous one kept to be restored by an undo.
    ///
    /// The cells are in positions relative to the size of the container, which don't change when
    /// it is resized, so the cells stay correct until they become narrower than the contact
    /// distance, or a grid becomes worth using, which happens only a few times.
    fn refresh_cells(&mut self) -> bool {
        if !CellList::enough_items(self.num_shapes()) {
            return false;
        }
        let extent = if self.cells.is_enabled() {
            self.extent
        } else {
            self.cell_extent()
        };
        let (nx, ny) = self.cell_divisions(extent);
        if self.cells.still_fits(nx, ny, self.num_shapes()) {
            return false;
        }
        std::mem::swap(&mut self.cells, &mut self.saved_cells.0);
        std::mem::swap(&mut self.extent, &mut self.saved_cells.1);
        self.rebuild_cells();
        true
    }

    /// Divide the container into a grid of cells at least as wide as the contact distance.
    ///
    /// The grid is periodic, which for a strip is the period of the strip. In a container with
    /// walls the wrapping only adds pairs of shapes on opposite sides, which can't be in contact,
    /// and a cluster, which has no size, uses the extent of its shapes, which may later spread
    /// beyond it with the same effect.
    fn rebuild_cells(&mut self) {
        self.extent = self.cell_extent();
        let (nx, ny) = self.cell_divisions(self.extent);
        let n = self.num_shapes();
        if !CellList::useful(nx, ny, n) {
            self.cells.disable();
            return;
        }
        let positions: Vec<(f64, f64)> = (0..n).map(|i| self.grid_position(i)).collect();
        self.cells.rebuild(nx, ny, &[n], |_, i| positions[i]);
    }

    /// The extent of the grid, in positions relative to the size of the container, which are
    /// fractions of its half width and half height
    fn cell_extent(&self) -> (f64, f64) {
        if !self.container.is_cluster() {
            return (2., 2.);
        }
        let inf = f64::INFINITY;
        let (lo, hi) = self.positions.iter().fold(
            (Vec2::new(inf, inf), Vec2::new(-inf, -inf)),
            |(lo, hi), p| {
                (
                    Vec2::new(lo.x.min(p.x), lo.y.min(p.y)),
                    Vec2::new(hi.x.max(p.x), hi.y.max(p.y)),
                )
            },
        );
        (
            hi.x - lo.x + 2. * self.contact / self.size.width,
            hi.y - lo.y + 2. * self.contact / self.size.height,
        )
    }

    /// The largest grid over an extent with cells at least as wide as the contact distance, or
    /// none for a small number of shapes
    fn cell_divisions(&self, extent: (f64, f64)) -> (usize, usize) {
        let n = self.num_shapes();
        if !CellList::enough_items(n) {
            return (0, 0);
        }
        let max = (2. * (n as f64).sqrt()).ceil() as usize;
        (
            CellList::divisions(0.5 * extent.0 * self.size.width, self.contact, max),
            CellList::divisions(0.5 * extent.1 * self.size.height, self.contact, max),
        )
    }

    /// The position of a shape as a fraction of the extent of the grid
    fn grid_position(&self, index: usize) -> (f64, f64) {
        let p = self.positions[index];
        (p.x / self.extent.0, p.y / self.extent.1)
    }

    fn update_cell(&mut self, index: usize) {
        if self.cells.is_enabled() {
            let (x, y) = self.grid_position(index);
            self.cells.update(0, index, x, y);
        }
    }

    /// The area of the container, or of the convex hull of a cluster
    pub fn area(&self) -> f64 {
        if self.container.is_cluster() {
            if let Some(area) = self.hull.get() {
                return area;
            }
            let mut circles = Vec::new();
            for placed in self.placed.iter() {
                S::hull_circles(placed, &mut circles);
            }
            let area = circles_hull_area(&circles);
            self.hull.set(Some(area));
            area
        } else {
            self.container.area(&self.size)
        }
    }

    /// Undo the last change of a shape or the size of the container, restoring the placed shapes
    /// from before it rather than placing them again. Only the last change can be undone.
    pub fn undo(&mut self) {
        match std::mem::replace(&mut self.undo, Undo::None) {
            Undo::None => {}
            Undo::Shape(index, position, angle, hull) => {
                self.positions[index] = position;
                self.angles[index] = angle;
                std::mem::swap(&mut self.placed[index], &mut self.saved[index]);
                self.hull.set(hull);
                self.update_cell(index);
            }
            Undo::Size(size, hull, rebuilt) => {
                self.size = size;
                std::mem::swap(&mut self.placed, &mut self.saved);
                self.hull.set(hull);
                if rebuilt {
                    std::mem::swap(&mut self.cells, &mut self.saved_cells.0);
                    std::mem::swap(&mut self.extent, &mut self.saved_cells.1);
                }
            }
        }
    }

    pub fn packing_fraction(&self) -> f64 {
        self.num_shapes() as f64 * self.shape.area() / self.area()
    }

    /// The periodic images along x which need checking, as multiples of the period
    fn images(&self) -> i64 {
        if self.container.is_periodic() {
            libm::ceil(self.contact / self.size.width) as i64 + 1
        } else {
            0
        }
    }

    /// Whether shape `j`, or any of its periodic images, overlaps shape `i`. For `i == j` this is
    /// whether the shape overlaps its own images.
    fn pair_overlaps(&self, i: usize, j: usize) -> bool {
        let images = self.images();
        (-images..=images).filter(|&k| i != j || k != 0).any(|k| {
            let shift = Vec2::new(k as f64 * self.size.width, 0.);
            let d = (self.position(j) + shift - self.position(i)).norm_sq();
            d < self.contact * self.contact
                && (d < self.core * self.core
                    || S::overlaps(&self.placed[i], &self.placed[j], shift))
        })
    }

    fn outside(&self, index: usize) -> bool {
        !self
            .container
            .contains::<S>(&self.size, &self.placed[index], self.shape.circumradius())
    }

    /// Whether a shape overlaps another shape, its own periodic images or the walls
    pub fn shape_overlaps(&self, index: usize) -> bool {
        if self.cells.is_enabled() {
            // The same checks, against only the shapes in the neighbouring cells
            return self.outside(index)
                || self.pair_overlaps(index, index)
                || self
                    .cells
                    .any_neighbour(0, index, |_, j| j != index && self.pair_overlaps(index, j));
        }
        self.outside(index) || (0..self.num_shapes()).any(|j| self.pair_overlaps(index, j))
    }

    pub fn any_overlap(&self) -> bool {
        let n = self.num_shapes();
        if self.cells.is_enabled() {
            return (0..n).any(|i| {
                self.outside(i)
                    || self.pair_overlaps(i, i)
                    || self
                        .cells
                        .any_neighbour(0, i, |_, j| j > i && self.pair_overlaps(i, j))
            });
        }
        (0..n).any(|i| self.outside(i) || (i..n).any(|j| self.pair_overlaps(i, j)))
    }

    pub fn structure(&self) -> Structure {
        let size = if self.container.is_cluster() {
            let side = self.area().sqrt();
            ContainerSize {
                width: side,
                height: side,
            }
        } else {
            self.size
        };
        Structure {
            shape: self.shape.name().to_string(),
            group: "finite".into(),
            sites: (0..self.num_shapes())
                .map(|i| {
                    let p = self.position(i);
                    SiteRecord {
                        letter: 'a',
                        x: p.x,
                        y: p.y,
                        angle: self.angles[i],
                    }
                })
                .collect(),
            lattice: Lattice::from_parameters(
                LatticeFamily::Rectangular,
                &[size.width, size.height],
            ),
            num_shapes: self.num_shapes(),
            packing_fraction: self.packing_fraction(),
            energy: None,
            container: Some(ContainerState {
                kind: self.container,
                width: size.width,
                height: size.height,
            }),
        }
    }

    /// Move a cluster so the centre of its shapes is at the origin
    fn centre(&mut self) {
        let n = self.num_shapes() as f64;
        let mean = self.positions.iter().fold(Vec2::ZERO, |a, &p| a + p) * (1. / n);
        for p in self.positions.iter_mut() {
            *p -= mean;
        }
        self.set_size(self.size);
        // Every position changed
        self.rebuild_cells();
    }

    /// Every placed shape, as (index, centre, orientation, placed shape), used for rendering
    pub fn placed_shapes(&self) -> Vec<(usize, Vec2, Mat2, S::Placed)> {
        (0..self.num_shapes())
            .map(|i| {
                (
                    i,
                    self.position(i),
                    Mat2::rotation(self.angles[i]),
                    self.placed[i].clone(),
                )
            })
            .collect()
    }
}

type Snapshot = (Vec<Vec2>, Vec<f64>, ContainerSize);

/// Isobaric Monte Carlo compression of shapes in a container
pub struct ConfinedOptimiser<S: HardShape> {
    confined: Confined<S>,
    config: MonteCarloConfig,
    schedule: PressureSchedule,
    rng: Rng,
    translate: StepSize,
    rotate: StepSize,
    resize: StepSize,
    step: u64,
    best: Snapshot,
    best_packing_fraction: f64,
}

impl<S: HardShape> ConfinedOptimiser<S> {
    pub fn new(confined: Confined<S>, config: MonteCarloConfig, seed: u64) -> Self {
        let best = (
            confined.positions.clone(),
            confined.angles.clone(),
            confined.size,
        );
        let best_packing_fraction = confined.packing_fraction();
        ConfinedOptimiser {
            confined,
            schedule: PressureSchedule::new(&config),
            config,
            rng: Rng::seed_from_u64(seed),
            translate: StepSize::new(0.2, 1e-12, 1.),
            rotate: StepSize::new(0.2, 1e-12, PI),
            resize: StepSize::new(0.02, 1e-12, 0.5),
            step: 0,
            best,
            best_packing_fraction,
        }
    }

    pub fn confined(&self) -> &Confined<S> {
        &self.confined
    }

    fn snapshot(&self) -> Snapshot {
        (
            self.confined.positions.clone(),
            self.confined.angles.clone(),
            self.confined.size,
        )
    }

    fn restore(&mut self, snapshot: &Snapshot) {
        self.confined.positions = snapshot.0.clone();
        self.confined.angles = snapshot.1.clone();
        self.confined.set_size(snapshot.2);
        // Every position changed
        self.confined.rebuild_cells();
    }

    fn update_best(&mut self) {
        let phi = self.confined.packing_fraction();
        if phi > self.best_packing_fraction {
            self.best_packing_fraction = phi;
            self.best = self.snapshot();
        }
    }

    fn pressure(&self) -> f64 {
        self.schedule.at(self.step)
    }

    /// Whether to accept a change in area, as a fraction of the area of a shape, at a fixed
    /// pressure or None for the pressure of the schedule
    fn accept_area(&mut self, change: f64, pressure: Option<f64>) -> bool {
        if change <= 0. {
            return true;
        }
        let pressure = pressure.unwrap_or_else(|| self.pressure());
        if pressure.is_infinite() {
            false
        } else {
            self.rng.uniform() < libm::exp(-pressure * change)
        }
    }

    fn move_shape(&mut self, index: usize, rotate: bool, pressure: Option<f64>) -> bool {
        let cluster = self.confined.container.is_cluster();
        let old_area = if cluster { self.confined.area() } else { 0. };
        let (old_position, old_angle) = (self.confined.position(index), self.confined.angle(index));
        let (position, angle) = if rotate {
            (
                old_position,
                old_angle + 2. * self.rotate.size * self.rng.centered(),
            )
        } else {
            let scale = 2. * self.translate.size * self.confined.shape.circumradius();
            let displacement = Vec2::new(self.rng.centered(), self.rng.centered()) * scale;
            (old_position + displacement, old_angle)
        };
        self.confined.set_shape(index, position, angle);
        if self.confined.shape_overlaps(index) {
            self.confined.undo();
            return false;
        }
        if cluster {
            // Moving a shape changes the area of the hull
            let change = (self.confined.area() - old_area) / self.confined.shape.area();
            if !self.accept_area(change, pressure) {
                self.confined.undo();
                return false;
            }
        }
        true
    }

    fn resize_container(&mut self, pressure: Option<f64>) -> bool {
        let pressure = pressure.unwrap_or_else(|| self.pressure());
        let old = self.confined.size;
        let delta = 2. * self.resize.size * self.rng.centered();
        let container = self.confined.container;
        let anisotropic = matches!(container, Container::Rectangle | Container::Strip);
        let new = match (anisotropic, self.rng.index(3)) {
            (true, 0) => ContainerSize {
                width: old.width * libm::exp(delta),
                height: old.height,
            },
            (true, 1) => ContainerSize {
                width: old.width,
                height: old.height * libm::exp(delta),
            },
            _ => ContainerSize {
                width: old.width * libm::exp(delta),
                height: old.height * libm::exp(delta),
            },
        };
        let area_old = self.confined.area();
        let area_new = if container.is_cluster() {
            // Scaling the positions of a cluster changes the area of its hull
            self.confined.set_size(new);
            self.confined.area()
        } else {
            container.area(&new)
        };
        let n = self.confined.num_shapes() as f64;
        let log_acceptance = if pressure.is_infinite() {
            if area_new <= area_old {
                0.
            } else {
                f64::NEG_INFINITY
            }
        } else {
            -pressure * (area_new - area_old) / self.confined.shape.area()
                + n * libm::log(area_new / area_old)
        };
        if log_acceptance < 0. && self.rng.uniform() >= libm::exp(log_acceptance) {
            if container.is_cluster() {
                self.confined.undo();
            }
            return false;
        }
        if !container.is_cluster() {
            self.confined.set_size(new);
        }
        if self.confined.any_overlap() {
            self.confined.undo();
            return false;
        }
        true
    }

    /// Perform a single move at a fixed pressure, or None for the pressure of the schedule
    fn step_once(&mut self, pressure: Option<f64>) {
        let n = self.confined.num_shapes();
        let interval = self.config.adapt_interval;
        let target = self.config.target_acceptance;
        let choice = self.rng.index(n + 1);
        if choice == n {
            let accepted = self.resize_container(pressure);
            self.resize.record(accepted, interval, target);
            if accepted {
                self.update_best();
            }
        } else {
            let rotate = self.rng.index(2) == 1;
            let accepted = self.move_shape(choice, rotate, pressure);
            if rotate {
                self.rotate.record(accepted, interval, target);
            } else {
                self.translate.record(accepted, interval, target);
            }
            if accepted && self.confined.container.is_cluster() {
                self.update_best();
            }
        }
    }

    pub fn run(&mut self, steps: u64) -> Progress {
        let end = (self.step + steps).min(self.config.steps);
        while self.step < end {
            self.step_once(None);
            self.step += 1;
        }
        self.progress()
    }

    pub fn progress(&self) -> Progress {
        Progress {
            steps: self.step,
            total_steps: self.config.steps,
            pressure: self.pressure(),
            packing_fraction: self.confined.packing_fraction(),
            best_packing_fraction: self.best_packing_fraction,
            best_energy: None,
        }
    }

    pub fn is_complete(&self) -> bool {
        self.step >= self.config.steps
    }

    /// Uniformly shrink the container until the shapes are in contact
    fn squeeze(&mut self) {
        let size = self.confined.size;
        let scaled = |s: f64| ContainerSize {
            width: size.width * s,
            height: size.height * s,
        };
        let mut high = 1.;
        let mut low = 0.99;
        loop {
            self.confined.set_size(scaled(low));
            if self.confined.any_overlap() || low < 1e-3 {
                break;
            }
            high = low;
            low *= 0.99;
        }
        for _ in 0..50 {
            let mid = 0.5 * (low + high);
            self.confined.set_size(scaled(mid));
            if self.confined.any_overlap() {
                low = mid;
            } else {
                high = mid;
            }
        }
        self.confined.set_size(scaled(high));
        self.update_best();
    }

    /// Compress the best structure found until the shapes are in contact, returning it
    pub fn finish(&mut self) -> Structure {
        let best = self.best.clone();
        self.restore(&best);
        if self.confined.container.is_cluster() {
            // Centre the cluster before squeezing, since moving shapes in contact afterwards can
            // make them overlap through rounding
            self.confined.centre();
            if self.confined.any_overlap() {
                self.restore(&best);
            } else {
                self.best = self.snapshot();
            }
        }
        self.squeeze();
        let rounds = 4;
        for _ in 0..rounds {
            self.translate.size = self.translate.size.min(1e-3);
            self.rotate.size = self.rotate.size.min(1e-3);
            self.resize.size = self.resize.size.min(1e-4);
            for _ in 0..self.config.squeeze_steps / rounds {
                self.step_once(Some(f64::INFINITY));
            }
            self.squeeze();
        }
        let best = self.best.clone();
        self.restore(&best);
        debug_assert!(!self.confined.any_overlap());
        self.confined.structure()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::hard_shape::{Molecule, Polygon};
    use crate::random::Rng;
    use crate::validate::validate_confined;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    fn optimise<S: HardShape>(shape: S, container: Container, n: usize, seed: u64) -> Structure {
        let mut rng = Rng::seed_from_u64(seed);
        let confined = Confined::new(shape, container, n, &mut rng).unwrap();
        let config = MonteCarloConfig {
            steps: 20_000,
            ..Default::default()
        };
        let mut optimiser = ConfinedOptimiser::new(confined, config, seed);
        optimiser.run(config.steps);
        optimiser.finish()
    }

    #[test]
    fn container_areas() {
        let size = ContainerSize {
            width: 2.,
            height: 3.,
        };
        assert_abs_diff_eq!(Container::Rectangle.area(&size), 6.);
        assert_abs_diff_eq!(
            Container::Circle.area(&ContainerSize {
                width: 2.,
                height: 2.
            }),
            PI
        );
    }

    #[test]
    fn walls_are_hard() {
        let mut rng = Rng::seed_from_u64(0);
        let mut confined = Confined::new(Molecule::disk(), Container::Circle, 1, &mut rng).unwrap();
        confined.set_size(ContainerSize {
            width: 2.,
            height: 2.,
        });
        confined.set_shape(0, Vec2::ZERO, 0.);
        assert!(!confined.shape_overlaps(0));
        confined.set_shape(0, Vec2::new(0.01, 0.), 0.);
        assert!(confined.shape_overlaps(0));

        let mut confined =
            Confined::new(Polygon::regular(4).unwrap(), Container::Square, 1, &mut rng).unwrap();
        let side = 2f64.sqrt();
        confined.set_size(ContainerSize {
            width: side,
            height: side,
        });
        confined.set_shape(0, Vec2::ZERO, PI / 4.);
        assert!(!confined.any_overlap());
        confined.set_shape(0, Vec2::ZERO, PI / 4. + 0.01);
        assert!(confined.any_overlap());
    }

    #[test]
    fn shapes_overlap_each_other() {
        let mut rng = Rng::seed_from_u64(0);
        let mut confined = Confined::new(Molecule::disk(), Container::Square, 2, &mut rng).unwrap();
        confined.set_size(ContainerSize {
            width: 10.,
            height: 10.,
        });
        confined.set_shape(0, Vec2::new(-1., 0.), 0.);
        confined.set_shape(1, Vec2::new(1., 0.), 0.);
        assert!(!confined.any_overlap());
        confined.set_shape(1, Vec2::new(0.99, 0.), 0.);
        assert!(confined.any_overlap());
        assert!(confined.shape_overlaps(1));
    }

    #[test]
    fn single_disk_fills_circle() {
        let structure = optimise(Molecule::disk(), Container::Circle, 1, 1);
        assert!(
            structure.packing_fraction > 0.999,
            "{}",
            structure.packing_fraction
        );
        validate_confined(&Molecule::disk(), &structure).unwrap();
    }

    #[test]
    fn two_squares_fill_rectangle() {
        let square = Polygon::regular(4).unwrap();
        let best = (0..4)
            .map(|seed| optimise(square.clone(), Container::Rectangle, 2, seed))
            .map(|s| s.packing_fraction)
            .fold(0., f64::max);
        assert!(best > 0.999, "{}", best);
    }

    #[test]
    fn results_are_valid_and_deterministic() {
        for (container, n) in [
            (Container::Square, 3),
            (Container::Circle, 4),
            (Container::Rectangle, 3),
            (Container::Strip, 3),
            (Container::Cluster, 3),
        ] {
            let pentagon = Polygon::regular(5).unwrap();
            let a = optimise(pentagon.clone(), container, n, 7);
            let b = optimise(pentagon.clone(), container, n, 7);
            assert_eq!(a, b);
            assert_eq!(a.num_shapes, n);
            assert_eq!(a.container.unwrap().kind, container);
            validate_confined(&pentagon, &a).unwrap();
            assert!(
                a.packing_fraction > 0.4 && a.packing_fraction < 1.,
                "{:?} {}",
                container,
                a.packing_fraction
            );
        }
    }

    #[test]
    fn validation_detects_problems() {
        let disk = Molecule::disk();
        let mut structure = optimise(disk.clone(), Container::Circle, 2, 3);
        validate_confined(&disk, &structure).unwrap();
        // Shrinking the container pushes the shapes through the wall
        let container = structure.container.as_mut().unwrap();
        container.width *= 0.95;
        container.height *= 0.95;
        structure.packing_fraction /= 0.95 * 0.95;
        assert!(validate_confined(&disk, &structure).is_err());
    }

    #[test]
    fn cell_list_kept_through_resizing() {
        let size = ContainerSize {
            width: 15.,
            height: 15.,
        };
        for container in [Container::Square, Container::Strip, Container::Cluster] {
            let mut confined = random(Molecule::disk(), container, 30, size, 0);
            assert!(confined.cells.is_enabled(), "{:?}", container);
            let before = confined.cells.clone();
            confined.set_size(ContainerSize {
                width: 14.9,
                height: 14.9,
            });
            assert_eq!(confined.cells, before, "{:?}", container);
            confined.set_size(ContainerSize {
                width: 7.,
                height: 7.,
            });
            assert_ne!(confined.cells, before, "{:?}", container);
            confined.undo();
            assert_eq!(confined.cells, before, "{:?}", container);
        }
    }

    #[test]
    fn cell_list_used_for_many_shapes() {
        let size = ContainerSize {
            width: 15.,
            height: 15.,
        };
        assert!(random(Molecule::disk(), Container::Square, 30, size, 0)
            .cells
            .is_enabled());
        assert!(!random(Molecule::disk(), Container::Square, 3, size, 0)
            .cells
            .is_enabled());
    }

    #[test]
    fn hull_of_disks() {
        let disks = |centres: &[(f64, f64)], r: f64| -> Vec<Outline> {
            centres
                .iter()
                .map(|&(x, y)| Outline::Circle(Vec2::new(x, y), r))
                .collect()
        };
        assert_abs_diff_eq!(hull_area(&disks(&[(0., 0.)], 1.)), PI, epsilon = 1e-12);
        assert_abs_diff_eq!(
            hull_area(&disks(&[(0., 0.), (2., 0.)], 1.)),
            PI + 4.,
            epsilon = 1e-12
        );
        let s3 = 3f64.sqrt();
        assert_abs_diff_eq!(
            hull_area(&disks(&[(0., 0.), (2., 0.), (1., s3)], 1.)),
            s3 + 6. + PI,
            epsilon = 1e-12
        );
        // A small disk inside a large one doesn't change the hull
        let mut outlines = disks(&[(0., 0.)], 2.);
        outlines.push(Outline::Circle(Vec2::new(0.5, 0.), 0.5));
        assert_abs_diff_eq!(hull_area(&outlines), 4. * PI, epsilon = 1e-12);
        // Two disks of radius 1 and 2 with centres 4 apart: a disk sector on each, joined by
        // tangents at an angle β = asin(1/4) to the line between the centres
        let outlines = [
            Outline::Circle(Vec2::ZERO, 2.),
            Outline::Circle(Vec2::new(4., 0.), 1.),
        ];
        let beta = (0.25f64).asin();
        let tangent = 4. * beta.cos();
        let expected = (PI + 2. * beta) * 4. / 2. + (PI - 2. * beta) / 2. + tangent * 3.;
        assert_abs_diff_eq!(hull_area(&outlines), expected, epsilon = 1e-12);
        // Two squares side by side make a rectangle
        let square = |x: f64| {
            Outline::Polygon(vec![
                Vec2::new(x, 0.),
                Vec2::new(x + 1., 0.),
                Vec2::new(x + 1., 1.),
                Vec2::new(x, 1.),
            ])
        };
        assert_abs_diff_eq!(hull_area(&[square(0.), square(1.)]), 2., epsilon = 1e-12);
        assert_abs_diff_eq!(hull_area(&[]), 0.);
    }

    fn strip(n: usize, width: f64, height: f64) -> Confined<Molecule> {
        let mut rng = Rng::seed_from_u64(0);
        let mut confined = Confined::new(Molecule::disk(), Container::Strip, n, &mut rng).unwrap();
        confined.set_size(ContainerSize { width, height });
        confined
    }

    #[test]
    fn strip_is_periodic() {
        let mut confined = strip(1, 2., 2.);
        confined.set_shape(0, Vec2::ZERO, 0.);
        assert!(!confined.any_overlap());
        assert_abs_diff_eq!(confined.packing_fraction(), PI / 4., epsilon = 1e-12);
        // The disk overlaps its own image when the period is less than its diameter
        confined.set_size(ContainerSize {
            width: 1.99,
            height: 2.,
        });
        assert!(confined.any_overlap());
        assert!(confined.shape_overlaps(0));
        // Positions wrap along x, but not y
        let mut confined = strip(1, 2., 2.);
        confined.set_shape(0, Vec2::new(1.5, 0.), 0.);
        assert_abs_diff_eq!(confined.position(0).x, -0.5, epsilon = 1e-12);
        assert!(!confined.any_overlap());
        confined.set_shape(0, Vec2::new(0., 0.01), 0.);
        assert!(confined.any_overlap());
    }

    #[test]
    fn strip_images_overlap() {
        // Two disks across the boundary of the period
        let mut confined = strip(2, 4., 2.);
        confined.set_shape(0, Vec2::new(-1.9, 0.), 0.);
        confined.set_shape(1, Vec2::new(1.9, 0.), 0.);
        assert!(confined.any_overlap());
        assert!(confined.shape_overlaps(0));
        confined.set_shape(1, Vec2::new(0.1, 0.), 0.);
        assert!(!confined.any_overlap());
        // The zig-zag arrangement of two disks per period
        let s3 = 3f64.sqrt();
        let mut confined = strip(2, 2., 2. + s3);
        confined.set_shape(0, Vec2::new(0., -0.5 * s3), 0.);
        confined.set_shape(1, Vec2::new(1., 0.5 * s3), 0.);
        assert!(!confined.any_overlap());
        assert_abs_diff_eq!(confined.packing_fraction(), PI / (2. + s3), epsilon = 1e-12);
        confined.set_size(ContainerSize {
            width: 1.99,
            height: 2. + s3,
        });
        assert!(confined.any_overlap());
    }

    #[test]
    fn cluster_has_no_walls() {
        let mut rng = Rng::seed_from_u64(0);
        let mut confined =
            Confined::new(Molecule::disk(), Container::Cluster, 2, &mut rng).unwrap();
        confined.set_shape(0, Vec2::new(-1., 0.), 0.);
        confined.set_shape(1, Vec2::new(1., 0.), 0.);
        assert!(!confined.any_overlap());
        assert_abs_diff_eq!(
            confined.packing_fraction(),
            2. * PI / (PI + 4.),
            epsilon = 1e-12
        );
        // Far from the origin is fine
        confined.set_shape(0, Vec2::new(1000., 0.), 0.);
        assert!(!confined.any_overlap());
        confined.set_shape(0, Vec2::new(-0.99, 0.), 0.);
        assert!(confined.any_overlap());
        // The saved size is a square of the same area as the hull
        confined.set_shape(0, Vec2::new(-1., 0.), 0.);
        let structure = confined.structure();
        assert_abs_diff_eq!(structure.container.unwrap().area(), PI + 4., epsilon = 1e-9);
        validate_confined(&Molecule::disk(), &structure).unwrap();
    }

    #[test]
    fn strip_and_cluster_of_disks() {
        let disk = Molecule::disk();
        let best = |container, n| {
            (0..4)
                .map(|seed| optimise(disk.clone(), container, n, seed))
                .map(|s| {
                    validate_confined(&disk, &s).unwrap();
                    s.packing_fraction
                })
                .fold(0., f64::max)
        };
        let strip = best(Container::Strip, 1);
        assert!(
            strip > PI / 4. - 1e-3 && strip <= PI / 4. + 1e-9,
            "{}",
            strip
        );
        let cluster = best(Container::Cluster, 2);
        let expected = 2. * PI / (PI + 4.);
        assert!(
            cluster > expected - 1e-3 && cluster <= expected + 1e-9,
            "{}",
            cluster
        );
    }

    #[test]
    fn validation_detects_strip_and_cluster_problems() {
        let disk = Molecule::disk();
        let mut structure = optimise(disk.clone(), Container::Strip, 2, 3);
        validate_confined(&disk, &structure).unwrap();
        // Shortening the period makes the shapes overlap their images
        let container = structure.container.as_mut().unwrap();
        container.width *= 0.9;
        structure.packing_fraction /= 0.9;
        assert!(validate_confined(&disk, &structure).is_err());

        let mut structure = optimise(disk.clone(), Container::Cluster, 3, 3);
        validate_confined(&disk, &structure).unwrap();
        // The area doesn't match the hull of the shapes
        let container = structure.container.as_mut().unwrap();
        container.width *= 0.99;
        structure.packing_fraction /= 0.99;
        assert!(validate_confined(&disk, &structure).is_err());
    }

    /// Shapes at random positions and orientations in a container of the given size
    fn random<S: HardShape>(
        shape: S,
        container: Container,
        n: usize,
        size: ContainerSize,
        seed: u64,
    ) -> Confined<S> {
        let mut rng = Rng::seed_from_u64(seed);
        let mut confined = Confined::new(shape, container, n, &mut rng).unwrap();
        confined.set_size(size);
        let size = confined.size();
        for i in 0..n {
            let position = Vec2::new(
                0.6 * size.width * rng.centered(),
                0.6 * size.height * rng.centered(),
            );
            confined.set_shape(i, position, 2. * PI * rng.uniform());
        }
        confined
    }

    /// The overlap checks agree with each other and with the independent validation
    fn check_overlaps<S: HardShape>(confined: &Confined<S>) -> Result<(), TestCaseError> {
        let overlap = confined.any_overlap();
        let valid = validate_confined(confined.shape(), &confined.structure());
        prop_assert_eq!(overlap, valid.is_err(), "{:?}", valid);
        let any_shape = (0..confined.num_shapes()).any(|i| confined.shape_overlaps(i));
        prop_assert_eq!(overlap, any_shape);
        Ok(())
    }

    fn containers() -> impl Strategy<Value = Container> {
        prop::sample::select(vec![
            Container::Square,
            Container::Rectangle,
            Container::Circle,
            Container::Strip,
            Container::Cluster,
        ])
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(200))]

        #[test]
        fn polygon_overlaps_match_validation(
            container in containers(),
            n in 1usize..5,
            width in 1.5f64..6.,
            height in 1.5f64..6.,
            seed in 0u64..10_000,
        ) {
            // A concave shape, so the overlaps of its pieces and their images are tested
            let shape = Polygon::from_radial("star", &[1., 0.5, 0.9, 0.4, 1., 0.6]).unwrap();
            let size = ContainerSize { width, height };
            check_overlaps(&random(shape, container, n, size, seed))?;
        }

        #[test]
        fn molecule_overlaps_match_validation(
            container in containers(),
            n in 1usize..5,
            width in 1.5f64..8.,
            height in 1.5f64..8.,
            seed in 0u64..10_000,
        ) {
            let trimer = Molecule::trimer(0.6, 120., 1.).unwrap();
            let size = ContainerSize { width, height };
            check_overlaps(&random(trimer, container, n, size, seed))?;
        }

        #[test]
        fn cell_list_matches_validation(
            container in containers(),
            n in 20usize..40,
            width in 10f64..20.,
            height in 10f64..20.,
            seed in 0u64..10_000,
        ) {
            // Many shapes in a large container, so the overlaps are found with the cell list
            let mut confined = random(Molecule::disk(), container, n, ContainerSize { width, height }, seed);
            prop_assume!(confined.cells.is_enabled());
            check_overlaps(&confined)?;
            // Moves, undoing them, and resizing all keep the cells up to date
            let mut rng = Rng::seed_from_u64(seed);
            for step in 0..20 {
                if step % 5 == 4 {
                    let size = confined.size();
                    let scale = 0.9 + 0.2 * rng.uniform();
                    confined.set_size(ContainerSize { width: size.width * scale, height: size.height * scale });
                } else {
                    let i = rng.index(n);
                    let p = confined.position(i) + Vec2::new(rng.centered(), rng.centered()) * 4.;
                    confined.set_shape(i, p, rng.uniform());
                }
                if rng.index(2) == 0 {
                    confined.undo();
                }
                check_overlaps(&confined)?;
            }
        }

        #[test]
        fn narrow_strips_match_validation(
            n in 1usize..4,
            period in 0.3f64..3.,
            height in 1.5f64..4.,
            seed in 0u64..10_000,
        ) {
            // A period shorter than the shapes, so shapes overlap images several periods away
            let shape = Polygon::regular(5).unwrap();
            let size = ContainerSize { width: period, height };
            check_overlaps(&random(shape, Container::Strip, n, size, seed))?;
        }
    }

    proptest! {
        #[test]
        fn circle_hull_matches_sampled_hull(
            circles in prop::collection::vec((-3f64..3., -3f64..3., 0f64..1.5), 1..8),
        ) {
            let circles: Vec<(Vec2, f64)> =
                circles.into_iter().map(|(x, y, r)| (Vec2::new(x, y), r)).collect();
            // Many points on each circle give a hull slightly smaller than the exact one
            let n = 4096;
            let points: Vec<Vec2> = circles
                .iter()
                .flat_map(|&(c, r)| {
                    (0..n).map(move |k| c + Vec2::from_angle(2. * PI * k as f64 / n as f64) * r)
                })
                .collect();
            let sampled = polygon_area(&convex_hull(&points));
            let exact = circle_hull_area(&circles);
            prop_assert!(exact >= sampled - 1e-9, "{} {}", exact, sampled);
            prop_assert!(exact <= sampled * (1. + 1e-5) + 1e-9, "{} {}", exact, sampled);
            // With equal radii it matches the area of the grown hull of the centres
            let equal: Vec<Outline> = circles.iter().map(|&(c, _)| Outline::Circle(c, 0.7)).collect();
            let equal_circles: Vec<(Vec2, f64)> = circles.iter().map(|&(c, _)| (c, 0.7)).collect();
            prop_assert!((hull_area(&equal) - circle_hull_area(&equal_circles)).abs() < 1e-9);
        }
    }
}
