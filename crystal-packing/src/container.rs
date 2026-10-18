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

use std::f64::consts::PI;

use anyhow::{bail, Error};
use serde::{Deserialize, Serialize};

use crate::engine::{SiteRecord, Structure};
use crate::geometry::{convex_hull, perimeter, polygon_area, wrap_angle, Mat2, Vec2};
use crate::hard_shape::{HardShape, Outline, OVERLAP_TOLERANCE};
use crate::lattice::Lattice;
use crate::monte_carlo::{MonteCarloConfig, Progress, StepSize};
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

/// The number of points used for the circles of molecules with atoms of different sizes
const CIRCLE_POINTS: usize = 64;

/// The area of the convex hull of a set of outlines.
///
/// The hull of circles of equal radius r is the hull of their centres grown by r, with an area of
/// A + P r + π r², from the area A and perimeter P of the hull of the centres, which is exact. With
/// circles of different sizes, each is shrunk by the smallest radius and replaced by a polygon
/// around it, slightly overestimating the area.
pub fn hull_area(outlines: &[Outline]) -> f64 {
    let r0 = outlines
        .iter()
        .map(|o| match o {
            Outline::Circle(_, r) => *r,
            Outline::Polygon(_) => 0.,
        })
        .fold(f64::INFINITY, f64::min);
    if !r0.is_finite() {
        return 0.;
    }
    let mut points = vec![];
    for outline in outlines {
        match outline {
            Outline::Polygon(p) => points.extend_from_slice(p),
            Outline::Circle(c, r) => {
                let rho = r - r0;
                if rho <= 1e-12 * r {
                    points.push(*c);
                } else {
                    // The vertices of a polygon containing the circle
                    let radius = rho / libm::cos(PI / CIRCLE_POINTS as f64);
                    points.extend((0..CIRCLE_POINTS).map(|k| {
                        *c + Vec2::from_angle(2. * PI * k as f64 / CIRCLE_POINTS as f64) * radius
                    }));
                }
            }
        }
    }
    let hull = convex_hull(&points);
    polygon_area(&hull) + perimeter(&hull) * r0 + PI * r0 * r0
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
                S::support(placed, Vec2::new(1., 0.)) <= w
                    && S::support(placed, Vec2::new(-1., 0.)) <= w
                    && S::support(placed, Vec2::new(0., 1.)) <= h
                    && S::support(placed, Vec2::new(0., -1.)) <= h
            }
            Container::Circle => S::farthest(placed, Vec2::ZERO) <= 0.5 * size.width + tolerance,
            Container::Strip => {
                let h = 0.5 * size.height + tolerance;
                S::support(placed, Vec2::new(0., 1.)) <= h
                    && S::support(placed, Vec2::new(0., -1.)) <= h
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
        };
        for (i, site) in structure.sites.iter().enumerate() {
            confined.set_shape(i, Vec2::new(site.x, site.y), site.angle);
        }
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
    }

    /// Move a shape to a position (in Cartesian coordinates) and orientation
    pub fn set_shape(&mut self, index: usize, position: Vec2, angle: f64) {
        let mut x = 2. * position.x / self.size.width;
        if self.container.is_periodic() {
            x -= 2. * libm::floor(0.5 * (x + 1.));
        }
        self.positions[index] = Vec2::new(x, 2. * position.y / self.size.height);
        self.angles[index] = wrap_angle(angle, 2. * PI);
        self.place(index);
    }

    /// Change the size of the container, scaling the positions of the shapes with it
    pub fn set_size(&mut self, size: ContainerSize) {
        self.size = self.container.constrain(size);
        for index in 0..self.positions.len() {
            self.place(index);
        }
    }

    /// The area of the container, or of the convex hull of a cluster
    pub fn area(&self) -> f64 {
        if self.container.is_cluster() {
            let outlines: Vec<Outline> = (0..self.num_shapes())
                .flat_map(|i| {
                    self.shape
                        .outlines(&Mat2::rotation(self.angles[i]), self.position(i))
                })
                .collect();
            hull_area(&outlines)
        } else {
            self.container.area(&self.size)
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

    fn image_overlaps(&self, i: usize, j: usize, shift: Vec2) -> bool {
        let d = (self.position(j) + shift - self.position(i)).norm_sq();
        if d >= self.contact * self.contact {
            return false;
        }
        d < self.core * self.core || S::overlaps(&self.placed[i], &self.placed[j], shift)
    }

    fn pair_overlaps(&self, i: usize, j: usize) -> bool {
        let images = self.images();
        (-images..=images)
            .any(|k| self.image_overlaps(i, j, Vec2::new(k as f64 * self.size.width, 0.)))
    }

    /// Whether a shape overlaps its own periodic images
    fn self_overlaps(&self, i: usize) -> bool {
        (1..=self.images())
            .any(|k| self.image_overlaps(i, i, Vec2::new(k as f64 * self.size.width, 0.)))
    }

    fn outside(&self, index: usize) -> bool {
        !self
            .container
            .contains::<S>(&self.size, &self.placed[index], self.shape.circumradius())
    }

    /// Whether a shape overlaps another shape or the walls of the container
    pub fn shape_overlaps(&self, index: usize) -> bool {
        self.outside(index)
            || self.self_overlaps(index)
            || (0..self.num_shapes()).any(|j| j != index && self.pair_overlaps(index, j))
    }

    pub fn any_overlap(&self) -> bool {
        let n = self.num_shapes();
        (0..n).any(|i| {
            self.outside(i) || self.self_overlaps(i) || (i + 1..n).any(|j| self.pair_overlaps(i, j))
        })
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
    }

    fn update_best(&mut self) {
        let phi = self.confined.packing_fraction();
        if phi > self.best_packing_fraction {
            self.best_packing_fraction = phi;
            self.best = self.snapshot();
        }
    }

    fn pressure(&self) -> f64 {
        let fraction = (self.step as f64 / self.config.steps.max(1) as f64).min(1.);
        self.config.pressure_start
            * libm::pow(
                self.config.pressure_end / self.config.pressure_start,
                fraction,
            )
    }

    /// Whether to accept a change in area, as a fraction of the area of a shape
    fn accept_area(&mut self, change: f64, pressure: f64) -> bool {
        if change <= 0. {
            true
        } else if pressure.is_infinite() {
            false
        } else {
            self.rng.uniform() < libm::exp(-pressure * change)
        }
    }

    fn move_shape(&mut self, index: usize, rotate: bool, pressure: f64) -> bool {
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
            self.confined.set_shape(index, old_position, old_angle);
            return false;
        }
        if cluster {
            // Moving a shape changes the area of the hull
            let change = (self.confined.area() - old_area) / self.confined.shape.area();
            if !self.accept_area(change, pressure) {
                self.confined.set_shape(index, old_position, old_angle);
                return false;
            }
        }
        true
    }

    fn resize_container(&mut self, pressure: f64) -> bool {
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
                self.confined.set_size(old);
            }
            return false;
        }
        self.confined.set_size(new);
        if self.confined.any_overlap() {
            self.confined.set_size(old);
            return false;
        }
        true
    }

    fn step_once(&mut self, pressure: f64) {
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
            let pressure = self.pressure();
            self.step_once(pressure);
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
                self.step_once(f64::INFINITY);
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
    use crate::validate::validate_confined;
    use approx::assert_abs_diff_eq;

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
        // A small disk inside a large one doesn't change the hull, up to the polygon around it
        let mut outlines = disks(&[(0., 0.)], 2.);
        outlines.push(Outline::Circle(Vec2::new(0.5, 0.), 0.5));
        let area = hull_area(&outlines);
        assert!(
            (4. * PI - 1e-9..4. * PI * 1.002).contains(&area),
            "{}",
            area
        );
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
}
