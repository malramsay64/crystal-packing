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

use std::f64::consts::PI;

use anyhow::{bail, Error};
use serde::{Deserialize, Serialize};

use crate::engine::{SiteRecord, Structure};
use crate::geometry::{wrap_angle, Mat2, Vec2};
use crate::hard_shape::{HardShape, OVERLAP_TOLERANCE};
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
            Container::Square | Container::Rectangle => size.width * size.height,
            Container::Circle => PI * size.width * size.width / 4.,
        }
    }

    /// Keep the size consistent with the container's shape
    fn constrain(&self, size: ContainerSize) -> ContainerSize {
        match self {
            Container::Rectangle => size,
            Container::Square | Container::Circle => ContainerSize {
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
        }
    }
}

/// A finite number of shapes inside a container.
///
/// The positions of the shapes are stored relative to the size of the container, as fractions of
/// its half width and half height, so changing the size of the container scales the positions.
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
        self.positions[index] = Vec2::new(
            2. * position.x / self.size.width,
            2. * position.y / self.size.height,
        );
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

    pub fn packing_fraction(&self) -> f64 {
        self.num_shapes() as f64 * self.shape.area() / self.container.area(&self.size)
    }

    fn pair_overlaps(&self, i: usize, j: usize) -> bool {
        let d = (self.position(j) - self.position(i)).norm_sq();
        if d >= self.contact * self.contact {
            return false;
        }
        d < self.core * self.core || S::overlaps(&self.placed[i], &self.placed[j], Vec2::ZERO)
    }

    fn outside(&self, index: usize) -> bool {
        !self
            .container
            .contains::<S>(&self.size, &self.placed[index], self.shape.circumradius())
    }

    /// Whether a shape overlaps another shape or the walls of the container
    pub fn shape_overlaps(&self, index: usize) -> bool {
        self.outside(index)
            || (0..self.num_shapes()).any(|j| j != index && self.pair_overlaps(index, j))
    }

    pub fn any_overlap(&self) -> bool {
        let n = self.num_shapes();
        (0..n).any(|i| self.outside(i) || (i + 1..n).any(|j| self.pair_overlaps(i, j)))
    }

    pub fn structure(&self) -> Structure {
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
                &[self.size.width, self.size.height],
            ),
            num_shapes: self.num_shapes(),
            packing_fraction: self.packing_fraction(),
            energy: None,
            container: Some(ContainerState {
                kind: self.container,
                width: self.size.width,
                height: self.size.height,
            }),
        }
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

    fn move_shape(&mut self, index: usize, rotate: bool) -> bool {
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
            false
        } else {
            true
        }
    }

    fn resize_container(&mut self, pressure: f64) -> bool {
        let old = self.confined.size;
        let delta = 2. * self.resize.size * self.rng.centered();
        let new = match (self.confined.container, self.rng.index(3)) {
            (Container::Rectangle, 0) => ContainerSize {
                width: old.width * libm::exp(delta),
                height: old.height,
            },
            (Container::Rectangle, 1) => ContainerSize {
                width: old.width,
                height: old.height * libm::exp(delta),
            },
            _ => ContainerSize {
                width: old.width * libm::exp(delta),
                height: old.height * libm::exp(delta),
            },
        };
        let container = self.confined.container;
        let (area_old, area_new) = (container.area(&old), container.area(&new));
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
        } else if self.rng.index(2) == 0 {
            let accepted = self.move_shape(choice, false);
            self.translate.record(accepted, interval, target);
        } else {
            let accepted = self.move_shape(choice, true);
            self.rotate.record(accepted, interval, target);
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
}
