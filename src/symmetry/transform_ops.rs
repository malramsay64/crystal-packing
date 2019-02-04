//
// symmetry_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]
use std::ops::{Add, Mul};

use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, DimName, Point, VectorN};

use crate::symmetry::Transform;

// Transform + Vector
//
impl<D> Add<VectorN<f64, D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn add(self, rhs: VectorN<f64, D>) -> Transform<D> {
        &self + &rhs
    }
}

impl<'a, D> Add<&'a VectorN<f64, D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn add(self, rhs: &'a VectorN<f64, D>) -> Transform<D> {
        &self + rhs
    }
}

impl<'a, D> Add<VectorN<f64, D>> for &'a Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn add(self, rhs: VectorN<f64, D>) -> Transform<D> {
        self + &rhs
    }
}

impl<'a, 'b, D> Add<&'b VectorN<f64, D>> for &'a Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn add(&'a self, rhs: &'b VectorN<f64, D>) -> Transform<D> {
        Transform {
            rotation: self.rotation,
            translation: self.translation + rhs,
        }
    }
}

// Transform * Vector
//
impl<D> Mul<VectorN<f64, D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: VectorN<f64, D>) -> VectorN<f64, D> {
        &self * &rhs
    }
}

impl<'a, D> Mul<&'a VectorN<f64, D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: &'a VectorN<f64, D>) -> VectorN<f64, D> {
        &self * rhs
    }
}

impl<'a, D> Mul<VectorN<f64, D>> for &'a Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: VectorN<f64, D>) -> VectorN<f64, D> {
        self * &rhs
    }
}

impl<'a, 'b, D> Mul<&'b VectorN<f64, D>> for &'a Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(&'a self, rhs: &'b VectorN<f64, D>) -> VectorN<f64, D> {
        // Multiplying by a vector doesn't apply the translational component, only the rotation.
        self.rotation * rhs
    }
}

// Transform * Point
//
impl<D> Mul<Point<f64, D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: Point<f64, D>) -> Point<f64, D> {
        &self * &rhs
    }
}

impl<'a, D> Mul<&'a Point<f64, D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: &'a Point<f64, D>) -> Point<f64, D> {
        &self * rhs
    }
}

impl<'a, D> Mul<Point<f64, D>> for &'a Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: Point<f64, D>) -> Point<f64, D> {
        self * &rhs
    }
}

impl<'a, 'b, D> Mul<&'b Point<f64, D>> for &'a Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(&'a self, rhs: &'b Point<f64, D>) -> Point<f64, D> {
        self.rotation * rhs + self.translation
    }
}

// Transform * Transform
//
impl<D> Mul<Transform<D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: Transform<D>) -> Transform<D> {
        &self * &rhs
    }
}

impl<'a, D> Mul<&'a Transform<D>> for Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(self, rhs: &'a Transform<D>) -> Transform<D> {
        &self * rhs
    }
}

impl<'a, D> Mul<Transform<D>> for &'a Transform<D>
where
    DefaultAllocator: Allocator<f64, D>,
    D: DimName,
{
    fn mul(self, rhs: Transform<D>) -> Transform<D> {
        self * &rhs
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<'a, 'b, D> Mul<&'b Transform<D>> for &'a Transform<D>
where
    D: DimName,
    DefaultAllocator: Allocator<f64, D>,
{
    fn mul(&'a self, rhs: &'b Transform<D>) -> Transform<D> {
        let shift = self.rotate(&rhs.translation);

        Transform {
            translation: self.translation + shift,
            rotation: self.rotation * rhs.rotation,
        }
    }
}
