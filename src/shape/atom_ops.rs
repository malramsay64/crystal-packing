//
// atom_ops.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

#![allow(clippy::op_ref)]

use std::ops::Mul;

use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, DimName};

use crate::shape::Atom;
use crate::symmetry::Transform;

impl<'a, 'b, D: DimName> Mul<&'b Transform<D>> for &'a Atom<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Atom<D>;

    fn mul(self, rhs: &Transform<D>) -> Self::Output {
        Atom {
            position: rhs * self.position.clone(),
            radius: self.radius,
        }
    }
}

impl<'a, D: DimName> Mul<Transform<D>> for &'a Atom<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Atom<D>;

    fn mul(self, rhs: Transform<D>) -> Self::Output {
        self * &rhs
    }
}

impl<'a, D: DimName> Mul<&'a Transform<D>> for Atom<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Atom<D>;

    fn mul(self, rhs: &Transform<D>) -> Self::Output {
        &self * rhs
    }
}

impl<D: DimName> Mul<Transform<D>> for Atom<D>
where
    DefaultAllocator: Allocator<f64, D>,
    DefaultAllocator: Allocator<f64, D, D>,
{
    type Output = Atom<D>;

    fn mul(self, rhs: Transform<D>) -> Self::Output {
        &self * &rhs
    }
}
