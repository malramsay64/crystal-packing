//
// shape_instance.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use std::fmt;
use std::fmt::{Debug, Display};
use std::ops::Mul;

use log::trace;
use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, U2};

use crate::shape::{Intersect, Shape};
use crate::symmetry::Transform;

/// Puts an abstract shape object in a physical space
///
/// This acts as a cache for computed values.
#[derive(PartialEq)]
pub struct ShapeInstance<I>
where
    I: Intersect<U2> + Debug + Display + Mul<Transform<U2>, Output = I>,
    for<'a> I: Mul<&'a Transform<U2>, Output = I>,
    for<'a, 'b> &'a I: Mul<&'b Transform<U2>, Output = I>,
    for<'a> &'a I: Mul<Transform<U2>, Output = I>,
    DefaultAllocator: Allocator<f64, U2>,
    DefaultAllocator: Allocator<f64, U2, U2>,
{
    pub items: Vec<I>,
}

impl<I> ShapeInstance<I>
where
    I: Intersect<U2> + Debug + Display + Mul<Transform<U2>, Output = I>,
    for<'a> I: Mul<&'a Transform<U2>, Output = I>,
    for<'a, 'b> &'a I: Mul<&'b Transform<U2>, Output = I>,
    for<'a> &'a I: Mul<Transform<U2>, Output = I>,
    DefaultAllocator: Allocator<f64, U2>,
    DefaultAllocator: Allocator<f64, U2, U2>,
{
    /// Create a ShapeInstance from a Shape and a Symmetry operation
    ///
    /// This takes the general shape typically centred around the origin, and transforms it into a
    /// position in the cell. In the simplest case, this is purely a transformation to put the
    /// shape in the appropriate coordinates. In other cases it performs both translations and
    /// rotations of the shape to the appropriate positions.
    ///
    pub fn from<S>(shape: &S, iso: &Transform<U2>) -> ShapeInstance<I>
    where
        S: Shape<U2, Component = I>,
    {
        trace!("Shape: {:?}, iso: {:?}", shape, iso);
        ShapeInstance {
            items: shape.iter().map(|p| p * iso).collect(),
        }
    }
}

impl<I> ShapeInstance<I>
where
    I: Intersect<U2> + Debug + Display + Mul<Transform<U2>, Output = I>,
    for<'a> I: Mul<&'a Transform<U2>, Output = I>,
    for<'a, 'b> &'a I: Mul<&'b Transform<U2>, Output = I>,
    for<'a> &'a I: Mul<Transform<U2>, Output = I>,
    DefaultAllocator: Allocator<f64, U2>,
    DefaultAllocator: Allocator<f64, U2, U2>,
{
    /// Check whether this shape intersects with another shape
    ///
    /// A ShapeInstance is considered to intersect with another when one of it's components
    /// intersects with a component of the other shape. For a square, there is an intersection
    /// when a line from one square crosses the other. Each component item of `self` is
    /// checked against `other`.
    ///
    pub fn intersects(&self, other: &ShapeInstance<I>) -> bool {
        // We want to compare every item of the current shape with every item of the other shape.
        for (index_a, item_a) in self.items.iter().enumerate() {
            for item_b in other.items.iter().skip(index_a) {
                if item_a.intersects(&item_b) {
                    return true;
                }
            }
        }
        false
    }
}

impl<I> Debug for ShapeInstance<I>
where
    I: Intersect<U2> + Debug + Display + Mul<Transform<U2>, Output = I>,
    for<'a> I: Mul<&'a Transform<U2>, Output = I>,
    for<'a, 'b> &'a I: Mul<&'b Transform<U2>, Output = I>,
    for<'a> &'a I: Mul<Transform<U2>, Output = I>,
    DefaultAllocator: Allocator<f64, U2>,
    DefaultAllocator: Allocator<f64, U2, U2>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ShapeInstance{{ {:?} }}", self.items)
    }
}

#[cfg(test)]
mod test {
    use std::f64::consts::PI;

    use approx::assert_abs_diff_eq;
    use itertools::izip;
    use nalgebra::Vector2;

    use crate::shape::{Line, LineShape};
    use crate::symmetry::Transform2;

    use super::*;

    #[test]
    fn lines() {
        let shape = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();
        let shape_i = ShapeInstance::from(&shape, &Transform2::identity());
        let expected_vec = vec![
            Line::new((0., 1.), (1., 0.)),
            Line::new((1., 0.), (0., -1.)),
            Line::new((0., -1.), (-1., 0.)),
            Line::new((-1., 0.), (0., 1.)),
        ];
        for (index, result, expected) in izip!(0.., shape_i.items.iter(), expected_vec.iter()) {
            println!("{}", index);
            assert_abs_diff_eq!(expected, result, epsilon = 1e-8);
        }
    }

    #[test]
    fn lines_translate() {
        let shape = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();
        let shape_i = ShapeInstance::from(&shape, &Transform2::new(na::Vector2::new(-1., 0.), 0.));
        let expected_vec = vec![
            Line::new((-1., 1.), (0., 0.)),
            Line::new((0., 0.), (-1., -1.)),
            Line::new((-1., -1.), (-2., 0.)),
            Line::new((-2., 0.), (-1., 1.)),
        ];
        for (index, result, expected) in izip!(0.., shape_i.items.iter(), expected_vec.iter()) {
            println!("{}", index);
            assert_abs_diff_eq!(*expected, *result, epsilon = 1e-8);
        }
    }

    #[test]
    fn lines_rotate_point() {
        let shape = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();
        let shape_i = ShapeInstance::from(&shape, &Transform2::new(na::Vector2::new(-1., 0.), PI));
        let expected_vec = vec![
            Line::new((-1., -1.), (-2., 0.)),
            Line::new((-2., 0.), (-1., 1.)),
            Line::new((-1., 1.), (0., 0.)),
            Line::new((0., 0.), (-1., -1.)),
        ];
        for (index, result, expected) in izip!(0.., shape_i.items.iter(), expected_vec.iter()) {
            println!("{}", index);
            assert_abs_diff_eq!(*expected, *result, epsilon = 1e-8);
        }
    }

    #[test]
    fn lines_translate_rotate() {
        let shape = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();
        let shape_i = ShapeInstance::from(&shape, &Transform2::new(na::Vector2::new(-1., 0.), PI));
        let expected_vec = vec![
            Line::new((-1., -1.), (-2., 0.)),
            Line::new((-2., 0.), (-1., 1.)),
            Line::new((-1., 1.), (0., 0.)),
            Line::new((0., 0.), (-1., -1.)),
        ];
        for (index, result, expected) in izip!(0.., shape_i.items.iter(), expected_vec.iter()) {
            println!("{}", index);
            assert_abs_diff_eq!(*expected, *result, epsilon = 1e-8);
        }
    }

    //
    // +-------------------|-------------------+
    // |                   |                   |
    // |          +        + 1       +         |
    // |           \       |       /           |
    // |             \     |     /             |
    // |               \   |   /               |
    // |                 \ | /                 |
    // |---------+---------+---------+---------|
    // |        -1       / | \       1         |
    // |               /   |   \               |
    // |             /     |     \             |
    // |           /       |       \           |
    // |         +         + -1      +         |
    // |                   |                   |
    // |                   |                   |
    // +-------------------|-------------------+
    //
    #[test]
    fn intersects() {
        let shape = LineShape::from_radial("Square", vec![1., 1., 1., 1.]).unwrap();
        let shape_i1 = ShapeInstance::from(&shape, &Transform2::new(Vector2::new(1., 0.), 0.));
        assert!(shape_i1.intersects(&shape_i1));

        let shape_i2 = ShapeInstance::from(&shape, &Transform2::new(Vector2::new(-1.001, 0.), 0.));
        assert!(!shape_i1.intersects(&shape_i2));

        let shape_i3 = ShapeInstance::from(&shape, &Transform2::new(Vector2::new(0., 0.), PI / 4.));
        assert!(shape_i1.intersects(&shape_i3));
        assert!(shape_i2.intersects(&shape_i3));
    }
}
