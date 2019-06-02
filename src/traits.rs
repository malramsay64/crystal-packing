//
// traits.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

pub trait Intersect {
    fn intersects(&self, other: &Self) -> bool;
}

pub trait Shape {
    type Component: Intersect + fmt::Debug + fmt::Display;

    fn area(&self) -> f64;
    fn enclosing_radius(&self) -> f64;
    fn get_items(&self) -> Vec<Self::Component>;
    fn rotational_symmetries(&self) -> u64 {
        1
    }
    fn iter(&self) -> slice::Iter<'_, Self::Component>;
}
