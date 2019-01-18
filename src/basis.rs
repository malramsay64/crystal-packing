//
// basis.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
extern crate rand;

use rand::Rng;

/// A Value which can be modified in many places
///
/// Rust's ownership rules preclude multiple mutable borrows from taking place. The SharedValue
/// struct is a wrapper around some of Rust's built in methods of handling this, providing a
/// nicer interface for that tasks I require.
///
pub struct SharedValue {
    value: std::rc::Rc<std::cell::RefCell<f64>>,
}

impl Clone for SharedValue {
    fn clone(&self) -> Self {
        return Self {
            value: std::rc::Rc::clone(&self.value),
        };
    }
}

impl SharedValue {
    pub fn get_value(&self) -> f64 {
        return *self.value.borrow();
    }
    pub fn set_value(&self, value: f64) {
        *self.value.borrow_mut() = value;
    }

    pub fn new(val: f64) -> Self {
        return Self {
            value: std::rc::Rc::new(std::cell::RefCell::new(val)),
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_shared_value() {
        let value = SharedValue::new(1.);
        assert_eq!(value.get_value(), 1.);
    }

    #[test]
    fn update_shared_value() {
        let value = SharedValue::new(1.);
        assert_eq!(value.get_value(), 1.);
        value.set_value(0.5);
        assert_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn cloned_shared_value() {
        let value1 = SharedValue::new(1.);
        let value2 = value1.clone();
        assert_eq!(value1.get_value(), 1.);
        assert_eq!(value2.get_value(), 1.);

        value2.set_value(0.5);
        assert_eq!(value1.get_value(), 0.5);
        assert_eq!(value2.get_value(), 0.5);
    }
}

pub trait Basis {
    fn set_value(&mut self, new_value: f64);
    fn get_value(&self) -> f64;
    fn reset_value(&self);
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64;
}

#[derive(Clone)]
pub struct StandardBasis {
    pub value: SharedValue,
    pub old: f64,
    pub min: f64,
    pub max: f64,
}

impl StandardBasis {
    pub fn new(value: f64, min: f64, max: f64) -> Self {
        Self {
            value: SharedValue::new(value),
            old: value,
            min,
            max,
        }
    }

    fn value_range(&self) -> f64 {
        return self.max - self.min;
    }
}

impl Basis for StandardBasis {
    fn get_value(&self) -> f64 {
        return self.value.get_value();
    }

    fn set_value(&mut self, new_value: f64) {
        self.old = self.get_value();
        self.value.set_value(match new_value {
            x if x < self.min => self.min,
            x if x > self.max => self.max,
            x => x,
        })
    }

    fn reset_value(&self) {
        self.value.set_value(self.old);
    }

    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        return self.get_value() + self.value_range() * rng.gen_range(-0.5, 0.5);
    }
}
