//
// basis.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
extern crate rand;

use approx;
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
        Self {
            value: std::rc::Rc::clone(&self.value),
        }
    }
}

impl SharedValue {
    pub fn get_value(&self) -> f64 {
        *self.value.borrow()
    }
    pub fn set_value(&self, value: f64) {
        *self.value.borrow_mut() = value;
    }

    pub fn new(val: f64) -> Self {
        Self {
            value: std::rc::Rc::new(std::cell::RefCell::new(val)),
        }
    }
}

#[cfg(test)]
mod shared_value_tests {
    use super::*;

    #[test]
    fn new() {
        let value = SharedValue::new(1.);
        assert_abs_diff_eq!(value.get_value(), 1.);
    }

    #[test]
    fn set_value() {
        let value = SharedValue::new(1.);
        assert_abs_diff_eq!(value.get_value(), 1.);
        value.set_value(0.5);
        assert_abs_diff_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn clone() {
        let value1 = SharedValue::new(1.);
        let value2 = value1.clone();
        assert_abs_diff_eq!(value1.get_value(), 1.);
        assert_abs_diff_eq!(value2.get_value(), 1.);

        value2.set_value(0.5);
        assert_abs_diff_eq!(value1.get_value(), 0.5);
        assert_abs_diff_eq!(value2.get_value(), 0.5);
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
    value: SharedValue,
    old: f64,
    min: f64,
    max: f64,
}

impl StandardBasis {
    pub fn new(value: &SharedValue, min: f64, max: f64) -> Self {
        Self {
            value: value.clone(),
            old: value.get_value(),
            min,
            max,
        }
    }

    fn value_range(&self) -> f64 {
        self.max - self.min
    }
}

impl Basis for StandardBasis {
    fn get_value(&self) -> f64 {
        self.value.get_value()
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
        self.get_value() + self.value_range() * rng.gen_range(-0.5, 0.5)
    }
}

#[cfg(test)]
mod standard_basis_tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn get_value() {
        let value = SharedValue::new(1.);
        let mut basis = StandardBasis::new(&value, 0., 1.);
        basis.set_value(0.5);
        assert_abs_diff_eq!(basis.get_value(), 0.5);
        assert_abs_diff_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn set_value_limits() {
        let value = SharedValue::new(1.);
        let mut basis = StandardBasis::new(&value, 0., 1.);

        // Over maximum value
        basis.set_value(1.1);
        assert_abs_diff_eq!(basis.get_value(), 1.);

        // Less than minimum value
        basis.set_value(-0.1);
        assert_abs_diff_eq!(basis.get_value(), 0.);
    }

    #[test]
    fn reset_value() {
        let value = SharedValue::new(1.);
        let mut basis = StandardBasis::new(&value, 0., 1.);
        basis.set_value(0.5);
        assert_abs_diff_eq!(basis.get_value(), 0.5);
        basis.reset_value();
        assert_abs_diff_eq!(basis.get_value(), 1.);
    }

    #[test]
    fn sample() {
        let value = SharedValue::new(1.);
        let basis = StandardBasis::new(&value, 0., 1.);
        let mut rng = thread_rng();
        for _ in 0..100 {
            let val = basis.sample(&mut rng);
            // Range of values which should be present
            assert!(0.5 <= val && val <= 1.5);
        }
    }

}
