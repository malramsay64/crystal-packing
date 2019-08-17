//
// basis.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

use rand::Rng;

use crate::traits::Basis;

/// A pointer to a value so it can be modified in many places
///
/// Rust's ownership rules preclude multiple mutable borrows from taking place. The basis interface
/// is designed around the assumption that this is possible. This provides an abstraction around
/// handling a value that is accessible in multiple locations.
///
/// This abstracts away the implementation details, allowing for a range of different methods to be
/// tested and implemented. The current implementation, based on using unsafe pointers, has the
/// best performance by a significant factor.
///
/// Using unsafe, does have a number of disadvantages, most notably there are no lifetime checks
/// for the pointer.
///
/// ```
/// use packing::SharedValue;
/// let mut x = 1.;
/// let shared_x = SharedValue::new(&mut x);
///
/// assert!(shared_x.get_value() == 1.);
/// assert!(shared_x.get_value() == x);
///
/// // Updating the SharedValue will update the linked variable
/// shared_x.set_value(2.);
/// assert!(shared_x.get_value() == x);
/// assert!(shared_x.get_value() == 2.);
/// assert!(x == 2.);
///
/// // When updating x, the SharedValue will also update
/// x = 10.;
/// assert!(shared_x.get_value() == x);
/// ```
///
#[derive(Debug)]
pub struct SharedValue {
    value: *mut f64,
}

impl Clone for SharedValue {
    /// Allow a value to be updated in another location
    ///
    /// This allows read/write access to the value being pointed to in an additional location.
    ///
    /// ```
    /// use packing::SharedValue;
    /// let mut x = 1.;
    /// let shared_x = SharedValue::new(&mut x);
    /// let shared_y = shared_x.clone();
    ///
    /// assert!(shared_x.get_value() == x);
    /// assert!(shared_y.get_value() == x);
    ///
    /// // Updating the value in one place will update all values
    /// x = 2.;
    /// assert!(shared_x.get_value() == x);
    /// assert!(shared_y.get_value() == x);
    /// ```
    ///
    fn clone(&self) -> Self {
        Self { value: self.value }
    }
}

impl SharedValue {
    /// Create a SharedValue allowing modification of the given value
    ///
    /// # Arguments
    ///
    /// * `val` - A mutable reference to a f64 value which can be updated in multiple locations
    ///
    /// # Remarks
    ///
    /// This provides a highly performant access to modifying the value of a variable in multiple
    /// locations.
    ///
    pub fn new(val: &mut f64) -> Self {
        Self {
            value: val as *mut f64,
        }
    }

    /// Get the value of the variable being shared
    pub fn get_value(&self) -> f64 {
        unsafe { *self.value }
    }

    /// This updates the value which is being shared
    ///
    /// # Arguments
    ///
    /// * `value` - The new value to assign to the shared value
    ///
    /// # Remarks
    ///
    /// This breaks the single mutability rules of rust, and is consequently unsafe to use in
    /// threaded code.
    ///
    /// # Example
    ///
    /// ```
    /// use packing::SharedValue;
    /// let mut x = 1.;
    /// let shared_x = SharedValue::new(&mut x);
    ///
    /// // Update the value of x through shared_x
    /// shared_x.set_value(2.);
    ///
    /// // The values of both x and shared_x will be updated
    /// assert!(x == 2.);
    /// assert!(shared_x.get_value() == 2.);
    /// ```
    ///
    pub fn set_value(&self, value: f64) {
        unsafe {
            *self.value = value;
        }
    }
}

#[cfg(test)]
mod shared_value_tests {
    // The floating point values should be the same, so using standard equality is valid.
    #![allow(clippy::float_cmp)]

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn new() {
        let value = SharedValue::new(&mut 1.);
        assert_eq!(value.get_value(), 1.);
    }

    #[test]
    fn new_from_var() {
        let mut x = 1.;
        let value = SharedValue::new(&mut x);
        assert_eq!(value.get_value(), x);
    }

    #[test]
    fn set_value() {
        let value = SharedValue::new(&mut 1.);
        assert_abs_diff_eq!(value.get_value(), 1.);
        value.set_value(0.5);
        assert_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn set_value_var() {
        // Setup
        let mut x = 1.;
        let value = SharedValue::new(&mut x);
        assert_eq!(value.get_value(), 1.);
        assert_eq!(value.get_value(), x);

        // Set from shared value
        value.set_value(0.5);
        assert_eq!(value.get_value(), x);
        assert_eq!(value.get_value(), 0.5);
        assert_eq!(x, 0.5);

        // Set from variable
        x = 2.;
        assert_eq!(value.get_value(), x);
        assert_eq!(value.get_value(), 2.);
        assert_eq!(x, 2.);
    }

    #[test]
    fn clone() {
        let value1 = SharedValue::new(&mut 1.);
        let value2 = value1.clone();
        assert_eq!(value1.get_value(), 1.);
        assert_eq!(value2.get_value(), 1.);

        value2.set_value(0.5);
        assert_eq!(value1.get_value(), 0.5);
        assert_eq!(value2.get_value(), 0.5);
    }

    #[test]
    fn clone_var() {
        // Setup
        let mut x = 1.;
        let value1 = SharedValue::new(&mut x);
        let value2 = value1.clone();
        assert_eq!(value1.get_value(), x);
        assert_eq!(value2.get_value(), x);

        // Set from value1
        value1.set_value(0.);
        assert_eq!(x, 0.);
        assert_eq!(value1.get_value(), x);
        assert_eq!(value2.get_value(), x);

        // Set from value2
        value2.set_value(0.5);
        assert_eq!(x, 0.5);
        assert_eq!(value1.get_value(), x);
        assert_eq!(value2.get_value(), x);

        // Set from x
        x = 3.;
        assert_eq!(x, 3.);
        assert_eq!(value1.get_value(), x);
        assert_eq!(value2.get_value(), x);
    }
}

#[derive(Clone, Debug)]
pub struct StandardBasis {
    value: SharedValue,
    old: f64,
    min: f64,
    max: f64,
}

impl StandardBasis {
    pub fn new(value: SharedValue, min: f64, max: f64) -> Self {
        Self {
            old: value.get_value(),
            value,
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

    fn sample<R: Rng + ?Sized>(&self, rng: &mut R, step_size: f64) -> f64 {
        self.get_value() + step_size * self.value_range() * rng.gen_range(-0.5, 0.5)
    }

    fn set_sampled<R: Rng + ?Sized>(&mut self, rng: &mut R, step_size: f64) {
        self.set_value(self.sample(rng, step_size));
    }
}

#[cfg(test)]
mod standard_basis_tests {
    use approx::assert_abs_diff_eq;
    use rand::thread_rng;

    use super::*;

    #[test]
    fn get_value() {
        let value = SharedValue::new(&mut 1.);
        let mut basis = StandardBasis::new(value.clone(), 0., 1.);
        basis.set_value(0.5);
        assert_abs_diff_eq!(basis.get_value(), 0.5);
        assert_abs_diff_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn set_value_limits() {
        let value = SharedValue::new(&mut 1.);
        let mut basis = StandardBasis::new(value, 0., 1.);

        // Over maximum value
        basis.set_value(1.1);
        assert_abs_diff_eq!(basis.get_value(), 1.);

        // Less than minimum value
        basis.set_value(-0.1);
        assert_abs_diff_eq!(basis.get_value(), 0.);
    }

    #[test]
    fn reset_value() {
        let value = SharedValue::new(&mut 1.);
        let mut basis = StandardBasis::new(value, 0., 1.);
        basis.set_value(0.5);
        assert_abs_diff_eq!(basis.get_value(), 0.5);
        basis.reset_value();
        assert_abs_diff_eq!(basis.get_value(), 1.);
    }

    #[test]
    fn sample() {
        let value = SharedValue::new(&mut 1.);
        let basis = StandardBasis::new(value, 0., 1.);
        let mut rng = thread_rng();
        for _ in 0..100 {
            let val = basis.sample(&mut rng, 1.);
            // Range of values which should be present
            assert!(0.5 <= val && val <= 1.5);
        }
    }

}
