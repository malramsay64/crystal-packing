use std::cell::UnsafeCell;
use std::fmt;

use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

struct F64Visitor;

impl<'de> Visitor<'de> for F64Visitor {
    type Value = f64;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("a float")
    }

    fn visit_f32<E>(self, value: f32) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(f64::from(value))
    }

    fn visit_f64<E>(self, value: f64) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(value)
    }
}

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
///
/// ```
/// use crystal_packing::SharedValue;
/// let x = SharedValue::new(1.);
/// let shared_x = &x;
///
/// assert!(shared_x.get_value() == 1.);
/// assert!(shared_x.get_value() == x.get_value());
///
/// // Updating the SharedValue will update the linked variable
/// shared_x.set_value(2.);
/// assert!(shared_x.get_value() == x.get_value());
/// assert!(shared_x.get_value() == 2.);
/// assert!(x.get_value() == 2.);
///
/// ```
///
#[derive(Debug)]
pub struct SharedValue {
    value: UnsafeCell<f64>,
}

impl Serialize for SharedValue {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_f64(self.get_value())
    }
}

impl<'de> Deserialize<'de> for SharedValue {
    fn deserialize<D>(deserializer: D) -> Result<SharedValue, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer
            .deserialize_f64(F64Visitor)
            .map(SharedValue::new)
    }
}

unsafe impl Send for SharedValue {}
unsafe impl Sync for SharedValue {}

impl SharedValue {
    /// Create a SharedValue allowing modification of the given value
    ///
    /// # Arguments
    ///
    /// * `val` - A reference to a f64 value which can be updated in multiple locations
    ///
    /// # Remarks
    ///
    /// This provides a highly performant access to modifying the value of a variable in multiple
    /// locations.
    ///
    /// modifying the value will result in a runtime memory fault. An alternative implementation
    /// which takes `&mut f64` would not suffer from the same issues, however this then has issues
    /// with mutability of lifetimes.
    ///
    ///
    #[allow(clippy::trivially_copy_pass_by_ref)]
    pub fn new(val: f64) -> SharedValue {
        SharedValue {
            value: UnsafeCell::new(val),
        }
    }

    /// Get the value of the variable being shared
    pub fn get_value(&self) -> f64 {
        unsafe { *self.value.get() }
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
    /// use crystal_packing::SharedValue;
    /// let x = SharedValue::new(1.);
    /// let shared_x = &x;
    ///
    /// // Update the value of x through shared_x
    /// shared_x.set_value(2.);
    ///
    /// // The values of both x and shared_x will be updated
    /// assert_eq!(x.get_value(), 2.);
    /// assert_eq!(shared_x.get_value(), 2.);
    /// ```
    ///
    pub fn set_value(&self, value: f64) {
        unsafe {
            self.value.get().write(value);
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
        let x = 1.;
        let value = SharedValue::new(x);
        assert_eq!(value.get_value(), 1.);
    }

    #[test]
    fn set_value() {
        let value = SharedValue::new(1.);
        assert_abs_diff_eq!(value.get_value(), 1.);
        value.set_value(0.5);
        assert_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn set_value_var() {
        // Setup
        let value = SharedValue::new(1.);
        assert_eq!(value.get_value(), 1.);
        assert_eq!(value.get_value(), 1.);

        // Set from shared value
        value.set_value(0.5);
        assert_eq!(value.get_value(), 0.5);
    }

    #[test]
    fn pointers() {
        let value1 = SharedValue::new(1.);
        let value2 = &value1;
        assert_eq!(value1.get_value(), 1.);
        assert_eq!(value2.get_value(), 1.);

        value2.set_value(0.5);
        assert_eq!(value1.get_value(), 0.5);
        assert_eq!(value2.get_value(), 0.5);
    }

    #[test]
    fn clone_var() {
        // Setup
        let value1 = SharedValue::new(1.0);
        let value2 = &value1;
        assert_eq!(value1.get_value(), 1.0);
        assert_eq!(value2.get_value(), 1.0);

        // Set from value1
        value1.set_value(0.);
        assert_eq!(value1.get_value(), 0.);
        assert_eq!(value2.get_value(), 0.);

        // Set from value2
        value2.set_value(0.5);
        assert_eq!(value1.get_value(), 0.5);
        assert_eq!(value2.get_value(), 0.5);
    }
}
