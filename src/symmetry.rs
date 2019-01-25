//
// symmetry.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//
//

use nalgebra::{IsometryMatrix2, Matrix2, Point2, Vector2};

/// Define the transformations of particle positions
///
/// These
#[derive(Debug, Clone)]
pub struct SymmetryTransform {
    pub isometry: IsometryMatrix2<f64>,
}

impl SymmetryTransform {
    fn parse_ops(ops: &str) -> (Vector2<f64>, f64) {
        let mut vec = Vector2::zeros();
        let mut sign = 1.;
        let mut constant = 0.;
        let mut operator: Option<char> = None;
        for c in ops.chars() {
            match c {
                'x' => {
                    vec[0] = sign;
                    sign = 1.;
                }
                'y' => {
                    vec[1] = sign;
                    sign = 1.;
                }
                '*' | '/' => {
                    operator = Some(c);
                }
                '-' => {
                    sign = -1.;
                }
                // This matches all digits from 0 to 9
                '0'...'9' => {
                    let val = c.to_string().parse::<u64>().unwrap() as f64;
                    // Is there an operator defined, i.e. is this the first digit
                    if let Some(op) = operator {
                        constant = match op {
                            '/' => sign * constant / val,
                            '*' => sign * constant * val,
                            _ => 0.,
                        };
                        operator = None;
                    } else {
                        constant = sign * val;
                    }
                    sign = 1.
                }
                // Default is do nothing (shouldn't encounter this at all)
                _ => {}
            };
        }

        (vec, constant)
    }

    pub fn new(sym_ops: &str) -> SymmetryTransform {
        let braces: &[_] = &['(', ')'];
        let operations: Vec<&str> = sym_ops
            // Remove braces from front and back
            .trim_matches(braces)
            // Split at the comma
            .split_terminator(',')
            .collect();
        let mut trans = Vector2::new(0., 0.);
        let mut rot: Matrix2<f64> = Matrix2::new(1., 0., 0., 1.);

        for (index, op) in operations.iter().enumerate() {
            let (r, t) = SymmetryTransform::parse_ops(op);
            rot.set_row(index, &r.transpose());
            trans[index] = t;
        }
        SymmetryTransform {
            isometry: IsometryMatrix2::from_parts(
                na::Translation2::from(trans),
                na::Rotation2::from_matrix_unchecked(rot),
            ),
        }
    }

    pub fn transform(&self, position: &Point2<f64>) -> Point2<f64> {
        self.isometry * position
    }

    pub fn rotate(&self, vect: &Vector2<f64>) -> Vector2<f64> {
        self.isometry * vect
    }
}

impl Default for SymmetryTransform {
    fn default() -> Self {
        Self {
            isometry: IsometryMatrix2::identity(),
        }
    }
}

#[cfg(test)]
mod symmetry_transform_tests {
    use super::*;

    fn create_identity() -> SymmetryTransform {
        SymmetryTransform {
            isometry: IsometryMatrix2::identity(),
        }
    }

    #[test]
    fn default() {
        let point = Point2::new(0.2, 0.2);
        let transform = SymmetryTransform::default();
        assert_eq!(transform.transform(&point), point);
    }

    #[test]
    fn identity_transform() {
        let identity = create_identity();
        let point = Point2::new(0.2, 0.2);
        assert_eq!(identity.transform(&point), point);

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(identity.rotate(&vec), vec);
    }

    #[test]
    fn transform() {
        let isometry = SymmetryTransform {
            isometry: IsometryMatrix2::new(Vector2::new(1., 1.), PI / 2.),
        };

        let point = Point2::new(0.2, 0.2);
        assert_eq!(isometry.transform(&point), Point2::new(0.8, 1.2));

        let vec = Vector2::new(0.2, 0.2);
        assert_eq!(isometry.rotate(&vec), Vector2::new(-0.2, 0.2));
    }

    #[test]
    fn parse_operation_default() {
        let input = String::from("(x, y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.1, 0.2));
    }

    #[test]
    fn parse_operation_xy() {
        let input = String::from("(-x, x+y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.1, 0.3));
    }

    #[test]
    fn parse_operation_consts() {
        let input = String::from("(x+1/2, -y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(0.6, -0.2));
    }

    #[test]
    fn parse_operation_neg_consts() {
        let input = String::from("(x-1/2, -y)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.4, -0.2));
    }

    #[test]
    fn parse_operation_zero_const() {
        let input = String::from("(-y, 0)");
        let st = SymmetryTransform::new(&input);
        let point = Point2::new(0.1, 0.2);
        assert_relative_eq!(st.transform(&point), Point2::new(-0.2, 0.));
    }
}
