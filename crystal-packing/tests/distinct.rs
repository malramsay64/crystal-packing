//
// distinct.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! The same crystal is found in many plane groups; check the search recognises them.
//!
//! ```text
//! cargo test --release -p crystal_packing --test distinct -- --ignored
//! ```

use crystal_packing::search::{Search, SearchConfig, ShapeSpec, Strategy};

fn distinct(shape: ShapeSpec) -> (usize, usize) {
    let config = SearchConfig {
        strategy: Strategy::Exhaustive,
        seed: 4,
        ..Default::default()
    };
    let mut search = Search::new(shape, config).unwrap();
    search.run();
    let near = search.near_best(0.02).len();
    let distinct = search.distinct_near_best(0.02, None).unwrap();
    for d in distinct.iter() {
        let equivalent: Vec<String> = d.equivalent.iter().map(|r| r.candidate.label()).collect();
        println!(
            "{} {:.6} = {}",
            d.result.candidate.label(),
            d.result.best.as_ref().unwrap().packing_fraction,
            equivalent.join(", ")
        );
    }
    (near, distinct.len())
}

#[test]
#[ignore]
fn pentagons() {
    // The double lattice packing is found in p2, p1g1, p2mg and p2gg
    let (near, distinct) = distinct(ShapeSpec::Polygon { sides: 5 });
    assert!(near >= 4, "{}", near);
    assert_eq!(distinct, 1);
}

#[test]
#[ignore]
fn disks() {
    // The hexagonal packing has the symmetry of p6mm, so is found in many groups
    let (near, distinct) = distinct(ShapeSpec::Disk);
    assert!(near >= 5, "{}", near);
    assert_eq!(distinct, 1);
}

#[test]
#[ignore]
fn heptagons() {
    // The double lattice in p2, p1g1, p2mg and p2gg, along with two different structures with
    // 3 fold symmetry within 2% of it
    let (near, distinct) = distinct(ShapeSpec::Polygon { sides: 7 });
    assert!(near >= 6, "{}", near);
    assert_eq!(distinct, 3);
}
