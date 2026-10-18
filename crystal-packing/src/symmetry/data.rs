//
// data.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Definitions of the plane groups from the International Tables for Crystallography, Volume A.
//!
//! Each group lists the coordinates of the general position (excluding centring translations)
//! and the Wyckoff positions as (letter, multiplicity, representative coordinates). The hexagonal
//! groups use the standard setting with a cell angle of 120°.

use super::LatticeFamily;
use LatticeFamily::*;

pub struct GroupDefinition {
    pub number: u8,
    pub name: &'static str,
    pub family: LatticeFamily,
    pub centred: bool,
    pub operations: &'static [&'static str],
    pub wyckoff: &'static [(char, usize, &'static str)],
}

/// Alternative short names for the groups, mapped to the names used here
pub const ALIASES: &[(&str, &str)] = &[
    ("pm", "p1m1"),
    ("pg", "p1g1"),
    ("cm", "c1m1"),
    ("pmm", "p2mm"),
    ("pmg", "p2mg"),
    ("pgg", "p2gg"),
    ("cmm", "c2mm"),
    ("p4m", "p4mm"),
    ("p4g", "p4gm"),
    ("p6m", "p6mm"),
];

const P6_OPERATIONS: [&str; 6] = ["x,y", "-y,x-y", "-x+y,-x", "-x,-y", "y,-x+y", "x-y,x"];

pub const GROUPS: &[GroupDefinition] = &[
    GroupDefinition {
        number: 1,
        name: "p1",
        family: Oblique,
        centred: false,
        operations: &["x,y"],
        wyckoff: &[('a', 1, "x,y")],
    },
    GroupDefinition {
        number: 2,
        name: "p2",
        family: Oblique,
        centred: false,
        operations: &["x,y", "-x,-y"],
        wyckoff: &[
            ('e', 2, "x,y"),
            ('d', 1, "1/2,1/2"),
            ('c', 1, "1/2,0"),
            ('b', 1, "0,1/2"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 3,
        name: "p1m1",
        family: Rectangular,
        centred: false,
        operations: &["x,y", "-x,y"],
        wyckoff: &[('c', 2, "x,y"), ('b', 1, "1/2,y"), ('a', 1, "0,y")],
    },
    GroupDefinition {
        number: 4,
        name: "p1g1",
        family: Rectangular,
        centred: false,
        operations: &["x,y", "-x,y+1/2"],
        wyckoff: &[('a', 2, "x,y")],
    },
    GroupDefinition {
        number: 5,
        name: "c1m1",
        family: Rectangular,
        centred: true,
        operations: &["x,y", "-x,y"],
        wyckoff: &[('b', 4, "x,y"), ('a', 2, "0,y")],
    },
    GroupDefinition {
        number: 6,
        name: "p2mm",
        family: Rectangular,
        centred: false,
        operations: &["x,y", "-x,-y", "-x,y", "x,-y"],
        wyckoff: &[
            ('i', 4, "x,y"),
            ('h', 2, "1/2,y"),
            ('g', 2, "0,y"),
            ('f', 2, "x,1/2"),
            ('e', 2, "x,0"),
            ('d', 1, "1/2,1/2"),
            ('c', 1, "1/2,0"),
            ('b', 1, "0,1/2"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 7,
        name: "p2mg",
        family: Rectangular,
        centred: false,
        operations: &["x,y", "-x,-y", "-x+1/2,y", "x+1/2,-y"],
        wyckoff: &[
            ('d', 4, "x,y"),
            ('c', 2, "1/4,y"),
            ('b', 2, "0,1/2"),
            ('a', 2, "0,0"),
        ],
    },
    GroupDefinition {
        number: 8,
        name: "p2gg",
        family: Rectangular,
        centred: false,
        operations: &["x,y", "-x,-y", "-x+1/2,y+1/2", "x+1/2,-y+1/2"],
        wyckoff: &[('c', 4, "x,y"), ('b', 2, "1/2,0"), ('a', 2, "0,0")],
    },
    GroupDefinition {
        number: 9,
        name: "c2mm",
        family: Rectangular,
        centred: true,
        operations: &["x,y", "-x,-y", "-x,y", "x,-y"],
        wyckoff: &[
            ('f', 8, "x,y"),
            ('e', 4, "0,y"),
            ('d', 4, "x,0"),
            ('c', 4, "1/4,1/4"),
            ('b', 2, "0,1/2"),
            ('a', 2, "0,0"),
        ],
    },
    GroupDefinition {
        number: 10,
        name: "p4",
        family: Square,
        centred: false,
        operations: &["x,y", "-x,-y", "-y,x", "y,-x"],
        wyckoff: &[
            ('d', 4, "x,y"),
            ('c', 2, "1/2,0"),
            ('b', 1, "1/2,1/2"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 11,
        name: "p4mm",
        family: Square,
        centred: false,
        operations: &[
            "x,y", "-x,-y", "-y,x", "y,-x", "-x,y", "x,-y", "y,x", "-y,-x",
        ],
        wyckoff: &[
            ('g', 8, "x,y"),
            ('f', 4, "x,x"),
            ('e', 4, "x,1/2"),
            ('d', 4, "x,0"),
            ('c', 2, "1/2,0"),
            ('b', 1, "1/2,1/2"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 12,
        name: "p4gm",
        family: Square,
        centred: false,
        operations: &[
            "x,y",
            "-x,-y",
            "-y,x",
            "y,-x",
            "-x+1/2,y+1/2",
            "x+1/2,-y+1/2",
            "y+1/2,x+1/2",
            "-y+1/2,-x+1/2",
        ],
        wyckoff: &[
            ('d', 8, "x,y"),
            ('c', 4, "x,x+1/2"),
            ('b', 2, "1/2,0"),
            ('a', 2, "0,0"),
        ],
    },
    GroupDefinition {
        number: 13,
        name: "p3",
        family: Hexagonal,
        centred: false,
        operations: &["x,y", "-y,x-y", "-x+y,-x"],
        wyckoff: &[
            ('d', 3, "x,y"),
            ('c', 1, "2/3,1/3"),
            ('b', 1, "1/3,2/3"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 14,
        name: "p3m1",
        family: Hexagonal,
        centred: false,
        operations: &["x,y", "-y,x-y", "-x+y,-x", "-y,-x", "-x+y,y", "x,x-y"],
        wyckoff: &[
            ('e', 6, "x,y"),
            ('d', 3, "x,-x"),
            ('c', 1, "2/3,1/3"),
            ('b', 1, "1/3,2/3"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 15,
        name: "p31m",
        family: Hexagonal,
        centred: false,
        operations: &["x,y", "-y,x-y", "-x+y,-x", "y,x", "x-y,-y", "-x,-x+y"],
        wyckoff: &[
            ('d', 6, "x,y"),
            ('c', 3, "x,0"),
            ('b', 2, "1/3,2/3"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 16,
        name: "p6",
        family: Hexagonal,
        centred: false,
        operations: &P6_OPERATIONS,
        wyckoff: &[
            ('d', 6, "x,y"),
            ('c', 3, "1/2,0"),
            ('b', 2, "1/3,2/3"),
            ('a', 1, "0,0"),
        ],
    },
    GroupDefinition {
        number: 17,
        name: "p6mm",
        family: Hexagonal,
        centred: false,
        operations: &[
            P6_OPERATIONS[0],
            P6_OPERATIONS[1],
            P6_OPERATIONS[2],
            P6_OPERATIONS[3],
            P6_OPERATIONS[4],
            P6_OPERATIONS[5],
            "-y,-x",
            "-x+y,y",
            "x,x-y",
            "y,x",
            "x-y,-y",
            "-x,-x+y",
        ],
        wyckoff: &[
            ('f', 12, "x,y"),
            ('e', 6, "x,-x"),
            ('d', 6, "x,0"),
            ('c', 3, "1/2,0"),
            ('b', 2, "1/3,2/3"),
            ('a', 1, "0,0"),
        ],
    },
];
