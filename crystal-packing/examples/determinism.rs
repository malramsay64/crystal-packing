//! Write a set of jobs and their outcomes as JSON, so other builds (e.g. WebAssembly) can check
//! they produce bit-identical results.

use crystal_packing::container::Container;
use crystal_packing::monte_carlo::MonteCarloConfig;
use crystal_packing::potential::LennardJones;
use crystal_packing::search::{run_job, Candidate, Job, ShapeSpec};

fn main() {
    let trimer = ShapeSpec::Trimer {
        radius: 0.637556,
        angle: 120.,
        distance: 1.,
    };
    let cases = [
        (ShapeSpec::Polygon { sides: 5 }, "p2", vec!['e'], None, None),
        (
            ShapeSpec::Polygon { sides: 3 },
            "p3m1",
            vec!['d'],
            None,
            None,
        ),
        (ShapeSpec::Disk, "p6mm", vec!['a'], None, None),
        (trimer.clone(), "p2gg", vec!['c'], None, None),
        (trimer, "p2", vec!['e'], Some(LennardJones::default()), None),
        (
            ShapeSpec::Polygon { sides: 5 },
            "finite",
            vec!['a'; 4],
            None,
            Some(Container::Circle),
        ),
        (
            ShapeSpec::Polygon { sides: 5 },
            "finite",
            vec!['a'; 3],
            None,
            Some(Container::Strip),
        ),
        (
            ShapeSpec::Disk,
            "finite",
            vec!['a'; 4],
            None,
            Some(Container::Cluster),
        ),
    ];
    let jobs: Vec<Job> = cases
        .into_iter()
        .enumerate()
        .map(|(i, (shape, group, sites, potential, container))| Job {
            shape,
            candidate: Candidate {
                group: group.into(),
                sites,
            },
            seed: 42 + i as u64,
            config: MonteCarloConfig {
                steps: 3000,
                ..Default::default()
            },
            start: None,
            potential,
            container,
        })
        .collect();
    let outcomes: Vec<_> = jobs.iter().map(run_job).collect();
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({ "jobs": jobs, "outcomes": outcomes }))
            .unwrap()
    );
}
