//! Write a set of jobs and their outcomes as JSON, so other builds (e.g. WebAssembly) can check
//! they produce bit-identical results.

use crystal_packing::monte_carlo::MonteCarloConfig;
use crystal_packing::search::{run_job, Candidate, Job, ShapeSpec};

fn main() {
    let cases = [
        (ShapeSpec::Polygon { sides: 5 }, "p2", vec!['e']),
        (ShapeSpec::Polygon { sides: 3 }, "p3m1", vec!['d']),
        (ShapeSpec::Disk, "p6mm", vec!['a']),
        (
            ShapeSpec::Trimer {
                radius: 0.637556,
                angle: 120.,
                distance: 1.,
            },
            "p2gg",
            vec!['c'],
        ),
    ];
    let jobs: Vec<Job> = cases
        .into_iter()
        .enumerate()
        .map(|(i, (shape, group, sites))| Job {
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
        })
        .collect();
    let outcomes: Vec<_> = jobs.iter().map(run_job).collect();
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({ "jobs": jobs, "outcomes": outcomes }))
            .unwrap()
    );
}
