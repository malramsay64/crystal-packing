use crystal_packing::engine::Packing;
use crystal_packing::hard_shape::{HardShape, Molecule, Polygon};
use crystal_packing::monte_carlo::{MonteCarloConfig, Optimiser};
use crystal_packing::symmetry::PlaneGroup;
use crystal_packing::Rng;

fn run<S: HardShape>(shape: S, group: &str, sites: &[char], steps: u64, reps: u64) {
    let g = PlaneGroup::from_name(group).unwrap();
    let t = std::time::Instant::now();
    let mut best: f64 = 0.;
    for seed in 0..reps {
        let mut rng = Rng::seed_from_u64(seed);
        let packing = Packing::new(shape.clone(), g, sites, &mut rng).unwrap();
        let mut opt = Optimiser::new(
            packing,
            MonteCarloConfig {
                steps,
                ..Default::default()
            },
            seed,
        );
        opt.run(steps);
        let s = opt.finish();
        best = best.max(s.packing_fraction);
    }
    println!(
        "{:10} {:6} {:?}: {:.5} in {:.2?}",
        shape.name(),
        group,
        sites,
        best,
        t.elapsed()
    );
}

fn main() {
    if std::env::var("PROFILE").is_ok() {
        profile();
        return;
    }
    let steps = std::env::args()
        .nth(1)
        .map(|s| s.parse().unwrap())
        .unwrap_or(20000);
    run(Molecule::disk(), "p1", &['a'], steps, 4);
    run(Polygon::regular(4).unwrap(), "p1", &['a'], steps, 4);
    run(Polygon::regular(3).unwrap(), "p2", &['e'], steps, 4);
    run(Polygon::regular(5).unwrap(), "p2", &['e'], steps, 4);
    run(Polygon::regular(7).unwrap(), "p2", &['e'], steps, 4);
    run(Polygon::regular(6).unwrap(), "p6", &['a'], steps, 4);
    run(Polygon::regular(5).unwrap(), "p2gg", &['c'], steps, 4);
}

#[allow(dead_code)]
pub fn profile() {
    let g = PlaneGroup::from_name("p6").unwrap();
    let mut rng = Rng::seed_from_u64(0);
    let packing = Packing::new(Polygon::regular(6).unwrap(), g, &['a'], &mut rng).unwrap();
    let mut opt = Optimiser::new(packing, MonteCarloConfig::default(), 0);
    let t = std::time::Instant::now();
    let p = opt.run(20000);
    println!("run {:?} {:?} steps {:?}", t.elapsed(), p, opt.step_sizes());
    let t = std::time::Instant::now();
    let s = opt.finish();
    println!("finish {:?} {}", t.elapsed(), s.packing_fraction);
}
