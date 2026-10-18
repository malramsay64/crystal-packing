# Improvement plan: crystal-packing optimiser (CLI + WASM, CPU only)

This plan builds on the existing Rust workspace
(`crystal-packing`, `crystal-packing-cli`, `crystal-packing-wasm`).
It keeps what works and replaces the parts that limit result quality and
speed. The 2D C reference (`malramsay64/2d-packing/src/reference`) is used
only as a data source and a comparison baseline.

## 1. Current state

### 1.1 Baseline measurements

These runs used the release build on 4 cores, with `--replications 8`.
φ_max is the known optimum.

| Shape / group | steps | φ found | wall time | φ_max |
|---|---|---|---|---|
| Pentagon / p2 | 2·10⁵ | 0.630 | 0.8 s | 0.9213 |
| Pentagon / p2 | 2·10⁶ | 0.630 | 1.6 s | 0.9213 |
| Square / p1 | 2·10⁶ | 0.767 | 0.7 s | 1.0 |
| Disk / p1 | 2·10⁶ | 0.788 | 0.4 s | 0.9069 |
| Triangle / p2 | 2·10⁶ | 0.449 | 1.4 s | 1.0 |

The run time grows far less than 10× when the step count grows 10×, and φ
does not change. The optimiser stalls well short of the optimum, so the
first problem is algorithmic quality, then speed.

### 1.2 Defects and limits found in the code

| Area | Issue |
|---|---|
| `optimisation.rs` | Returns the **final** state, not the best one seen. With kT > 0 the run can end worse than its peak. |
| `optimisation.rs` | Step adaptation `step_ratio *= inner / (rejections + 1)` has its fixed point at about **100 % rejection**. Steps grow until nearly every move fails, so the chain freezes. The comment targets 75 % rejection. |
| `optimisation.rs` | The temperature acts on raw φ differences. No pressure/area term is used, and the CLI chains three runs at kT = 0, kT, kT = 0. |
| `cell.rs` | The length can only shrink (`max = initial`), the ratio stays in [0.1, 1] and the angle in [π/4, π/2]. There are no area-preserving shear or aspect moves, and the cell is never reduced, so the cell distorts and jams. |
| `state/packed.rs` | The periodic shell count (1, 2 or 3) comes from a ratio/angle heuristic and is not guaranteed to cover all neighbours. Every step runs a full O(M²·images) check. Each image calls `shape.transform`, which allocates a new `Vec` and clones the name `String`. |
| `line_shape.rs` / `line2.rs` | Overlap is detected only by edge crossings. Parallel or collinear overlapping edges return false (`u_b == 0`), containment is never detected, and the cost is O(n·m) per pair. |
| `basis/shared_value.rs` | `UnsafeCell` with `unsafe impl Sync`, and `Basis` holds references into the state. This is unsound across threads and makes undo and the best-state snapshot awkward. |
| `wallpaper.rs` | Only 7 groups, hand-entered. Only the general position is supported: `degrees_of_freedom()` always returns all true, there are no special Wyckoff sites and no isopointal enumeration. Hexagonal and tetragonal families are never exercised, and the user picks the group by hand. |
| `site.rs` / CLI | Every replica starts from the identical initial state, with x = y = −0.5 + 0.5/mult and θ = 0. Only the seed differs. |
| `Ord for PackedState` | `partial_cmp().unwrap()` panics on any invalid state, and this ordering drives `.max()` in the CLI. |
| `crystal-packing-wasm` | A `seed` framework counter template with no packing functionality. |
| Hygiene | Edition 2018. Dated dependencies: structopt, clap 2, nalgebra 0.27, criterion 0.3. The workspace release profile is `opt-level = 's'`, which slows the native CLI. A doctest in `ops_macros.rs` fails. |

### 1.3 Keep

- The workspace layout.
- The `Transform2` type and its `from_operations` symmetry parser.
- The general shape of the `Shape` / `State` traits.
- `MolecularShape2` (disks and trimers) and the LJ `PotentialState`.
- SVG output, the criterion benchmarks and the proptest setup.

## 2. Target architecture

```
crystal-packing/        core: no IO, no global state, no unsafe, builds for wasm32
  data/wallpaper.json   17 groups + Wyckoff sites (generated, embedded via include_str!)
  src/geometry/         convex polygons, SAT, disks, convex decomposition
  src/symmetry/         wallpaper groups, Wyckoff sites, isopointal enumeration
  src/cell.rs           h-matrix cell, Gauss reduction, image enumeration
  src/state/            parameter vector, image cache, overlap engine, undo log
  src/optimise/         floppy-box MC, step control, squeeze, search scheduler
crystal-packing-cli/    clap 4, rayon, JSON/CSV/SVG output
crystal-packing-wasm/   wasm-bindgen API + Web Worker pool (replaces the seed template)
```

- **State without shared mutation.** Replace `SharedValue`/`Basis` with a
  plain `Vec<f64>` parameter vector plus a static list of
  degree-of-freedom descriptors (kind, bounds, owning site). Moves write to
  the vector and push an undo record. Snapshots of the best state are
  `Vec<f64>` copies. `State: Send + Sync` then follows automatically, with
  no `unsafe`.
- **Resumable optimiser API.** `Optimiser::new(spec, seed)`,
  `run(n_steps) -> Progress`, `best() -> Structure`. The CLI loops over it.
  The browser runs it in time slices inside a worker so it can report
  progress and cancel.
- **Total ordering.** Results are ranked by `f64::total_cmp` on a
  validated score. No `unwrap` on `partial_cmp`.

### Libraries

| Purpose | Crate |
|---|---|
| Linear algebra | `nalgebra` (upgrade to current; it is already used throughout, so switching is not worth the churn) |
| Deterministic trig | `libm` in the hot path (bit-identical native vs wasm) |
| RNG | `rand` 0.9 + `rand_pcg` (already in use); `getrandom` (`wasm_js`) only for unseeded runs |
| Data and serialisation | `serde`, `serde_json`; `thiserror` in the library, `anyhow` in the binaries |
| Combinatorics | `itertools` (already in use) |
| CLI | `clap` 4 (derive; replaces structopt/paw), `rayon`, `indicatif`, `log` + `env_logger` |
| Output | `svg` (already in use) |
| WASM | `wasm-bindgen`, `js-sys`, `web-sys`, `serde-wasm-bindgen`, `console_error_panic_hook`; built with `wasm-pack` and `wasm-opt` |
| Testing | `proptest`, `approx`, `insta`, `wasm-bindgen-test`; `parry2d-f64` as a dev-only overlap oracle |
| Benchmarks | `criterion` 0.5 |

## 3. Algorithmic changes

### 3.1 Symmetry data and candidate enumeration

- Generate `data/wallpaper.json` for all 17 groups from the Bilbao data
  (the 2d-packing scraper or `reference/groups.c`). Each Wyckoff site
  records:
  - its letter and multiplicity;
  - its symmetry operations, as `from_operations` strings;
  - its free coordinates and site symmetry (rotation order, mirror
    directions).
- Correct degrees of freedom per site: fixed, line, or general for (x, y),
  and θ free, discrete (mirror-constrained) or fixed by the site rotation
  order.
- **Isopointal enumeration.** Take multisets of Wyckoff sites compatible
  with the shape's symmetry. The shape must contain the site symmetry, with
  rotation orders dividing and mirrors aligned. Cap the total shapes per
  cell, deduplicate, and enumerate all 17 groups automatically.
  `--group` becomes an optional filter.

### 3.2 Exact, fast overlap

- `ConvexPolygon` with cached edge normals uses the **Separating Axis
  Theorem**: O(n+m) per pair, and it correctly handles containment and
  collinear edges. `LineShape::from_radial` stays as a constructor that
  produces a polygon.
- Non-convex shapes are split into convex pieces once at construction.
  Disks and trimers keep the exact circle path in `MolecularShape2`.
- Each pair test runs from cheap to exact: circumradius distance →
  inradius overlap (an immediate overlap) → per-piece bounding circle → SAT.

### 3.3 Overlap engine

- An image cache holds transformed vertices for every symmetry image of
  every site, in flat arrays, rebuilt in place with no per-step allocation.
- **Asymmetric-unit checking.** One representative per occupied site is
  checked against **all** images within the cutoff, excluding only itself
  untranslated. This is complete, and costs k·M tests instead of M².
- Lattice translations are enumerated exactly from `|n·a + m·b| ≤ 2·R_max`
  on the **reduced** cell. This replaces the shell heuristic.
- **Incremental checks.** A site move updates only that orbit's images and
  checks only pairs involving it. A cell move rechecks everything, with
  early exit.

### 3.4 Cell

- Store the cell as an upper-triangular `h` matrix constrained by the
  family:
  - oblique: a, b, γ free;
  - rectangular: a, b free, γ = 90°;
  - square: a = b, γ = 90°;
  - hexagonal: a = b, γ = 60°.
- **Gauss lattice reduction** runs when the shear passes a threshold. Site
  fractional coordinates are remapped by the same unimodular matrix. This
  removes the hard angle/ratio bounds and the shrink-only length.

### 3.5 Optimiser: floppy-box isobaric Monte Carlo

- Move types:
  - site translation;
  - site rotation;
  - discrete mirror-orientation flips;
  - cell log-area scaling, shear and aspect moves.
- Cell-move acceptance: `min(1, exp(−βP·ΔA + N·ln(A'/A)))`. Any overlap
  rejects the move.
- A pressure ramp replaces the kT ramp: low βP to explore, rising to about
  10⁶ to compress.
- **Per-move-type step control** targets 30–40 % acceptance during
  equilibration, with clamped step bounds. Steps are then frozen for the
  production phase. This fixes the 100 %-rejection fixed point.
- The best-seen snapshot is tracked and returned.
- A **final squeeze** bisects the uniform cell scale to first contact, then
  runs coordinate line searches on each free parameter. This takes φ from
  MC noise to contact.
- The LJ `PotentialState` reuses the same move set and step control at
  constant temperature or pressure.

### 3.6 Search across candidates

- Diverse starts: random fractional positions and orientations, with a
  dilute initial cell.
- An **upper bound** prunes candidates: `φ ≤ N·A_shape / A_min`, where
  A_min comes from the inradius exclusion.
- **Successive halving:** every candidate gets a short run, and the top 1/η
  get η× more budget, repeated until one round remains.
- Optional **parallel tempering in pressure** across the replicas of the
  top candidates.

## 4. Computational changes

- No allocation, `String` clone or dynamic dispatch in the hot loop.
  Shapes are stored once, and images are vertex buffers.
- The release profile uses `opt-level = 3`, `lto = "fat"` and
  `codegen-units = 1` for the CLI. A separate `wasm-release` profile uses
  `opt-level = 3`, `panic = "abort"` and `+simd128`, followed by `wasm-opt -O3`.
- Parallelism happens only at the candidate and replica level:
  - CLI: `rayon` over the work queue;
  - browser: `navigator.hardwareConcurrency` Web Workers pulling work
    items, which avoids the COOP/COEP requirement of `wasm-bindgen-rayon`.
- Each work item gets its own seed, `hash(seed, group, sites, replica)`, so
  results do not depend on thread scheduling.
- CLI outputs:
  - ranked JSON of all candidates;
  - SVGs of the top N;
  - a CSV summary;
  - resuming from a saved structure (the existing `--start-config`, now
    implemented).

## 5. Tests

### 5.1 Unit and property tests (`cargo test`)

- **Geometry:**
  - SAT agrees with the `parry2d` oracle on random convex polygons;
  - the overlap test is symmetric and invariant under translation and
    rotation;
  - containment and collinear edges are covered, including regressions for
    the current `Line2` misses;
  - convex decomposition preserves area.
- **Symmetry data:**
  - there are exactly 17 groups;
  - each group's operations are closed modulo lattice translations;
  - each Wyckoff multiplicity equals its orbit size;
  - each special position is invariant under its site symmetry;
  - the isopointal count for small cases matches the C reference.
- **Cell:**
  - reduction preserves |det h| and the lattice;
  - remapped sites have identical Cartesian positions;
  - the image cutoff matches a brute-force 7×7 enumeration.
- **Overlap engine:**
  - the asymmetric-unit check equals a naive all-pairs check over a 5×5
    supercell;
  - the incremental check equals a full recheck after random move
    sequences;
  - move followed by undo restores the parameter vector exactly.
- **Optimiser:**
  - the step controller converges to its target acceptance;
  - the best snapshot is never worse than the final state;
  - the same seed gives the same result;
  - a hard-disk equation-of-state check at fixed low βP holds.

### 5.2 Known-answer regression suite (release mode, CI)

Every result is re-validated by an independent naive supercell overlap
check that shares no code path with the optimiser.

| Shape | Target φ (tolerance 1e-3) | Current |
|---|---|---|
| Disk | 0.9069 | 0.788 |
| Square | 1.0 | 0.767 |
| Triangle | 1.0 | 0.449 |
| Hexagon | 1.0 | — |
| Pentagon | 0.9213 (p2) | 0.630 |
| Heptagon | 0.8926 (p2) | — |
| Trimer (default parameters) | `insta` snapshot | — |

### 5.3 WASM

- `wasm-bindgen-test` runs in headless Chromium. A seeded run must give an
  identical trajectory hash to the native build.
- A smoke test: a worker-pool run completes and returns a valid structure.

### 5.4 Benchmarks (criterion; the existing `benches/intersection.rs` is extended)

- Pair overlap time in ns for polygons with 3–64 sides and for disks,
  comparing the old `LineShape` with SAT.
- State check time, full and incremental.
- MC steps per second, and time to reach within 1e-3 of the regression
  targets.
- Targets:
  - ≥ 10× steps/s on polygons;
  - the regression suite reached in under 10 s per shape on 4 cores;
  - wasm within 2× of native single-thread.

### 5.5 CI

Replace the deprecated `actions-rs` steps with `dtolnay/rust-toolchain`.
CI runs:

- `fmt --check`, `clippy -D warnings`;
- `cargo test --workspace`, including doctests;
- the regression suite;
- `wasm-pack test --headless --chrome`;
- a benchmark smoke run.

## 6. Phases

Each phase lands as its own PR, with that phase's tests passing.

0. **Hygiene.** Move to edition 2021, upgrade the dependencies (clap 4,
   nalgebra, rand 0.9, criterion 0.5), fix the failing doctest, set up the
   profiles, update CI and add the baseline benchmark numbers.
1. **State refactor.** Replace `SharedValue`/`Basis` with a parameter
   vector and undo log. Fix best-state tracking and `Ord`. Behaviour is
   otherwise unchanged and the benchmarks are recorded.
2. **Geometry.** SAT `ConvexPolygon`, convex decomposition and the oracle
   tests.
3. **Symmetry data.** All 17 groups, Wyckoff sites and isopointal
   enumeration.
4. **Cell and overlap engine.** The h-matrix cell, reduction, exact images
   and incremental checks.
5. **Optimiser.** Floppy-box MC, step control and the squeeze. The
   regression suite passes.
6. **Search and CLI.** Pruning, successive halving and the new outputs.
7. **WASM.** Real bindings, the worker pool, a minimal viewer, and the
   determinism test.

## 7. Status

All phases are implemented on the `optimiser-rewrite` branch. The original implementation has been
removed, with the Lennard-Jones potential ported to the new engine (`potential.rs`).

### Results

| Measure | Before | After |
|---|---|---|
| Pentagon φ (full 17 group search) | 0.630 (p2 only) | 0.921311 in 2.0 s |
| Square / disk / triangle φ | 0.767 / 0.788 / 0.449 | 1.0 / 0.906900 / 1.0 |
| Pentagon state overlap check | 1.8 µs | 0.19 µs (full), 0.19 µs (one site) |
| 10,000 MC steps, pentagons in p2 | 55 ms | 6.2 ms |
| Native vs wasm results | — | bit-identical (checked in CI) |

### Differences from the plan

- `nalgebra` was replaced by small `Vec2`/`Mat2` types, and `rand`/`rand_pcg` by an in-crate
  xoshiro256++ generator. This keeps the dependencies small and makes the random number stream
  part of the code, which the cross-platform determinism relies on.
- Errors use `anyhow` throughout rather than `thiserror`.
- Values cross the wasm boundary as JSON strings rather than through `serde-wasm-bindgen`.
- Wasm is tested with node (determinism against native) and playwright (a full search in
  headless Chromium) instead of `wasm-bindgen-test`.
- Step sizes keep adapting throughout the run instead of being frozen for a production phase,
  since the aim is optimisation rather than sampling.
- The trimer regression test checks a threshold rather than an `insta` snapshot.

### Not yet done

- Pruning candidates with an analytical upper bound on φ, and parallel tempering in pressure;
  successive halving has been sufficient for the regression shapes.
- Candidates are de-duplicated by origin shifts only, not by changes of basis (e.g. p2 sites
  ab, ac and ad are all searched).
- Concave polygons are split into triangles; merging them into larger convex pieces
  (Hertel–Mehlhorn) would speed up their overlap checks.
- Per pair, the separating axis test is slower than the old line-crossing test for polygons with
  many sides (1.1 µs vs 0.08 µs at 64 sides when overlapping), though the old test missed
  containment and collinear overlaps. A rotating calipers or GJK test would help shapes with many
  vertices.
- Mirror images of chiral shapes (the "flips" of the original C code) aren't supported.
- The MC speed-up is 9× rather than the targeted 10×; the overall time to solution is far better
  since the previous optimiser never reached the optimum.
