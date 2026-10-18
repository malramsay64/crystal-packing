# Packing

[![DOI](https://zenodo.org/badge/165979598.svg)](https://zenodo.org/badge/latestdoi/165979598)

Find the densest periodic packings of 2D shapes, searching over the 17 plane (wallpaper) groups
and the ways shapes can occupy their Wyckoff positions. Runs as a command line tool or in the
browser using WebAssembly.

| Crate | Purpose |
|---|---|
| `crystal-packing` | The library: symmetry, overlap checks, optimisation and search |
| `crystal-packing-cli` | The `packing` command line tool |
| `crystal-packing-wasm` | WebAssembly bindings and a browser front end |

## Usage

```bash
cargo install --path crystal-packing-cli

# Search every plane group for the densest packing of regular pentagons
packing search --shape polygon --sides 5 --output pentagons

# Optimise every candidate of all 17 groups, reporting the best structure of each group
packing survey --shape polygon --sides 5 --output pentagon-survey

# Independent shapes in a periodic cell without any symmetry, for 1 to 4 shapes
packing unconstrained --shape polygon --sides 5 --shapes 1,2,3,4

# Or in a container with hard walls: square, rectangle (free aspect ratio) or circle
packing unconstrained --shape disk --shapes 2,3,7 --container circle

# Periodic in one direction (a strip between walls), or in none (a cluster)
packing unconstrained --shape disk --shapes 1,2 --container strip
packing unconstrained --shape disk --shapes 2,3,4 --container cluster

# Optimise a single group and choice of Wyckoff sites
packing optimise --shape polygon --sides 7 --group p2 --sites e

# Molecules of three disks, minimising the Lennard-Jones lattice energy
packing search --shape trimer --potential lj

# Render a saved structure
packing render pentagons/best.json
```

Shapes are regular polygons (`--sides`), disks, trimers (`--radius --angle --distance`), radial
shapes (`--radii 1,0.8,1,0.8`) or polygons from a JSON list of vertices (`--vertices file.json`),
which can be concave. Every command lists the structures within `--within` percent (default 2%)
of the best, best first. Candidates which found the same crystal (e.g. the pentagon double lattice
in p2, p1g1, p2mg and p2gg) are grouped into one distinct structure, listing where else it was
found; `--show-equivalent` lists them separately and `--max-results` limits how many are shown. A search writes `results.json`, `results.csv` (every candidate), and the structures near the
best to `near-best/`. A survey also writes `survey.csv` and
the best structure of each group to `groups/`.

For the browser front end see [crystal-packing-wasm](crystal-packing-wasm/README.md). It can run a
search, a survey of all groups or an unconstrained search, showing every optimiser as it runs and
every structure near the best. Shapes can be drawn by hand, and the optimisation of a single group
and set of sites can be animated.

## How it works

- **Symmetry**: the plane groups are defined by their general position and the representative
  coordinates of each Wyckoff position from the International Tables. The orbits, site symmetries
  and the orientations a shape can take on a special position are derived from these.
- **Candidates**: every multiset of Wyckoff positions compatible with the shape's symmetry (up to
  `--max-sites` sites and `--max-shapes` shapes per cell) is a candidate, with candidates related
  by a shift of the origin removed.
- **Overlaps**: polygons use the separating axis theorem on convex pieces. Using the symmetry,
  only the first image of each site needs checking against the other shapes, and moving a site
  only rechecks the pairs involving it. Periodic images are found exactly from the reduced
  lattice basis.
- **Optimisation**: an isobaric ("floppy box") Monte Carlo compression with a rising pressure and
  adaptive step sizes, followed by a squeeze to contact. Every result is checked by an
  independent validator.
- **Unconstrained**: N independent shapes in a fully flexible periodic cell, which is the general
  position of p1 occupied N times, so the same engine and optimiser are used. The budget scales
  with the number of shapes (`--steps-per-shape`, 2500 by default), and from independent shapes,
  2 to 32 pentagons recover the double lattice packing.
- **Containers**: a finite number of shapes in a square, rectangle or circle with hard walls,
  compressed around them, reproducing known optima such as 7 disks in a circle.
- **Fewer periodic directions**: a strip is periodic along its length with walls either side,
  with both the period and the width compressed, so two disks per period find the zig-zag at
  π/(2 + √3). A cluster has no walls or periodic boundaries, and the packing fraction uses the area
  of the convex hull of the shapes. This is exact: A + Pr + πr² from the hull of the centres for
  polygons and equal disks, and from the support function of the hull for atoms of different sizes. Moving a shape on the edge of a cluster changes the area, so these moves
  are accepted with the same pressure as the compression.
- **Distinct structures**: structures are compared by the distances between the vertices of each
  shape and its neighbours, which doesn't depend on the group, cell, origin or orientation.
- **Search**: successive halving. Every candidate gets a short optimisation, then the best third
  get three times as many steps, until one candidate remains, which is optimised until it stops
  improving. Every candidate within `--within` of the best then gets the same budget as the best,
  so the structures reported near the best are fairly compared. A survey instead optimises every
  candidate each round. Jobs are independent with their own seeds, so results are reproducible
  and identical natively and in WebAssembly.

## Results

A full search over all 17 plane groups, with default settings on 4 cores:

| Shape | Packing fraction | Known optimum | Time |
|---|---|---|---|
| Disk | 0.906900 | π/√12 = 0.906900 | 1.3 s |
| Triangle | 1.000000 | 1 | 3.5 s |
| Square | 1.000000 | 1 | 5.3 s |
| Pentagon | 0.921311 | (5−√5)/3 = 0.921311 | 1.8 s |
| Hexagon | 1.000000 | 1 | 3.0 s |
| Heptagon | 0.892691 | 0.892691 | 2.6 s |
| Trimer | 0.946284 | — | 1.0 s |

A survey optimising every candidate of every group takes longer, about 22 s for pentagons.

For comparison, the previous optimiser reached 0.630 for pentagons, 0.767 for squares, 0.788 for
disks and 0.449 for triangles, and didn't improve with more steps.

## Testing

The minimum supported Rust version is 1.88. CI checks formatting, clippy (native and wasm) and the
documentation with warnings as errors, runs every test suite, and audits the dependencies.

- Unit and property tests compare every overlap check (convex, concave and reflected polygons,
  molecules, the symmetry-reduced crystal checks, containers, strips and clusters) with the brute
  force checks of the independent validator, check the
  symmetry data against the International Tables, and test the search scheduling with fake results.
- `tests/regression.rs`: the known densest packings above.
- `tests/adversarial.rs`: shapes chosen to stress the optimisation, most with known answers since
  they tile the plane: every triangle and quadrilateral (including concave ones) in p2, the
  L-tromino and plus pentomino, a 20:1 rectangle, shapes with their origin at a vertex or outside
  the shape, a 64-gon, a pentagram and a trimer of heavily overlapping atoms.
- `tests/unconstrained.rs`: independent shapes find the known packings without symmetry,
  including systems of 12 pentagons and 24 disks.
- `tests/confined.rs`: known optimal packings of disks, squares and triangles in containers,
  strips and clusters.
- `tests/distinct.rs`: full searches find a single distinct structure near the best for pentagons
  and disks, and three for heptagons.
- `tests/known_geometries.rs`: known densest arrangements are constructed by hand (the hexagonal
  disk packing, square lattice, triangle tiling and honeycomb, disks and squares in containers,
  the zig-zag of disks in a strip and the triangle and rhombus clusters of disks)
  and the structures the search finds are compared with them, not just their packing fractions.
- `crystal-packing-cli/tests`: run the `packing` binary.
- `crystal-packing-wasm/tests`: the WebAssembly build matches the native build exactly, and the
  web front end runs a search with live views of every worker, filters the results near the best,
  animates an optimisation, rejects and accepts drawn shapes, and runs an unconstrained search.

```bash
cargo test --workspace
# The slower suites, in release mode
cargo test --release --workspace -- --include-ignored
cargo bench -p crystal_packing --bench engine
```
