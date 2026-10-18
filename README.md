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

# Optimise a single group and choice of Wyckoff sites
packing optimise --shape polygon --sides 7 --group p2 --sites e

# Molecules of three disks, minimising the Lennard-Jones lattice energy
packing search --shape trimer --potential lj

# Render a saved structure
packing render pentagons/best.json
```

Shapes are regular polygons (`--sides`), disks, trimers (`--radius --angle --distance`), radial
shapes (`--radii 1,0.8,1,0.8`) or polygons from a JSON list of vertices (`--vertices file.json`),
which can be concave. A search writes `results.json`, `results.csv`, and JSON and SVG files of
the best structures to the output directory.

For the browser front end see [crystal-packing-wasm](crystal-packing-wasm/README.md).

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
- **Search**: successive halving. Every candidate gets a short optimisation, then the best third
  get three times as many steps, until one candidate remains. Jobs are independent with their own
  seeds, so results are reproducible and identical natively and in WebAssembly.

## Results

A full search over all 17 plane groups, with default settings on 4 cores:

| Shape | Packing fraction | Known optimum | Time |
|---|---|---|---|
| Disk | 0.906900 | π/√12 = 0.906900 | 0.7 s |
| Triangle | 1.000000 | 1 | 2.1 s |
| Square | 1.000000 | 1 | 1.2 s |
| Pentagon | 0.921311 | (5−√5)/3 = 0.921311 | 2.0 s |
| Hexagon | 1.000000 | 1 | 1.5 s |
| Heptagon | 0.892691 | 0.892691 | 3.4 s |
| Trimer | 0.946284 | — | 1.0 s |

For comparison, the previous optimiser reached 0.630 for pentagons, 0.767 for squares, 0.788 for
disks and 0.449 for triangles, and didn't improve with more steps.

## Development

```bash
cargo test --workspace
# The known packings above, in release mode
cargo test --release -p crystal_packing --test regression -- --ignored
cargo bench -p crystal_packing --bench engine
```
