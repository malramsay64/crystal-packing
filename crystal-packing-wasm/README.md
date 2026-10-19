# crystal-packing-wasm

WebAssembly bindings and a browser front end for the crystal packing search.

There are two pages in `www`:

- `index.html` runs a shape in every structure at once. Choose a disk, a regular polygon, a trimer
  or draw your own shape, then run it in all 17 plane groups, in a strip periodic in one direction
  and as a cluster in a square. Every optimisation is shown as it runs, only the first row of plane
  groups is shown until expanded, and the distinct structures within 2% of the best are listed at
  the end with the strip and the cluster. The styling follows the
  [Serene](https://github.com/isunjn/serene) theme for Zola, in light and dark mode, using the
  theme's own CSS variables and session theme preference so it fits in a Serene blog. The theme's
  defaults are wrapped in `:where()`, so the blog's stylesheet takes precedence when it's loaded.
- `explorer.html` has every option of the search.

The explorer can search for the densest packing, survey every plane group, or pack independent shapes
without symmetry, in a periodic cell, a strip periodic in one direction, a cluster with no boundaries, or a square, rectangular or circular container. Every worker's current optimisation is shown as it runs, and the results list
the distinct structures within a chosen percentage of the best (2% by default), best first, with
the other candidates which found the same structure, updating as the
percentage or the maximum number of results is changed. Shapes can be drawn by hand: click to add
vertices, drag to move them and double click to remove them, with the shape checked as it is drawn.
The optimisation of a single group and set of sites can be animated, with a plot of the packing
fraction as it is compressed.

The search runs on the page's main thread and hands each round's jobs to a pool of web workers,
one per core, which run the Monte Carlo optimisation in slices so progress can be reported.
All values cross the JavaScript boundary as JSON, using the same format as the command line tool.

## Building

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-bindgen-cli --version <version of wasm-bindgen in Cargo.lock>
./build.sh
python3 -m http.server --directory www
```

then open <http://localhost:8000>. Module workers are used, so a recent browser is required.

## Tests

`./build.sh node` additionally builds bindings for node, which are used by
`tests/determinism.mjs` to check the WebAssembly build gives results identical to the native
build for the same seed. `tests/browser.mjs` uses playwright to run a search and watch an optimisation in headless
Chromium with the explorer, checking the structures near the best and the animation, and
`tests/dashboard.mjs` runs pentagons in every structure with `index.html` (set `SCREENSHOTS` to a
directory to save screenshots in light and dark mode). Playwright is pinned in
`package-lock.json`:

```bash
./build.sh node
npm ci && npx playwright install chromium
npm run test:browser
```
