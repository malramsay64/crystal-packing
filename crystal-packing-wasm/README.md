# crystal-packing-wasm

WebAssembly bindings and a browser front end for the crystal packing search.

The page can search for the densest packing, survey every plane group, or pack independent shapes
without symmetry. Every worker's current optimisation is shown as it runs, and the results list
every structure within a chosen percentage of the best (2% by default), best first, updating as the
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
Chromium, checking the structures near the best and the animation.
