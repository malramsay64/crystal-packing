# crystal-packing-wasm

WebAssembly bindings and a browser front end for the crystal packing search.

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
build for the same seed. `tests/browser.mjs` runs a small search in headless Chromium using
playwright.
