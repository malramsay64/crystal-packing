#!/usr/bin/env bash
# Build the WebAssembly module and JavaScript bindings.
#
# Requires the wasm32-unknown-unknown target and a wasm-bindgen-cli matching the version of the
# wasm-bindgen crate in Cargo.lock. wasm-opt (from binaryen) is used when available.
#
#   ./build.sh          builds www/pkg for the browser
#   ./build.sh node     also builds node-pkg, used by the determinism test
set -euo pipefail
cd "$(dirname "$0")"

cargo build --release --target wasm32-unknown-unknown -p crystal-packing-wasm
WASM=../target/wasm32-unknown-unknown/release/crystal_packing_wasm.wasm

if command -v wasm-opt >/dev/null; then
    wasm-opt -O3 --enable-bulk-memory --enable-nontrapping-float-to-int "$WASM" -o "$WASM"
fi

wasm-bindgen --target web --out-dir www/pkg "$WASM"
if [[ "${1:-}" == "node" ]]; then
    wasm-bindgen --target nodejs --out-dir node-pkg "$WASM"
fi
