// Check the WebAssembly build gives identical results to the native build.
//
// Usage: node tests/determinism.mjs reference.json
// where reference.json is the output of `cargo run --example determinism -p crystal_packing`.
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const { runJob } = require("../node-pkg/crystal_packing_wasm.js");

// Serialise with sorted keys, so only the values are compared
function canonical(value) {
  if (Array.isArray(value)) return `[${value.map(canonical).join(",")}]`;
  if (value && typeof value === "object") {
    return `{${Object.keys(value)
      .sort()
      .map((k) => `${JSON.stringify(k)}:${canonical(value[k])}`)
      .join(",")}}`;
  }
  return JSON.stringify(value);
}

const reference = JSON.parse(readFileSync(process.argv[2], "utf8"));
let failures = 0;
reference.jobs.forEach((job, i) => {
  const outcome = JSON.parse(runJob(JSON.stringify(job)));
  const expected = reference.outcomes[i];
  const label = `${job.candidate.group}:${job.candidate.sites.join("")}`;
  if (canonical(outcome) === canonical(expected)) {
    console.log(`ok   ${label} phi=${outcome.structure.packing_fraction}`);
  } else {
    failures += 1;
    console.log(`FAIL ${label}\n  wasm:   ${JSON.stringify(outcome)}\n  native: ${JSON.stringify(expected)}`);
  }
});
if (failures > 0) {
  process.exit(1);
}
