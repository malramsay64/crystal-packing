// A web worker running packing jobs. Each message contains a job, which is run in slices so
// progress can be reported back to the page.
import init, { WasmJob } from "./pkg/crystal_packing_wasm.js";

const ready = init();
const SLICE = 2000;

self.onmessage = async (event) => {
  await ready;
  const { id, job } = event.data;
  let outcome;
  try {
    const runner = new WasmJob(JSON.stringify(job));
    while (!runner.isComplete()) {
      const progress = JSON.parse(runner.run(SLICE));
      self.postMessage({ id, type: "progress", progress });
    }
    outcome = JSON.parse(runner.finish());
    runner.free();
  } catch (error) {
    outcome = { candidate: job.candidate, structure: null, error: String(error) };
  }
  self.postMessage({ id, type: "done", outcome });
};
