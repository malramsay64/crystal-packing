// A web worker running packing jobs. Each message contains a job, which is run in slices so
// progress can be reported back to the page. When watching, the current structure is sent with
// the progress after each slice so the optimisation can be animated.
import init, { WasmJob } from "./pkg/crystal_packing_wasm.js";

const ready = init();

self.onmessage = async (event) => {
  await ready;
  const { id, job, watch } = event.data;
  const slice = watch ? Math.max(100, Math.round(job.config.steps / 100)) : 2000;
  let outcome;
  try {
    const runner = new WasmJob(JSON.stringify(job));
    while (!runner.isComplete()) {
      const progress = JSON.parse(runner.run(slice));
      const message = { id, type: "progress", progress };
      if (watch) message.structure = JSON.parse(runner.current());
      self.postMessage(message);
    }
    outcome = JSON.parse(runner.finish());
    runner.free();
  } catch (error) {
    outcome = { candidate: job.candidate, structure: null, error: String(error) };
  }
  self.postMessage({ id, type: "done", outcome });
};
