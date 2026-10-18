// Drive a search over the plane groups, running the jobs of each round on a pool of workers.
import init, { WasmSearch, renderSvg, version } from "./pkg/crystal_packing_wasm.js";

const form = document.querySelector("#controls");
const status = document.querySelector("#status");
const ranking = document.querySelector("#ranking tbody");
const image = document.querySelector("#structure");
const stop = document.querySelector("#stop");
const progress = document.querySelector("#progress");

let workers = [];
let cancelled = false;

function shapeSpec(data) {
  switch (data.get("shape")) {
    case "disk":
      return { type: "disk" };
    case "trimer":
      return {
        type: "trimer",
        radius: Number(data.get("radius")),
        angle: Number(data.get("angle")),
        distance: Number(data.get("distance")),
      };
    default:
      return { type: "polygon", sides: Number(data.get("sides")) };
  }
}

function createPool(size) {
  return Array.from({ length: size }, () => new Worker("./worker.js", { type: "module" }));
}

// Run the jobs on the pool, resolving with the outcomes in the order of the jobs
function runJobs(jobs, onProgress) {
  return new Promise((resolve) => {
    const outcomes = new Array(jobs.length);
    let next = 0;
    let finished = 0;
    const dispatch = (worker) => {
      if (cancelled || next >= jobs.length) return;
      const id = next++;
      worker.postMessage({ id, job: jobs[id] });
    };
    for (const worker of workers) {
      worker.onmessage = (event) => {
        const { id, type } = event.data;
        if (type === "progress") return;
        outcomes[id] = event.data.outcome;
        finished += 1;
        onProgress(finished, jobs.length);
        if (finished === jobs.length || cancelled) {
          resolve(outcomes.filter(Boolean));
        } else {
          dispatch(worker);
        }
      };
      dispatch(worker);
    }
  });
}

function showRanking(results) {
  ranking.replaceChildren(
    ...results.slice(0, 15).map((result, index) => {
      const row = document.createElement("tr");
      const phi = result.best ? result.best.packing_fraction.toFixed(5) : "—";
      const shapes = result.best ? result.best.num_shapes : "—";
      for (const value of [index + 1, `${result.candidate.group}:${result.candidate.sites.join("")}`, phi, shapes]) {
        const cell = document.createElement("td");
        cell.textContent = value;
        row.append(cell);
      }
      return row;
    }),
  );
}

async function search(event) {
  event.preventDefault();
  const data = new FormData(form);
  const shape = shapeSpec(data);
  const config = {
    max_sites: Number(data.get("max_sites")),
    replicas: Number(data.get("replicas")),
    steps: Number(data.get("steps")),
    seed: Number(data.get("seed")),
  };
  const groups = String(data.get("groups")).split(",").map((g) => g.trim()).filter(Boolean);
  if (groups.length) config.groups = groups;

  let searcher;
  try {
    searcher = new WasmSearch(JSON.stringify(shape), JSON.stringify(config));
  } catch (error) {
    status.textContent = `Error: ${error}`;
    return;
  }
  cancelled = false;
  workers = createPool(navigator.hardwareConcurrency || 4);
  form.querySelector("button[type=submit]").disabled = true;
  stop.disabled = false;
  const start = performance.now();

  while (!searcher.isFinished() && !cancelled) {
    const jobs = JSON.parse(searcher.jobs());
    const round = searcher.round() + 1;
    status.textContent = `Round ${round}: ${searcher.numActive()} of ${searcher.numCandidates()} candidates`;
    const outcomes = await runJobs(jobs, (done, total) => {
      progress.value = done / total;
    });
    if (cancelled) break;
    searcher.submit(JSON.stringify(outcomes));
    showRanking(JSON.parse(searcher.ranking()));
    const best = JSON.parse(searcher.best());
    if (best) image.innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(best), 1);
  }

  const seconds = ((performance.now() - start) / 1000).toFixed(1);
  const best = JSON.parse(searcher.best());
  status.textContent = cancelled
    ? "Stopped"
    : `Best packing fraction ${best.packing_fraction.toFixed(6)} in ${best.group} (${best.sites
        .map((s) => s.letter)
        .join("")}) after ${seconds} s`;
  document.body.dataset.result = JSON.stringify(best);
  for (const worker of workers) worker.terminate();
  searcher.free();
  form.querySelector("button[type=submit]").disabled = false;
  stop.disabled = true;
}

stop.addEventListener("click", () => {
  cancelled = true;
  for (const worker of workers) worker.terminate();
});

await init();
document.querySelector("#version").textContent = version();
form.addEventListener("submit", search);
document.body.dataset.ready = "true";
