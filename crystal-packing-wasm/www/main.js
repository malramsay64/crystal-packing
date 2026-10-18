// Drive searches over the plane groups, running the jobs of each round on a pool of workers, and
// animate the optimisation of a single candidate.
import init, { WasmSearch, renderSvg, version } from "./pkg/crystal_packing_wasm.js";

const form = document.querySelector("#controls");
const status = document.querySelector("#status");
const progressBar = document.querySelector("#progress");
const buttons = {
  search: document.querySelector("#search"),
  survey: document.querySelector("#survey"),
  watch: document.querySelector("#watch"),
  stop: document.querySelector("#stop"),
};

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

function searchConfig(data, strategy) {
  const config = {
    strategy,
    max_sites: Number(data.get("max_sites")),
    replicas: Number(data.get("replicas")),
    steps: Number(data.get("steps")),
    seed: Number(data.get("seed")),
    within: Number(data.get("within")) / 100,
  };
  const groups = String(data.get("groups"))
    .split(",")
    .map((g) => g.trim())
    .filter(Boolean);
  if (groups.length) config.groups = groups;
  if (data.get("potential") === "lj") config.potential = {};
  return config;
}

function label(candidate) {
  return `${candidate.group}:${candidate.sites.join("")}`;
}

function formatEnergy(structure) {
  return structure?.energy != null ? structure.energy.toFixed(4) : "—";
}

function cell(row, value) {
  const td = document.createElement("td");
  td.textContent = value;
  row.append(td);
}

function setRunning(running) {
  buttons.search.disabled = running;
  buttons.survey.disabled = running;
  buttons.watch.disabled = running;
  buttons.stop.disabled = !running;
}

function startPool(size) {
  workers = Array.from({ length: size }, () => new Worker("./worker.js", { type: "module" }));
}

function stopPool() {
  for (const worker of workers) worker.terminate();
  workers = [];
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
        if (type !== "done") return;
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
  const body = document.querySelector("#ranking tbody");
  body.replaceChildren(
    ...results.slice(0, 15).map((result, index) => {
      const row = document.createElement("tr");
      cell(row, index + 1);
      cell(row, label(result.candidate));
      cell(row, result.best ? result.best.packing_fraction.toFixed(5) : "—");
      cell(row, formatEnergy(result.best));
      cell(row, result.best ? result.best.num_shapes : "—");
      return row;
    }),
  );
}

function score(structure) {
  return structure.energy != null ? -structure.energy : structure.packing_fraction;
}

function percentBelow(structure, best) {
  return (100 * (score(best) - score(structure))) / Math.abs(score(best));
}

function showNearBest(shape, near, best, within) {
  document.querySelector("#near-section").hidden = false;
  document.querySelector("#near-title").textContent =
    `${near.length} structure${near.length === 1 ? "" : "s"} within ${within}% of the best`;
  document.querySelector("#near-best").replaceChildren(
    ...near.map((result) => {
      const card = document.createElement("div");
      card.className = "card";
      card.dataset.phi = result.best.packing_fraction;
      card.innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(result.best), 1);
      const caption = document.createElement("p");
      const energy = result.best.energy != null ? `, energy ${formatEnergy(result.best)}` : "";
      caption.textContent = `${label(result.candidate)}: φ ${result.best.packing_fraction.toFixed(5)}${energy}, ${percentBelow(result.best, best).toFixed(2)}% below the best`;
      card.append(caption);
      return card;
    }),
  );
}

function showGroups(perGroup, best) {
  document.querySelector("#groups-section").hidden = false;
  document.querySelector("#groups tbody").replaceChildren(
    ...perGroup.map(([group, result]) => {
      const row = document.createElement("tr");
      const structure = result?.best;
      cell(row, group);
      cell(row, result ? result.candidate.sites.join("") : "—");
      cell(row, structure ? structure.packing_fraction.toFixed(5) : "—");
      cell(row, formatEnergy(structure));
      cell(row, structure ? percentBelow(structure, best).toFixed(2) : "—");
      return row;
    }),
  );
}

async function search(strategy) {
  const data = new FormData(form);
  const shape = shapeSpec(data);
  const within = Number(data.get("within"));
  let searcher;
  try {
    searcher = new WasmSearch(JSON.stringify(shape), JSON.stringify(searchConfig(data, strategy)));
  } catch (error) {
    status.textContent = `Error: ${error}`;
    return;
  }
  cancelled = false;
  delete document.body.dataset.result;
  setRunning(true);
  startPool(navigator.hardwareConcurrency || 4);
  document.querySelector("#results-section").hidden = false;
  document.querySelector("#groups-section").hidden = true;
  document.querySelector("#near-section").hidden = true;
  const start = performance.now();

  while (!searcher.isFinished() && !cancelled) {
    const jobs = JSON.parse(searcher.jobs());
    status.textContent = `Round ${searcher.round() + 1} (${searcher.phase()}): ${searcher.numActive()} of ${searcher.numCandidates()} candidates`;
    const outcomes = await runJobs(jobs, (done, total) => {
      progressBar.value = done / total;
    });
    if (cancelled) break;
    searcher.submit(JSON.stringify(outcomes));
    showRanking(JSON.parse(searcher.ranking()));
    const best = JSON.parse(searcher.best());
    if (best) {
      document.querySelector("#structure").innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(best), 1);
    }
  }
  stopPool();
  setRunning(false);

  const best = JSON.parse(searcher.best());
  if (cancelled || !best) {
    status.textContent = "Stopped";
  } else {
    const seconds = ((performance.now() - start) / 1000).toFixed(1);
    const summary = best.energy != null ? `Lowest energy ${best.energy.toFixed(6)}, packing fraction` : "Best packing fraction";
    status.textContent = `${summary} ${best.packing_fraction.toFixed(6)} in ${best.group} (${best.sites.map((s) => s.letter).join("")}) after ${seconds} s`;
    showNearBest(shape, JSON.parse(searcher.nearBest(within)), best, within);
    if (strategy === "exhaustive") showGroups(JSON.parse(searcher.bestPerGroup()), best);
    document.body.dataset.result = JSON.stringify(best);
  }
  searcher.free();
}

// Draw the packing fraction against the number of steps
function drawChart(history, totalSteps) {
  const chart = document.querySelector("#watch-chart");
  const values = history.map((h) => h.phi);
  const low = Math.min(...values);
  const high = Math.max(...values, low + 1e-6);
  const x = (steps) => 30 + (360 * steps) / Math.max(totalSteps, 1);
  const y = (phi) => 180 - (160 * (phi - low)) / (high - low);
  const points = history.map((h) => `${x(h.steps).toFixed(1)},${y(h.phi).toFixed(1)}`).join(" ");
  chart.innerHTML = `
    <line x1="30" y1="180" x2="390" y2="180" stroke="currentColor" opacity="0.3" />
    <line x1="30" y1="20" x2="30" y2="180" stroke="currentColor" opacity="0.3" />
    <text x="0" y="24">${high.toFixed(3)}</text>
    <text x="0" y="184">${low.toFixed(3)}</text>
    <text x="330" y="196">steps</text>
    <polyline points="${points}" />`;
}

async function watch() {
  const data = new FormData(form);
  const shape = shapeSpec(data);
  const config = searchConfig(data, "halving");
  const candidate = {
    group: String(data.get("watch_group")).trim(),
    sites: [...String(data.get("watch_sites")).trim()],
  };
  const job = {
    shape,
    candidate,
    seed: config.seed,
    config: {
      steps: config.steps,
      pressure_start: 0.1,
      pressure_end: 1e6,
      target_acceptance: 0.3,
      adapt_interval: 100,
      squeeze_steps: 2000,
    },
    start: null,
    potential: config.potential ?? null,
  };
  cancelled = false;
  delete document.body.dataset.watchResult;
  setRunning(true);
  startPool(1);
  const view = document.querySelector("#watch-view");
  const watchStatus = document.querySelector("#watch-status");
  document.querySelector("#watch-section").hidden = false;
  document.querySelector("#watch-label").textContent = label(candidate);
  view.dataset.frames = "0";
  view.replaceChildren();
  const history = [];

  const outcome = await new Promise((resolve) => {
    workers[0].onmessage = (event) => {
      const { type } = event.data;
      if (type === "progress") {
        const { progress, structure } = event.data;
        history.push({ steps: progress.steps, phi: progress.packing_fraction });
        drawChart(history, progress.total_steps);
        if (structure) {
          view.innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(structure), 1);
          view.dataset.frames = String(Number(view.dataset.frames) + 1);
        }
        progressBar.value = progress.steps / progress.total_steps;
        watchStatus.textContent = `Step ${progress.steps} of ${progress.total_steps}, pressure ${progress.pressure.toPrecision(3)}, φ ${progress.packing_fraction.toFixed(4)} (best ${progress.best_packing_fraction.toFixed(4)})`;
      } else if (type === "done") {
        resolve(event.data.outcome);
      }
    };
    workers[0].postMessage({ id: 0, job, watch: true });
  });
  stopPool();
  setRunning(false);
  if (outcome.structure) {
    view.innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(outcome.structure), 1);
    watchStatus.textContent = `Finished: φ ${outcome.structure.packing_fraction.toFixed(6)} after compressing to contact`;
    document.body.dataset.watchResult = JSON.stringify(outcome.structure);
  } else {
    watchStatus.textContent = `Error: ${outcome.error}`;
  }
}

buttons.search.addEventListener("click", () => search("halving"));
buttons.survey.addEventListener("click", () => search("exhaustive"));
buttons.watch.addEventListener("click", watch);
buttons.stop.addEventListener("click", () => {
  cancelled = true;
  stopPool();
  setRunning(false);
});

await init();
document.querySelector("#version").textContent = version();
document.body.dataset.ready = "true";
