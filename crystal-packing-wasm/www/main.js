// Drive searches over the plane groups, running the jobs of each round on a pool of workers and
// showing every optimiser as it runs, animate the optimisation of a single candidate, and draw
// arbitrary shapes to pack.
import init, { WasmSearch, renderSvg, shapeInfo, version } from "./pkg/crystal_packing_wasm.js";

const form = document.querySelector("#controls");
const status = document.querySelector("#status");
const progressBar = document.querySelector("#progress");
const buttons = {
  search: document.querySelector("#search"),
  survey: document.querySelector("#survey"),
  watch: document.querySelector("#watch"),
  unconstrained: document.querySelector("#unconstrained"),
  stop: document.querySelector("#stop"),
};
const runButtons = [buttons.search, buttons.survey, buttons.watch, buttons.unconstrained];

let workers = [];
let cancelled = false;
let running = false;
// The shape and results of the latest search, kept so the results can be filtered again
let lastSearch = null;

// --- Drawing a shape -------------------------------------------------------------------------

const drawArea = document.querySelector("#draw-area");
const drawStatus = document.querySelector("#draw-status");
let drawnVertices = [];
let drawnValid = false;

/// The vertices moved so the area centroid is at the origin
function centred(vertices) {
  let area = 0;
  let cx = 0;
  let cy = 0;
  vertices.forEach(([x0, y0], i) => {
    const [x1, y1] = vertices[(i + 1) % vertices.length];
    const cross = x0 * y1 - x1 * y0;
    area += cross / 2;
    cx += ((x0 + x1) * cross) / 6;
    cy += ((y0 + y1) * cross) / 6;
  });
  if (Math.abs(area) < 1e-12) return vertices;
  return vertices.map(([x, y]) => [x - cx / area, y - cy / area]);
}

function drawnShape() {
  return { type: "vertices", name: "Drawn", vertices: centred(drawnVertices) };
}

function validateDrawing() {
  if (drawnVertices.length < 3) {
    drawnValid = false;
    drawStatus.textContent = "Add at least three vertices";
    return;
  }
  try {
    const info = JSON.parse(shapeInfo(JSON.stringify(drawnShape())));
    drawnValid = true;
    drawStatus.textContent = `Valid polygon, ${drawnVertices.length} vertices, area ${info.area.toFixed(3)}, ${info.convex ? "convex" : "concave"}, ${info.rotations}-fold rotation, ${info.mirrors} mirror${info.mirrors === 1 ? "" : "s"}`;
  } catch (error) {
    drawnValid = false;
    drawStatus.textContent = `Invalid shape: ${error.message ?? error}`;
  }
}

function renderDrawing() {
  const grid = [];
  for (let i = -20; i <= 20; i += 1) {
    const v = i / 10;
    grid.push(`<line class="${i === 0 ? "axis" : "grid"}" x1="${v}" y1="-2" x2="${v}" y2="2" />`);
    grid.push(`<line class="${i === 0 ? "axis" : "grid"}" x1="-2" y1="${v}" x2="2" y2="${v}" />`);
  }
  const points = drawnVertices.map(([x, y]) => `${x},${-y}`).join(" ");
  const vertices = drawnVertices
    .map(([x, y], i) => `<circle data-index="${i}" cx="${x}" cy="${-y}" r="0.05" />`)
    .join("");
  const shape = drawnVertices.length >= 2 ? `<polygon class="${drawnValid ? "" : "invalid"}" points="${points}" />` : "";
  drawArea.innerHTML = grid.join("") + shape + vertices;
}

function updateDrawing() {
  validateDrawing();
  renderDrawing();
  updateButtons();
}

/// The position of a pointer event in the drawing coordinates, with y upwards
function drawingPoint(event) {
  const matrix = drawArea.getScreenCTM().inverse();
  const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix);
  let [x, y] = [point.x, -point.y];
  if (document.querySelector("input[name=snap]").checked) {
    [x, y] = [Math.round(x * 10) / 10, Math.round(y * 10) / 10];
  }
  return [x, y];
}

let dragging = null;
let moved = false;
drawArea.addEventListener("pointerdown", (event) => {
  const index = event.target.dataset?.index;
  if (index !== undefined) {
    dragging = Number(index);
    moved = false;
    drawArea.setPointerCapture(event.pointerId);
  }
});
drawArea.addEventListener("pointermove", (event) => {
  if (dragging === null) return;
  moved = true;
  drawnVertices[dragging] = drawingPoint(event);
  updateDrawing();
});
drawArea.addEventListener("pointerup", (event) => {
  if (dragging !== null) {
    dragging = null;
    return;
  }
  if (event.target.dataset?.index === undefined) {
    drawnVertices.push(drawingPoint(event));
    updateDrawing();
  }
});
drawArea.addEventListener("dblclick", (event) => {
  const index = event.target.dataset?.index;
  if (index !== undefined && !moved) {
    drawnVertices.splice(Number(index), 1);
    updateDrawing();
  }
});
document.querySelector("#draw-undo").addEventListener("click", () => {
  drawnVertices.pop();
  updateDrawing();
});
document.querySelector("#draw-clear").addEventListener("click", () => {
  drawnVertices = [];
  updateDrawing();
});

// --- Shapes and configuration ----------------------------------------------------------------

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
    case "drawn":
      return drawnShape();
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

function updateButtons() {
  const drawing = form.elements.shape.value === "drawn";
  document.querySelector("#draw-section").hidden = !drawing;
  for (const button of runButtons) button.disabled = running || (drawing && !drawnValid);
  buttons.stop.disabled = !running;
}

form.elements.shape.addEventListener("change", updateButtons);

function setRunning(state) {
  running = state;
  updateButtons();
}

// --- Display helpers -------------------------------------------------------------------------

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

function score(structure) {
  return structure.energy != null ? -structure.energy : structure.packing_fraction;
}

function percentBelow(structure, best) {
  return (100 * (score(best) - score(structure))) / Math.abs(score(best));
}

// --- Workers ---------------------------------------------------------------------------------

function startPool(size) {
  workers = Array.from({ length: size }, () => new Worker("./worker.js", { type: "module" }));
}

function stopPool() {
  for (const worker of workers) worker.terminate();
  workers = [];
}

/// A tile for each worker showing the job it is running as it is optimised
function createTiles() {
  const live = document.querySelector("#live");
  document.querySelector("#live-section").hidden = false;
  live.replaceChildren(
    ...workers.map((_, i) => {
      const tile = document.createElement("div");
      tile.className = "tile";
      tile.dataset.frames = "0";
      tile.innerHTML = `<strong>Worker ${i + 1}</strong><div class="view"></div><p class="caption muted">Waiting</p>`;
      return tile;
    }),
  );
  return [...live.querySelectorAll(".tile")];
}

/// Draw the latest structure of each tile once per animation frame
function tileRenderer(shape, tiles) {
  const pending = new Map();
  let scheduled = false;
  const draw = () => {
    scheduled = false;
    for (const [index, { job, progress, structure }] of pending) {
      const tile = tiles[index];
      tile.querySelector(".view").innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(structure), 1);
      tile.querySelector(".caption").textContent = `${label(job.candidate)}: φ ${progress.packing_fraction.toFixed(4)}, step ${progress.steps} of ${progress.total_steps}`;
      tile.dataset.frames = String(Number(tile.dataset.frames) + 1);
    }
    pending.clear();
  };
  return (index, update) => {
    pending.set(index, update);
    if (!scheduled) {
      scheduled = true;
      requestAnimationFrame(draw);
    }
  };
}

// Run the jobs on the pool, resolving with the outcomes in the order of the jobs
function runJobs(jobs, onProgress, onFrame) {
  return new Promise((resolve) => {
    const outcomes = new Array(jobs.length);
    let next = 0;
    let finished = 0;
    const dispatch = (worker) => {
      if (cancelled || next >= jobs.length) return;
      const id = next++;
      worker.postMessage({ id, job: jobs[id], frames: 20 });
    };
    workers.forEach((worker, index) => {
      worker.onmessage = (event) => {
        const { id, type } = event.data;
        if (type === "progress") {
          if (event.data.structure) {
            onFrame(index, { job: jobs[id], progress: event.data.progress, structure: event.data.structure });
          }
          return;
        }
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
    });
  });
}

// --- Results ---------------------------------------------------------------------------------

/// Show the results within the chosen percentage of the best, best first. By default the
/// candidates which found the same crystal are grouped into one distinct structure.
function showResults() {
  if (!lastSearch) return;
  const { shape, ranking, searcher } = lastSearch;
  const data = new FormData(form);
  const within = Number(data.get("within"));
  const limit = data.get("max_results") ? Number(data.get("max_results")) : Infinity;
  const showEquivalent = data.get("show_equivalent") !== null;
  const found = ranking.filter((r) => r.best);
  if (!found.length) return;
  const best = found[0].best;
  const threshold = score(best) - (within / 100) * Math.abs(score(best));
  const near = found.filter((r) => score(r.best) >= threshold - 1e-12);
  const distinct = JSON.parse(searcher.distinctNearBest(within));
  const rows = showEquivalent
    ? near.map((result) => ({ result, equivalent: [] }))
    : distinct;
  const shown = rows.slice(0, limit);

  const title = document.querySelector("#results-title");
  title.dataset.count = String(near.length);
  title.dataset.distinct = String(distinct.length);
  title.dataset.candidates = String(ranking.length);
  const limited = shown.length < rows.length ? `, showing the best ${shown.length}` : "";
  title.textContent = `${near.length} of ${ranking.length} results within ${within}% of the best, ${distinct.length} distinct structure${distinct.length === 1 ? "" : "s"}${limited}`;
  document.querySelector("#results tbody").replaceChildren(
    ...shown.map(({ result, equivalent }, index) => {
      const row = document.createElement("tr");
      const also = equivalent.map((r) => label(r.candidate)).join(" ");
      row.dataset.phi = String(result.best.packing_fraction);
      row.dataset.label = label(result.candidate);
      row.dataset.equivalent = also;
      cell(row, index + 1);
      cell(row, label(result.candidate));
      cell(row, result.best.packing_fraction.toFixed(5));
      cell(row, formatEnergy(result.best));
      cell(row, percentBelow(result.best, best).toFixed(2));
      cell(row, result.best.num_shapes);
      cell(row, also);
      return row;
    }),
  );
  document.querySelector("#near-section").hidden = false;
  document.querySelector("#near-best").replaceChildren(
    ...shown.map(({ result, equivalent }) => {
      const card = document.createElement("div");
      card.className = "card";
      card.dataset.phi = result.best.packing_fraction;
      card.innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(result.best), 1);
      const caption = document.createElement("p");
      const energy = result.best.energy != null ? `, energy ${formatEnergy(result.best)}` : "";
      const also = equivalent.length ? `. Also found in ${equivalent.map((r) => label(r.candidate)).join(", ")}` : "";
      caption.textContent = `${label(result.candidate)}: φ ${result.best.packing_fraction.toFixed(5)}${energy}, ${percentBelow(result.best, best).toFixed(2)}% below the best${also}`;
      card.append(caption);
      return card;
    }),
  );
}

for (const name of ["within", "max_results", "show_equivalent"]) {
  form.elements[name].addEventListener("change", showResults);
  form.elements[name].addEventListener("input", showResults);
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

// --- Searching -------------------------------------------------------------------------------

async function search(strategy, extra = {}) {
  const data = new FormData(form);
  const shape = shapeSpec(data);
  let searcher;
  try {
    const config = { ...searchConfig(data, strategy), ...extra };
    searcher = new WasmSearch(JSON.stringify(shape), JSON.stringify(config));
  } catch (error) {
    status.textContent = `Error: ${error.message ?? error}`;
    return;
  }
  cancelled = false;
  delete document.body.dataset.result;
  setRunning(true);
  startPool(navigator.hardwareConcurrency || 4);
  const tiles = createTiles();
  const onFrame = tileRenderer(shape, tiles);
  document.querySelector("#results-section").hidden = false;
  document.querySelector("#groups-section").hidden = true;
  const start = performance.now();

  while (!searcher.isFinished() && !cancelled) {
    const jobs = JSON.parse(searcher.jobs());
    status.textContent = `Round ${searcher.round() + 1} (${searcher.phase()}): ${searcher.numActive()} of ${searcher.numCandidates()} candidates`;
    const outcomes = await runJobs(
      jobs,
      (done, total) => {
        progressBar.value = done / total;
      },
      onFrame,
    );
    if (cancelled) break;
    searcher.submit(JSON.stringify(outcomes));
    if (lastSearch && lastSearch.searcher !== searcher) lastSearch.searcher.free();
    lastSearch = { shape, searcher, ranking: JSON.parse(searcher.ranking()) };
    showResults();
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
    if (strategy === "exhaustive" && !extra.unconstrained) {
      showGroups(JSON.parse(searcher.bestPerGroup()), best);
    }
    document.body.dataset.result = JSON.stringify(best);
  }
  // The search is kept so the results can be filtered again, unless it produced no results
  if (!lastSearch || lastSearch.searcher !== searcher) searcher.free();
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
    workers[0].postMessage({ id: 0, job, frames: 100 });
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

function unconstrained() {
  const counts = String(new FormData(form).get("unconstrained_shapes"))
    .split(",")
    .map((n) => Number(n.trim()))
    .filter((n) => Number.isInteger(n) && n > 0);
  // Larger systems need a proportionally larger budget
  search("exhaustive", { unconstrained: counts, steps_per_shape: 2500 });
}

buttons.search.addEventListener("click", () => search("halving"));
buttons.survey.addEventListener("click", () => search("exhaustive"));
buttons.watch.addEventListener("click", watch);
buttons.unconstrained.addEventListener("click", unconstrained);
buttons.stop.addEventListener("click", () => {
  cancelled = true;
  stopPool();
  setRunning(false);
});

await init();
document.querySelector("#version").textContent = version();
updateDrawing();
document.body.dataset.ready = "true";
