// Run a shape in every structure at once: a survey of all 17 plane groups, then independent
// shapes in a strip and in a square, showing every optimisation as it runs and listing the
// densest structures at the end.
import init, { WasmSearch, renderSvg, version } from "./pkg/crystal_packing_wasm.js";
import { centred, drawingPad } from "./draw.js";

const SHAPES = [
  { key: "disk", name: "Disk", note: "circle", spec: { type: "disk" } },
  { key: "triangle", name: "Triangle", note: "3 sides", spec: { type: "polygon", sides: 3 } },
  { key: "square", name: "Square", note: "4 sides", spec: { type: "polygon", sides: 4 } },
  { key: "pentagon", name: "Pentagon", note: "5 sides", spec: { type: "polygon", sides: 5 } },
  { key: "hexagon", name: "Hexagon", note: "6 sides", spec: { type: "polygon", sides: 6 } },
  { key: "heptagon", name: "Heptagon", note: "7 sides", spec: { type: "polygon", sides: 7 } },
  { key: "trimer", name: "Trimer", note: "3 disks", spec: { type: "trimer", radius: 0.637556, angle: 120, distance: 1 } },
  { key: "drawn", name: "Draw", note: "your own" },
];

// The plane groups by their names in the library, the short names shown and the lattice
const GROUPS = [
  ["p1", "p1", "oblique"],
  ["p2", "p2", "oblique"],
  ["p1m1", "pm", "rectangular"],
  ["p1g1", "pg", "rectangular"],
  ["c1m1", "cm", "centred"],
  ["p2mm", "p2mm", "rectangular"],
  ["p2mg", "p2mg", "rectangular"],
  ["p2gg", "p2gg", "rectangular"],
  ["c2mm", "c2mm", "centred"],
  ["p4", "p4", "square"],
  ["p4mm", "p4mm", "square"],
  ["p4gm", "p4gm", "square"],
  ["p3", "p3", "hexagonal"],
  ["p3m1", "p3m1", "hexagonal"],
  ["p31m", "p31m", "hexagonal"],
  ["p6", "p6", "hexagonal"],
  ["p6mm", "p6mm", "hexagonal"],
];
const SHORT_NAMES = Object.fromEntries(GROUPS.map(([name, short]) => [name, short]));

// The settings of each search. The confined searches pack independent shapes, with a budget
// which scales with the number of shapes.
const WITHIN = 2;
const SURVEY = { strategy: "exhaustive", max_sites: 1, replicas: 2, steps: 6000, seed: 0, within: WITHIN / 100 };
const CONFINED = [
  { key: "strip", title: "strip", note: "periodic in 1D", counts: [2], container: "strip" },
  { key: "cluster", title: "cluster", note: "in a square", counts: [3, 4], container: "square" },
];

const DEFAULT_DRAWING = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0], [0, 0], [0, 0.5], [-0.5, 0.5]];

const statusLine = document.querySelector("#status");
const progressBar = document.querySelector("#progress div");
const runButton = document.querySelector("#run");
const expandButton = document.querySelector("#expand");
const groupGrid = document.querySelector("#group-tiles");
const bestList = document.querySelector("#best");
const drawSection = document.querySelector("#draw-section");
const drawStatus = document.querySelector("#draw-status");

let selected = "pentagon";
let running = false;
let cancelled = false;
let expanded = false;
let workers = [];
// The finished searches of the latest run, kept for the list of the best structures
let searches = { survey: null, confined: {} };

// --- Theme ---------------------------------------------------------------------------------

document.querySelector("#theme-toggle").addEventListener("click", () => {
  const dark = document.body.classList.toggle("dark");
  try {
    sessionStorage.setItem("theme", dark ? "dark" : "light");
  } catch {
    // The preference only lasts for the page without storage
  }
});

// --- Shapes --------------------------------------------------------------------------------

const pad = drawingPad({
  area: document.querySelector("#draw-area"),
  status: drawStatus,
  initial: DEFAULT_DRAWING,
  onChange: () => {
    drawStatus.classList.toggle("invalid", !pad.valid());
    renderShapes();
    updateControls();
  },
});
document.querySelector("#draw-undo").addEventListener("click", pad.undo);
document.querySelector("#draw-clear").addEventListener("click", pad.clear);

function shapeSpec(key) {
  return key === "drawn" ? pad.shape() : SHAPES.find((s) => s.key === key).spec;
}

/// An icon of the shape, in a box from -12 to 12
function icon(shape) {
  const circle = (x, y, r) => `<circle cx="${x}" cy="${y}" r="${r}" />`;
  const polygon = (vertices) => `<polygon points="${vertices.map(([x, y]) => `${x.toFixed(2)},${y.toFixed(2)}`).join(" ")}" />`;
  if (shape.key === "disk") return circle(0, 0, 10);
  if (shape.key === "trimer") return circle(0, -2.5, 5) + circle(-6, 3.5, 5) + circle(6, 3.5, 5);
  if (shape.key === "drawn") {
    const vertices = pad.vertices().length >= 3 ? centred(pad.vertices()) : [];
    const size = Math.max(1e-6, ...vertices.map(([x, y]) => Math.hypot(x, y)));
    return polygon(vertices.map(([x, y]) => [(10 * x) / size, (-10 * y) / size]));
  }
  const n = shape.spec.sides;
  return polygon(Array.from({ length: n }, (_, k) => {
    const angle = (2 * Math.PI * k) / n - Math.PI / 2;
    return [10 * Math.cos(angle), 10 * Math.sin(angle)];
  }));
}

function renderShapes() {
  document.querySelector("#shapes").replaceChildren(
    ...SHAPES.map((shape) => {
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.shape = shape.key;
      button.setAttribute("aria-pressed", String(shape.key === selected));
      button.disabled = running;
      button.innerHTML = `<svg viewBox="-12 -12 24 24" aria-hidden="true">${icon(shape)}</svg><span>${shape.name}</span><small>${shape.note}</small>`;
      button.addEventListener("click", () => {
        if (selected === shape.key) return;
        selected = shape.key;
        renderShapes();
        resetTiles();
        updateControls();
      });
      return button;
    }),
  );
}

function updateControls() {
  const drawing = selected === "drawn";
  drawSection.hidden = !drawing;
  runButton.textContent = running ? "Stop" : "Run all";
  runButton.disabled = !running && drawing && !pad.valid();
  for (const button of document.querySelectorAll("#shapes button")) button.disabled = running;
}

// --- Tiles ---------------------------------------------------------------------------------

function createTile(parent, key, title, note) {
  const tile = document.createElement("article");
  tile.className = "tile";
  tile.dataset.key = key;
  tile.innerHTML = `
    <div class="tile-head"><span><span class="group"></span> <span class="note"></span></span><span class="state"></span></div>
    <div class="view"></div>
    <div class="bar"><div></div></div>
    <div class="phi">—</div>
    <div class="detail">&nbsp;</div>`;
  tile.querySelector(".group").textContent = title;
  tile.querySelector(".note").textContent = note;
  parent.append(tile);
  return {
    el: tile,
    view: tile.querySelector(".view"),
    phi: tile.querySelector(".phi"),
    bar: tile.querySelector(".bar div"),
    detail: tile.querySelector(".detail"),
    state: tile.querySelector(".state"),
    best: null,
    bestLabel: "",
    frames: 0,
  };
}

const tiles = {};
for (const [name, short, note] of GROUPS) tiles[name] = createTile(groupGrid, name, short, note);
for (const c of CONFINED) tiles[c.key] = createTile(document.querySelector("#confined-tiles"), c.key, c.title, c.note);

function setState(tile, state, text) {
  tile.el.dataset.state = state;
  tile.state.textContent = text;
}

function resetTiles() {
  for (const tile of Object.values(tiles)) {
    tile.best = null;
    tile.bestLabel = "";
    tile.frames = 0;
    tile.el.dataset.frames = "0";
    delete tile.el.dataset.phi;
    tile.view.replaceChildren();
    tile.phi.textContent = "—";
    tile.bar.style.width = "0";
    tile.detail.innerHTML = "&nbsp;";
    setState(tile, "idle", "");
  }
  for (const search of [searches.survey, ...Object.values(searches.confined)]) search?.free();
  searches = { survey: null, confined: {} };
  renderBest();
  statusLine.textContent = "17 plane groups, a strip and a cluster in a square";
  progressBar.style.width = "0";
}

function draw(tile, shape, structure) {
  // Neighbouring cells are only drawn for periodic structures
  tile.view.innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(structure), structure.container ? 0 : 1);
}

function label(candidate) {
  if (candidate.group === "finite") {
    return `${candidate.sites.length} shapes`;
  }
  return `${SHORT_NAMES[candidate.group] ?? candidate.group}:${candidate.sites.join("")}`;
}

/// Show the first row of plane groups unless every group is expanded
function applyExpansion() {
  const columns = getComputedStyle(groupGrid).gridTemplateColumns.split(" ").filter(Boolean).length || 1;
  GROUPS.forEach(([name], index) => {
    tiles[name].el.hidden = !expanded && index >= columns;
  });
  const more = GROUPS.length - Math.min(columns, GROUPS.length);
  expandButton.setAttribute("aria-expanded", String(expanded));
  expandButton.querySelector("span").textContent = expanded ? "Show the first row only" : `Show all 17 plane groups (${more} more)`;
  expandButton.parentElement.hidden = more === 0;
}
expandButton.addEventListener("click", () => {
  expanded = !expanded;
  applyExpansion();
});
new ResizeObserver(applyExpansion).observe(groupGrid);

// --- Workers -------------------------------------------------------------------------------

function startPool(size) {
  workers = Array.from({ length: size }, () => new Worker("./worker.js", { type: "module" }));
}

function stopPool() {
  for (const worker of workers) worker.terminate();
  workers = [];
}

/// Draw the latest frame of each tile once per animation frame
function frameRenderer(shape) {
  const pending = new Map();
  let scheduled = false;
  const render = () => {
    scheduled = false;
    for (const [tile, { job, progress, structure }] of pending) {
      // A frame can arrive after its search has finished
      if (tile.el.dataset.state === "done") continue;
      draw(tile, shape, structure);
      const phi = Math.max(tile.best?.packing_fraction ?? 0, progress.best_packing_fraction ?? progress.packing_fraction);
      tile.phi.textContent = phi.toFixed(4);
      tile.bar.style.width = `${(100 * progress.steps) / progress.total_steps}%`;
      tile.detail.textContent = `${label(job.candidate)} · step ${progress.steps} of ${progress.total_steps}`;
      tile.frames += 1;
      tile.el.dataset.frames = String(tile.frames);
      setState(tile, "running", "Running");
    }
    pending.clear();
  };
  return (tile, update) => {
    pending.set(tile, update);
    if (!scheduled) {
      scheduled = true;
      requestAnimationFrame(render);
    }
  };
}

/// Run the jobs on the pool, resolving with the outcomes in the order of the jobs
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
    for (const worker of workers) {
      worker.onmessage = (event) => {
        const { id, type } = event.data;
        if (type === "progress") {
          if (event.data.structure) onFrame(jobs[id], event.data.progress, event.data.structure);
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
    }
  });
}

/// Run a search to completion, calling `tileFor` to find where each job is shown and keeping
/// the best structure of each tile after every round
async function runSearch(shape, config, title, tileFor) {
  const searcher = new WasmSearch(JSON.stringify(shape), JSON.stringify(config));
  const frame = frameRenderer(shape);
  while (!searcher.isFinished() && !cancelled) {
    const jobs = JSON.parse(searcher.jobs());
    statusLine.textContent = `${title}: round ${searcher.round() + 1}, ${searcher.numActive()} of ${searcher.numCandidates()} candidates`;
    progressBar.style.width = "0";
    const outcomes = await runJobs(
      jobs,
      (done, total) => {
        progressBar.style.width = `${(100 * done) / total}%`;
      },
      (job, progress, structure) => frame(tileFor(job.candidate), { job, progress, structure }),
    );
    if (cancelled) break;
    searcher.submit(JSON.stringify(outcomes));
    for (const result of JSON.parse(searcher.ranking())) {
      const tile = result.best && tileFor(result.candidate);
      if (tile && (!tile.best || result.best.packing_fraction > tile.best.packing_fraction)) {
        tile.best = result.best;
        tile.bestLabel = label(result.candidate);
      }
    }
  }
  return searcher;
}

/// Show the best structure found by each tile once its search has finished
function finishTile(tile, shape) {
  if (tile.best) {
    draw(tile, shape, tile.best);
    tile.phi.textContent = tile.best.packing_fraction.toFixed(4);
    tile.el.dataset.phi = String(tile.best.packing_fraction);
    tile.bar.style.width = "100%";
    tile.detail.textContent = `best ${tile.bestLabel}`;
    setState(tile, "done", "Done");
  } else {
    setState(tile, "none", cancelled ? "Stopped" : "No candidates");
  }
}

// --- Best structures -----------------------------------------------------------------------

function renderBest() {
  const items = [];
  const shape = searches.shape;
  if (searches.survey) {
    for (const { result, equivalent } of JSON.parse(searches.survey.distinctNearBest(WITHIN))) {
      items.push({ label: label(result.candidate), structure: result.best, also: equivalent.map((r) => label(r.candidate)), crystal: true });
    }
  }
  for (const c of CONFINED) {
    const best = searches.confined[c.key] && JSON.parse(searches.confined[c.key].best());
    if (best) items.push({ label: `${c.title} · ${best.num_shapes} shapes`, structure: best, also: [], crystal: false });
  }
  items.sort((a, b) => b.structure.packing_fraction - a.structure.packing_fraction);
  const top = items[0]?.structure.packing_fraction;
  document.querySelector("#best-empty").hidden = items.length > 0;
  bestList.replaceChildren(
    ...items.map((item, index) => {
      const phi = item.structure.packing_fraction;
      const below = (100 * (top - phi)) / top;
      const li = document.createElement("li");
      li.dataset.phi = String(phi);
      li.dataset.label = item.label;
      li.classList.toggle("near", below <= WITHIN + 1e-9);
      li.innerHTML = `
        <span class="rank">${index + 1}</span>
        <div class="view"></div>
        <div><div class="label"></div><div class="also"></div></div>
        <div class="phi">${phi.toFixed(5)}<span class="below">${index === 0 ? "best" : `${below.toFixed(2)}% below`}</span></div>`;
      li.querySelector(".label").textContent = item.label;
      li.querySelector(".also").textContent = item.also.length
        ? `Also found in ${item.also.join(", ")}`
        : item.crystal ? "Crystal" : "Independent shapes";
      li.querySelector(".view").innerHTML = renderSvg(JSON.stringify(shape), JSON.stringify(item.structure), item.structure.container ? 0 : 1);
      return li;
    }),
  );
}

// --- Running -------------------------------------------------------------------------------

async function runAll() {
  const shape = shapeSpec(selected);
  resetTiles();
  searches.shape = shape;
  cancelled = false;
  running = true;
  delete document.body.dataset.result;
  document.body.dataset.running = "true";
  updateControls();
  startPool(navigator.hardwareConcurrency || 4);
  const start = performance.now();
  for (const tile of Object.values(tiles)) setState(tile, "queued", "Queued");

  try {
    const survey = await runSearch(shape, SURVEY, "Plane groups", (candidate) => tiles[candidate.group]);
    searches.survey = survey;
    for (const [name] of GROUPS) finishTile(tiles[name], shape);
    if (!cancelled) renderBest();

    for (const c of CONFINED) {
      if (cancelled) break;
      const config = { ...SURVEY, unconstrained: c.counts, steps_per_shape: 2500, container: c.container };
      searches.confined[c.key] = await runSearch(shape, config, `The ${c.title}`, () => tiles[c.key]);
      finishTile(tiles[c.key], shape);
      if (!cancelled) renderBest();
    }
  } catch (error) {
    statusLine.textContent = `Error: ${error.message ?? error}`;
    cancelled = true;
  }

  stopPool();
  running = false;
  delete document.body.dataset.running;
  updateControls();
  for (const tile of Object.values(tiles)) {
    if (tile.el.dataset.state === "queued" || tile.el.dataset.state === "running") finishTile(tile, shape);
  }
  progressBar.style.width = cancelled ? "0" : "100%";
  if (!cancelled) {
    const best = JSON.parse(searches.survey.best());
    const seconds = ((performance.now() - start) / 1000).toFixed(1);
    statusLine.textContent = `Finished in ${seconds} s: the densest crystal has a packing fraction of ${best.packing_fraction.toFixed(6)} in ${SHORT_NAMES[best.group] ?? best.group}`;
    document.body.dataset.result = JSON.stringify(best);
  } else if (!statusLine.textContent.startsWith("Error")) {
    statusLine.textContent = "Stopped";
  }
}

runButton.addEventListener("click", () => {
  if (running) {
    cancelled = true;
    stopPool();
  } else {
    runAll();
  }
});

await init();
document.querySelector("#version").textContent = version();
pad.update();
resetTiles();
applyExpansion();
document.body.dataset.ready = "true";
