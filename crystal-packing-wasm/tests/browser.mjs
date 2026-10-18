// Run a small search in headless Chromium using the web front end and its worker pool.
//
// Usage: node tests/browser.mjs   (after ./build.sh)
// Requires playwright; set CHROMIUM_PATH to use a specific browser executable.
import { createServer } from "node:http";
import { readFile } from "node:fs/promises";
import { createRequire } from "node:module";
import { extname, join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const require = createRequire(import.meta.url);
const { chromium } = require("playwright");

const root = join(dirname(fileURLToPath(import.meta.url)), "..", "www");
const types = { ".html": "text/html", ".js": "text/javascript", ".wasm": "application/wasm" };
const server = createServer(async (request, response) => {
  const path = join(root, request.url === "/" ? "index.html" : request.url);
  try {
    const body = await readFile(path);
    response.writeHead(200, { "Content-Type": types[extname(path)] ?? "application/octet-stream" });
    response.end(body);
  } catch {
    response.writeHead(404);
    response.end();
  }
});
await new Promise((resolve) => server.listen(0, resolve));
const url = `http://localhost:${server.address().port}/`;

const browser = await chromium.launch({ executablePath: process.env.CHROMIUM_PATH || undefined });
let failed = false;
function check(condition, message) {
  if (!condition) {
    console.error(`FAIL: ${message}`);
    failed = true;
  }
}
try {
  const page = await browser.newPage();
  page.on("pageerror", (error) => {
    console.error(`Page error: ${error}`);
    failed = true;
  });
  page.setDefaultTimeout(120_000);
  await page.goto(url);
  await page.waitForSelector("body[data-ready=true]");

  // Results within a configurable percentage of the best, best first
  async function results() {
    return page.$$eval("#results tbody tr", (rows) => rows.map((row) => Number(row.dataset.phi)));
  }
  async function resultCount() {
    return Number(await page.getAttribute("#results-title", "data-count"));
  }

  // A search, with every optimiser shown as it runs
  await page.selectOption("select[name=shape]", "polygon");
  await page.fill("input[name=sides]", "5");
  await page.fill("input[name=groups]", "");
  await page.fill("input[name=steps]", "4000");
  await page.fill("input[name=replicas]", "2");
  await page.fill("input[name=within]", "2");
  await page.click("#search");
  await page.waitForFunction(
    () => {
      const tiles = [...document.querySelectorAll("#live .tile")];
      return tiles.length > 0 && tiles.every((tile) => Number(tile.dataset.frames) > 0);
    },
    null,
    { timeout: 60_000 },
  );
  const tiles = await page.$$eval("#live .tile", (tiles) => tiles.length);
  const workers = await page.evaluate(() => navigator.hardwareConcurrency || 4);
  console.log(`${tiles} live optimisers for ${workers} workers`);
  check(tiles === workers, "Expected a live view of every worker");
  await page.waitForSelector("body[data-result]", { timeout: 300_000 });
  const best = JSON.parse(await page.getAttribute("body", "data-result"));
  console.log(await page.textContent("#status"));
  check(Math.abs(best.packing_fraction - 0.921311) < 1e-3, "Expected the pentagon double lattice");

  let phis = await results();
  console.log(`${phis.length} results within 2%: ${phis.map((p) => p.toFixed(5)).join(", ")}`);
  check(phis.length === (await resultCount()), "Every result within 2% is listed");
  check(phis.length >= 1, "Expected at least one result");
  check(phis.every((p, i) => i === 0 || p <= phis[i - 1]), "Results are ordered best first");
  check(phis.every((p) => p >= 0.98 * best.packing_fraction - 1e-9), "Results are within 2%");
  check(
    (await page.$$("#near-best .card")).length === phis.length,
    "Every result within 2% is drawn",
  );

  // Changing the tolerance and limit updates the results without searching again
  await page.fill("input[name=within]", "100");
  await page.dispatchEvent("input[name=within]", "change");
  phis = await results();
  const candidates = Number(await page.getAttribute("#results-title", "data-candidates"));
  check(phis.length === candidates, `Expected all ${candidates} candidates within 100%, found ${phis.length}`);
  check(phis.every((p, i) => i === 0 || p <= phis[i - 1]), "All results are ordered best first");
  await page.fill("input[name=max_results]", "3");
  await page.dispatchEvent("input[name=max_results]", "change");
  check((await results()).length === 3, "The number of results can be limited");
  check((await page.textContent("#results-title")).includes("showing the best 3"), "The limit is reported");
  await page.fill("input[name=max_results]", "");
  await page.fill("input[name=within]", "2");
  await page.dispatchEvent("input[name=within]", "change");

  // Watch a single optimisation, which should animate the structure and plot its progress
  await page.fill("input[name=watch_group]", "p2");
  await page.fill("input[name=watch_sites]", "e");
  await page.fill("input[name=steps]", "20000");
  await page.click("#watch");
  await page.waitForSelector("body[data-watch-result]", { timeout: 120_000 });
  const watched = JSON.parse(await page.getAttribute("body", "data-watch-result"));
  const frames = Number(await page.getAttribute("#watch-view", "data-frames"));
  const points = await page.$$eval("#watch-chart polyline", (lines) =>
    lines.length ? lines[0].getAttribute("points").trim().split(/\s+/).length : 0,
  );
  console.log(`Watched ${frames} frames, ${points} points, packing fraction ${watched.packing_fraction}`);
  check(frames >= 5, "Expected the structure to be drawn as it is optimised");
  check(points >= 5, "Expected the packing fraction to be plotted");
  check(watched.packing_fraction > 0.85, "Expected the pentagons to compress");

  // Draw a shape by clicking, which is checked as it is drawn
  async function draw(points) {
    await page.click("#draw-clear");
    // Convert from the drawing coordinates (y upwards) to the page using the SVG's transform
    const screen = await page.evaluate((points) => {
      const area = document.querySelector("#draw-area");
      const matrix = area.getScreenCTM();
      return points.map(([x, y]) => {
        const p = new DOMPoint(x, -y).matrixTransform(matrix);
        return [p.x, p.y];
      });
    }, points);
    for (const [x, y] of screen) {
      await page.mouse.click(x, y);
    }
  }
  await page.selectOption("select[name=shape]", "drawn");
  await draw([[0, 0], [1.5, 1.5], [1.5, 0], [0, 1.2]]);
  const invalid = await page.textContent("#draw-status");
  console.log(`Bow tie: ${invalid}`);
  check(invalid.includes("not simple"), "A self-intersecting shape is rejected");
  check(await page.isDisabled("#search"), "Searching is disabled for an invalid shape");

  await draw([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]);
  const valid = await page.textContent("#draw-status");
  console.log(`Square: ${valid}`);
  check(valid.includes("4 vertices") && valid.includes("area 1.000") && valid.includes("4-fold"), "The drawn square is valid");
  await page.fill("input[name=groups]", "p1,p2");
  await page.fill("input[name=steps]", "3000");
  await page.click("#search");
  await page.waitForSelector("body[data-result]", { timeout: 120_000 });
  const drawn = JSON.parse(await page.getAttribute("body", "data-result"));
  console.log(`Drawn square: ${drawn.packing_fraction}`);
  check(drawn.packing_fraction > 0.99, "Expected the drawn square to fill space");

  // Independent shapes without symmetry in a periodic cell
  await page.fill("input[name=unconstrained_shapes]", "1,2");
  await page.click("#unconstrained");
  await page.waitForSelector("body[data-result]", { timeout: 120_000 });
  const free = JSON.parse(await page.getAttribute("body", "data-result"));
  const counts = await page.$$eval("#results tbody tr", (rows) => rows.map((row) => row.dataset.label));
  console.log(`Unconstrained: ${free.packing_fraction} ${counts.join(", ")}`);
  check(free.group === "p1" && free.packing_fraction > 0.99, "Expected free squares to fill space");
} finally {
  await browser.close();
  server.close();
}
process.exit(failed ? 1 : 0);
