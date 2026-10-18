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
  await page.goto(url);
  await page.waitForSelector("body[data-ready=true]");

  // A search, which should list every structure close to the best
  await page.fill("input[name=sides]", "4");
  await page.fill("input[name=groups]", "p1,p2,p4,p4mm");
  await page.fill("input[name=steps]", "3000");
  await page.fill("input[name=replicas]", "2");
  await page.fill("input[name=within]", "2");
  await page.click("#search");
  await page.waitForSelector("body[data-result]", { timeout: 120_000 });
  const result = JSON.parse(await page.getAttribute("body", "data-result"));
  console.log(await page.textContent("#status"));
  check(result.packing_fraction > 0.99, "Expected squares to fill space");
  check(await page.$("#structure svg"), "Expected the best structure to be rendered");
  const cards = await page.$$eval("#near-best .card", (cards) =>
    cards.map((card) => ({ phi: Number(card.dataset.phi), svg: !!card.querySelector("svg") })),
  );
  console.log(`${cards.length} structures within 2% of the best`);
  check(cards.length >= 1, "Expected at least one structure near the best");
  check(cards.every((card) => card.svg), "Every structure near the best is drawn");
  check(
    cards.every((card) => card.phi >= 0.98 * result.packing_fraction - 1e-9),
    "Structures near the best are within 2%",
  );

  // Watch a single optimisation, which should animate the structure and plot its progress
  await page.selectOption("select[name=shape]", "polygon");
  await page.fill("input[name=sides]", "5");
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
  check(await page.$("#watch-view svg"), "Expected the final structure to be drawn");
} finally {
  await browser.close();
  server.close();
}
process.exit(failed ? 1 : 0);
