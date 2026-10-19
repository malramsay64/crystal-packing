// Run every structure of a shape in headless Chromium using the dashboard (index.html).
//
// Usage: node tests/dashboard.mjs   (after ./build.sh)
// Requires playwright; set CHROMIUM_PATH to use a specific browser executable and SCREENSHOTS to
// a directory to save screenshots of the finished page in light and dark mode.
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
  const page = await browser.newPage({ viewport: { width: 1000, height: 900 } });
  page.on("pageerror", (error) => {
    console.error(`Page error: ${error}`);
    failed = true;
  });
  page.setDefaultTimeout(120_000);
  await page.goto(url);
  await page.waitForSelector("body[data-ready=true]");

  // Only the first row of plane groups is shown until expanded
  const visible = () => page.$$eval("#group-tiles .tile", (tiles) => tiles.filter((t) => !t.hidden).length);
  const firstRow = await visible();
  check(firstRow > 0 && firstRow < 17, `Expected only the first row of groups, found ${firstRow}`);
  await page.click("#expand");
  check((await visible()) === 17, "Expected every plane group when expanded");
  await page.click("#expand");
  check((await visible()) === firstRow, "Expected the first row again when collapsed");

  // Pentagons in every structure, shown as they run
  await page.click("#shapes button[data-shape=pentagon]");
  await page.click("#run");
  await page.waitForFunction(() => [...document.querySelectorAll("#group-tiles .tile")].some((t) => Number(t.dataset.frames) > 0));
  await page.waitForSelector("body[data-result]", { timeout: 600_000 });
  console.log(await page.textContent("#status"));
  const best = JSON.parse(await page.getAttribute("body", "data-result"));
  check(Math.abs(best.packing_fraction - 0.921311) < 1e-3, "Expected the pentagon double lattice");
  const groups = await page.$$eval("#group-tiles .tile", (tiles) => tiles.map((t) => t.dataset.state));
  check(groups.every((state) => state === "done"), `Expected every group to finish: ${groups.join(" ")}`);
  const confined = await page.$$eval("#confined-tiles .tile", (tiles) => tiles.map((t) => [t.dataset.key, Number(t.dataset.phi)]));
  console.log(`Confined: ${confined.map(([k, p]) => `${k} ${p.toFixed(4)}`).join(", ")}`);
  check(confined.length === 2 && confined.every(([, p]) => p > 0.5), "Expected the strip and cluster to pack");
  const rows = await page.$$eval("#best li", (rows) => rows.map((r) => [r.dataset.label, Number(r.dataset.phi)]));
  console.log(`Best: ${rows.map(([l, p]) => `${l} ${p.toFixed(5)}`).join(", ")}`);
  check(rows.length >= 3, "Expected the crystals, strip and cluster to be listed");
  check(rows.every(([, p], i) => i === 0 || p <= rows[i - 1][1]), "Expected the best structures in order");
  check(Math.abs(rows[0][1] - best.packing_fraction) < 1e-9, "Expected the densest crystal first");

  if (process.env.SCREENSHOTS) {
    await page.click("#expand");
    await page.screenshot({ path: join(process.env.SCREENSHOTS, "light.png"), fullPage: true });
    await page.click("#theme-toggle");
    await page.screenshot({ path: join(process.env.SCREENSHOTS, "dark.png"), fullPage: true });
    await page.click("#theme-toggle");
  }

  // A drawn shape is checked as it is drawn and can't be run while invalid
  await page.click("#shapes button[data-shape=drawn]");
  check(await page.isVisible("#draw-area"), "Expected the drawing area");
  check(!(await page.isDisabled("#run")), "Expected the default drawing to be valid");
  await page.click("#draw-clear");
  check(await page.isDisabled("#run"), "Expected running to be disabled without a shape");
} finally {
  await browser.close();
  server.close();
}
process.exit(failed ? 1 : 0);
