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
try {
  const page = await browser.newPage();
  page.on("pageerror", (error) => {
    console.error(`Page error: ${error}`);
    failed = true;
  });
  await page.goto(url);
  await page.waitForSelector("body[data-ready=true]");
  await page.fill("input[name=sides]", "4");
  await page.fill("input[name=groups]", "p1,p2,p4");
  await page.fill("input[name=steps]", "3000");
  await page.fill("input[name=replicas]", "2");
  await page.click("button[type=submit]");
  await page.waitForSelector("body[data-result]", { timeout: 120_000 });
  const result = JSON.parse(await page.getAttribute("body", "data-result"));
  const status = await page.textContent("#status");
  const svg = await page.$("#structure svg");
  console.log(status);
  if (!(result.packing_fraction > 0.99) || !svg) {
    console.error("Expected squares to fill space and the structure to be rendered");
    failed = true;
  }
} finally {
  await browser.close();
  server.close();
}
process.exit(failed ? 1 : 0);
