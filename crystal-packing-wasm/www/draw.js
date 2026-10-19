// Draw a polygon by hand: click to add a vertex, drag a vertex to move it, and double click a
// vertex to remove it. The shape is checked by the library as it is drawn.
import { shapeInfo } from "./pkg/crystal_packing_wasm.js";

/// The vertices moved so the area centroid is at the origin
export function centred(vertices) {
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

/// Attach a drawing pad to the SVG `area`, reporting whether the shape is valid in `status`.
/// `snap` returns whether to snap vertices to a grid of 0.1 and `onChange` is called after every
/// change to the drawing.
export function drawingPad({ area, status, snap = () => true, onChange = () => {}, initial = [] }) {
  let vertices = initial.map(([x, y]) => [x, y]);
  let valid = false;

  const shape = () => ({ type: "vertices", name: "Drawn", vertices: centred(vertices) });

  function validate() {
    if (vertices.length < 3) {
      valid = false;
      status.textContent = "Add at least three vertices";
      return;
    }
    try {
      const info = JSON.parse(shapeInfo(JSON.stringify(shape())));
      valid = true;
      status.textContent = `Valid polygon, ${vertices.length} vertices, area ${info.area.toFixed(3)}, ${info.convex ? "convex" : "concave"}, ${info.rotations}-fold rotation, ${info.mirrors} mirror${info.mirrors === 1 ? "" : "s"}`;
    } catch (error) {
      valid = false;
      status.textContent = `Invalid shape: ${error.message ?? error}`;
    }
  }

  function render() {
    const grid = [];
    for (let i = -20; i <= 20; i += 1) {
      const v = i / 10;
      const kind = i === 0 ? "axis" : i % 5 === 0 ? "grid major" : "grid";
      grid.push(`<line class="${kind}" x1="${v}" y1="-2" x2="${v}" y2="2" />`);
      grid.push(`<line class="${kind}" x1="-2" y1="${v}" x2="2" y2="${v}" />`);
    }
    const points = vertices.map(([x, y]) => `${x},${-y}`).join(" ");
    const handles = vertices
      .map(([x, y], i) => `<circle data-index="${i}" cx="${x}" cy="${-y}" r="0.05" />`)
      .join("");
    const polygon = vertices.length >= 2 ? `<polygon class="${valid ? "" : "invalid"}" points="${points}" />` : "";
    area.innerHTML = grid.join("") + polygon + handles;
  }

  function update() {
    validate();
    render();
    onChange();
  }

  /// The position of a pointer event in the drawing coordinates, with y upwards
  function point(event) {
    const matrix = area.getScreenCTM().inverse();
    const p = new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix);
    let [x, y] = [p.x, -p.y];
    if (snap()) [x, y] = [Math.round(x * 10) / 10, Math.round(y * 10) / 10];
    return [x, y];
  }

  let dragging = null;
  let moved = false;
  area.addEventListener("pointerdown", (event) => {
    const index = event.target.dataset?.index;
    if (index !== undefined) {
      dragging = Number(index);
      moved = false;
      area.setPointerCapture(event.pointerId);
    }
  });
  area.addEventListener("pointermove", (event) => {
    if (dragging === null) return;
    moved = true;
    vertices[dragging] = point(event);
    update();
  });
  area.addEventListener("pointerup", (event) => {
    if (dragging !== null) {
      dragging = null;
      return;
    }
    if (event.target.dataset?.index === undefined) {
      vertices.push(point(event));
      update();
    }
  });
  area.addEventListener("dblclick", (event) => {
    const index = event.target.dataset?.index;
    if (index !== undefined && !moved) {
      vertices.splice(Number(index), 1);
      update();
    }
  });

  return {
    shape,
    vertices: () => vertices,
    valid: () => valid,
    update,
    undo() {
      vertices.pop();
      update();
    },
    clear() {
      vertices = [];
      update();
    },
  };
}
