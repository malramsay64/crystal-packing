//
// lib.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! WebAssembly bindings for the crystal packing library.
//!
//! Values cross the boundary as JSON strings, using the same serialisation as the command line
//! output. A search runs on the main thread, handing out [`Job`]s which are run by a pool of web
//! workers using [`WasmJob`], so the page stays responsive and all cores are used.

use wasm_bindgen::prelude::*;

use crystal_packing::engine::Structure;
use crystal_packing::render::render_structure;
use crystal_packing::search::{run_job, Job, JobRunner, Outcome, Search, SearchConfig, ShapeSpec};

fn to_js<E: std::fmt::Display>(error: E) -> JsError {
    JsError::new(&error.to_string())
}

#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
}

/// The version of the packing library compiled into this module.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// A search over the plane groups for the densest packing of a shape
#[wasm_bindgen]
pub struct WasmSearch {
    search: Search,
}

#[wasm_bindgen]
impl WasmSearch {
    /// Create a search from JSON descriptions of the shape and the search configuration.
    ///
    /// Fields missing from the configuration take their default values.
    #[wasm_bindgen(constructor)]
    pub fn new(shape: &str, config: &str) -> Result<WasmSearch, JsError> {
        let shape: ShapeSpec = serde_json::from_str(shape).map_err(to_js)?;
        let mut merged = serde_json::to_value(SearchConfig::default()).map_err(to_js)?;
        let overrides: serde_json::Value = serde_json::from_str(config).map_err(to_js)?;
        if let (Some(base), Some(overrides)) = (merged.as_object_mut(), overrides.as_object()) {
            for (key, value) in overrides {
                base.insert(key.clone(), value.clone());
            }
        }
        let config: SearchConfig = serde_json::from_value(merged).map_err(to_js)?;
        Ok(WasmSearch {
            search: Search::new(shape, config).map_err(to_js)?,
        })
    }

    #[wasm_bindgen(js_name = isFinished)]
    pub fn is_finished(&self) -> bool {
        self.search.is_finished()
    }

    pub fn round(&self) -> usize {
        self.search.round()
    }

    #[wasm_bindgen(js_name = numCandidates)]
    pub fn num_candidates(&self) -> usize {
        self.search.results.len()
    }

    #[wasm_bindgen(js_name = numActive)]
    pub fn num_active(&self) -> usize {
        self.search.num_active()
    }

    /// The jobs of the current round as a JSON array
    pub fn jobs(&self) -> Result<String, JsError> {
        serde_json::to_string(&self.search.jobs()).map_err(to_js)
    }

    /// Submit a JSON array of the outcomes of the current round's jobs
    pub fn submit(&mut self, outcomes: &str) -> Result<(), JsError> {
        let outcomes: Vec<Outcome> = serde_json::from_str(outcomes).map_err(to_js)?;
        self.search.submit(outcomes);
        Ok(())
    }

    /// The candidates ranked by their best packing fraction, as JSON
    pub fn ranking(&self) -> Result<String, JsError> {
        serde_json::to_string(&self.search.ranking()).map_err(to_js)
    }

    /// The best structure found so far as JSON, or null
    pub fn best(&self) -> Result<String, JsError> {
        serde_json::to_string(&self.search.best()).map_err(to_js)
    }

    /// The candidates within `within_percent` of the best score, best first, as JSON
    #[wasm_bindgen(js_name = nearBest)]
    pub fn near_best(&self, within_percent: f64) -> Result<String, JsError> {
        serde_json::to_string(&self.search.near_best(within_percent / 100.)).map_err(to_js)
    }

    /// The distinct structures within `within_percent` of the best score, best first, as JSON
    /// objects with the `result` and the `equivalent` results which found the same crystal
    #[wasm_bindgen(js_name = distinctNearBest)]
    pub fn distinct_near_best(&self, within_percent: f64) -> Result<String, JsError> {
        let distinct = self
            .search
            .distinct_near_best(within_percent / 100., None)
            .map_err(to_js)?;
        serde_json::to_string(&distinct).map_err(to_js)
    }

    /// The best candidate of every group searched, as JSON pairs of [group, result or null]
    #[wasm_bindgen(js_name = bestPerGroup)]
    pub fn best_per_group(&self) -> Result<String, JsError> {
        serde_json::to_string(&self.search.best_per_group()).map_err(to_js)
    }

    /// The current phase of the search
    pub fn phase(&self) -> String {
        format!("{:?}", self.search.phase())
    }

    /// Run the remaining search on the current thread, returning the best structure as JSON
    #[wasm_bindgen(js_name = runToCompletion)]
    pub fn run_to_completion(&mut self) -> Result<String, JsError> {
        serde_json::to_string(&self.search.run()).map_err(to_js)
    }
}

/// A single job, run in slices so progress can be reported from a web worker
#[wasm_bindgen]
pub struct WasmJob {
    runner: Option<JobRunner>,
}

#[wasm_bindgen]
impl WasmJob {
    #[wasm_bindgen(constructor)]
    pub fn new(job: &str) -> Result<WasmJob, JsError> {
        let job: Job = serde_json::from_str(job).map_err(to_js)?;
        Ok(WasmJob {
            runner: Some(JobRunner::new(&job)),
        })
    }

    /// Run up to `steps` steps, returning the progress as JSON
    pub fn run(&mut self, steps: u32) -> Result<String, JsError> {
        let runner = self
            .runner
            .as_mut()
            .ok_or_else(|| JsError::new("Job already finished"))?;
        serde_json::to_string(&runner.run(steps as u64)).map_err(to_js)
    }

    /// The structure as it currently is, as JSON
    pub fn current(&self) -> Result<String, JsError> {
        let runner = self
            .runner
            .as_ref()
            .ok_or_else(|| JsError::new("Job already finished"))?;
        serde_json::to_string(&runner.current()).map_err(to_js)
    }

    #[wasm_bindgen(js_name = isComplete)]
    pub fn is_complete(&self) -> bool {
        self.runner.as_ref().is_none_or(|r| r.is_complete())
    }

    /// Finish the job, returning the outcome as JSON
    pub fn finish(&mut self) -> Result<String, JsError> {
        let runner = self
            .runner
            .take()
            .ok_or_else(|| JsError::new("Job already finished"))?;
        serde_json::to_string(&runner.finish()).map_err(to_js)
    }
}

/// Run a job to completion, taking and returning JSON
#[wasm_bindgen(js_name = runJob)]
pub fn run_job_json(job: &str) -> Result<String, JsError> {
    let job: Job = serde_json::from_str(job).map_err(to_js)?;
    serde_json::to_string(&run_job(&job)).map_err(to_js)
}

/// Check a shape is valid, returning a JSON summary of its properties, or an error explaining
/// why it is invalid.
#[wasm_bindgen(js_name = shapeInfo)]
pub fn shape_info(shape: &str) -> Result<String, JsError> {
    let spec: ShapeSpec = serde_json::from_str(shape).map_err(to_js)?;
    let shape = spec.build().map_err(to_js)?;
    let symmetry = shape.symmetry();
    serde_json::to_string(&serde_json::json!({
        "name": shape.name(),
        "area": shape.area(),
        "circumradius": shape.circumradius(),
        "convex": shape.is_convex(),
        "rotations": symmetry.rotations,
        "mirrors": symmetry.mirrors.len(),
    }))
    .map_err(to_js)
}

/// Render a structure of a shape as an SVG image
#[wasm_bindgen(js_name = renderSvg)]
pub fn render_svg(shape: &str, structure: &str, cells: i32) -> Result<String, JsError> {
    let shape: ShapeSpec = serde_json::from_str(shape).map_err(to_js)?;
    let structure: Structure = serde_json::from_str(structure).map_err(to_js)?;
    Ok(render_structure(&shape, &structure, cells)
        .map_err(to_js)?
        .to_string())
}
