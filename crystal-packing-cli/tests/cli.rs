//
// cli.rs
// Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
// Distributed under terms of the MIT license.
//

//! Tests running the `packing` binary.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn packing(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_packing"))
        .args(args)
        .output()
        .expect("Failed to run packing")
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).to_string()
}

/// A fresh, empty directory for the output of a test
fn output_dir(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&dir);
    dir
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "Command failed\nstdout:\n{}\nstderr:\n{}",
        stdout(output),
        String::from_utf8_lossy(&output.stderr)
    );
}

const FAST: &[&str] = &["--steps", "400", "--replicas", "1", "--seed", "1"];

#[test]
fn survey_reports_every_group() {
    let dir = output_dir("survey");
    let mut args = vec!["survey", "--shape", "polygon", "--sides", "4", "-o"];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let text = stdout(&output);
    for group in [
        "p1", "p2", "p1m1", "p1g1", "c1m1", "p2mm", "p2mg", "p2gg", "c2mm", "p4", "p4mm", "p4gm",
        "p3", "p3m1", "p31m", "p6", "p6mm",
    ] {
        assert!(
            text.lines()
                .any(|l| l.split_whitespace().next() == Some(group)),
            "{} missing from\n{}",
            group,
            text
        );
    }
    let csv = fs::read_to_string(dir.join("survey.csv")).unwrap();
    assert_eq!(csv.lines().count(), 18, "{}", csv);
    assert!(csv.starts_with("group,"));
    // A square can't occupy any site of the hexagonal groups with 3 fold symmetry, but can
    // occupy the general position of every group, so every group has a structure
    assert!(dir.join("groups").join("p4.svg").exists());
    assert!(dir.join("groups").join("p6mm.json").exists());
}

#[test]
fn survey_of_selected_groups() {
    let dir = output_dir("survey-groups");
    let mut args = vec!["survey", "--shape", "disk", "--groups", "p1,p6mm", "-o"];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let csv = fs::read_to_string(dir.join("survey.csv")).unwrap();
    assert_eq!(csv.lines().count(), 3, "{}", csv);
}

#[test]
fn search_lists_results_near_the_best() {
    let dir = output_dir("search-near");
    let mut args = vec![
        "search",
        "--shape",
        "polygon",
        "--sides",
        "4",
        "--groups",
        "p1,p2,p4,p4mm",
        "--within",
        "2",
        "-o",
    ];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let text = stdout(&output);
    assert!(text.contains("within 2% of the best"), "{}", text);

    let near: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("near-best.json")).unwrap()).unwrap();
    let near = near.as_array().unwrap();
    assert!(!near.is_empty());
    let best = near[0]["structure"]["packing_fraction"].as_f64().unwrap();
    for entry in near {
        let phi = entry["structure"]["packing_fraction"].as_f64().unwrap();
        assert!(phi >= 0.98 * best - 1e-12, "{} vs {}", phi, best);
    }
    let images = fs::read_dir(dir.join("near-best")).unwrap().count();
    // A JSON and an SVG file for each structure
    assert_eq!(images, 2 * near.len());
}

/// The packing fractions listed in the results table, which follows the line starting with
/// `title` and ends at the first blank line
fn table_fractions(text: &str, title: &str) -> Vec<f64> {
    text.lines()
        .skip_while(|l| !l.contains(title))
        .skip(2)
        .take_while(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().nth(2).unwrap().parse().unwrap())
        .collect()
}

/// The number of results, and distinct structures, from the line introducing the table
fn counts(text: &str, title: &str) -> (usize, usize) {
    let line = text.lines().find(|l| l.contains(title)).unwrap();
    let results: usize = line.split_whitespace().next().unwrap().parse().unwrap();
    let distinct: usize = line
        .split(", ")
        .find(|part| part.contains("distinct"))
        .and_then(|part| part.split_whitespace().next())
        .unwrap()
        .parse()
        .unwrap();
    (results, distinct)
}

#[test]
fn search_shows_distinct_structures_near_the_best_in_order() {
    let dir = output_dir("search-order");
    let mut args = vec!["search", "--shape", "polygon", "--sides", "5", "-o"];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let text = stdout(&output);
    // By default the distinct structures within 2% of the best are listed
    let (results, distinct) = counts(&text, "within 2% of the best");
    assert!(distinct <= results);
    let fractions = table_fractions(&text, "within 2% of the best");
    assert_eq!(fractions.len(), distinct, "{}", text);
    assert!(
        fractions.windows(2).all(|w| w[0] >= w[1]),
        "Not ordered:\n{}",
        text
    );
    assert!(fractions.iter().all(|&f| f >= 0.98 * fractions[0] - 1e-9));
    // The CSV has every candidate, not only those near the best
    let csv = fs::read_to_string(dir.join("results.csv")).unwrap();
    assert!(csv.lines().count() - 1 > results);
    // Every equivalent candidate is listed with its structure
    let near: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("near-best.json")).unwrap()).unwrap();
    let listed: usize = near
        .as_array()
        .unwrap()
        .iter()
        .map(|d| 1 + d["equivalent"].as_array().unwrap().len())
        .sum();
    assert_eq!(listed, results);
}

#[test]
fn equivalent_structures_can_be_listed_separately() {
    let dir = output_dir("search-equivalent");
    let mut args = vec![
        "search",
        "--shape",
        "disk",
        "--show-equivalent",
        "--within",
        "5",
        "-o",
    ];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let text = stdout(&output);
    let (results, _) = counts(&text, "within 5% of the best");
    assert_eq!(
        table_fractions(&text, "within 5% of the best").len(),
        results,
        "{}",
        text
    );
}

#[test]
fn within_is_configurable() {
    let dir = output_dir("search-within");
    let mut args = vec![
        "search", "--shape", "polygon", "--sides", "5", "--within", "100", "-o",
    ];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let text = stdout(&output);
    let csv = fs::read_to_string(dir.join("results.csv")).unwrap();
    // Within 100% of the best includes every candidate
    let (results, _) = counts(&text, "within 100% of the best");
    assert_eq!(results, csv.lines().count() - 1, "{}", text);
    let fractions = table_fractions(&text, "within 100% of the best");
    assert!(fractions.windows(2).all(|w| w[0] >= w[1]));
}

#[test]
fn number_of_results_can_be_limited() {
    let dir = output_dir("search-limit");
    let mut args = vec![
        "search",
        "--shape",
        "polygon",
        "--sides",
        "5",
        "--within",
        "100",
        "--max-results",
        "3",
        "-o",
    ];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let text = stdout(&output);
    assert!(text.contains("showing the best 3"), "{}", text);
    assert_eq!(table_fractions(&text, "within 100% of the best").len(), 3);
    assert_eq!(fs::read_dir(dir.join("near-best")).unwrap().count(), 6);
}

#[test]
fn unconstrained_shapes() {
    let dir = output_dir("unconstrained");
    let mut args = vec![
        "unconstrained",
        "--shape",
        "polygon",
        "--sides",
        "4",
        "--shapes",
        "1,2,3",
        "-o",
    ];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    let text = stdout(&output);
    let csv = fs::read_to_string(dir.join("unconstrained.csv")).unwrap();
    assert!(csv.starts_with("shapes,"));
    assert_eq!(csv.lines().count(), 4, "{}", csv);
    for n in 1..=3 {
        assert!(dir.join(format!("shapes-{}.svg", n)).exists());
    }
    assert!(text.contains("shapes"), "{}", text);
}

#[test]
fn optimise_and_render() {
    let dir = output_dir("optimise");
    let mut args = vec![
        "optimise", "--shape", "polygon", "--sides", "5", "--group", "p2", "--sites", "e", "-o",
    ];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    assert_success(&packing(&args));
    let best = dir.join("best.json");
    assert!(best.exists());
    let rendered = dir.join("rendered.svg");
    let output = packing(&[
        "render",
        best.to_str().unwrap(),
        "-o",
        rendered.to_str().unwrap(),
    ]);
    assert_success(&output);
    assert!(fs::read_to_string(rendered).unwrap().contains("<polygon"));
}

#[test]
fn rejects_invalid_shape() {
    let dir = output_dir("invalid");
    fs::create_dir_all(&dir).unwrap();
    let vertices = dir.join("bow-tie.json");
    fs::write(&vertices, "[[0, 0], [2, 1], [2, 0], [0, 1.5]]").unwrap();
    let output = packing(&[
        "search",
        "--shape",
        "vertices",
        "--vertices",
        vertices.to_str().unwrap(),
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    let error = String::from_utf8_lossy(&output.stderr);
    assert!(error.contains("not simple"), "{}", error);
}

#[test]
fn lennard_jones_reports_energy() {
    let dir = output_dir("lj");
    let mut args = vec![
        "search",
        "--shape",
        "trimer",
        "--potential",
        "lj",
        "--groups",
        "p1,p2",
        "-o",
    ];
    args.push(dir.to_str().unwrap());
    args.extend(FAST);
    let output = packing(&args);
    assert_success(&output);
    assert!(
        stdout(&output).contains("Lowest energy"),
        "{}",
        stdout(&output)
    );
}
