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
    let best = near[0]["packing_fraction"].as_f64().unwrap();
    for entry in near {
        let phi = entry["packing_fraction"].as_f64().unwrap();
        assert!(phi >= 0.98 * best - 1e-12, "{} vs {}", phi, best);
    }
    let images = fs::read_dir(dir.join("near-best")).unwrap().count();
    // A JSON and an SVG file for each structure
    assert_eq!(images, 2 * near.len());
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
