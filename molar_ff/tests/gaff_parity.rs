//! GAFF atom-typing parity harness.
//!
//! Loads every molecule in `tests/data/gaff_ref/references.json`, runs
//! `apply_ff(FFType::Gaff)` on the corresponding SDF, and compares the assigned
//! `type_name` against antechamber's reference type (atom index `i` is aligned between
//! the SDF and the reference, so this is a direct per-atom comparison).
//!
//! - `gaff_parity_report` always passes and prints the accuracy / confusion report
//!   (run with `--nocapture`). It is the development driver.
//! - `gaff_parity_threshold` asserts the accuracy target and is expected to fail until
//!   the typing port is complete.

use std::collections::BTreeMap;

use molar::prelude::*;
use molar_ff::{ApplyFF, FFType};

const REFS: &str = "tests/data/gaff_ref/references.json";
const REFS_GAFF2: &str = "tests/data/gaff_ref/references_gaff2.json";
const SDF_DIR: &str = "tests/data/gaff_ref/sdf";

/// Per-atom accuracy required for the strict test.
const TARGET: f64 = 0.995;

#[derive(serde::Deserialize)]
struct RefFile {
    molecules: Vec<RefMol>,
}

#[derive(serde::Deserialize)]
struct RefMol {
    name: String,
    atoms: Vec<RefAtom>,
}

#[derive(serde::Deserialize)]
struct RefAtom {
    gaff: String,
    #[allow(dead_code)]
    el: String,
}

#[derive(Default)]
struct Stats {
    mols: usize,
    mols_perfect: usize,
    load_err: usize,
    apply_err: usize,
    atoms: usize,
    matched: usize,
    /// (expected, got) -> count, for mismatches only.
    confusion: BTreeMap<(String, String), usize>,
    /// (molecule, local index, element-ish, expected, got) samples.
    samples: Vec<(String, usize, String, String)>,
}

fn load_refs(path: &str) -> RefFile {
    let s = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    serde_json::from_str(&s).unwrap_or_else(|e| panic!("cannot parse {path}: {e}"))
}

fn run_parity(refs_path: &str, sdf_dir: &str, ff: FFType) -> Stats {
    let refs = load_refs(refs_path);
    let mut st = Stats::default();

    for mol in &refs.molecules {
        st.mols += 1;
        let path = format!("{sdf_dir}/{}.sdf", mol.name);

        let mut sys = match System::from_file(&path) {
            Ok(s) => s,
            Err(_) => {
                st.load_err += 1;
                continue;
            }
        };

        if sys.apply_ff(ff).is_err() {
            st.apply_err += 1;
            continue;
        }

        let got: Vec<String> =
            sys.iter_atoms().map(|a| a.get_type_name().unwrap_or("").to_string()).collect();

        if got.len() != mol.atoms.len() {
            // atom-count mismatch between reader and reference: count all as wrong
            st.samples.push((
                mol.name.clone(),
                usize::MAX,
                format!("natom {} != {}", got.len(), mol.atoms.len()),
                String::new(),
            ));
            st.atoms += mol.atoms.len();
            continue;
        }

        let mut mol_ok = true;
        for (i, exp) in mol.atoms.iter().enumerate() {
            st.atoms += 1;
            if got[i] == exp.gaff {
                st.matched += 1;
            } else {
                mol_ok = false;
                *st.confusion
                    .entry((exp.gaff.clone(), got[i].clone()))
                    .or_default() += 1;
                if st.samples.len() < 40 {
                    st.samples.push((
                        mol.name.clone(),
                        i,
                        exp.gaff.clone(),
                        got[i].clone(),
                    ));
                }
            }
        }
        if mol_ok {
            st.mols_perfect += 1;
        }
    }
    st
}

fn print_report(st: &Stats) {
    let acc = if st.atoms > 0 {
        st.matched as f64 / st.atoms as f64
    } else {
        0.0
    };
    println!("\n=== GAFF parity report ===");
    println!("molecules      : {}", st.mols);
    println!("  perfect      : {}", st.mols_perfect);
    println!("  load errors  : {}", st.load_err);
    println!("  apply errors : {}", st.apply_err);
    println!(
        "atoms          : {} matched / {} = {:.2}%",
        st.matched,
        st.atoms,
        acc * 100.0
    );

    // Top mismatches (expected -> got) by frequency.
    let mut conf: Vec<_> = st.confusion.iter().collect();
    conf.sort_by(|a, b| b.1.cmp(a.1));
    println!("\ntop mismatches (expected -> got : count):");
    for ((exp, got), c) in conf.iter().take(25) {
        println!("  {exp:>4} -> {:<4} : {c}", if got.is_empty() { "∅" } else { got });
    }

    println!("\nsample mismatches (molecule / atom / expected / got):");
    for (name, i, exp, got) in st.samples.iter().take(40) {
        if *i == usize::MAX {
            println!("  {name}: {exp}");
        } else {
            println!("  {name} #{i}: {exp} -> {}", if got.is_empty() { "∅" } else { got });
        }
    }
    println!("==========================\n");
}

#[test]
fn gaff_parity_report() {
    let st = run_parity(REFS, SDF_DIR, FFType::Gaff);
    print_report(&st);
}

#[test]
fn gaff_parity_threshold() {
    let st = run_parity(REFS, SDF_DIR, FFType::Gaff);
    print_report(&st);
    let acc = st.matched as f64 / st.atoms.max(1) as f64;
    assert_eq!(st.load_err, 0, "{} molecules failed to load", st.load_err);
    assert_eq!(st.apply_err, 0, "{} molecules failed apply_ff", st.apply_err);
    assert!(
        acc >= TARGET,
        "per-atom GAFF accuracy {:.2}% < target {:.2}%",
        acc * 100.0,
        TARGET * 100.0
    );
}

#[test]
fn gaff2_parity_report() {
    let st = run_parity(REFS_GAFF2, SDF_DIR, FFType::Gaff2);
    print_report(&st);
}

#[test]
fn gaff2_parity_threshold() {
    let st = run_parity(REFS_GAFF2, SDF_DIR, FFType::Gaff2);
    print_report(&st);
    let acc = st.matched as f64 / st.atoms.max(1) as f64;
    assert_eq!(st.load_err, 0, "{} molecules failed to load", st.load_err);
    assert_eq!(st.apply_err, 0, "{} molecules failed apply_ff", st.apply_err);
    assert!(
        acc >= TARGET,
        "per-atom GAFF2 accuracy {:.2}% < target {:.2}%",
        acc * 100.0,
        TARGET * 100.0
    );
}

/// The blanket impl must resolve both on `System` and on a bound selection, and give
/// identical results when the selection covers the whole system.
#[test]
fn apply_ff_resolves_on_system_and_selection() {
    let refs = load_refs(REFS);
    let name = &refs.molecules[0].name;
    let path = format!("{SDF_DIR}/{name}.sdf");

    let mut a = System::from_file(&path).expect("load");
    a.apply_ff(FFType::Gaff).expect("apply on System");
    let types_sys: Vec<String> = a.iter_atoms().map(|x| x.get_type_name().unwrap_or("").to_string()).collect();

    let mut b = System::from_file(&path).expect("load");
    b.select_all_bound_mut()
        .apply_ff(FFType::Gaff)
        .expect("apply on selection");
    let types_sel: Vec<String> = b.iter_atoms().map(|x| x.get_type_name().unwrap_or("").to_string()).collect();

    assert_eq!(types_sys, types_sel, "System and whole-selection typing must agree");
}

/// Fixtures drawn from AmberTools' own `antechamber` test suite (sustiva/tp/fpph),
/// converted from their `.ac` output so the SDF carries antechamber's own Kekulé bond
/// orders. Because the input bond orders are exactly what antechamber typed against,
/// these must reproduce its GAFF types atom-for-atom.
#[test]
fn antechamber_suite_parity() {
    const SUITE_REFS: &str = "tests/data/gaff_ref/antechamber_suite/references.json";
    const SUITE_DIR: &str = "tests/data/gaff_ref/antechamber_suite";
    let st = run_parity(SUITE_REFS, SUITE_DIR, FFType::Gaff);
    print_report(&st);
    assert_eq!(st.load_err, 0, "a fixture failed to load");
    assert_eq!(st.apply_err, 0, "a fixture failed apply_ff");
    assert_eq!(st.matched, st.atoms, "antechamber-suite fixtures must type exactly");
}
