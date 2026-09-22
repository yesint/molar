//! Bond-order / formal-charge perception parity harness (delivery step 4 of the
//! bond-perception plan).
//!
//! For every molecule in the GAFF reference corpus (`tests/data/gaff_ref/sdf/*.sdf`, which
//! carries real Kekulé bond orders and `M CHG` formal charges), the harness:
//!
//! 1. loads the molecule as the ground truth;
//! 2. erases every bond order and formal charge, keeping only connectivity;
//! 3. re-derives them with [`assign_bond_orders`], constrained to the molecule's known net
//!    charge;
//! 4. compares the result with the truth.
//!
//! Comparison canonicalizes aromaticity on both sides with [`perceive`] (so the arbitrary
//! Kekulé resonance form does not count as a difference). Two metrics are reported:
//!
//! - **strict**: bond order equal at every bond and formal charge equal at every atom;
//! - **pattern**: resonance-tolerant — the two structures have the same *multiset* of per-atom
//!   signatures `(element, sorted incident bond orders, formal charge)` and the same net
//!   charge. Because the signatures are matched as an unordered multiset, relabeling a
//!   symmetric group (which oxygen of a carboxylate carries the double bond, which nitrogen of
//!   a guanidinium carries the charge) does not count as a difference — only a genuinely
//!   different distribution of orders and charges does.
//!
//! `bond_order_parity_report` always passes and prints the report (run with `--nocapture`); it
//! is the development driver. `bond_order_parity_threshold` asserts the current target.

use std::collections::BTreeMap;

use molar::prelude::*;

const SDF_DIR: &str = "tests/data/gaff_ref/sdf";

/// Resonance-tolerant per-molecule target.
const TARGET_PATTERN: f64 = 0.90;

#[derive(Default)]
struct Stats {
    mols: usize,
    load_err: usize,
    solve_err: usize,
    /// solve-error kind -> count
    solve_err_kind: BTreeMap<String, usize>,
    bonds: usize,
    bonds_match: usize,
    fc_atoms: usize,
    fc_match: usize,
    strict_perfect: usize,
    pattern_perfect: usize,
    /// molecule -> short reason, for pattern failures.
    samples: Vec<(String, String)>,
}

fn order_key(o: BondOrder) -> u8 {
    match o {
        BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}

/// A resonance-invariant fingerprint of a structure: the sorted multiset of per-atom
/// signatures `(element, sorted incident bond orders, formal charge)`. Two assignments with
/// equal fingerprints differ at most by relabeling chemically equivalent atoms.
fn fingerprint(top: &Topology) -> Vec<(u8, Vec<u8>, i32)> {
    let n = top.atoms.len();
    let mut orders = vec![Vec::new(); n];
    for b in top.bonds.iter() {
        let [i, j] = b.pair();
        let k = order_key(b.order());
        orders[i].push(k);
        orders[j].push(k);
    }
    let mut sig: Vec<(u8, Vec<u8>, i32)> = top
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            orders[i].sort_unstable();
            (
                a.get_atomic_number(),
                std::mem::take(&mut orders[i]),
                a.get_formal_charge().unwrap_or(0),
            )
        })
        .collect();
    sig.sort();
    sig
}

fn charges(top: &Topology) -> Vec<i32> {
    top.atoms.iter().map(|a| a.get_formal_charge().unwrap_or(0)).collect()
}

fn run() -> Stats {
    let mut st = Stats::default();
    let mut names: Vec<String> = std::fs::read_dir(SDF_DIR)
        .unwrap_or_else(|e| panic!("cannot read {SDF_DIR}: {e}"))
        .filter_map(|e| e.ok())
        .filter_map(|e| e.file_name().into_string().ok())
        .filter(|n| n.ends_with(".sdf"))
        .collect();
    names.sort();

    for name in &names {
        st.mols += 1;
        let path = format!("{SDF_DIR}/{name}");

        let sys = match System::from_file(&path) {
            Ok(s) => s,
            Err(_) => {
                st.load_err += 1;
                continue;
            }
        };

        // Ground truth (aromatized for comparison).
        let mut truth = sys.topology().clone();
        let net: i32 = charges(&truth).iter().sum();

        // Strip orders and charges, keeping connectivity.
        let mut stripped = truth.clone();
        for b in 0..stripped.bonds.len() {
            stripped.bonds.set_order(b, BondOrder::Unspecified);
        }
        for a in 0..stripped.atoms.len() {
            stripped.atoms.get_mut(a).unwrap().set_formal_charge(0);
        }

        // Re-derive, constrained to the known net charge (fall back to unconstrained when the
        // molecule is more than one bonded fragment, which the v1 constraint cannot span).
        let opts_constrained = BondOrderOptions {
            total_charge: Some(net),
            ..BondOrderOptions::default()
        };
        let assignment = match assign_bond_orders(&stripped, None, &opts_constrained) {
            Ok(a) => Ok(a),
            Err(BondPerceptionError::TotalChargeWithMultipleComponents) => {
                assign_bond_orders(&stripped, None, &BondOrderOptions::default())
            }
            Err(e) => Err(e),
        };
        let assignment = match assignment {
            Ok(a) => a,
            Err(e) => {
                st.solve_err += 1;
                let kind = format!("{e:?}");
                let kind = kind.split_whitespace().next().unwrap_or("?").to_string();
                *st.solve_err_kind.entry(kind).or_default() += 1;
                if st.samples.len() < 60 {
                    st.samples.push((name.clone(), format!("solve error: {e}")));
                }
                continue;
            }
        };
        assignment.apply_to(&mut stripped).expect("apply on the same topology");

        // Canonicalize aromaticity on both sides.
        perceive(&mut truth);
        perceive(&mut stripped);

        // Strict per-bond / per-atom comparison.
        let mut strict_ok = true;
        for b in 0..truth.bonds.len() {
            st.bonds += 1;
            if truth.bonds.get(b).unwrap().order() == stripped.bonds.get(b).unwrap().order() {
                st.bonds_match += 1;
            } else {
                strict_ok = false;
            }
        }
        let truth_fc = charges(&truth);
        let solved_fc = charges(&stripped);
        for a in 0..truth.atoms.len() {
            st.fc_atoms += 1;
            if truth_fc[a] == solved_fc[a] {
                st.fc_match += 1;
            } else {
                strict_ok = false;
            }
        }
        if strict_ok {
            st.strict_perfect += 1;
        }

        // Resonance-tolerant comparison.
        let pattern_ok = fingerprint(&truth) == fingerprint(&stripped)
            && truth_fc.iter().sum::<i32>() == solved_fc.iter().sum::<i32>();
        if pattern_ok {
            st.pattern_perfect += 1;
        } else if st.samples.len() < 60 {
            st.samples.push((name.clone(), format!("fingerprint differs (net truth {net})")));
        }
    }
    st
}

fn print_report(st: &Stats) {
    let solved = st.mols - st.load_err - st.solve_err;
    let bacc = st.bonds_match as f64 / st.bonds.max(1) as f64;
    let facc = st.fc_match as f64 / st.fc_atoms.max(1) as f64;
    println!("\n=== bond-order parity report ===");
    println!("molecules        : {}", st.mols);
    println!("  load errors    : {}", st.load_err);
    println!("  solve errors   : {}  {:?}", st.solve_err, st.solve_err_kind);
    println!("  solved         : {solved}");
    println!(
        "bond order (strict, after aromatization): {}/{} = {:.2}%",
        st.bonds_match,
        st.bonds,
        bacc * 100.0
    );
    println!(
        "formal charge (strict, per atom)        : {}/{} = {:.2}%",
        st.fc_match,
        st.fc_atoms,
        facc * 100.0
    );
    println!(
        "molecules perfect (strict)   : {}/{} = {:.2}%",
        st.strict_perfect,
        st.mols,
        st.strict_perfect as f64 / st.mols.max(1) as f64 * 100.0
    );
    println!(
        "molecules perfect (pattern)  : {}/{} = {:.2}%",
        st.pattern_perfect,
        st.mols,
        st.pattern_perfect as f64 / st.mols.max(1) as f64 * 100.0
    );
    println!("\nsample failures (molecule : reason):");
    for (name, reason) in st.samples.iter().take(60) {
        println!("  {name}: {reason}");
    }
    println!("================================\n");
}

#[test]
fn bond_order_parity_report() {
    let st = run();
    print_report(&st);
}

/// The molecules that failed to solve in the first baseline. Fast to iterate on: it prints the
/// per-molecule outcome and time so search-robustness changes can be checked without the full
/// 597-molecule run.
#[test]
fn bond_order_hard_cases() {
    const HARD: &[&str] = &[
        "BleomycinA2",
        "Calcein",
        "Curare",
        "Itraconazole",
        "Ritonavir",
        "Vinblastine",
    ];
    for name in HARD {
        let path = format!("{SDF_DIR}/{name}.sdf");
        let sys = System::from_file(&path).expect("load");
        let mut stripped = sys.topology().clone();
        let net: i32 = charges(&stripped).iter().sum();
        for b in 0..stripped.bonds.len() {
            stripped.bonds.set_order(b, BondOrder::Unspecified);
        }
        for a in 0..stripped.atoms.len() {
            stripped.atoms.get_mut(a).unwrap().set_formal_charge(0);
        }
        let opts = BondOrderOptions {
            total_charge: Some(net),
            ..BondOrderOptions::default()
        };
        let t = std::time::Instant::now();
        let res = match assign_bond_orders(&stripped, None, &opts) {
            Err(BondPerceptionError::TotalChargeWithMultipleComponents) => {
                assign_bond_orders(&stripped, None, &BondOrderOptions::default())
            }
            other => other,
        };
        let dt = t.elapsed();
        match res {
            Ok(a) => println!(
                "{name:<14} OK in {:?}  warnings={:?}",
                dt,
                a.warnings()
            ),
            Err(e) => println!("{name:<14} ERR in {:?}: {e}", dt),
        }
    }
}

#[test]
fn bond_order_parity_threshold() {
    let st = run();
    print_report(&st);
    let pattern = st.pattern_perfect as f64 / st.mols.max(1) as f64;
    assert!(
        pattern >= TARGET_PATTERN,
        "resonance-tolerant per-molecule parity {:.2}% < target {:.2}%",
        pattern * 100.0,
        TARGET_PATTERN * 100.0
    );
}

/// Solve a molecule from connectivity only (orders and charges erased), constrained to its
/// known net charge, and return truth + solved topologies, both aromatized for comparison.
fn solve_pair(name: &str) -> (Topology, Topology) {
    let path = format!("{SDF_DIR}/{name}.sdf");
    let sys = System::from_file(&path).expect("load");
    let mut truth = sys.topology().clone();
    let net: i32 = charges(&truth).iter().sum();
    let mut stripped = truth.clone();
    for b in 0..stripped.bonds.len() {
        stripped.bonds.set_order(b, BondOrder::Unspecified);
    }
    for a in 0..stripped.atoms.len() {
        stripped.atoms.get_mut(a).unwrap().set_formal_charge(0);
    }
    let opts = BondOrderOptions {
        total_charge: Some(net),
        ..BondOrderOptions::default()
    };
    let assignment = match assign_bond_orders(&stripped, None, &opts) {
        Err(BondPerceptionError::TotalChargeWithMultipleComponents) => {
            assign_bond_orders(&stripped, None, &BondOrderOptions::default())
        }
        other => other,
    }
    .expect("solve");
    assignment.apply_to(&mut stripped).expect("apply");
    perceive(&mut truth);
    perceive(&mut stripped);
    (truth, stripped)
}

/// Prints, for each genuinely-failing molecule, the per-atom signatures present in the truth
/// structure but not the solved one (and vice versa). This names the functional group the
/// solver assigns differently, to guide the later functional-group work (plan step 5).
#[test]
fn bond_order_failure_diff() {
    const FAILING: &[&str] = &[
        "Omeprazole",
        "Lansoprazole",
        "Pantoprazole",
        "Mesoridazine",
        "Zidovudine",
        "Bremazocine",
    ];
    // multiset difference a \ b
    fn diff(
        a: &[(u8, Vec<u8>, i32)],
        b: &[(u8, Vec<u8>, i32)],
    ) -> Vec<(u8, Vec<u8>, i32)> {
        let mut rest = b.to_vec();
        let mut out = Vec::new();
        for x in a {
            if let Some(pos) = rest.iter().position(|y| y == x) {
                rest.remove(pos);
            } else {
                out.push(x.clone());
            }
        }
        out
    }
    for name in FAILING {
        let (truth, solved) = solve_pair(name);
        let tf = fingerprint(&truth);
        let sf = fingerprint(&solved);
        println!("{name}:");
        println!("  truth-only : {:?}", diff(&tf, &sf));
        println!("  solved-only: {:?}", diff(&sf, &tf));
    }
}
