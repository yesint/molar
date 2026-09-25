//! Aromatic input that has lost a ring N-H must be refused by both typing and charging, with an
//! error that points at the nitrogen. The fixtures are real pipeline failures: PDB ligands written
//! as aromatic (order-4) mol blocks without the N-H as an explicit atom.

use molar::prelude::*;
use molar_ff::{ApplyCharges, ApplyFF, ChargeModel, FFError, FFType};

const DIR: &str = "tests/data/invalid_aromatic";

/// Tryptophan with the indole N-H missing (3zte): 9 atoms of the indole need a double bond.
#[test]
fn missing_indole_nh_is_rejected_by_typing_and_charging() {
    let path = format!("{DIR}/trp_no_indole_nh.sdf");

    let mut sys = System::from_file(&path).unwrap();
    let err = sys.apply_ff(FFType::Gaff).unwrap_err();
    let FFError::InvalidAromatic(KekulizeError::NonKekulizable(systems)) = &err else {
        panic!("expected InvalidAromatic(NonKekulizable), got {err:?}");
    };
    assert_eq!(systems.len(), 1);
    assert_eq!(systems[0].n_demanding, 9);
    assert_eq!(
        systems[0].h_candidates,
        vec![5],
        "the indole N is atom 5 (0-based)"
    );
    assert!(
        err.to_string().contains("hydrogen is probably missing"),
        "{err}"
    );

    let mut sys = System::from_file(&path).unwrap();
    assert!(sys.apply_charges(ChargeModel::Espaloma).is_err());
}

/// The same molecule with the N-H restored types and charges, and the N gets `na`.
#[test]
fn restored_indole_nh_types_and_charges() {
    let mut sys = System::from_file(format!("{DIR}/trp.sdf")).unwrap();
    sys.apply_ff(FFType::Gaff).unwrap();
    assert_eq!(sys.get_atom(5).unwrap().get_type_name(), Some("na"));
    sys.apply_charges(ChargeModel::Espaloma).unwrap();
}

/// AIQ (1q63) lost two N-H: one on an imidazole, one on a quinazolinone. Both systems are
/// reported in one error, not only the first.
#[test]
fn every_broken_aromatic_system_is_reported() {
    let mut sys = System::from_file(format!("{DIR}/aiq_two_nh_missing.sdf")).unwrap();
    let err = sys.apply_charges(ChargeModel::Espaloma).unwrap_err();
    let molar_ff::ChargeError::Kekulize(KekulizeError::NonKekulizable(systems)) = &err else {
        panic!("expected Kekulize(NonKekulizable), got {err:?}");
    };
    let mut sizes: Vec<usize> = systems.iter().map(|s| s.n_demanding).collect();
    sizes.sort_unstable();
    assert_eq!(sizes, vec![5, 9]);
    assert!(
        systems.iter().all(|s| s.h_candidates.len() == 2),
        "{systems:?}"
    );
}

/// GDP (7f0w) lost the guanine N1-H. N1 and N3 are both valid homes for it; N7 is a bare
/// aromatic N too, but an H there leaves the six-ring unpaired, so it must not be suggested.
#[test]
fn only_nitrogens_that_fix_the_system_are_suggested() {
    let mut sys = System::from_file(format!("{DIR}/gdp_no_n1h.sdf")).unwrap();
    let err = sys.apply_charges(ChargeModel::Espaloma).unwrap_err();
    let molar_ff::ChargeError::Kekulize(KekulizeError::NonKekulizable(systems)) = &err else {
        panic!("expected Kekulize(NonKekulizable), got {err:?}");
    };
    assert_eq!(
        systems[0].h_candidates,
        vec![2, 27],
        "N3 and N1, not N7 (atom 5)"
    );
}

/// Strip every hydrogen and bond order, then re-infer them from heavy-atom geometry. The
/// fixtures list their hydrogens after the heavy atoms, so heavy-atom indices are unchanged.
fn reinferred_hydrogens(file: &str) -> Vec<u8> {
    reinferred_hydrogens_at(&format!("{DIR}/{file}"))
}

fn reinferred_hydrogens_at(path: &str) -> Vec<u8> {
    let mut sys = System::from_file(path).unwrap();
    let hs: Vec<usize> = sys
        .iter_atoms()
        .enumerate()
        .filter(|(_, a)| a.get_atomic_number() == 1)
        .map(|(i, _)| i)
        .collect();
    sys.remove(hs.into_iter()).unwrap();
    let mut top = sys.topology().clone();
    for b in 0..top.bonds.len() {
        top.bonds.set_order(b, BondOrder::Unspecified);
    }
    let sys = System::new(top, sys.state().clone()).unwrap();
    let opts = BondOrderOptions {
        hydrogens: HydrogenPolicy::InferFromGeometry,
        ..BondOrderOptions::default()
    };
    sys.assign_bond_orders(&opts)
        .unwrap()
        .implicit_hydrogens()
        .to_vec()
}

/// Re-inferring every hydrogen from heavy-atom geometry puts the lost N-H back on the indole
/// nitrogen. Before the sp2-carbon tie-break, the solver saturated the planar CG instead: the
/// same cost otherwise, one implicit H either way.
#[test]
fn stripped_trp_gets_its_indole_nh_from_geometry() {
    // Heavy atoms: N-alpha, CA, CB, CG, CD1, NE1, CE2, CZ2, CH2, CZ3, CE3, CD2, C, O, OXT.
    assert_eq!(
        reinferred_hydrogens("trp_no_indole_nh.sdf"),
        vec![2, 1, 2, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1],
        "C11H12N2O2 with the N-H on NE1"
    );
}

/// Guanine (GDP): the lactam N1-H, not N3-H, and the C6=O stays a carbonyl — the enol would
/// make the six-ring aromatic, but its 0.123 nm C-O bond rules the O-H out.
#[test]
fn stripped_guanine_takes_the_lactam_n1_h() {
    let ih = reinferred_hydrogens("gdp_no_n1h.sdf");
    assert_eq!(ih[27], 1, "N1 (file atom 28) carries the H");
    assert_eq!(ih[2], 0, "N3 does not");
    assert_eq!(ih[26], 0, "O6 stays a carbonyl");
    assert_eq!(ih[0], 2, "the exocyclic amine stays NH2");
}

/// 2-amino-6-chloropurine (6GU): the H goes on the five-ring (N9 here), keeping both rings
/// aromatic; N1-H or N3-H would break the six-ring.
#[test]
fn stripped_purine_keeps_both_rings_aromatic() {
    let ih = reinferred_hydrogens("6gu_no_nh.sdf");
    // Heavy atoms: N1 (amine), C2, N3, C4, Cl5, C6, N7, C8, N9, C10, N11 (file numbering).
    assert_eq!(ih[0], 2, "amine NH2");
    assert_eq!(ih[2] + ih[10], 0, "no H on a six-ring N");
    assert_eq!(ih[6] + ih[8], 1, "one H on a five-ring N");
}

/// AIQ: 3H-quinazolin-4-one (the N-H beside the C=O) plus an imidazole N-H.
#[test]
fn stripped_quinazolinone_takes_the_lactam_n3_h() {
    let ih = reinferred_hydrogens("aiq_two_nh_missing.sdf");
    assert_eq!(ih[15], 1, "N16, beside C17=O18, carries the H");
    assert_eq!(ih[12], 0, "N13 does not");
    assert_eq!(ih[17], 0, "O18 stays a carbonyl");
    assert_eq!(ih[7] + ih[10], 1, "one imidazole N-H");
}

/// Cytosine keeps its amino-oxo form. The lactam rule alone would pick imino-oxo (a ring N-H
/// beside the C=O plus an exocyclic C=NH), so the amino-over-imino rule has to come first.
#[test]
fn stripped_cytidine_keeps_amino_oxo() {
    let path = "tests/data/gaff_ref/sdf/2,3_Dideoxycytidine.sdf";
    let sys = System::from_file(path).unwrap();
    let z: Vec<u8> = sys.iter_atoms().map(|a| a.get_atomic_number()).collect();
    let mut truth = vec![0u8; z.len()];
    for [i, j] in sys.topology().bonds.iter_pairs() {
        if z[j] == 1 {
            truth[i] += 1;
        }
        if z[i] == 1 {
            truth[j] += 1;
        }
    }
    let heavy: Vec<u8> = (0..z.len())
        .filter(|&i| z[i] != 1)
        .map(|i| truth[i])
        .collect();
    assert_eq!(reinferred_hydrogens_at(path), heavy);
}
