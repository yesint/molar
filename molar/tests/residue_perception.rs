//! Residue-template perception through the full solver (delivery step 6).
//!
//! Builds a hydrogen-free di-glycine backbone and perceives it with residue templates on. The
//! templates pin each backbone carbonyl, the peptide link and the rest are solved from geometry,
//! and the missing hydrogens are inferred. This is the tractable end-to-end check; perceiving a
//! whole protein in one connected fragment is a solver-scaling task tracked separately.

use molar::prelude::*;

#[test]
fn diglycine_backbone_perceived_with_residue_templates() {
    // Two glycine residues, heavy atoms only: N-CA-C(=O)-N-CA-C(=O).
    let names = ["N", "CA", "C", "O", "N", "CA", "C", "O"];
    let resindex = [0usize, 0, 0, 0, 1, 1, 1, 1];
    let z = [7u8, 6, 6, 8, 7, 6, 6, 8];
    // Extended backbone with ~0.15 nm links (single by geometry) and short C=O carbonyls.
    let xyz = [
        [0.00, 0.0, 0.0],
        [0.15, 0.0, 0.0],
        [0.30, 0.0, 0.0],
        [0.30, 0.123, 0.0], // O3 carbonyl
        [0.45, 0.0, 0.0],
        [0.60, 0.0, 0.0],
        [0.75, 0.0, 0.0],
        [0.75, 0.123, 0.0], // O7 carbonyl
    ];
    let bonds = [
        [0, 1],
        [1, 2],
        [2, 3], // C2=O3
        [2, 4], // peptide bond
        [4, 5],
        [5, 6],
        [6, 7], // C6=O7
    ];

    let mut top = Topology::default();
    for k in 0..8 {
        top.atoms.push(
            &Atom::new()
                .with_atomic_number(z[k])
                .with_name(names[k])
                .with_resname("GLY")
                .with_resindex(resindex[k]),
        );
    }
    for &[i, j] in &bonds {
        top.bonds.push(&Bond::new(i, j));
    }
    let mut state = State::default();
    state.coords = xyz.iter().map(|p| Pos::new(p[0], p[1], p[2])).collect();
    let sys = System::new(top, state).unwrap();

    let opts = BondOrderOptions {
        hydrogens: HydrogenPolicy::InferFromGeometry,
        total_charge: None,
        ..BondOrderOptions::default()
    };
    let assignment = sys.assign_bond_orders(&opts).expect("perceive di-glycine");

    let mut solved = sys.topology().clone();
    assignment.apply_to(&mut solved).unwrap();
    let order = |b: usize| solved.bonds.get(b).unwrap().order();

    // Residue templates pin both backbone carbonyls double.
    assert_eq!(order(2), BondOrder::Double, "res0 C=O");
    assert_eq!(order(6), BondOrder::Double, "res1 C=O");
    // The peptide bond and the N-CA / CA-C links stay single.
    assert_eq!(order(3), BondOrder::Single, "peptide C-N");
    assert_eq!(order(0), BondOrder::Single, "N-CA");
    assert_eq!(order(1), BondOrder::Single, "CA-C");

    // Charges are neutral, and the inferred hydrogens complete the backbone: the amine N takes
    // two, the amide N one, each CA two, each terminal carbonyl carbon one.
    let ih = assignment.implicit_hydrogens();
    assert!(solved.atoms.iter().all(|a| a.get_formal_charge().unwrap_or(0) == 0));
    assert_eq!(ih[0], 2, "N-terminal amine");
    assert_eq!(ih[4], 1, "amide nitrogen");
    assert_eq!(ih[1], 2, "CA");
    assert_eq!(ih[2], 0, "internal carbonyl carbon (=O, C-CA, C-N: valence full)");
    assert_eq!(ih[6], 1, "C-terminal carbonyl carbon (aldehyde end)");
}
