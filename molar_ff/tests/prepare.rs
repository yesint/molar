//! Opt-in force-field preparation wrapper (delivery step 8).

use molar::prelude::*;
use molar_ff::{ApplyFF, FFType, PrepareForFF, PrepareOptions};

const ACETIC_ACID: &str = "tests/data/gaff_ref/sdf/Acetic_acid.sdf";

fn type_names(sys: &System) -> Vec<String> {
    sys.iter_atoms()
        .map(|a| a.get_type_name().unwrap_or("").to_string())
        .collect()
}

#[test]
fn prepare_enables_typing_of_an_order_less_input() {
    // Baseline: the SDF carries bond orders, so it types directly.
    let mut base = System::from_file(ACETIC_ACID).expect("load");
    base.apply_ff(FFType::Gaff).expect("baseline typing");
    let expected = type_names(&base);

    // Mimic an order-less input (PDB/GRO): keep connectivity, erase every bond order.
    let fresh = System::from_file(ACETIC_ACID).expect("load");
    let mut top = fresh.topology().clone();
    for b in 0..top.bonds.len() {
        top.bonds.set_order(b, BondOrder::Unspecified);
    }
    let mut sys = System::new(top, fresh.state().clone()).unwrap();

    // Strict typing rejects it.
    assert!(sys.apply_ff(FFType::Gaff).is_err(), "order-less input must be rejected");

    // Preparation perceives the orders; typing then succeeds and matches the baseline.
    sys.prepare_for_ff(&PrepareOptions::default()).expect("prepare");
    sys.apply_ff(FFType::Gaff).expect("typing after preparation");
    assert_eq!(type_names(&sys), expected, "prepared typing matches the SDF baseline");
}

#[test]
fn prepare_builds_a_full_molecule_from_bare_heavy_atoms() {
    // Six bare carbons on a hexagon, no bonds and no hydrogens — the hardest case.
    let r = 0.139;
    let mut top = Topology::default();
    let mut state = State::default();
    for k in 0..6 {
        let t = std::f64::consts::PI / 3.0 * k as f64;
        top.atoms.push(&Atom::new().with_atomic_number(6));
        state.coords.push(Pos::new(r * t.cos() as f32, r * t.sin() as f32, 0.0));
    }
    let mut sys = System::new(top, state).unwrap();

    let options = PrepareOptions {
        bond_orders: BondOrderOptions {
            hydrogens: HydrogenPolicy::InferFromGeometry,
            input_orders: InputOrders::PreserveKnown,
            ..BondOrderOptions::default()
        },
        add_hydrogens: Some(HydrogenOptions::default()),
        ..PrepareOptions::default()
    };
    sys.prepare_for_ff(&options)
        .expect("prepare perceives connectivity, orders and adds hydrogens");

    // Benzene: six carbons and six added hydrogens, all typed by GAFF.
    assert_eq!(sys.len(), 12, "C6H6");
    sys.apply_ff(FFType::Gaff).expect("typing the prepared benzene");
    assert!(
        sys.iter_atoms().all(|a| !a.get_type_name().unwrap_or("").is_empty()),
        "every atom, including the added hydrogens, is typed"
    );
    // The ring carbons take the aromatic-carbon GAFF type.
    for a in sys.topology().atoms.iter().take(6) {
        assert_eq!(a.get_type_name(), Some("ca"), "aromatic ring carbon");
    }
}
