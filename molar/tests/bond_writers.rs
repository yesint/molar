//! Shared selection bond remapping for every writer that emits connectivity.
use molar::prelude::*;
use std::{
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};
struct Output(PathBuf);
impl Output {
    fn new(ext: &str) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        Self(std::env::temp_dir().join(format!(
            "molar-bonds-{}-{}.{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed),
            ext
        )))
    }
}
impl Drop for Output {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}
fn system() -> System {
    let atoms = (0..6).map(|i| {
        Atom::new()
            .with_name(&format!("C{i}"))
            .with_resname("MOL")
            .with_resid(1)
            .with_chain('A')
            .with_atomic_number(6)
            .with_mass(12.0)
    });
    let top = Topology {
        atoms: atoms.collect(),
        bonds: [
            Bond::with_order(0, 1, BondOrder::Single),
            Bond::with_order(1, 3, BondOrder::Double),
            Bond::with_order(3, 5, BondOrder::Triple),
            Bond::with_order(2, 5, BondOrder::Aromatic),
            Bond::with_order(4, 5, BondOrder::Single),
        ]
        .into_iter()
        .collect(),
        ..Default::default()
    };
    let state = State {
        coords: (0..6).map(|i| Pos::new(i as Float, 1.0, 2.0)).collect(),
        pbox: Some(PeriodicBox::from_matrix(Matrix3f::identity() * 10.0).unwrap()),
        ..Default::default()
    };
    System::new(top, state).unwrap()
}
#[test]
fn bond_writers_roundtrip_selected_connectivity() {
    let system = system();
    let cases = [
        (
            vec![1, 3, 5],
            vec![(0, 1, BondOrder::Double), (1, 2, BondOrder::Triple)],
        ),
        (
            vec![3, 4, 5],
            vec![(0, 2, BondOrder::Triple), (1, 2, BondOrder::Single)],
        ),
        (vec![5], vec![]),
        (
            (0..6).collect(),
            vec![
                (0, 1, BondOrder::Single),
                (1, 3, BondOrder::Double),
                (3, 5, BondOrder::Triple),
                (2, 5, BondOrder::Aromatic),
                (4, 5, BondOrder::Single),
            ],
        ),
    ];
    for ext in ["sdf", "mol", "cif", "mmcif", "data"] {
        for (indices, expected) in &cases {
            let selection = system.select_bound(indices.clone()).unwrap();
            let output = Output::new(ext);
            {
                FileHandler::create(&output.0)
                    .unwrap()
                    .write(&selection)
                    .unwrap();
            }
            let (top, state) = FileHandler::open(&output.0).unwrap().read().unwrap();
            assert_eq!(top.len(), indices.len(), "{ext}");
            assert_eq!(top.bonds.len(), expected.len(), "{ext} {indices:?}");
            let mut actual: Vec<_> = top
                .bonds
                .iter()
                .map(|b| (b.i1(), b.i2(), b.order()))
                .collect();
            let mut expected: Vec<_> = expected
                .iter()
                .map(|&(i, j, order)| {
                    (
                        i,
                        j,
                        if ext == "data" {
                            BondOrder::Unspecified
                        } else {
                            order
                        },
                    )
                })
                .collect();
            actual.sort_by_key(|&(i, j, _)| (i, j));
            expected.sort_by_key(|&(i, j, _)| (i, j));
            assert_eq!(actual, expected, "{ext} {indices:?}");
            for (pos, &index) in state.coords.iter().zip(indices) {
                assert!(
                    (pos - system.get_pos(index).unwrap()).norm() < 1e-4,
                    "{ext}"
                );
            }
        }
    }
}
#[test]
fn shared_bond_mapping_retains_orders_and_checks_endpoints() {
    let system = system();
    let selection = system.select_bound(vec![1, 3, 5]).unwrap();
    let bonds = selection.bonds_for_write().unwrap();
    assert_eq!(
        bonds,
        vec![
            Bond::with_order(0, 1, BondOrder::Double),
            Bond::with_order(1, 2, BondOrder::Triple)
        ]
    );
    drop(selection);
    let mut top = Topology {
        atoms: system.iter_atoms().map(|atom| Atom::from(&atom)).collect(),
        ..Default::default()
    };
    top.bonds.push(&Bond::new(0, 999));
    assert!(matches!(
        top.bonds_for_write(),
        Err(FileFormatError::InvalidBondEndpoints(0, 999))
    ));
}
