use molar::prelude::*;
use std::{
    io::Cursor,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

const UNWRAPPED: &str = include_str!("lammps/unwrapped.data");
const WRAPPED: &str = include_str!("lammps/wrapped_images.data");
const SMALL: &str = "polymer\n\n3 atoms\n2 atom types\n1 bonds\n2 bond types\n\n-2 8 xlo xhi\n-3 7 ylo yhi\n-4 6 zlo zhi\n1 2 3 xy xz yz\n\nAtoms # molecular\n\n30 7 2 0 0 0 1 -1 2\n10 7 1 1 1 1 0 0 0\n20 0 1 2 2 2 0 0 0\n\nBonds\n\n5 2 30 10\n\nMasses\n\n2 2\n1 1\n";
fn reader(text: &str) -> FileHandler {
    FileHandler::from_reader("data", Cursor::new(text.as_bytes().to_vec())).unwrap()
}
struct Output(PathBuf);
impl Output {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        Self(std::env::temp_dir().join(format!(
            "molar-lammps-{}-{}.data",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        )))
    }
}
impl Drop for Output {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}
fn assert_pos(a: &Pos, b: &Pos) {
    assert!((a - b).norm() < 1e-4, "{a:?} != {b:?}");
}

#[test]
fn lammps_examples_are_equivalent() {
    let (top_a, state_a) = reader(UNWRAPPED).read().unwrap();
    let (top_b, state_b) = reader(WRAPPED).read().unwrap();
    for top in [&top_a, &top_b] {
        assert_eq!(top.atoms.len(), 420);
        assert_eq!(top.bonds.len(), 400);
        assert_eq!(top.molecules.len(), 20);
        for (i, m) in top.molecules.iter().enumerate() {
            assert_eq!(*m, [21 * i, 21 * i + 20]);
        }
        assert!(
            top.bonds
                .iter()
                .all(|b| b.order() == BondOrder::Unspecified)
        );
        assert!(
            top.atoms
                .iter()
                .all(|a| a.get_mass() == 1.0 && a.get_atomic_number() == 0)
        );
    }
    assert!(state_b.velocities.is_empty());
    for (a, b) in state_a.coords.iter().zip(&state_b.coords) {
        assert_pos(a, b);
    }
    assert_pos(&state_a.coords[0], &Pos::new(13.38118, 10.89952, 5.74699));
    assert_eq!(
        state_a.pbox.unwrap().get_matrix(),
        Matrix3f::identity() * 40.0
    );
    let pairs_a: Vec<_> = top_a.bonds.iter().map(|b| b.pair()).collect();
    let pairs_b: Vec<_> = top_b.bonds.iter().map(|b| b.pair()).collect();
    assert_eq!(pairs_a, pairs_b);
    assert_eq!(
        top_a
            .atoms
            .iter()
            .map(|a| a.get_type_id())
            .collect::<Vec<_>>(),
        top_b
            .atoms
            .iter()
            .map(|a| a.get_type_id())
            .collect::<Vec<_>>()
    );
}
#[test]
fn lammps_sparse_ids_triclinic_and_scales() {
    let mut h = FileHandler::from_lammps_reader(
        Cursor::new(SMALL),
        LammpsOptions {
            length_scale: 0.1,
            mass_scale: 3.0,
        },
    )
    .unwrap();
    let (top, state) = h.read().unwrap();
    assert_eq!(
        top.atoms.iter().map(|a| a.get_resid()).collect::<Vec<_>>(),
        vec![0, 7, 7]
    );
    assert_eq!(top.molecules, vec![[1, 2]]);
    assert_eq!(top.bonds.get(0).unwrap().pair(), [2, 1]);
    assert_eq!(top.bonds.get(0).unwrap().order(), BondOrder::Unspecified);
    assert_eq!(top.atoms.get(2).unwrap().get_mass(), 6.0);
    assert_pos(&state.coords[2], &Pos::new(1.5, -0.1, 2.4));
    assert!(matches!(h.read().unwrap_err().kind(), FileFormatError::Eof));
}
#[test]
fn lammps_separate_reads_and_iteration() {
    let mut h = reader(SMALL);
    assert_eq!(h.read_topology().unwrap().len(), 3);
    assert_eq!(h.read_state().unwrap().len(), 3);
    assert!(matches!(
        h.read_state().unwrap_err().kind(),
        FileFormatError::Eof
    ));
    let mut h = reader(SMALL);
    assert_eq!(h.read_state().unwrap().len(), 3);
    assert_eq!(h.read_topology().unwrap().len(), 3);
    assert!(matches!(
        h.read_topology().unwrap_err().kind(),
        FileFormatError::Eof
    ));
    assert_eq!(reader(SMALL).into_iter().count(), 1);
}
#[test]
fn lammps_roundtrip_and_single_write() {
    let (top, state) = reader(WRAPPED).read().unwrap();
    let system = System::new(top, state).unwrap();
    let output = Output::new();
    let mut w = FileHandler::create_lammps(
        &output.0,
        LammpsOptions {
            length_scale: 0.25,
            mass_scale: 2.0,
        },
    )
    .unwrap();
    w.write(&system).unwrap();
    assert!(w.write(&system).is_err());
    let (top, state) = FileHandler::open_lammps(
        &output.0,
        LammpsOptions {
            length_scale: 0.25,
            mass_scale: 2.0,
        },
    )
    .unwrap()
    .read()
    .unwrap();
    assert_eq!(top.bonds.len(), 400);
    for (a, b) in state.coords.iter().zip(system.iter_pos()) {
        assert_pos(a, b);
    }
    assert!(top.atoms.iter().all(|a| a.get_mass() == 1.0));
    let text = std::fs::read_to_string(&output.0).unwrap();
    assert!(text.contains("1 bond types"));
    assert!(text.contains("Atoms # molecular"));
}
#[test]
fn lammps_selection_remaps_and_filters_bonds() {
    let (top, state) = reader(UNWRAPPED).read().unwrap();
    let system = System::new(top, state).unwrap();
    let selection = system.select_bound("resid 2").unwrap();
    let output = Output::new();
    FileHandler::create(&output.0)
        .unwrap()
        .write(&selection)
        .unwrap();
    let (top, state) = FileHandler::open(&output.0).unwrap().read().unwrap();
    assert_eq!(top.atoms.len(), 21);
    assert_eq!(top.bonds.len(), 20);
    assert_eq!(top.bonds.get(0).unwrap().pair(), [0, 1]);
    for (a, b) in state.coords.iter().zip(selection.iter_pos()) {
        assert_pos(a, b);
    }
}
#[test]
fn lammps_separate_writes() {
    let (top, state) = reader(SMALL).read().unwrap();
    let output = Output::new();
    let mut w = FileHandler::create(&output.0).unwrap();
    w.write_topology(&top).unwrap();
    w.write_state(&state).unwrap();
    let (read_top, read_state) = FileHandler::open(&output.0).unwrap().read().unwrap();
    assert_eq!(
        read_top.bonds.get(0).unwrap().pair(),
        top.bonds.get(0).unwrap().pair()
    );
    for (a, b) in read_state.coords.iter().zip(&state.coords) {
        assert_pos(a, b);
    }
}
#[test]
fn lammps_rejects_invalid_data_with_location() {
    let cases = [
        SMALL.replace("30 7 2", "10 7 2"),        // duplicate ID
        SMALL.replace("30 7 2", "30 7 3"),        // type range
        SMALL.replace("5 2 30 10", "5 2 30 999"), // absent endpoint
        SMALL.replace("5 2 30 10", "5 2 30 30"),  // self bond
        SMALL.replace("molecular", "full"),
        SMALL.replace("30 7 2 0", "30 7 2 NaN"),
        SMALL.replace("-2 8 xlo xhi", "8 -2 xlo xhi"),
        SMALL.replace("30 7 2 0 0 0 1 -1 2", "30 7 2 0 0 0"),
        SMALL.replace("1 1\n", ""), // missing mass
        SMALL.replace("1 2 3 xy xz yz", "1 0 0 avec"),
        SMALL.replace("Bonds\n", "Custom Section\n"),
        SMALL.replace("3 atoms", "4 atoms"),
        SMALL.replace("2 2\n1 1", "2 -2\n1 1"),
        SMALL.replace("20 0 1", "20 -1 1"),
    ];
    for text in cases {
        let error = reader(&text).read().unwrap_err();
        assert!(
            matches!(
                error.kind(),
                FileFormatError::Lammps(LammpsHandlerError::Parse {
                    section: _,
                    line: _,
                    message: _
                })
            ),
            "{error:?}"
        );
    }
    assert!(reader("").read().is_err());
}
#[test]
fn lammps_scales_are_checked_before_creation() {
    let output = Output::new();
    for scale in [0.0, -1.0, Float::INFINITY, Float::NAN] {
        assert!(
            FileHandler::create_lammps(
                &output.0,
                LammpsOptions {
                    length_scale: scale,
                    ..Default::default()
                }
            )
            .is_err()
        );
        assert!(!output.0.exists());
        assert!(
            FileHandler::from_lammps_reader(
                Cursor::new(SMALL),
                LammpsOptions {
                    mass_scale: scale,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
}

#[test]
fn lammps_comments_styles_and_discarded_sections() {
    for comment in ["", " # id mol type xu yu zu", " # bond", " # angle"] {
        let data = UNWRAPPED.replace(" # id mol type xu yu zu", comment);
        assert_eq!(reader(&data).read().unwrap().0.len(), 420);
    }
    let data = SMALL.replace("2 bond types", "2 bond types\n1 angles\n1 angle types")
        + "\nPair Coeffs # lj/cut\n\n1 1 1\n2 2 2 # ignored\n\nBond Coeffs # fene\n\n1 30 1.5 1 1\n2 30 1.5 1 1\n\nAngles\n\n1 1 30 10 20\n\nAngle Coeffs\n\n1 20 180\n";
    assert_eq!(reader(&data).read().unwrap().0.bonds.len(), 1);
}

#[test]
fn lammps_writer_rejects_invalid_properties() {
    let (top, state) = reader(SMALL).read().unwrap();
    for case in 0..5 {
        let mut top = top.clone();
        let mut state = state.clone();
        match case {
            0 => state.pbox = None,
            1 => state.coords[0].x = Float::NAN,
            2 => top.atoms.get_mut(0).unwrap().set_mass(4.0),
            3 => top.atoms.get_mut(0).unwrap().set_resid(-1),
            4 => {
                state.pbox = Some(
                    PeriodicBox::from_matrix(Matrix3f::new(
                        10.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 10.0,
                    ))
                    .unwrap(),
                )
            }
            _ => unreachable!(),
        }
        let output = Output::new();
        let system = System::new(top, state).unwrap();
        assert!(
            FileHandler::create(&output.0)
                .unwrap()
                .write(&system)
                .is_err()
        );
        assert_eq!(std::fs::metadata(&output.0).unwrap().len(), 0);
    }
}
