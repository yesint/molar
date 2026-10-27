//! PDB alternate positions must not become duplicate atoms in one conformation.
use molar::prelude::*;
use std::io::Cursor;

fn reader(text: String) -> FileHandler {
    FileHandler::from_reader("pdb", Cursor::new(text.into_bytes())).unwrap()
}

fn atom(record: &str, serial: usize, name: &str, alt: char, x: f32, occupancy: f32) -> String {
    format!(
        "{record:<6}{serial:5} {name:<4}{alt}CYS A 917    {x:8.3}{:8.3}{:8.3}{occupancy:6.2}{:6.2}           C\n",
        0.0, 0.0, 0.0
    )
}

#[test]
fn pdb_keeps_blank_and_a_sites_regardless_of_record_type_or_occupancy() {
    let mut text = String::new();
    for record in ["ATOM", "HETATM"] {
        for (i, alt) in [' ', 'A', 'B', 'C', '1'].into_iter().enumerate() {
            text += &atom(
                record,
                text.lines().count() + 1,
                "CB",
                alt,
                i as f32,
                if alt == 'A' { 0.25 } else { 0.75 },
            );
        }
    }
    let (top, state) = reader(text).read().unwrap();
    assert_eq!(top.len(), 4);
    assert_eq!(state.coords.len(), 4);
    for (i, p) in state.coords.iter().enumerate() {
        assert!((p.x - (i % 2) as Float * 0.1).abs() < 1e-6);
    }
    assert_eq!(top.atoms.get(1).unwrap().get_occupancy(), 0.25);
    assert_eq!(top.atoms.get(3).unwrap().get_occupancy(), 0.25);
}

#[test]
fn pdb_conect_skips_discarded_serials_and_remaps_retained_sites() {
    let text = atom("ATOM", 10, "CA", ' ', 0.0, 1.0)
        + &atom("ATOM", 20, "CB", 'B', 1.0, 0.5)
        + &atom("ATOM", 30, "CB", 'A', 1.1, 0.5)
        + &atom("HETATM", 40, "SG", 'A', 2.0, 0.5)
        + &atom("HETATM", 50, "SG", 'B', 2.1, 0.5)
        + "CONECT   10   20   30\nCONECT   20   50\nCONECT   30   40   50\nCONECT   50   40\nEND\n";
    let (top, state) = reader(text).read().unwrap();
    assert_eq!(top.len(), 3);
    assert_eq!(state.coords.len(), 3);
    assert_eq!(
        top.bonds.iter_pairs().collect::<Vec<_>>(),
        vec![[0, 1], [1, 2]]
    );
}

#[test]
fn pdb_altloc_filter_applies_to_every_trajectory_model() {
    let mut text = String::new();
    for model in 1..=2 {
        text += &format!("MODEL     {model:4}\n");
        text += &atom("ATOM", 1, "CA", ' ', model as f32, 1.0);
        text += &atom("ATOM", 2, "CB", 'B', 8.0, 0.5);
        text += &atom("ATOM", 3, "CB", 'A', model as f32 + 1.0, 0.5);
        text += "ENDMDL\n";
    }
    let mut fh = reader(text);
    assert_eq!(fh.read_topology().unwrap().len(), 2);
    for model in 1..=2 {
        let state = fh.read_state().unwrap();
        assert_eq!(state.coords.len(), 2);
        assert!((state.coords[0].x - model as Float * 0.1).abs() < 1e-6);
        assert!((state.coords[1].x - (model + 1) as Float * 0.1).abs() < 1e-6);
    }
    assert!(matches!(
        fh.read_state().unwrap_err().kind(),
        FileFormatError::Eof
    ));
}

#[test]
fn pdb_model_with_only_rejected_alternates_reports_an_error() {
    let rejected = atom("ATOM", 1, "CB", 'B', 1.0, 1.0);
    let error = reader(rejected.clone()).read().unwrap_err();
    assert!(!matches!(error.kind(), FileFormatError::Eof));
    assert!(matches!(error.kind(), FileFormatError::Pdb(cause)
        if cause.to_string().contains("blank or A alternate location")));
    // An empty filtered model must not masquerade as EOF after a valid frame.
    let mut fh = reader(format!(
        "MODEL        1\n{}ENDMDL\nMODEL        2\n{rejected}ENDMDL\n",
        atom("ATOM", 1, "CB", 'A', 0.0, 1.0)
    ));
    fh.read().unwrap();
    assert!(!matches!(
        fh.read().unwrap_err().kind(),
        FileFormatError::Eof
    ));
}

#[test]
fn jak2_cys917_has_one_copy_of_each_atom_in_conformer_a() {
    let (top, state) = reader(include_str!("cys917_altloc.pdb").into())
        .read()
        .unwrap();
    let names: Vec<_> = top.atoms.iter().map(|a| a.get_name().to_string()).collect();
    assert_eq!(
        names,
        ["N", "CA", "C", "O", "CB", "SG", "HA", "HB1", "HB2", "HG"]
    );
    assert_eq!(state.coords.len(), names.len());
    // These are the A coordinates; the screenshot selected the rejected CB B site.
    assert!((state.coords[4] - Pos::new(4.8470, -1.2578, 0.7830)).norm() < 1e-6);
    assert!((state.coords[5] - Pos::new(4.9253, -1.1753, 0.9193)).norm() < 1e-6);
}
