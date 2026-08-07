//! SMILES writer + kekulizer parity harness over the 597-molecule GAFF corpus.
//!
//! Lives in `molar_ff` rather than `molar` because that is where the corpus is
//! (`tests/data/gaff_ref/`, excluded from the published package) and `molar_ff` already depends
//! on `molar`. Modelled on `gaff_parity.rs`: a `*_report` test that always passes and prints
//! statistics (the development driver, run with `--nocapture`) plus a strict threshold test.
//!
//! Two independent families of check, both dependency-free:
//!
//! **A. Kekulizer round trip.** Record the SDF's Kekulé orders, aromatize with `perceive` (which
//! destroys them by design), then `kekulize` and require the recovered structure to be *a* valid
//! Kekulé form — per-atom summed bond order preserved. Deliberately not order-for-order
//! identity: resonance forms are not unique and returning a different valid one is correct.
//!
//! **B. Writer graph fidelity.** Rebuild the heavy-atom graph from the emitted string and check
//! it against the molecule: same atoms (via `atom_order`), same bonds with the same orders, same
//! fragment count, same total charge, and — with the OpenSMILES valence table re-derived here
//! independently — the same total hydrogen count.
//!
//! The scanner in `parse` is **not** a SMILES reader: it handles only the narrow token subset
//! this writer emits (no lowercase aromatics, no stereo, no isotopes). A real reader is separate,
//! deferred work.

use std::collections::HashMap;

use molar::prelude::*;

const REFS: &str = "tests/data/gaff_ref/references.json";
const SDF_DIR: &str = "tests/data/gaff_ref/sdf";

#[derive(serde::Deserialize)]
struct RefFile {
    molecules: Vec<RefMol>,
}

#[derive(serde::Deserialize)]
struct RefMol {
    name: String,
}

// ---------------------------------------------------------------------------
// A minimal scanner for the subset molar's writer emits
// ---------------------------------------------------------------------------

struct PAtom {
    z: u8,
    /// Hydrogen count stated in a bracket; `None` when the token is bare and the count is implied.
    h: Option<u8>,
    fc: i32,
}

struct Parsed {
    atoms: Vec<PAtom>,
    /// `(atom token index, atom token index, order)`.
    bonds: Vec<(usize, usize, u8)>,
    fragments: usize,
}

/// The organic-subset default valences, written out again here on purpose: re-deriving them
/// independently of `molar/src/smiles/valence.rs` is what makes the hydrogen check meaningful
/// rather than circular.
fn implied_h(z: u8, bond_order_sum: u32) -> u8 {
    let valences: &[u8] = match z {
        5 => &[3],
        6 => &[4],
        7 => &[3, 5],
        8 => &[2],
        15 => &[3, 5],
        16 => &[2, 4, 6],
        9 | 17 | 35 | 53 => &[1],
        _ => return 0, // not in the subset: such an atom is always bracketed
    };
    valences
        .iter()
        .copied()
        .find(|&v| u32::from(v) >= bond_order_sum)
        .map_or(0, |v| (u32::from(v) - bond_order_sum) as u8)
}

fn symbol_to_z(sym: &str) -> Option<u8> {
    Some(match sym {
        "H" => 1,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Na" => 11,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Br" => 35,
        "I" => 53,
        _ => return None,
    })
}

fn parse(s: &str) -> Result<Parsed, String> {
    let b = s.as_bytes();
    let mut i = 0;
    let mut atoms: Vec<PAtom> = Vec::new();
    let mut bonds: Vec<(usize, usize, u8)> = Vec::new();
    let mut fragments = 1;

    let mut prev: Option<usize> = None;
    let mut branch: Vec<usize> = Vec::new();
    let mut pending: Option<u8> = None; // bond order stated before the next atom / ring digit
    let mut open: HashMap<u16, (usize, u8)> = HashMap::new();

    while i < b.len() {
        let c = b[i] as char;
        match c {
            '(' => {
                branch.push(prev.ok_or("branch opened with no preceding atom")?);
                i += 1;
            }
            ')' => {
                prev = Some(branch.pop().ok_or("unbalanced ')'")?);
                i += 1;
            }
            '=' => {
                pending = Some(2);
                i += 1;
            }
            '#' => {
                pending = Some(3);
                i += 1;
            }
            '.' => {
                prev = None;
                fragments += 1;
                i += 1;
            }
            '%' | '0'..='9' => {
                let (num, adv) = if c == '%' {
                    let t = s.get(i + 1..i + 3).ok_or("truncated %nn ring number")?;
                    (t.parse::<u16>().map_err(|e| e.to_string())?, 3)
                } else {
                    (u16::from(b[i] - b'0'), 1)
                };
                let a = prev.ok_or("ring closure with no preceding atom")?;
                match open.remove(&num) {
                    // Second endpoint: the order was stated at the opening.
                    Some((other, order)) => bonds.push((other, a, order)),
                    None => {
                        open.insert(num, (a, pending.unwrap_or(1)));
                    }
                }
                pending = None;
                i += adv;
            }
            '[' => {
                let end = s[i..].find(']').ok_or("unterminated bracket atom")? + i;
                let inner = &s[i + 1..end];
                let (atom, consumed) = parse_bracket(inner)?;
                if consumed != inner.len() {
                    return Err(format!("unparsed bracket remainder in [{inner}]"));
                }
                link(&mut atoms, &mut bonds, &mut prev, &mut pending, atom);
                i = end + 1;
            }
            _ => {
                // A bare organic-subset symbol. Try two letters first so `Cl`/`Br` win over
                // `C`/`B`; no concatenation of bare tokens can spell another proper-case symbol,
                // so this cannot mis-split (`CC` -> C, C; `CI` -> C, I).
                let two = s.get(i..i + 2).and_then(symbol_to_z).map(|z| (z, 2));
                let one = s.get(i..i + 1).and_then(symbol_to_z).map(|z| (z, 1));
                let (z, adv) = two
                    .or(one)
                    .ok_or_else(|| format!("unrecognised token at byte {i} of {s:?}"))?;
                link(&mut atoms, &mut bonds, &mut prev, &mut pending, PAtom { z, h: None, fc: 0 });
                i += adv;
            }
        }
    }

    if !branch.is_empty() {
        return Err("unbalanced '('".into());
    }
    if !open.is_empty() {
        return Err(format!("{} ring closure(s) never closed", open.len()));
    }
    Ok(Parsed { atoms, bonds, fragments })
}

/// Append an atom and, unless it starts a fragment, the bond joining it to the previous atom.
fn link(
    atoms: &mut Vec<PAtom>,
    bonds: &mut Vec<(usize, usize, u8)>,
    prev: &mut Option<usize>,
    pending: &mut Option<u8>,
    atom: PAtom,
) {
    let idx = atoms.len();
    atoms.push(atom);
    if let Some(p) = *prev {
        bonds.push((p, idx, pending.unwrap_or(1)));
    }
    *pending = None;
    *prev = Some(idx);
}

/// `<symbol>[H[n]][+|-[n]]` — the only bracket forms this writer produces.
fn parse_bracket(inner: &str) -> Result<(PAtom, usize), String> {
    let b = inner.as_bytes();
    let mut i = 0;
    if b.is_empty() {
        return Err("empty bracket atom".into());
    }

    // Element symbol: an uppercase letter, optionally followed by a lowercase one.
    let mut end = i + 1;
    if b.get(end).is_some_and(|c| c.is_ascii_lowercase()) {
        end += 1;
    }
    let sym = &inner[i..end];
    let z = symbol_to_z(sym).ok_or_else(|| format!("unknown bracket element {sym:?}"))?;
    i = end;

    // Hydrogen count. Bare `H` means one; `H` plus digits means that many.
    let mut h = 0u8;
    if b.get(i) == Some(&b'H') {
        i += 1;
        let ds = i;
        while b.get(i).is_some_and(|c| c.is_ascii_digit()) {
            i += 1;
        }
        h = if i == ds { 1 } else { inner[ds..i].parse().map_err(|_| "bad H count")? };
    }

    // Charge.
    let mut fc = 0i32;
    if let Some(&sign) = b.get(i) {
        if sign == b'+' || sign == b'-' {
            i += 1;
            let ds = i;
            while b.get(i).is_some_and(|c| c.is_ascii_digit()) {
                i += 1;
            }
            let mag: i32 = if i == ds { 1 } else { inner[ds..i].parse().map_err(|_| "bad charge")? };
            fc = if sign == b'+' { mag } else { -mag };
        }
    }

    Ok((PAtom { z, h: Some(h), fc }, i))
}

// ---------------------------------------------------------------------------
// Checks
// ---------------------------------------------------------------------------

/// Compare the emitted string's graph against the molecule it came from.
fn check_writer(sys: &System) -> Result<(), String> {
    let smi = sys.to_smiles().map_err(|e| e.to_string())?;
    let parsed = parse(smi.as_str())?;
    let order = smi.atom_order();

    // -- atoms: one token per written atom, same element, in the mapping's order
    if parsed.atoms.len() != order.len() {
        return Err(format!(
            "{} atom tokens but atom_order has {}",
            parsed.atoms.len(),
            order.len()
        ));
    }
    let z_of: Vec<u8> = sys.iter_atoms().map(|a| a.get_atomic_number()).collect();
    let fc_of: Vec<i32> = sys.iter_atoms().map(|a| a.get_formal_charge().unwrap_or(0)).collect();
    for (k, pa) in parsed.atoms.iter().enumerate() {
        if pa.z != z_of[order[k]] {
            return Err(format!(
                "token {k} is Z={} but atom {} is Z={}",
                pa.z, order[k], z_of[order[k]]
            ));
        }
    }

    // -- heavy atoms: none may be dropped or invented
    let heavy_mol = z_of.iter().filter(|&&z| z != 1).count();
    let heavy_str = parsed.atoms.iter().filter(|a| a.z != 1).count();
    if heavy_mol != heavy_str {
        return Err(format!("{heavy_str} heavy atom tokens for {heavy_mol} heavy atoms"));
    }

    // -- bonds: every written bond is a real bond of the same order, and none is missing
    let written: Vec<bool> = {
        let mut w = vec![false; z_of.len()];
        for &g in order {
            w[g] = true;
        }
        w
    };
    let mut real: HashMap<(usize, usize), u8> = HashMap::new();
    for bond in sys.iter_bonds() {
        let [i, j] = bond.pair();
        if !written[i] || !written[j] {
            continue; // a folded hydrogen
        }
        let code = match bond.order() {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => 1,
        };
        real.insert((i.min(j), i.max(j)), code);
    }
    if parsed.bonds.len() != real.len() {
        return Err(format!("{} bonds written for {} real bonds", parsed.bonds.len(), real.len()));
    }
    for &(a, b2, code) in &parsed.bonds {
        let (g1, g2) = (order[a], order[b2]);
        let key = (g1.min(g2), g1.max(g2));
        match real.get(&key) {
            None => return Err(format!("written bond {g1}-{g2} does not exist")),
            Some(&want) if want != code => {
                return Err(format!("bond {g1}-{g2} written as order {code}, real order {want}"))
            }
            Some(_) => {}
        }
    }

    // -- fragments
    let n = z_of.len();
    let mut top_bonds = BondStorage::default();
    for bond in sys.iter_bonds() {
        let [i, j] = bond.pair();
        top_bonds.push(&Bond::new(i, j));
    }
    let adj = BondAdjacency::build(n, top_bonds.iter_pairs());
    let frags = connected_components(&adj).iter().max().map_or(0, |m| m + 1);
    if parsed.fragments != frags {
        return Err(format!("{} fragments written, molecule has {frags}", parsed.fragments));
    }

    // -- total charge
    let charge_str: i32 = parsed.atoms.iter().map(|a| a.fc).sum();
    let charge_mol: i32 = fc_of.iter().sum();
    if charge_str != charge_mol {
        return Err(format!("written charge {charge_str}, molecule charge {charge_mol}"));
    }

    // -- hydrogens: nothing lost in the folding
    let mut sum: u32 = 0;
    let mut bond_sum = vec![0u32; parsed.atoms.len()];
    for &(a, b2, code) in &parsed.bonds {
        bond_sum[a] += u32::from(code);
        bond_sum[b2] += u32::from(code);
    }
    for (k, pa) in parsed.atoms.iter().enumerate() {
        sum += u32::from(pa.h.unwrap_or_else(|| implied_h(pa.z, bond_sum[k])));
        if pa.z == 1 {
            sum += 1; // a hydrogen written as its own atom
        }
    }
    let h_mol = z_of.iter().filter(|&&z| z == 1).count() as u32;
    if sum != h_mol {
        return Err(format!("string accounts for {sum} hydrogens, molecule has {h_mol}"));
    }

    Ok(())
}

/// Total bond order incident on each atom — what a Kekulé structure must preserve.
fn valences(sys: &System, orders: &[BondOrder]) -> Vec<u32> {
    let mut v = vec![0u32; sys.len()];
    for (bi, bond) in sys.iter_bonds().enumerate() {
        let w = match orders[bi] {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 0, // must not survive
            _ => 1,
        };
        let [i, j] = bond.pair();
        v[i] += w;
        v[j] += w;
    }
    v
}

/// Aromatize with `perceive`, then require `kekulize` to recover a valid Kekulé form.
fn check_kekulize(sys: &mut System) -> Result<bool, String> {
    let before_orders: Vec<BondOrder> = sys.iter_bonds().map(|b| b.order()).collect();
    let before = valences(sys, &before_orders);

    sys.perceive();
    let aromatized: Vec<BondOrder> = sys.iter_bonds().map(|b| b.order()).collect();
    if !aromatized.contains(&BondOrder::Aromatic) {
        return Ok(false); // nothing aromatic in this molecule; not a round-trip case
    }

    let n = sys.len();
    let z: Vec<u8> = sys.iter_atoms().map(|a| a.get_atomic_number()).collect();
    let fc: Vec<i32> = sys.iter_atoms().map(|a| a.get_formal_charge().unwrap_or(0)).collect();
    let mut pairs = BondStorage::default();
    for bond in sys.iter_bonds() {
        let [i, j] = bond.pair();
        pairs.push(&Bond::new(i, j));
    }
    let adj = BondAdjacency::build(n, pairs.iter_pairs());

    let recovered = kekulize(&z, &fc, &aromatized, &adj).map_err(|e| e.to_string())?;
    if recovered.contains(&BondOrder::Aromatic) {
        return Err("an aromatic bond survived kekulization".into());
    }
    let after = valences(sys, &recovered);
    if after != before {
        let bad = (0..n).find(|&i| after[i] != before[i]).unwrap();
        return Err(format!(
            "atom {bad} (Z={}) valence {} before, {} after",
            z[bad], before[bad], after[bad]
        ));
    }
    Ok(true)
}

#[derive(Default)]
struct Stats {
    mols: usize,
    load_err: usize,
    writer_ok: usize,
    writer_fail: Vec<(String, String)>,
    kek_cases: usize,
    kek_ok: usize,
    kek_fail: Vec<(String, String)>,
}

fn run() -> Stats {
    let refs: RefFile = {
        let s = std::fs::read_to_string(REFS).unwrap_or_else(|e| panic!("cannot read {REFS}: {e}"));
        serde_json::from_str(&s).unwrap_or_else(|e| panic!("cannot parse {REFS}: {e}"))
    };
    let mut st = Stats::default();

    for mol in &refs.molecules {
        st.mols += 1;
        let path = format!("{SDF_DIR}/{}.sdf", mol.name);
        let Ok(mut sys) = System::from_file(&path) else {
            st.load_err += 1;
            continue;
        };

        // Writer first, on the untouched Kekulé structure from the file.
        match check_writer(&sys) {
            Ok(()) => st.writer_ok += 1,
            Err(e) => st.writer_fail.push((mol.name.clone(), e)),
        }

        // Then the kekulizer round trip, which mutates the system.
        match check_kekulize(&mut sys) {
            Ok(true) => {
                st.kek_cases += 1;
                st.kek_ok += 1;
            }
            Ok(false) => {}
            Err(e) => {
                st.kek_cases += 1;
                st.kek_fail.push((mol.name.clone(), e));
            }
        }
    }
    st
}

fn report(st: &Stats) {
    println!("\n=== SMILES parity over {} corpus molecules ===", st.mols);
    println!("  load failures         : {}", st.load_err);
    println!(
        "  writer graph fidelity : {}/{} ({:.2}%)",
        st.writer_ok,
        st.mols - st.load_err,
        100.0 * st.writer_ok as f64 / (st.mols - st.load_err).max(1) as f64
    );
    println!(
        "  kekulize round trip   : {}/{} aromatic molecules ({:.2}%)",
        st.kek_ok,
        st.kek_cases,
        100.0 * st.kek_ok as f64 / st.kek_cases.max(1) as f64
    );
    for (label, fails) in [("writer", &st.writer_fail), ("kekulize", &st.kek_fail)] {
        for (name, why) in fails.iter().take(15) {
            println!("  [{label}] {name}: {why}");
        }
        if fails.len() > 15 {
            println!("  [{label}] ... and {} more", fails.len() - 15);
        }
    }
    println!();
}

/// Development driver: always passes, prints the numbers. `-- --nocapture`.
#[test]
fn smiles_parity_report() {
    report(&run());
}

/// Dump every corpus molecule's SMILES for the offline RDKit cross-check, which is the real
/// correctness gate — the in-tree checks above verify the string encodes molar's own graph, but
/// only RDKit can confirm it means the same molecule to another toolkit.
///
/// ```sh
/// cargo test -p molar_ff --test smiles_parity -- --ignored dump
/// python tests/data/gaff_ref/check_smiles.py
/// ```
#[test]
#[ignore = "development harness: writes a dump for check_smiles.py"]
fn smiles_dump_for_rdkit() {
    let refs: RefFile =
        serde_json::from_str(&std::fs::read_to_string(REFS).unwrap()).unwrap();
    let mut out = Vec::new();
    for mol in &refs.molecules {
        let Ok(sys) = System::from_file(format!("{SDF_DIR}/{}.sdf", mol.name)) else { continue };
        let smi = match sys.to_smiles() {
            Ok(s) => s.into_string(),
            Err(e) => format!("!ERROR: {e}"),
        };
        out.push(format!("  {{\"name\": {:?}, \"smiles\": {:?}}}", mol.name, smi));
    }
    let path = "tests/data/gaff_ref/smiles_molar.json";
    std::fs::write(path, format!("[\n{}\n]\n", out.join(",\n"))).unwrap();
    println!("wrote {} molecules to {path}", out.len());
}

/// The gate: every corpus molecule must write a string that encodes its graph exactly, and every
/// aromatic one must survive the aromatize/kekulize round trip.
#[test]
fn smiles_parity_threshold() {
    let st = run();
    report(&st);
    assert_eq!(st.load_err, 0, "every corpus SDF must load");
    assert!(st.writer_fail.is_empty(), "{} writer failures", st.writer_fail.len());
    assert!(st.kek_fail.is_empty(), "{} kekulize failures", st.kek_fail.len());
    assert!(st.kek_cases > 100, "corpus should exercise many aromatic molecules");
}
