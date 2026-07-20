//! Post-typing conjugation parity split.
//!
//! After rule matching, atoms in conjugated systems are all given the "c/e/g" (or
//! nitrogen / phosphorus) base type. These routines 2-colour each conjugated system
//! along the Kekulé single/double bonds — the sign stays across single-like bonds and
//! flips across double/triple bonds — and rename the atoms that end up with sign −1
//! (`cc→cd`, `ce→cf`, `cg→ch`, `nc→nd`, `ne→nf`, `pc→pd`, `pe→pf`; and `cp→cq`).

use super::LocalBond;

/// Types that participate in the conjugated-system parity colouring.
fn is_atadjust(t: &str) -> bool {
    matches!(t, "cc" | "ce" | "cg" | "pc" | "pe" | "nc" | "ne")
}

/// Map a sign −1 base type to its partner.
fn flip_atadjust(t: &str) -> Option<&'static str> {
    Some(match t {
        "cc" => "cd",
        "ce" => "cf",
        "cg" => "ch",
        "pc" => "pd",
        "pe" => "pf",
        "nc" => "nd",
        "ne" => "nf",
        _ => return None,
    })
}

/// Run both conjugation adjustments in place, in order: the general c/e/g/n/p split
/// first, then the biphenyl `cp`/`cq` split.
pub(crate) fn adjust(types: &mut [String], bonds: &[LocalBond]) {
    atadjust(types, bonds);
    cpadjust(types, bonds);
}

fn atadjust(types: &mut [String], bonds: &[LocalBond]) {
    let n = types.len();
    let mut index1 = vec![0i32; n]; // colour (+1 / −1), 0 = uncoloured
    let mut index2 = vec![false; n]; // participates in the colouring

    let mut seeded = false;
    let mut num = 0usize;
    for i in 0..n {
        if is_atadjust(&types[i]) {
            index2[i] = true;
            if !seeded {
                index1[i] = 1;
                seeded = true;
            }
            num += 1;
        }
    }
    if num == 0 {
        return;
    }

    // `num - 1` propagation sweeps (one per participating atom, less one).
    for _ in 0..(num - 1) {
        let mut flag = false;
        for b in bonds {
            let (bi, bj) = (b.i, b.j);
            if !(index2[bi] && index2[bj]) {
                continue;
            }
            if !flag && index1[bi] == 0 && index1[bj] == 0 {
                index1[bi] = 1;
            }
            if index1[bi] == 0 && index1[bj] != 0 {
                flag = true;
                index1[bi] = if b.order == 1 { index1[bj] } else { -index1[bj] };
            }
            if index1[bj] == 0 && index1[bi] != 0 {
                flag = true;
                index1[bj] = if b.order == 1 { index1[bi] } else { -index1[bi] };
            }
        }
    }

    for i in 0..n {
        if index1[i] == -1 {
            if let Some(t) = flip_atadjust(&types[i]) {
                types[i] = t.to_string();
            }
        }
    }
}

fn cpadjust(types: &mut [String], bonds: &[LocalBond]) {
    let n = types.len();
    let mut index1 = vec![0i32; n];
    let mut index2 = vec![false; n];

    let mut seeded = false;
    let mut num = 0usize;
    for i in 0..n {
        if types[i] == "cp" {
            index2[i] = true;
            if !seeded {
                index1[i] = 1;
                seeded = true;
            }
            num += 1;
        }
    }
    if num == 0 {
        return;
    }

    for _ in 0..(num - 1) {
        for b in bonds {
            let (bi, bj) = (b.i, b.j);
            if !(index2[bi] && index2[bj]) {
                continue;
            }
            if index1[bi] == 0 && index1[bj] != 0 {
                index1[bi] = if b.order == 1 { index1[bj] } else { -index1[bj] };
            }
            if index1[bj] == 0 && index1[bi] != 0 {
                index1[bj] = if b.order == 1 { index1[bi] } else { -index1[bi] };
            }
        }
    }

    for i in 0..n {
        if index1[i] == -1 && types[i] == "cp" {
            types[i] = "cq".to_string();
        }
    }
}
