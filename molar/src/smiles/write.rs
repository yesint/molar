//! SMILES emission: spanning-tree traversal, ring-closure bookkeeping and token formatting.
//!
//! The molecule arrives here already reduced to a local `0..n` subgraph with Kekulé bond orders
//! and hydrogens folded into per-atom counts (see the parent module). This file is only about
//! turning that into a string.
//!
//! Both passes over the graph consume neighbours in the order given by a single **rank** array,
//! which is the one place traversal order is decided. v1 ranks by local index (i.e. ascending
//! global atom index), which makes the output deterministic; canonical output later replaces the
//! ranks and nothing else here changes.

use super::{MolView, SmilesError};
use crate::bond::BondOrder;
use crate::periodic_table::ELEMENT_NAME;
use crate::smiles::valence::implied_h;

/// Highest ring-closure number the notation can express (`%99`).
const MAX_RING_NUM: u16 = 99;

/// A pending emission step. The traversal is an explicit stack rather than recursion because a
/// selection can be a whole protein — tens of thousands of atoms deep — which would overflow the
/// stack on a recursive walk.
enum Work {
    /// Emit atom `a`, arrived over bond `via`, then queue its children.
    Visit { a: usize, via: Option<usize> },
    /// Emit a literal — a branch parenthesis or a fragment separator.
    Text(&'static str),
}

/// Emit the whole molecule: every bonded fragment, joined with `.`.
///
/// Returns the string together with the global atom index of each atom written, in written
/// order (suppressed hydrogens excluded).
pub(super) fn write(mol: &MolView<'_>) -> Result<(String, Vec<usize>), SmilesError> {
    let n = mol.adj.n_atoms();
    let n_bonds = mol.adj.n_bonds();

    // Neighbour lists in rank order, computed once and shared by both passes so they agree on
    // the tree they are walking. `sort_by_key` is stable, so equal ranks keep the adjacency's
    // ascending-bond-index order.
    let mut nbrs: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (a, list) in nbrs.iter_mut().enumerate() {
        if mol.suppressed[a] {
            continue;
        }
        list.extend(
            mol.adj
                .neighbors(a)
                .iter()
                .filter(|nb| !mol.suppressed[nb.atom()])
                .map(|nb| (nb.atom(), nb.bond())),
        );
        list.sort_by_key(|&(nb, _)| mol.rank[nb]);
    }

    // Fragments, ordered by their lowest-ranked written atom.
    let labels = crate::perception::connected_components(mol.adj);
    let n_frags = labels.iter().max().map_or(0, |m| m + 1);
    let mut starts: Vec<(u32, usize)> = Vec::new(); // (rank of start atom, start atom)
    for f in 0..n_frags {
        let start = (0..n)
            .filter(|&a| labels[a] == f && !mol.suppressed[a])
            .min_by_key(|&a| mol.rank[a]);
        // A fragment with nothing to write cannot happen: a hydrogen is only suppressed into a
        // heavy neighbour, which is in the same fragment and is itself written.
        if let Some(start) = start {
            starts.push((mol.rank[start], start));
        }
    }
    starts.sort_unstable();

    let mut out = String::new();
    let mut atom_order = Vec::new();
    let mut emitted = vec![false; n];
    let mut ring_open: Vec<Option<u16>> = vec![None; n_bonds];
    let mut num_in_use = [false; MAX_RING_NUM as usize + 1];

    for (i, &(_, start)) in starts.iter().enumerate() {
        if i > 0 {
            out.push('.');
        }
        let (is_tree, is_ring) = spanning_tree(&nbrs, start, n, n_bonds);
        emit_fragment(
            mol,
            &nbrs,
            start,
            &is_tree,
            &is_ring,
            &mut emitted,
            &mut ring_open,
            &mut num_in_use,
            &mut out,
            &mut atom_order,
        )?;
    }

    Ok((out, atom_order))
}

/// Classify the fragment's bonds into a spanning tree plus the ring-closure bonds left over.
///
/// This has to be a separate pass: a ring-closure number appears at **both** endpoints, and the
/// first endpoint is emitted long before the traversal discovers the bond closes a ring.
///
/// The tree need not be a depth-first tree for the output to be valid — SMILES ring-closure
/// bonds may join any two atoms, not only an ancestor and a descendant — but a true DFS is used
/// anyway because it yields the conventional, more readable nesting.
fn spanning_tree(
    nbrs: &[Vec<(usize, usize)>],
    start: usize,
    n: usize,
    n_bonds: usize,
) -> (Vec<bool>, Vec<bool>) {
    let mut is_tree = vec![false; n_bonds];
    let mut is_ring = vec![false; n_bonds];
    let mut visited = vec![false; n];

    // (atom, index of the next neighbour to consider) — an explicit DFS.
    let mut stack = vec![(start, 0usize)];
    visited[start] = true;
    while let Some(&(a, ni)) = stack.last() {
        if ni >= nbrs[a].len() {
            stack.pop();
            continue;
        }
        stack.last_mut().unwrap().1 += 1;

        let (nb, bond) = nbrs[a][ni];
        if is_tree[bond] || is_ring[bond] {
            continue; // already classified, including the bond back to the parent
        }
        if visited[nb] {
            is_ring[bond] = true;
        } else {
            visited[nb] = true;
            is_tree[bond] = true;
            stack.push((nb, 0));
        }
    }
    (is_tree, is_ring)
}

/// Walk one fragment's spanning tree, appending to `out`.
#[allow(clippy::too_many_arguments)]
fn emit_fragment(
    mol: &MolView<'_>,
    nbrs: &[Vec<(usize, usize)>],
    start: usize,
    is_tree: &[bool],
    is_ring: &[bool],
    emitted: &mut [bool],
    ring_open: &mut [Option<u16>],
    num_in_use: &mut [bool; MAX_RING_NUM as usize + 1],
    out: &mut String,
    atom_order: &mut Vec<usize>,
) -> Result<(), SmilesError> {
    let mut work = vec![Work::Visit { a: start, via: None }];

    while let Some(w) = work.pop() {
        let Work::Visit { a, via } = w else {
            if let Work::Text(t) = w {
                out.push_str(t);
            }
            continue;
        };

        if let Some(b) = via {
            out.push_str(bond_symbol(mol.orders[b]));
        }
        out.push_str(&atom_token(mol, a)?);
        atom_order.push(mol.global[a]);
        emitted[a] = true;

        // Ring-closure numbers come directly after the atom token, before any branch.
        for &(_, bond) in &nbrs[a] {
            if !is_ring[bond] {
                continue;
            }
            match ring_open[bond].take() {
                // Second endpoint: repeat the number and hand it back for reuse. The bond order
                // was already stated at the opening, so it is not repeated here.
                Some(num) => {
                    push_ring_num(out, num);
                    num_in_use[num as usize] = false;
                }
                // First endpoint: allocate, and state the bond order here.
                None => {
                    let num = (1..=MAX_RING_NUM)
                        .find(|&k| !num_in_use[k as usize])
                        .ok_or(SmilesError::TooManyRingClosures)?;
                    num_in_use[num as usize] = true;
                    ring_open[bond] = Some(num);
                    out.push_str(bond_symbol(mol.orders[bond]));
                    push_ring_num(out, num);
                }
            }
        }

        // Tree children, minus the parent (already emitted). All but the last are parenthesised;
        // the last continues the main chain.
        let children: Vec<(usize, usize)> = nbrs[a]
            .iter()
            .copied()
            .filter(|&(nb, bond)| is_tree[bond] && !emitted[nb])
            .collect();
        if let Some(&(last, last_bond)) = children.last() {
            // Pushed in reverse: the stack is LIFO, so the first child pops first.
            work.push(Work::Visit { a: last, via: Some(last_bond) });
            for &(nb, bond) in children[..children.len() - 1].iter().rev() {
                work.push(Work::Text(")"));
                work.push(Work::Visit { a: nb, via: Some(bond) });
                work.push(Work::Text("("));
            }
        }
    }
    Ok(())
}

/// `1`–`9` are written bare; `10`–`99` need the `%` prefix.
fn push_ring_num(out: &mut String, num: u16) {
    if num < 10 {
        out.push((b'0' + num as u8) as char);
    } else {
        out.push('%');
        out.push_str(&num.to_string());
    }
}

/// The bond symbol between two atoms. Single is the default bond in a Kekulé (all-uppercase)
/// SMILES and so is never written, which is why `-` never appears in this writer's output.
fn bond_symbol(o: BondOrder) -> &'static str {
    match o {
        BondOrder::Single => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        // Rejected before the writer runs; kekulization guarantees neither survives.
        BondOrder::Unspecified => unreachable!("unspecified orders are rejected up front"),
        BondOrder::Aromatic => unreachable!("kekulization resolves every aromatic bond"),
    }
}

/// One atom: the bare organic-subset symbol when its implied hydrogen count happens to be right
/// and it carries no charge, otherwise the explicit bracket form.
fn atom_token(mol: &MolView<'_>, a: usize) -> Result<String, SmilesError> {
    let z = mol.z[a];
    let symbol = ELEMENT_NAME
        .get(z as usize)
        .filter(|_| z != 0)
        .ok_or(SmilesError::UnknownElement { global: mol.global[a], z })?;

    let h = mol.h_count[a];
    let fc = mol.fc[a];

    // Summed order of the bonds actually written — suppressed hydrogens are counted in `h`.
    let bond_order_sum: u32 = mol
        .adj
        .neighbors(a)
        .iter()
        .filter(|nb| !mol.suppressed[nb.atom()])
        .map(|nb| match mol.orders[nb.bond()] {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => 1,
        })
        .sum();

    if fc == 0 && implied_h(z, bond_order_sum) == Some(h) {
        return Ok((*symbol).to_string());
    }

    // Bracket form: no implicit hydrogens, so the count is stated outright and is always exact.
    let mut s = String::with_capacity(8);
    s.push('[');
    s.push_str(symbol);
    match h {
        0 => {}
        1 => s.push('H'),
        k => {
            s.push('H');
            s.push_str(&k.to_string());
        }
    }
    match fc {
        0 => {}
        1 => s.push('+'),
        -1 => s.push('-'),
        k if k > 0 => {
            s.push('+');
            s.push_str(&k.to_string());
        }
        k => {
            s.push('-');
            s.push_str(&(-k).to_string());
        }
    }
    s.push(']');
    Ok(s)
}
