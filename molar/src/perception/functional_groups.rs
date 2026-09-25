//! Functional-group constraints (delivery step 5 of the bond-perception plan).
//!
//! The general valence search finds a chemically valid assignment, but for a few groups it can
//! land on an equivalent resonance form that is not the conventional one (azide is the clearest
//! case). This pass recognizes such groups and pins their canonical bond orders. The pinned
//! orders enter [`assign_bond_orders`](super::assign_bond_orders) as fixed bonds; formal charges
//! then follow from the valence model, so a template only needs to fix **orders**, never
//! charges.
//!
//! Every pattern here is matched on connectivity that does not depend on explicit hydrogen
//! (terminal oxygens, a terminal-nitrogen chain), so it applies under both hydrogen policies.
//! Choosing *which* terminal atom carries the double bond can depend on hydrogen, though: with
//! hydrogens removed, a carboxylic acid's O-H oxygen is terminal too. Where that matters the
//! choice follows bond length when coordinates are given. Groups the
//! general search already assigns canonically (nitrile, amide, aromatic rings) are left to it.
//! The set is deliberately small and extensible; residue templates (step 6) will build on the
//! same "fix orders, let charges follow" mechanism.

use crate::prelude::*;

/// Canonical bond orders for the recognized functional groups, indexed by bond. `None` means no
/// template applies to that bond.
///
/// `coords` (one position per atom, or `None`) decides between otherwise equivalent terminal
/// oxygens by bond length; without it the lowest index wins.
pub(super) fn functional_group_orders(
    z: &[u8],
    adj: &BondAdjacency,
    coords: Option<&[Pos]>,
) -> Vec<Option<BondOrder>> {
    let mut orders: Vec<Option<BondOrder>> = vec![None; adj.n_bonds()];
    for a in 0..adj.n_atoms() {
        match z[a] {
            7 => {
                nitro(a, z, adj, &mut orders);
                n_oxide(a, z, adj, &mut orders);
                azide(a, z, adj, &mut orders);
            }
            6 => carboxyl(a, z, adj, coords, &mut orders),
            16 => sulfoxide(a, z, adj, &mut orders),
            _ => {}
        }
    }
    orders
}

fn degree(adj: &BondAdjacency, a: usize) -> usize {
    adj.neighbors(a).len()
}

/// Claim a bond for a template order. The first template to reach a bond wins; because atoms are
/// scanned in ascending order the choice is deterministic.
fn claim(orders: &mut [Option<BondOrder>], bond: usize, o: BondOrder) {
    if orders[bond].is_none() {
        orders[bond] = Some(o);
    }
}

/// R-NO2: a nitrogen with three neighbors, two of them terminal oxygens. Pin one N=O and one
/// N-O; the valence model then makes the nitrogen +1 and the single-bonded oxygen -1.
fn nitro(a: usize, z: &[u8], adj: &BondAdjacency, orders: &mut [Option<BondOrder>]) {
    if degree(adj, a) != 3 {
        return;
    }
    let mut terminal_o: Vec<(usize, usize)> = Vec::new();
    let mut others = 0;
    for nb in adj.neighbors(a) {
        if z[nb.atom()] == 8 && degree(adj, nb.atom()) == 1 {
            terminal_o.push((nb.atom(), nb.bond()));
        } else {
            others += 1;
        }
    }
    if terminal_o.len() == 2 && others == 1 {
        terminal_o.sort_by_key(|&(atom, _)| atom);
        claim(orders, terminal_o[0].1, BondOrder::Double);
        claim(orders, terminal_o[1].1, BondOrder::Single);
    }
}

/// R3N+-O-: a nitrogen with four neighbors and a terminal oxygen (an amine oxide, or the
/// N-oxide of a pyridine). Pin the N-O bond single so the valence model makes the nitrogen +1
/// and the oxygen -1, rather than a neutral N=O the four-coordinate nitrogen cannot support.
fn n_oxide(a: usize, z: &[u8], adj: &BondAdjacency, orders: &mut [Option<BondOrder>]) {
    if degree(adj, a) != 4 {
        return;
    }
    for nb in adj.neighbors(a) {
        if z[nb.atom()] == 8 && degree(adj, nb.atom()) == 1 {
            claim(orders, nb.bond(), BondOrder::Single);
        }
    }
}

/// R-N=N+=N-: the middle nitrogen of a three-nitrogen chain whose far end is terminal. Pin both
/// N-N bonds double; the valence model then makes the ends 0 / -1 and the middle +1.
fn azide(a: usize, z: &[u8], adj: &BondAdjacency, orders: &mut [Option<BondOrder>]) {
    let neighbors = adj.neighbors(a);
    if neighbors.len() != 2 || !neighbors.iter().all(|nb| z[nb.atom()] == 7) {
        return;
    }
    if neighbors.iter().any(|nb| degree(adj, nb.atom()) == 1) {
        for nb in neighbors {
            claim(orders, nb.bond(), BondOrder::Double);
        }
    }
}

/// A sulfoxide/sulfinyl sulfur: three neighbors, exactly one a terminal oxygen (the other two
/// carbon or the like). Pin the S-O bond single, giving the charge-separated S+-O- form the
/// reference structures use for a sulfoxide (as opposed to a sulfone, whose four-coordinate
/// sulfur keeps its S=O double bonds). The valence model then makes the sulfur +1 and the
/// oxygen -1.
fn sulfoxide(a: usize, z: &[u8], adj: &BondAdjacency, orders: &mut [Option<BondOrder>]) {
    if degree(adj, a) != 3 {
        return;
    }
    let terminal_o: Vec<usize> = adj
        .neighbors(a)
        .iter()
        .filter(|nb| z[nb.atom()] == 8 && degree(adj, nb.atom()) == 1)
        .map(|nb| nb.bond())
        .collect();
    if terminal_o.len() == 1 {
        claim(orders, terminal_o[0], BondOrder::Single);
    }
}

/// A carboxyl/carboxylate/ester carbon: bonded to exactly two oxygens, at least one terminal.
/// Pin a terminal oxygen as the C=O and the other as a single bond; the charge (carboxylate -1,
/// or neutral for an acid/ester whose second oxygen keeps a substituent or hydrogen) then
/// follows from the valence model. Between two terminal oxygens the shorter C-O bond is the
/// double one when `coords` is given — with hydrogens stripped that separates C=O (~0.12 nm)
/// from C-OH (~0.135 nm) — and the lowest index otherwise.
fn carboxyl(
    a: usize,
    z: &[u8],
    adj: &BondAdjacency,
    coords: Option<&[Pos]>,
    orders: &mut [Option<BondOrder>],
) {
    let oxygens: Vec<(usize, usize, usize)> = adj
        .neighbors(a)
        .iter()
        .filter(|nb| z[nb.atom()] == 8)
        .map(|nb| (nb.atom(), nb.bond(), degree(adj, nb.atom())))
        .collect();
    // Exactly two oxygens and at least one non-oxygen neighbor (the R group). This excludes
    // carbon dioxide (O=C=O, two oxygens and nothing else), which is not a carboxyl.
    if oxygens.len() != 2 || !adj.neighbors(a).iter().any(|nb| z[nb.atom()] != 8) {
        return;
    }
    let mut terminal: Vec<&(usize, usize, usize)> =
        oxygens.iter().filter(|&&(_, _, d)| d == 1).collect();
    if terminal.is_empty() {
        return;
    }
    terminal.sort_by_key(|&&(atom, _, _)| atom);
    if let Some(c) = coords {
        // Stable sort: equal lengths keep the lowest-index order.
        terminal.sort_by(|&&(p, _, _), &&(q, _, _)| {
            (c[p] - c[a]).norm().total_cmp(&(c[q] - c[a]).norm())
        });
    }
    let double_bond = terminal[0].1;
    claim(orders, double_bond, BondOrder::Double);
    for &(_, bond, _) in &oxygens {
        if bond != double_bond {
            claim(orders, bond, BondOrder::Single);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn adjacency(n: usize, bonds: &[[usize; 2]]) -> BondAdjacency {
        BondAdjacency::build(n, bonds.iter().copied())
    }

    #[test]
    fn nitro_group_is_pinned() {
        // C0-N1(-O2)(-O3): nitromethane's heavy atoms (methyl carbon standing in for R).
        let z = [6u8, 7, 8, 8];
        let bonds = [[0, 1], [1, 2], [1, 3]];
        let adj = adjacency(4, &bonds);
        let o = functional_group_orders(&z, &adj, None);
        assert_eq!(o[0], None, "C-N left to the solver");
        assert_eq!(o[1], Some(BondOrder::Double), "N=O on the lower-index oxygen");
        assert_eq!(o[2], Some(BondOrder::Single), "N-O on the other oxygen");
    }

    #[test]
    fn azide_chain_is_pinned_double_double() {
        // C0-N1=N2=N3
        let z = [6u8, 7, 7, 7];
        let bonds = [[0, 1], [1, 2], [2, 3]];
        let adj = adjacency(4, &bonds);
        let o = functional_group_orders(&z, &adj, None);
        assert_eq!(o[0], None, "R-N stays single via the solver");
        assert_eq!(o[1], Some(BondOrder::Double));
        assert_eq!(o[2], Some(BondOrder::Double));
    }

    #[test]
    fn carboxyl_pins_one_double_one_single() {
        // C0(-O1)(-O2)-C3, both oxygens terminal (carboxylate).
        let z = [6u8, 8, 8, 6];
        let bonds = [[0, 1], [0, 2], [0, 3]];
        let adj = adjacency(4, &bonds);
        let o = functional_group_orders(&z, &adj, None);
        assert_eq!(o[0], Some(BondOrder::Double), "=O on the lower-index oxygen");
        assert_eq!(o[1], Some(BondOrder::Single));
        assert_eq!(o[2], None, "C-C left to the solver");
    }

    /// Hydrogen-stripped carboxylic acid: both oxygens are terminal. The short bond is the C=O
    /// even when it is on the higher-index oxygen.
    #[test]
    fn carboxyl_double_bond_follows_geometry() {
        // C0(-O1 0.135 nm)(-O2 0.120 nm)-C3
        let z = [6u8, 8, 8, 6];
        let bonds = [[0, 1], [0, 2], [0, 3]];
        let adj = adjacency(4, &bonds);
        let coords = [
            Pos::new(0.0, 0.0, 0.0),
            Pos::new(0.135, 0.0, 0.0),
            Pos::new(-0.060, 0.104, 0.0),
            Pos::new(-0.075, -0.130, 0.0),
        ];
        let o = functional_group_orders(&z, &adj, Some(&coords));
        assert_eq!(o[0], Some(BondOrder::Single), "long C-O is the C-OH");
        assert_eq!(o[1], Some(BondOrder::Double), "short C-O is the C=O");
    }

    #[test]
    fn n_oxide_bond_is_pinned_single() {
        // (C0)(C1)(C2)N3-O4: a tertiary amine oxide. The N-O bond must be single.
        let z = [6u8, 6, 6, 7, 8];
        let bonds = [[0, 3], [1, 3], [2, 3], [3, 4]];
        let adj = adjacency(5, &bonds);
        let o = functional_group_orders(&z, &adj, None);
        assert_eq!(o[3], Some(BondOrder::Single), "N-O of the oxide is single");
    }

    #[test]
    fn sulfoxide_bond_is_pinned_single() {
        // (C0)(C1)S2-O3: a dialkyl sulfoxide. The S-O bond is single (S+, O-).
        let z = [6u8, 6, 16, 8];
        let bonds = [[0, 2], [1, 2], [2, 3]];
        let adj = adjacency(4, &bonds);
        let o = functional_group_orders(&z, &adj, None);
        assert_eq!(o[2], Some(BondOrder::Single), "S-O of the sulfoxide is single");
    }

    #[test]
    fn sulfone_sulfur_is_not_a_sulfoxide() {
        // (C0)(C1)S2(-O3)(-O4): a sulfone (four-coordinate S) must NOT be pinned by the
        // sulfoxide rule; its S=O double bonds are left to the general search.
        let z = [6u8, 6, 16, 8, 8];
        let bonds = [[0, 2], [1, 2], [2, 3], [2, 4]];
        let adj = adjacency(5, &bonds);
        let o = functional_group_orders(&z, &adj, None);
        assert_eq!(o[2], None);
        assert_eq!(o[3], None);
    }

    #[test]
    fn amine_nitrogen_is_not_mistaken_for_a_group() {
        // A plain amine N (C-N with two hydrogens) must match nothing.
        let z = [6u8, 7, 1, 1];
        let bonds = [[0, 1], [1, 2], [1, 3]];
        let adj = adjacency(4, &bonds);
        assert!(
            functional_group_orders(&z, &adj, None)
                .iter()
                .all(|o| o.is_none())
        );
    }
}
