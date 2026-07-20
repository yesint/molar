//! The GAFF rule-matching engine.
//!
//! It comprises an atomic-property matcher, a wildatom matcher, the recursive
//! chemical-environment engine (path enumeration plus a distinctness check across
//! branches), and the top-level per-atom typing loop. All rule data is fully structured
//! at compile time (see [`super::tables`]); nothing is parsed from strings here.
//!
//! Bond-code simplification (verified for Kekulé input): only bond orders 1/2/3 occur,
//! so a required-order test of ±1 ⇔ order 1, ±2 ⇔ order 2, ±3 ⇔ order 3; codes 4/9
//! (aromatic / delocalised) never match. `AB`/`DL` counts are always 0.

use std::collections::HashMap;

use crate::{FFError, FFType};

use super::props::Props;
use super::tables::{Atom, Bead, Pred, Prop, Rule, WILDATOMS};
use super::LocalBond;

/// A matched path: `id` = chain index it satisfies, `at` = the path atoms (excluding the
/// central atom).
struct Schain {
    id: usize,
    at: Vec<usize>,
}

/// Everything the matcher needs about the molecule, in the local index space.
pub(crate) struct Ctx<'a> {
    z: &'a [u8],
    con: &'a [Vec<usize>],
    props: &'a Props,
    rg: &'a [[u16; 11]],
    ar: &'a [[u16; 6]],
    nr: &'a [bool],
    ewd: &'a [i8],
    /// The rule set to match against (GAFF or GAFF2), in definition order.
    rules: &'static [Rule],
    /// Which force field is being applied (for error reporting).
    ff: FFType,
    /// Bond order between an unordered atom pair.
    bond_order: HashMap<(usize, usize), u8>,
}

impl<'a> Ctx<'a> {
    pub(crate) fn new(
        z: &'a [u8],
        con: &'a [Vec<usize>],
        bonds: &'a [LocalBond],
        props: &'a Props,
        rg: &'a [[u16; 11]],
        ar: &'a [[u16; 6]],
        nr: &'a [bool],
        ewd: &'a [i8],
        rules: &'static [Rule],
        ff: FFType,
    ) -> Self {
        let mut bond_order = HashMap::new();
        for b in bonds {
            bond_order.insert((b.i.min(b.j), b.i.max(b.j)), b.order);
        }
        Ctx { z, con, props, rg, ar, nr, ewd, rules, ff, bond_order }
    }

    /// True iff atoms `a` and `b` are joined by a bond of the given order.
    fn bond_is(&self, a: usize, b: usize, order: u8) -> bool {
        if order == 0 {
            return false;
        }
        self.bond_order
            .get(&(a.min(b), a.max(b)))
            .is_some_and(|&o| o == order)
    }

    // ---- atomic-property matching ----------------------------------------------

    /// Evaluate a structured atomic-property constraint (the `[...]` field) against atom
    /// `atmid`: an AND over groups, each an OR over its [`Pred`]s. `pre` is the immediate
    /// predecessor in a chem-env chain (`None` at top level).
    fn apcheck(&self, atmid: usize, pre: Option<usize>, prop: &[&[Pred]]) -> bool {
        prop.iter()
            .all(|group| group.iter().any(|pred| self.pred_ok(atmid, pre, pred)))
    }

    /// Evaluate one [`Pred`] against atom `atmid` (with predecessor `pre`).
    fn pred_ok(&self, atmid: usize, pre: Option<usize>, pred: &Pred) -> bool {
        // `> 0` when the predicate has no count, `== count` otherwise.
        let cnt_ok = |val: u16| pred.n.map_or(val > 0, |k| val == k as u16);
        // Quote flag on a bond token: `'` requires the bond to `pre` to be that order,
        // `''` requires it not to be; either fails when there is no predecessor.
        let bond_quote = |order: u8| match pred.q {
            0 => true,
            1 => pre.is_some_and(|p| self.bond_is(atmid, p, order)),
            2 => pre.is_some_and(|p| !self.bond_is(atmid, p, order)),
            _ => false,
        };
        match pred.p {
            Prop::Ring(k) => cnt_ok(self.rg[atmid][k as usize]),
            Prop::Arom(k) => cnt_ok(self.ar[atmid][k as usize]),
            Prop::NonRing => cnt_ok(self.nr[atmid] as u16),
            Prop::Sb => cnt_ok(self.props.sb[atmid]) && bond_quote(1),
            Prop::Db => cnt_ok(self.props.db[atmid]) && bond_quote(2),
            Prop::Tb => cnt_ok(self.props.tb[atmid]) && bond_quote(3),
            Prop::Dl => cnt_ok(0),
            Prop::Ab => cnt_ok(0),
        }
    }

    // ---- wildatom matching -----------------------------------------------------

    /// Does atom `a` satisfy wildatom index `w` (its `(z, connum)` in the wildatom's
    /// pair list; a pair `connum` of `0` matches any coordination)?
    fn wild_ok(&self, w: u8, a: usize) -> bool {
        for &(anum, cnum) in WILDATOMS[w as usize].pairs {
            if self.z[a] == anum && (cnum == 0 || self.props.connum[a] as u8 == cnum) {
                return true;
            }
        }
        false
    }

    // ---- chemical-environment engine ------------------------------------------

    /// DFS from the central atom, recording every path that satisfies some chain.
    #[allow(clippy::too_many_arguments)]
    fn cematch(
        &self,
        caid: usize,
        chains: &[&[Bead]],
        maxchain: usize,
        path: &mut Vec<usize>,
        startnum: usize,
        cesindex: &mut [u32],
        schains: &mut Vec<Schain>,
    ) {
        path.push(startnum);
        let selectnum = path.len(); // path length including the central atom
        for (k, ch) in chains.iter().enumerate() {
            if selectnum - 1 == ch.len() && self.match_chain(caid, path, ch) {
                cesindex[k] += 1;
                schains.push(Schain { id: k, at: path[1..].to_vec() });
            }
        }
        if selectnum <= maxchain {
            let deg = self.con[startnum].len().min(6);
            for idx in 0..deg {
                let nb = self.con[startnum][idx];
                if path.contains(&nb) {
                    continue;
                }
                self.cematch(caid, chains, maxchain, path, nb, cesindex, schains);
            }
        }
        path.pop();
    }

    /// Match the current DFS path (excluding the central atom) against one chain.
    fn match_chain(&self, caid: usize, path: &[usize], ch: &[Bead]) -> bool {
        for (b, bead) in ch.iter().enumerate() {
            let a = path[b + 1];
            // 1. required connum (`None` = any).
            if let Some(nn) = bead.n {
                if self.props.connum[a] != nn as usize {
                    return false;
                }
            }
            // 2. the required atom.
            match bead.atom {
                Atom::Z(z) => {
                    if self.z[a] != z {
                        return false;
                    }
                }
                Atom::Wild(w) => {
                    if !self.wild_ok(w, a) {
                        return false;
                    }
                }
                Atom::Ew => {
                    if self.ewd[a] != 1 {
                        return false;
                    }
                }
            }
            // 3. bead atomic property.
            if !bead.prop.is_empty() {
                let pred = if b == 0 { caid } else { path[b] };
                if !self.apcheck(a, Some(pred), bead.prop) {
                    return false;
                }
            }
        }
        true
    }

    /// Is there an assignment of one matched path to each chain slot such that all
    /// branches are simultaneously and distinctly satisfied?
    fn dccheck(
        &self,
        slot: usize,
        chain_count: usize,
        schains: &[Schain],
        sci: &mut [usize],
        chains: &[&[Bead]],
    ) -> bool {
        for (i, sc) in schains.iter().enumerate() {
            if sc.id != slot {
                continue;
            }
            sci[slot] = i;
            let done = if slot + 1 == chain_count {
                self.chain_check(sci, schains, chains, chain_count)
            } else {
                self.dccheck(slot + 1, chain_count, schains, sci, chains)
            };
            if done {
                return true;
            }
        }
        false
    }

    fn chain_check(
        &self,
        sci: &[usize],
        schains: &[Schain],
        chains: &[&[Bead]],
        chain_count: usize,
    ) -> bool {
        for i in 0..chain_count {
            for j in (i + 1)..chain_count {
                let si = sci[i];
                let sj = sci[j];
                if si == sj {
                    return false;
                }
                let a = &schains[si].at;
                let b = &schains[sj].at;
                let min = a.len().min(b.len());
                let mut differ = false;
                for k in 0..min {
                    if a[k] != b[k] {
                        differ = true;
                        break;
                    }
                }
                if !differ {
                    return false; // one path is a prefix of the other
                }
                for k in 0..min {
                    let ci = chains[i][k].cesname;
                    let cj = chains[j][k].cesname;
                    if a[k] == b[k] && ci != cj {
                        return false;
                    }
                    if a[k] != b[k] && ci == cj {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Decide whether atom `atomno` satisfies the structured chem-env constraint `env`
    /// (a set of required neighbour chains).
    fn jatspecial(&self, atomno: usize, env: &[&[Bead]]) -> bool {
        if env.is_empty() {
            return false;
        }
        let maxchain = env.iter().map(|c| c.len()).max().unwrap_or(0);
        let mut cesindex = vec![0u32; env.len()];
        let mut schains: Vec<Schain> = Vec::new();
        let mut path: Vec<usize> = Vec::new();
        self.cematch(atomno, env, maxchain, &mut path, atomno, &mut cesindex, &mut schains);
        if cesindex.iter().any(|&c| c == 0) {
            return false;
        }
        let mut sci = vec![0usize; env.len()];
        self.dccheck(0, env.len(), &schains, &mut sci, env)
    }

    // ---- top-level typing loop -------------------------------------------------

    /// Type a single atom `i` by testing the rules in DEF order (first match wins).
    fn type_atom(&self, i: usize) -> Option<&'static str> {
        for rule in self.rules.iter() {
            if let Some(name) = self.try_rule(i, rule) {
                return Some(name);
            }
        }
        None
    }

    /// Evaluate one rule against atom `i`. `Some(name)` on match, `None` to try the next
    /// rule. Every present field (scalar, `prop`, `env`) must be satisfied.
    fn try_rule(&self, i: usize, rule: &'static Rule) -> Option<&'static str> {
        if let Some(v) = rule.z {
            if v != self.z[i] {
                return None;
            }
        }
        if let Some(v) = rule.connum {
            if v as usize != self.props.connum[i] {
                return None;
            }
        }
        if let Some(v) = rule.nh {
            if v as usize != self.props.nh[i] {
                return None;
            }
        }
        if let Some(v) = rule.ew {
            // Electron-withdrawing neighbours of the first-bonded atom (H typing).
            let first = *self.con[i].first().unwrap_or(&i);
            if v as usize != self.props.ewd_neigh[first] {
                return None;
            }
        }
        if !rule.prop.is_empty() && !self.apcheck(i, None, rule.prop) {
            return None;
        }
        if !rule.env.is_empty() && !self.jatspecial(i, rule.env) {
            return None;
        }
        Some(rule.name)
    }
}

/// Assign a raw GAFF type name to every atom (before the conjugation parity split).
pub(crate) fn jat(ctx: &Ctx) -> Result<Vec<String>, FFError> {
    let n = ctx.z.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        match ctx.type_atom(i) {
            Some(name) => out.push(name.to_string()),
            None => {
                return Err(FFError::UntypedAtom { ff: ctx.ff, local: i, z: ctx.z[i] })
            }
        }
    }
    Ok(out)
}
