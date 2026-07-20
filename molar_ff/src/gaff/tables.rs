//! Compile-time GAFF atom-type tables.
//!
//! The rule set from `ATOMTYPE_GFF.DEF` is held directly as the `const` [`RULES`] array
//! (first match wins, in file order) — no runtime string parsing. Every field is fully
//! structured:
//!
//! - the scalar fields ([`Rule::z`], [`Rule::connum`], [`Rule::nh`], [`Rule::ew`]) are
//!   `Option<u8>` (`None` = unconstrained);
//! - [`Rule::prop`] (the atomic-property `[...]` field) is `&[&[Pred]]`: an AND of
//!   OR-groups of [`Pred`]s (empty = unconstrained);
//! - [`Rule::env`] (the chemical-environment `(...)` field) is `&[&[Bead]]`: the required
//!   neighbour chains, each a path of [`Bead`]s from the central atom outward (empty =
//!   unconstrained).
//!
//! Regenerate with `tools/gen_tables.py` after editing `ATOMTYPE_GFF.DEF`.

/// One GAFF atom-type rule. A `None` scalar or an empty `prop`/`env` imposes no
/// constraint; a rule matches when every present field is satisfied.
pub struct Rule {
    /// The atom-type name emitted on a match.
    pub name: &'static str,
    /// Atomic number.
    pub z: Option<u8>,
    /// Number of attached atoms (coordination).
    pub connum: Option<u8>,
    /// Number of attached hydrogens.
    pub nh: Option<u8>,
    /// Number of electron-withdrawing neighbours of the bonded heavy atom (H typing).
    pub ew: Option<u8>,
    /// Atomic-property constraint: an AND of OR-groups (empty = none).
    pub prop: &'static [&'static [Pred]],
    /// Chemical-environment constraint: the required neighbour chains (empty = none).
    pub env: &'static [&'static [Bead]],
}

/// One atomic-property predicate, e.g. `RG6`, `1RG6`, `AR2`, `sb'`, `0DL`.
pub struct Pred {
    /// Required count: `None` means "> 0", `Some(k)` means "== k".
    pub n: Option<u8>,
    /// The property being tested.
    pub p: Prop,
    /// Bond-flag quote code: `0` bare, `1` (`'`) bond order to the predecessor, `2` (`''`)
    /// *not* that bond order to the predecessor.
    pub q: u8,
}

/// An atomic property that a [`Pred`] can test. Some variants (`NonRing`, `Ab`) belong to
/// the full GAFF property vocabulary but are not exercised by the current rule corpus.
#[allow(dead_code)]
pub enum Prop {
    /// Ring membership: `Ring(0)` = any ring (`RG`); `Ring(3..=10)` = an n-membered ring.
    Ring(u8),
    /// Aromaticity class `AR1..=AR5`.
    Arom(u8),
    /// Non-ring atom (`NR`).
    NonRing,
    /// Single bond (`SB`/`sb`).
    Sb,
    /// Double bond (`DB`/`db`).
    Db,
    /// Triple bond (`TB`/`tb`).
    Tb,
    /// Delocalised bond (`DL`).
    Dl,
    /// Aromatic bond (`AB`).
    Ab,
}

/// One bead of a chemical-environment chain: a required neighbour at some depth.
pub struct Bead {
    /// Which atom this bead stands for.
    pub atom: Atom,
    /// Required coordination number (`None` = any).
    pub n: Option<u8>,
    /// Bead atomic-property constraint (same shape as [`Rule::prop`]; empty = none).
    pub prop: &'static [&'static [Pred]],
    /// Unique bead id (an incrementing per-token counter) used to tell branches apart
    /// when checking that matched paths are distinct.
    pub cesname: u32,
}

/// The atom a [`Bead`] requires. `Ew` belongs to the full GAFF vocabulary but is not
/// exercised by the current rule corpus.
#[allow(dead_code)]
pub enum Atom {
    /// A specific element by atomic number.
    Z(u8),
    /// A wildcard atom `XX`/`XA`/`XB`/`XC`/`XD`, indexed `0..=4` into [`WILDATOMS`].
    Wild(u8),
    /// The `EW` pseudo-atom (an electron-withdrawing neighbour).
    Ew,
}

/// A wildcard atom: a name and the `(atomic number, connum)` pairs it stands for
/// (`connum == 0` matches any coordination). The `name` is documentation; matching is by
/// index into [`WILDATOMS`].
#[allow(dead_code)]
pub struct WildAtom {
    pub name: &'static str,
    pub pairs: &'static [(u8, u8)],
}

/// Wildcard-atom definitions (`XX`, `XA`, `XB`, `XC`, `XD`).
pub const WILDATOMS: &[WildAtom] = &[
    WildAtom { name: "XX", pairs: &[(6, 0), (7, 0), (8, 0), (16, 0), (15, 0)] },
    WildAtom { name: "XA", pairs: &[(8, 0), (16, 0)] },
    WildAtom { name: "XB", pairs: &[(7, 0), (15, 0)] },
    WildAtom { name: "XC", pairs: &[(9, 0), (17, 0), (35, 0), (53, 0)] },
    WildAtom { name: "XD", pairs: &[(16, 0), (15, 0)] },
];

/// Element symbols indexed by atomic number (index 0 unused).
#[allow(dead_code)]
const SYMBOLS: [&str; 119] = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
    "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// The element symbol for atomic number `z` (`""` if out of range / `z == 0`).
#[allow(dead_code)]
pub fn element_symbol(z: u8) -> &'static str {
    SYMBOLS.get(z as usize).copied().unwrap_or("")
}

/// All GAFF atom-type rules, in definition order (first match wins).
pub const RULES: &[Rule] = &[
    Rule { name: "cx", z: Some(6), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[] },
    Rule { name: "cy", z: Some(6), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[] },
    Rule { name: "c3", z: Some(6), connum: Some(4), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "c", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: Some(2), p: Prop::Dl, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(1), n: Some(1), prop: &[], cesname: 1 }]] },
    Rule { name: "c", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: Some(1), p: Prop::Db, q: 0 }], &[Pred { n: Some(0), p: Prop::Dl, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(1), n: Some(1), prop: &[], cesname: 1 }]] },
    Rule { name: "c", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: Some(3), p: Prop::Sb, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(1), n: Some(1), prop: &[], cesname: 1 }]] },
    Rule { name: "cz", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(7), n: Some(3), prop: &[], cesname: 1 }], &[Bead { atom: Atom::Z(7), n: Some(3), prop: &[], cesname: 2 }], &[Bead { atom: Atom::Z(7), n: Some(3), prop: &[], cesname: 3 }]] },
    Rule { name: "cp", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }], &[Pred { n: Some(1), p: Prop::Ring(6), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(0), n: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }]], cesname: 1 }], &[Bead { atom: Atom::Wild(0), n: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }]], cesname: 2 }], &[Bead { atom: Atom::Wild(0), n: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }]], cesname: 3 }]] },
    Rule { name: "ca", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }]], env: &[] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[] },
    Rule { name: "cc", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[] },
    Rule { name: "ce", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "ce", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "ce", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "ce", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "ce", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "cu", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[] },
    Rule { name: "cv", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[] },
    Rule { name: "c2", z: Some(6), connum: Some(3), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "cg", z: Some(6), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Tb, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "cg", z: Some(6), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Tb, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "cg", z: Some(6), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Tb, q: 0 }]], env: &[&[Bead { atom: Atom::Z(7), n: Some(1), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "cg", z: Some(6), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Tb, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "c1", z: Some(6), connum: Some(2), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "c1", z: Some(6), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "hn", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(7), n: None, prop: &[], cesname: 1 }]] },
    Rule { name: "ho", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(8), n: None, prop: &[], cesname: 1 }]] },
    Rule { name: "hs", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(16), n: None, prop: &[], cesname: 1 }]] },
    Rule { name: "hp", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(15), n: None, prop: &[], cesname: 1 }]] },
    Rule { name: "hx", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: None, prop: &[], cesname: 1 }, Bead { atom: Atom::Z(7), n: Some(4), prop: &[], cesname: 2 }]] },
    Rule { name: "hw", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(8), n: None, prop: &[], cesname: 1 }, Bead { atom: Atom::Z(1), n: Some(1), prop: &[], cesname: 2 }]] },
    Rule { name: "h3", z: Some(1), connum: Some(1), nh: None, ew: Some(3), prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(4), prop: &[], cesname: 1 }]] },
    Rule { name: "h2", z: Some(1), connum: Some(1), nh: None, ew: Some(2), prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(4), prop: &[], cesname: 1 }]] },
    Rule { name: "h1", z: Some(1), connum: Some(1), nh: None, ew: Some(1), prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(4), prop: &[], cesname: 1 }]] },
    Rule { name: "hc", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(4), prop: &[], cesname: 1 }]] },
    Rule { name: "h5", z: Some(1), connum: Some(1), nh: None, ew: Some(2), prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }]] },
    Rule { name: "h4", z: Some(1), connum: Some(1), nh: None, ew: Some(1), prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }]] },
    Rule { name: "ha", z: Some(1), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "f", z: Some(9), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "cl", z: Some(17), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "br", z: Some(35), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "i", z: Some(53), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "pc", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "pb", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }]], env: &[] },
    Rule { name: "pe", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pe", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pe", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(1), n: Some(1), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pe", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "pe", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "pe", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "p2", z: Some(15), connum: Some(2), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "p2", z: Some(15), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "px", z: Some(15), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "px", z: Some(15), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "px", z: Some(15), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "px", z: Some(15), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "p4", z: Some(15), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(1), n: Some(1), prop: &[], cesname: 1 }]] },
    Rule { name: "p3", z: Some(15), connum: Some(3), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "py", z: Some(15), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "py", z: Some(15), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "py", z: Some(15), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "py", z: Some(15), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "p5", z: Some(15), connum: Some(4), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "p5", z: Some(15), connum: Some(5), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "p5", z: Some(15), connum: Some(6), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "ni", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(1), n: Some(1), prop: &[], cesname: 2 }]] },
    Rule { name: "nj", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(1), n: Some(1), prop: &[], cesname: 2 }]] },
    Rule { name: "n", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(1), n: Some(1), prop: &[], cesname: 2 }]] },
    Rule { name: "nk", z: Some(7), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[] },
    Rule { name: "nl", z: Some(7), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[] },
    Rule { name: "n4", z: Some(7), connum: Some(4), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "no", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(8), n: Some(1), prop: &[], cesname: 1 }], &[Bead { atom: Atom::Z(8), n: Some(1), prop: &[], cesname: 2 }]] },
    Rule { name: "na", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }, Pred { n: None, p: Prop::Arom(2), q: 0 }, Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[] },
    Rule { name: "nm", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(0), n: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }, Pred { n: None, p: Prop::Arom(2), q: 0 }, Pred { n: None, p: Prop::Arom(3), q: 0 }]], cesname: 1 }]] },
    Rule { name: "nm", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nm", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nm", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(15), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nn", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(0), n: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }, Pred { n: None, p: Prop::Arom(2), q: 0 }, Pred { n: None, p: Prop::Arom(3), q: 0 }]], cesname: 1 }]] },
    Rule { name: "nn", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nn", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[&[Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nn", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[&[Bead { atom: Atom::Z(15), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nh", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Wild(0), n: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }, Pred { n: None, p: Prop::Arom(2), q: 0 }, Pred { n: None, p: Prop::Arom(3), q: 0 }]], cesname: 1 }]] },
    Rule { name: "nh", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nh", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nh", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[], env: &[&[Bead { atom: Atom::Z(15), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "np", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[] },
    Rule { name: "nq", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[] },
    Rule { name: "n3", z: Some(7), connum: Some(3), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "nb", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Arom(1), q: 0 }]], env: &[] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(2), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(3), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Z(6), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 1 }, Bead { atom: Atom::Wild(2), n: Some(2), prop: &[], cesname: 2 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }], &[Pred { n: None, p: Prop::Arom(3), q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], env: &[&[Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], cesname: 1 }, Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], cesname: 2 }, Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], cesname: 3 }]] },
    Rule { name: "nc", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], env: &[&[Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], cesname: 1 }], &[Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], cesname: 2 }, Bead { atom: Atom::Z(7), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Ring(5), q: 0 }]], cesname: 3 }]] },
    Rule { name: "ne", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "ne", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "ne", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(1), n: Some(1), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "ne", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "ne", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "ne", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Sb, q: 0 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "n1", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: Some(2), p: Prop::Db, q: 0 }]], env: &[] },
    Rule { name: "n1", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Tb, q: 0 }], &[Pred { n: None, p: Prop::Sb, q: 0 }]], env: &[] },
    Rule { name: "n2", z: Some(7), connum: Some(2), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "n1", z: Some(7), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "o", z: Some(8), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "oh", z: Some(8), connum: Some(2), nh: Some(1), ew: None, prop: &[], env: &[] },
    Rule { name: "oh", z: Some(8), connum: Some(2), nh: Some(2), ew: None, prop: &[], env: &[] },
    Rule { name: "oh", z: Some(8), connum: Some(3), nh: Some(1), ew: None, prop: &[], env: &[] },
    Rule { name: "oh", z: Some(8), connum: Some(3), nh: Some(2), ew: None, prop: &[], env: &[] },
    Rule { name: "oh", z: Some(8), connum: Some(3), nh: Some(3), ew: None, prop: &[], env: &[] },
    Rule { name: "op", z: Some(8), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[] },
    Rule { name: "oq", z: Some(8), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[] },
    Rule { name: "os", z: Some(8), connum: Some(2), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "os", z: Some(8), connum: Some(3), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "os", z: Some(8), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "s", z: Some(16), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "s2", z: Some(16), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[] },
    Rule { name: "s2", z: Some(16), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Tb, q: 0 }]], env: &[] },
    Rule { name: "sh", z: Some(16), connum: Some(2), nh: Some(1), ew: None, prop: &[], env: &[] },
    Rule { name: "sh", z: Some(16), connum: Some(2), nh: Some(2), ew: None, prop: &[], env: &[] },
    Rule { name: "sp", z: Some(16), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(3), q: 0 }]], env: &[] },
    Rule { name: "sq", z: Some(16), connum: Some(2), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Ring(4), q: 0 }]], env: &[] },
    Rule { name: "ss", z: Some(16), connum: Some(2), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "sx", z: Some(16), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "sx", z: Some(16), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "sx", z: Some(16), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "sx", z: Some(16), connum: Some(3), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "s4", z: Some(16), connum: Some(3), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "sy", z: Some(16), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(2), n: Some(2), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "sy", z: Some(16), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Z(6), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }]], cesname: 1 }]] },
    Rule { name: "sy", z: Some(16), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(3), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "sy", z: Some(16), connum: Some(4), nh: None, ew: None, prop: &[&[Pred { n: None, p: Prop::Db, q: 0 }]], env: &[&[Bead { atom: Atom::Wild(4), n: Some(4), prop: &[&[Pred { n: None, p: Prop::Sb, q: 1 }], &[Pred { n: None, p: Prop::Db, q: 0 }]], cesname: 1 }]] },
    Rule { name: "s6", z: Some(16), connum: Some(4), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "s6", z: Some(16), connum: Some(5), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "s6", z: Some(16), connum: Some(6), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "He", z: Some(2), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Li", z: Some(3), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Be", z: Some(4), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "B", z: Some(5), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ne", z: Some(10), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Na", z: Some(11), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Mg", z: Some(12), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Al", z: Some(13), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Si", z: Some(14), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ar", z: Some(18), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "K", z: Some(19), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ca", z: Some(20), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Sc", z: Some(21), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ti", z: Some(22), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "V", z: Some(23), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Cr", z: Some(24), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Mn", z: Some(25), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Fe", z: Some(26), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Co", z: Some(27), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ni", z: Some(28), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Cu", z: Some(29), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Zn", z: Some(30), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ga", z: Some(31), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ge", z: Some(32), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "As", z: Some(33), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Se", z: Some(34), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Kr", z: Some(36), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Rb", z: Some(37), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Sr", z: Some(38), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Y", z: Some(39), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Zr", z: Some(40), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Nb", z: Some(41), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Mo", z: Some(42), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Tc", z: Some(43), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ru", z: Some(44), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Rh", z: Some(45), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Pd", z: Some(46), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ag", z: Some(47), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Cd", z: Some(48), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "In", z: Some(49), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Sn", z: Some(50), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Sb", z: Some(51), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Te", z: Some(52), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Xe", z: Some(54), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Cs", z: Some(55), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ba", z: Some(56), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "La", z: Some(57), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ce", z: Some(58), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Pr", z: Some(59), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Nd", z: Some(60), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Pm", z: Some(61), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Sm", z: Some(62), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Eu", z: Some(63), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Gd", z: Some(64), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Tb", z: Some(65), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Dy", z: Some(66), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ho", z: Some(67), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Er", z: Some(68), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Tm", z: Some(69), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Yb", z: Some(70), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Lu", z: Some(71), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Hf", z: Some(72), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ta", z: Some(73), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "W", z: Some(74), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Re", z: Some(75), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Os", z: Some(76), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ir", z: Some(77), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Pt", z: Some(78), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Au", z: Some(79), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Hg", z: Some(80), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Tl", z: Some(81), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Pb", z: Some(82), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Bi", z: Some(83), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Po", z: Some(84), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "At", z: Some(85), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Rn", z: Some(86), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Fr", z: Some(87), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ra", z: Some(88), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ac", z: Some(89), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Th", z: Some(90), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Pa", z: Some(91), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "U", z: Some(92), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Np", z: Some(93), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Pu", z: Some(94), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Am", z: Some(95), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Cm", z: Some(96), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Bk", z: Some(97), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Cf", z: Some(98), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Es", z: Some(99), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Fm", z: Some(100), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Md", z: Some(101), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "No", z: Some(102), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Lr", z: Some(103), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Rf", z: Some(104), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Db", z: Some(105), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Sg", z: Some(106), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Bh", z: Some(107), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Hs", z: Some(108), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Mt", z: Some(109), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "Ds", z: Some(103), connum: None, nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "LP", z: Some(0), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "lp", z: Some(0), connum: Some(1), nh: None, ew: None, prop: &[], env: &[] },
    Rule { name: "DU", z: None, connum: None, nh: None, ew: None, prop: &[], env: &[] },
];
