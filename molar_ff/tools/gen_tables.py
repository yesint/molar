#!/usr/bin/env python3
"""Generate the structured GAFF rule tables (`src/gaff/tables.rs`).

Each ATD rule in `ATOMTYPE_GFF.DEF` becomes a fully-structured `Rule` value: the scalar
fields (atomic number / coordination / attached-H / EW-neighbour count) are `Option<u8>`,
the atomic-property field (`f8`, `[...]`) becomes `&[&[Pred]]` (an AND of OR-groups), and
the chemical-environment field (`f9`, `(...)`) becomes `&[&[Bead]]` (the required
neighbour chains). No rule text is parsed at run time.

The f8/f9 grammar reproduced here is the one the matcher used to interpret at run time:
its structure and, crucially, the per-chain `cesname` bead-id counter are byte-for-byte
equivalent to that behaviour.
"""

DEF = "/home/semen/work/Projects/molar/molar_ff/src/gaff/ATOMTYPE_GFF.DEF"
OUT = "/home/semen/work/Projects/molar/molar_ff/src/gaff/tables.rs"

# Element symbols indexed by atomic number (index 0 unused). Kept identical to the
# `SYMBOLS` table emitted into the generated file.
SYMBOLS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
    "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
SYM2Z = {s: i for i, s in enumerate(SYMBOLS) if s}

# Wildcard-atom names in table order (index becomes `Atom::Wild(idx)`).
WILD_NAMES = ["XX", "XA", "XB", "XC", "XD"]

RING_MAP = {"RG": 0, "RG3": 3, "RG4": 4, "RG5": 5, "RG6": 6,
            "RG7": 7, "RG8": 8, "RG9": 9, "RG10": 10}
AROM_MAP = {"AR1": 1, "AR2": 2, "AR3": 3, "AR4": 4, "AR5": 5}


def scal(x):
    return "None" if x == "*" else f"Some({int(x)})"


def is_alpha(ch):
    return ch.isascii() and ch.isalpha()


# ---- f8 atomic-property grammar --------------------------------------------------

def parse_unit(unit):
    """One property unit (e.g. `RG6`, `1RG6`, `AR2`, `sb'`, `0DL`) -> (n, Prop, q).

    `n` is the leading integer count (`None` if absent), `q` is the trailing-quote code
    (0 bare, 1 `'`, 2 `''`) which only bond-flag tokens carry.
    """
    i = 0
    while i < len(unit) and unit[i].isdigit():
        i += 1
    n = int(unit[:i]) if i > 0 else None
    name = unit[i:]

    if name in RING_MAP:
        return (n, f"Prop::Ring({RING_MAP[name]})", 0)
    if name in AROM_MAP:
        return (n, f"Prop::Arom({AROM_MAP[name]})", 0)
    if name == "NR":
        return (n, "Prop::NonRing", 0)

    if len(name) < 2:
        raise SystemExit(f"unrecognised property unit {unit!r}")
    two = name[:2]
    if two in ("SB", "sb"):
        variant = "Prop::Sb"
    elif two in ("DB", "db"):
        variant = "Prop::Db"
    elif two in ("TB", "tb"):
        variant = "Prop::Tb"
    elif two == "DL":
        variant = "Prop::Dl"
    elif two == "AB":
        variant = "Prop::Ab"
    else:
        raise SystemExit(f"unrecognised property unit {unit!r}")

    q = 0
    if len(name) > 2 and name[2] == "'":
        q = 1
        if len(name) > 3 and name[3] == "'":
            q = 2
    return (n, variant, q)


def parse_prop(s):
    """`[...]` field -> list of AND-groups, each a list of OR-alternative units."""
    if not s or s == "*":
        return []
    constraints, units, cur = [], [], ""
    for ch in s:
        if ch == '[':
            continue
        elif ch == ']':
            units.append(cur); cur = ""
            constraints.append(units); units = []
            break
        elif ch == '.':
            units.append(cur); cur = ""
        elif ch == ',':
            units.append(cur); cur = ""
            constraints.append(units); units = []
        else:
            cur += ch
    return [[parse_unit(u) for u in group] for group in constraints]


# ---- f9 chemical-environment grammar ---------------------------------------------

def parse_cenv(keyword):
    """`(...)` field -> list of chains, each a list of beads (from the central atom out).

    Reproduces the stateful walk: two-letter tokens, `[..]` bead properties, `<..>` named
    atoms (skipped), depth (`layer`), chain emission on `,`/`)` (unless right after `)`),
    and the incrementing per-token `cesname` id used to tell branches apart.
    """
    if not keyword or keyword == "*":
        return []
    kw = keyword
    n = len(kw)

    def get(i):
        return kw[i] if 0 <= i < n else '\0'

    def getm(i):
        return get(i - 1) if i != 0 else '\0'

    SZ = 64
    atname = [""] * SZ
    atconnum = [0] * SZ
    apindex = [False] * SZ
    ap = [""] * SZ
    cesname = [0] * SZ

    chains = []
    layer = 0
    index0 = False
    tmpapindex = False
    tmpap = ""
    cesname_index = False
    cea_id = 1

    def make_beed(j):
        return {
            "atname": atname[j],
            "atconnum": atconnum[j],
            "ap": ap[j] if apindex[j] else None,
            "cesname": cesname[j],
        }

    for i in range(n):
        c = kw[i]

        # (A) skip the first letter of a two-letter token (outside [..]/<..>).
        if (not tmpapindex) and (not cesname_index) and is_alpha(c) and is_alpha(get(i + 1)):
            continue

        if c == '(':
            layer += 1
        if c == ')':
            layer = max(0, layer - 1)

        # atomic property [ ... ]
        if (not tmpapindex) and c == '[':
            tmpapindex = True
            tmpap = "["
            continue
        if tmpapindex and c == ']':
            apindex[layer] = True
            tmpap += "]"
            ap[layer] = tmpap
            tmpapindex = False
            continue
        if tmpapindex:
            tmpap += c
            continue

        # named atom < ... > (unused by GAFF, but still parsed and skipped)
        if (not cesname_index) and c == '<':
            cesname_index = True
            continue
        if cesname_index and c == '>':
            cesname_index = False
            continue
        if cesname_index:
            continue

        # record a chain on ',' or ')' unless the previous char was ')'
        if c == ',' and getm(i) != ')':
            chains.append([make_beed(j + 1) for j in range(layer)])
        if c == ')' and getm(i) != ')':
            chains.append([make_beed(j + 1) for j in range(layer + 1)])

        # (B) skip first letter of a two-letter token (as in step (A) above).
        if is_alpha(c) and is_alpha(get(i + 1)):
            continue

        # atom-name character.
        if is_alpha(c):
            index0 = True
            atname[layer] = (getm(i) + c) if is_alpha(getm(i)) else c
            ap[layer] = ""
            apindex[layer] = False
            cesname[layer] = cea_id
            cea_id += 1

        # connum digit (or reset to 0 right after the atom name).
        if c.isdigit():
            atconnum[layer] = int(c)
        elif index0:
            atconnum[layer] = 0
            index0 = False

    return chains


# ---- Rust emission ----------------------------------------------------------------

def emit_pred(pred):
    n, variant, q = pred
    ns = "None" if n is None else f"Some({n})"
    return f"Pred {{ n: {ns}, p: {variant}, q: {q} }}"


def emit_prop(groups):
    if not groups:
        return "&[]"
    gs = ["&[" + ", ".join(emit_pred(p) for p in g) + "]" for g in groups]
    return "&[" + ", ".join(gs) + "]"


def emit_atom(atname):
    if atname == "EW":
        return "Atom::Ew"
    if atname in WILD_NAMES:
        return f"Atom::Wild({WILD_NAMES.index(atname)})"
    z = SYM2Z.get(atname)
    if z is None:
        raise SystemExit(f"unknown atom token {atname!r}")
    return f"Atom::Z({z})"


def emit_bead(b):
    ns = "None" if b["atconnum"] == 0 else f"Some({b['atconnum']})"
    prop = emit_prop(parse_prop(b["ap"])) if b["ap"] else "&[]"
    return (f"Bead {{ atom: {emit_atom(b['atname'])}, n: {ns}, "
            f"prop: {prop}, cesname: {b['cesname']} }}")


def emit_env(chains):
    if not chains:
        return "&[]"
    cs = ["&[" + ", ".join(emit_bead(b) for b in ch) + "]" for ch in chains]
    return "&[" + ", ".join(cs) + "]"


rows = []
for line in open(DEF):
    toks = line.split()
    if not toks or toks[0] != "ATD":
        continue
    name = toks[1]
    vals = []
    for t in toks[2:]:
        if t == "&":
            break
        vals.append(t)
    while len(vals) < 7:
        vals.append("*")
    _f3, f4, f5, f6, f7, f8, f9 = vals[:7]
    prop = emit_prop(parse_prop(f8))
    env = emit_env(parse_cenv(f9))
    rows.append(
        f'    Rule {{ name: "{name}", z: {scal(f4)}, connum: {scal(f5)}, '
        f'nh: {scal(f6)}, ew: {scal(f7)}, prop: {prop}, env: {env} }},'
    )

SYMBOLS_RS = (
    '    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",\n'
    '    "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",\n'
    '    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",\n'
    '    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",\n'
    '    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",\n'
    '    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",\n'
    '    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",\n'
    '    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",\n'
    '    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",'
)

hdr = '''//! Compile-time GAFF atom-type tables.
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
''' + SYMBOLS_RS + '''
];

/// The element symbol for atomic number `z` (`""` if out of range / `z == 0`).
#[allow(dead_code)]
pub fn element_symbol(z: u8) -> &'static str {
    SYMBOLS.get(z as usize).copied().unwrap_or("")
}

/// All GAFF atom-type rules, in definition order (first match wins).
pub const RULES: &[Rule] = &[
'''

with open(OUT, "w") as f:
    f.write(hdr)
    f.write("\n".join(rows))
    f.write("\n];\n")

print(f"wrote {OUT}: {len(rows)} rules")
