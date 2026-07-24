//! PyMOL-style secondary-structure assignment (`dss`).
//!
//! A clean-room Rust reimplementation of the algorithm behind PyMOL's `dss`
//! command (`SelectorAssignSS`, Open-Source PyMOL, © Schrödinger LLC, BSD-style
//! license — algorithm reimplemented from the source, no code copied). It is a
//! fast hydrogen-bond-pattern + φ/ψ-geometry method, *not* DSSP/STRIDE. Compared
//! with [`crate::dssp::Dssp`] it produces cleaner, shorter elements (minimum
//! length 3, every terminal strand residue must be genuinely H-bond-paired),
//! which matches author/STRIDE assignments better — DSSP tends to over-extend
//! β-strands.
//!
//! Output is 3-state (`H`/`E`/`~` → [`SS::AlphaHelix`]/[`SS::BetaSheet`]/
//! [`SS::Coil`]), like PyMOL. Distances are in **nm** (PyMOL's Å constants /10).
//!
//! # References
//! - PyMOL `SelectorAssignSS` (layer3/Selector.cpp) and the backbone H-bond test
//!   `ObjectMoleculeGetCheckHBond`/`TestHBond` (layer2/ObjectMolecule2.cpp).

use crate::prelude::*;
use std::collections::BTreeMap;

//──────────────────────────────────────────────────────────────────────────────
// Constants (PyMOL `dss` H-bond criteria; distances scaled Å → nm)
//──────────────────────────────────────────────────────────────────────────────

const MAX_ANGLE: Float = 63.0; // degrees
const MAX_DIST_AT_MAX_ANGLE: Float = 0.32; // nm (3.2 Å)
const MAX_DIST_AT_ZERO: Float = 0.40; // nm (4.0 Å)
const POWER_A: Float = 1.6;
const POWER_B: Float = 5.0;
const H_BOND_LEN: Float = 0.1; // nm (virtual amide H, 1 Å from N)

// φ/ψ classification windows (degrees), PyMOL `ss_*` setting defaults.
const HELIX_PHI_TARGET: Float = -57.0;
const HELIX_PSI_TARGET: Float = -48.0;
const HELIX_INCLUDE: Float = 55.0;
const HELIX_EXCLUDE: Float = 85.0;
const STRAND_PHI_TARGET: Float = -129.0;
const STRAND_PSI_TARGET: Float = 124.0;
const STRAND_INCLUDE: Float = 40.0;
const STRAND_PHI_EXCLUDE: Float = 100.0;
const STRAND_PSI_EXCLUDE: Float = 90.0;

// Per-residue flag bits (mirroring PyMOL's cSS* masks).
const HELIX3: u32 = 0x0001;
const HELIX4: u32 = 0x0002;
const HELIX5: u32 = 0x0004;
const HELIX_HB: u32 = HELIX3 | HELIX4 | HELIX5;
const PHIPSI_HELIX: u32 = 0x0010;
const PHIPSI_NOT_HELIX: u32 = 0x0020;
const PHIPSI_STRAND: u32 = 0x0040;
const PHIPSI_NOT_STRAND: u32 = 0x0080;
const ANTI_SINGLE: u32 = 0x0100;
const ANTI_DOUBLE: u32 = 0x0200;
const ANTI_BULGE: u32 = 0x0400;
const ANTI_SKIP: u32 = 0x0800;
const PARA_SINGLE: u32 = 0x1000;
const PARA_DOUBLE: u32 = 0x2000;
const PARA_SKIP: u32 = 0x4000;

/// Padding (blank residues) inserted at chain breaks and both ends, so the
/// ±2…±5 neighbor indexing in the rules never crosses a break or runs off the
/// array (PyMOL's `cSSBreakSize`).
const PAD: usize = 5;

//──────────────────────────────────────────────────────────────────────────────
// Public API
//──────────────────────────────────────────────────────────────────────────────

/// PyMOL-style secondary structure assignment for a selection.
pub struct Dss {
    ss: Vec<SS>,
}

impl Dss {
    /// Assign secondary structure to `sel` using PyMOL's `dss` algorithm.
    /// Result is ordered by ascending residue index, one entry per distinct
    /// `resindex` present in the selection (same convention as [`Dssp::ss`]).
    pub fn new(sel: &impl ParticleIterProvider) -> Self {
        Dss { ss: assign(sel) }
    }

    /// Per-residue secondary structure codes, ordered by residue index.
    pub fn ss(&self) -> &[SS] {
        &self.ss
    }

    /// Compact string of codes, one character per residue.
    pub fn ss_string(&self) -> String {
        self.ss.iter().map(|s| s.to_char()).collect()
    }

    pub fn len(&self) -> usize {
        self.ss.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ss.is_empty()
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Internal residue record
//──────────────────────────────────────────────────────────────────────────────

#[derive(Clone)]
struct R {
    real: bool,
    resindex: usize,
    n: Pos,
    ca: Pos,
    c: Pos,
    o: Pos,
    flags: u32,
    ss: u8, // b'L' / b'H' / b'S' / b'h'
    acc: Vec<usize>, // residues whose N donates to THIS residue's O (this = acceptor)
    don: Vec<usize>, // residues whose O accepts from THIS residue's N (this = donor)
}

impl R {
    fn blank() -> Self {
        R {
            real: false,
            resindex: 0,
            n: Pos::origin(),
            ca: Pos::origin(),
            c: Pos::origin(),
            o: Pos::origin(),
            flags: 0,
            ss: b'L',
            acc: Vec::new(),
            don: Vec::new(),
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Driver
//──────────────────────────────────────────────────────────────────────────────

fn assign(sel: &impl ParticleIterProvider) -> Vec<SS> {
    // Group backbone atoms by residue (BTreeMap keeps ascending resindex order).
    struct Bb {
        n: Option<Pos>,
        ca: Option<Pos>,
        c: Option<Pos>,
        o: Option<Pos>,
        chain: char,
    }
    let mut by_res: BTreeMap<usize, Bb> = BTreeMap::new();
    for p in sel.iter_particle() {
        let e = by_res.entry(p.atom.get_resindex()).or_insert(Bb {
            n: None,
            ca: None,
            c: None,
            o: None,
            chain: p.atom.get_chain(),
        });
        match p.atom.name() {
            "N" => e.n = Some(*p.pos),
            "CA" => e.ca = Some(*p.pos),
            "C" => e.c = Some(*p.pos),
            "O" | "OT1" | "OXT" => {
                if e.o.is_none() {
                    e.o = Some(*p.pos);
                }
            }
            _ => {}
        }
    }

    let all_resindices: Vec<usize> = by_res.keys().copied().collect();

    // Build the padded residue array: PAD blanks, then runs of contiguous,
    // same-chain, complete residues separated by PAD blanks, then PAD blanks.
    let mut res: Vec<R> = vec![R::blank(); PAD];
    let mut prev: Option<(char, usize)> = None;
    for &ri in &all_resindices {
        let bb = &by_res[&ri];
        let complete = bb.n.is_some() && bb.ca.is_some() && bb.c.is_some() && bb.o.is_some();
        if !complete {
            // An incomplete residue acts as a chain break.
            if prev.is_some() {
                res.extend(std::iter::repeat_with(R::blank).take(PAD));
                prev = None;
            }
            continue;
        }
        let break_here = matches!(prev, Some((pc, pr)) if bb.chain != pc || ri != pr + 1);
        if break_here {
            res.extend(std::iter::repeat_with(R::blank).take(PAD));
        }
        res.push(R {
            real: true,
            resindex: ri,
            n: bb.n.unwrap(),
            ca: bb.ca.unwrap(),
            c: bb.c.unwrap(),
            o: bb.o.unwrap(),
            flags: 0,
            ss: b'L',
            acc: Vec::new(),
            don: Vec::new(),
        });
        prev = Some((bb.chain, ri));
    }
    res.extend(std::iter::repeat_with(R::blank).take(PAD));
    let n = res.len();

    if n > 2 * PAD {
        detect_hbonds(&mut res);
        classify_phi_psi(&mut res);
        pattern_flags(&mut res);
        assign_helices(&mut res);
        assign_strands(&mut res);
        cleanup(&mut res);
    }

    // Map computed assignments back onto every distinct resindex (Coil default).
    let mut ss_by_ri: BTreeMap<usize, SS> = BTreeMap::new();
    for r in &res {
        if r.real {
            ss_by_ri.insert(
                r.resindex,
                match r.ss {
                    b'H' => SS::AlphaHelix,
                    b'S' => SS::BetaSheet,
                    _ => SS::Coil,
                },
            );
        }
    }
    all_resindices
        .iter()
        .map(|ri| ss_by_ri.get(ri).copied().unwrap_or(SS::Coil))
        .collect()
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 1: backbone hydrogen bonds (PyMOL angle-dependent N···O test)
//──────────────────────────────────────────────────────────────────────────────

fn detect_hbonds(res: &mut [R]) {
    let n = res.len();
    let mut bonds: Vec<(usize, usize)> = Vec::new(); // (acceptor i, donor j)
    for i in PAD..n - PAD {
        if !res[i].real {
            continue;
        }
        for j in PAD..n - PAD {
            if i == j || !res[j].real {
                continue;
            }
            // Exclude residues close in sequence (PyMOL skips ≤5 covalent bonds).
            if (i as isize - j as isize).abs() <= 2 {
                continue;
            }
            // donor j's amide H direction needs the previous residue's C.
            let c_prev = if res[j - 1].real { Some(res[j - 1].c) } else { None };
            if is_hbond(&res[j].n, &res[j].ca, c_prev.as_ref(), &res[i].o, &res[i].c) {
                bonds.push((i, j));
            }
        }
    }
    for (i, j) in bonds {
        if res[i].acc.len() < 5 {
            res[i].acc.push(j);
        }
        if res[j].don.len() < 5 {
            res[j].don.push(i);
        }
    }
}

/// PyMOL backbone H-bond test for donor amide N (`n`, neighbors `ca`/`c_prev`)
/// and acceptor carbonyl O (`o`, carbon `c_acc`). Heavy-atom only; the amide H is
/// inferred as the in-plane bisector pointing away from CA and C_prev.
fn is_hbond(n: &Pos, ca: &Pos, c_prev: Option<&Pos>, o: &Pos, c_acc: &Pos) -> bool {
    let don_to_acc = o - n;
    let dist = don_to_acc.norm();
    if dist > MAX_DIST_AT_ZERO {
        return false; // cheap reject (the cutoff never exceeds this)
    }

    // Virtual amide H direction (unit).
    let h_dir = match c_prev {
        Some(cp) => {
            let u_ca = (ca - n).normalize();
            let u_cp = (cp - n).normalize();
            let s = u_ca + u_cp;
            if s.norm() < 1e-6 {
                (n - ca).normalize()
            } else {
                -s.normalize()
            }
        }
        None => (n - ca).normalize(),
    };
    let h = n + h_dir * H_BOND_LEN;

    // Acceptor cone test: H must be on the carbon side of O (cone_dangle = 0).
    let plane = (o - c_acc).normalize();
    let h_to_acc = (o - h).normalize();
    if h_to_acc.dot(&plane) > 0.0 {
        return false;
    }

    // Donor angle = angle at N between N→H and N→O.
    let dangle = h_dir.dot(&don_to_acc.normalize());
    let angle = if dangle < 1.0 && dangle > 0.0 {
        dangle.clamp(-1.0, 1.0).acos().to_degrees()
    } else if dangle > 0.0 {
        0.0
    } else {
        90.0
    };
    if angle > MAX_ANGLE {
        return false;
    }

    // Angle-dependent distance cutoff: 0.40 nm at 0°, 0.32 nm at 63°.
    let factor_a = 0.5 / MAX_ANGLE.powf(POWER_A);
    let factor_b = 0.5 / MAX_ANGLE.powf(POWER_B);
    let curve = angle.powf(POWER_A) * factor_a + angle.powf(POWER_B) * factor_b;
    let cutoff = MAX_DIST_AT_MAX_ANGLE * curve + MAX_DIST_AT_ZERO * (1.0 - curve);
    dist <= cutoff
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 2: φ/ψ geometry classification
//──────────────────────────────────────────────────────────────────────────────

fn classify_phi_psi(res: &mut [R]) {
    let n = res.len();
    for a in PAD..n - PAD {
        if !(res[a].real && res[a - 1].real && res[a + 1].real) {
            continue;
        }
        let phi = dihedral(&res[a - 1].c, &res[a].n, &res[a].ca, &res[a].c);
        let psi = dihedral(&res[a].n, &res[a].ca, &res[a].c, &res[a + 1].n);

        let h_phi = circ_delta(phi, HELIX_PHI_TARGET);
        let h_psi = circ_delta(psi, HELIX_PSI_TARGET);
        let s_phi = circ_delta(phi, STRAND_PHI_TARGET);
        let s_psi = circ_delta(psi, STRAND_PSI_TARGET);

        let f = &mut res[a].flags;
        if h_psi > HELIX_EXCLUDE || h_phi > HELIX_EXCLUDE {
            *f |= PHIPSI_NOT_HELIX;
        } else if h_psi < HELIX_INCLUDE && h_phi < HELIX_INCLUDE {
            *f |= PHIPSI_HELIX;
        }
        if s_psi > STRAND_PSI_EXCLUDE || s_phi > STRAND_PHI_EXCLUDE {
            *f |= PHIPSI_NOT_STRAND;
        } else if s_psi < STRAND_INCLUDE && s_phi < STRAND_INCLUDE {
            *f |= PHIPSI_STRAND;
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 3: H-bond pattern → flags
//──────────────────────────────────────────────────────────────────────────────

fn pattern_flags(res: &mut [R]) {
    let n = res.len();
    let mut add: Vec<(usize, u32)> = Vec::new();
    for a in PAD..n - PAD {
        if !res[a].real {
            continue;
        }
        // Helix i+3/4/5 H-bonds (this residue acceptor) and i-3/4/5 (donor).
        for &acc in &res[a].acc {
            if acc == a + 3 {
                add.push((a, HELIX3));
            } else if acc == a + 4 {
                add.push((a, HELIX4));
            } else if acc == a + 5 {
                add.push((a, HELIX5));
            }
        }
        for &don in &res[a].don {
            if don + 3 == a {
                add.push((a, HELIX3));
            } else if don + 4 == a {
                add.push((a, HELIX4));
            } else if don + 5 == a {
                add.push((a, HELIX5));
            }
        }

        // Antiparallel double H-bond (reciprocal acceptance a ↔ r2).
        for &r2 in &res[a].acc {
            if res[r2].real && res[r2].acc.contains(&a) {
                add.push((a, ANTI_DOUBLE));
                add.push((r2, ANTI_DOUBLE));
            }
        }
        // Antiparallel β-bulge (reciprocal partner offset by +1).
        for &acc in &res[a].acc {
            let r2 = acc + 1;
            if r2 < n && res[r2].real && res[r2].acc.contains(&a) {
                add.push((a, ANTI_DOUBLE));
                add.push((r2, ANTI_BULGE));
                add.push((r2 - 1, ANTI_BULGE));
            }
        }
        // Antiparallel ladder (single/double): a accepts from r2, a+2 from r2-2.
        if res[a + 1].real && res[a + 2].real {
            for &acc in &res[a].acc {
                if acc < 2 {
                    continue;
                }
                let r2 = acc - 2;
                if res[r2].real && res[r2].acc.contains(&(a + 2)) {
                    add.push((a, ANTI_SINGLE));
                    add.push((a + 1, ANTI_SKIP));
                    add.push((a + 2, ANTI_SINGLE));
                    add.push((r2, ANTI_SINGLE));
                    add.push((r2 + 1, ANTI_SKIP));
                    add.push((r2 + 2, ANTI_SINGLE));
                }
            }
        }
        // Parallel ladder: a accepts from r2, r2 accepts from a+2.
        if res[a + 1].real && res[a + 2].real {
            for &acc in &res[a].acc {
                let r2 = acc;
                if res[r2].real && res[r2].acc.contains(&(a + 2)) {
                    add.push((a, PARA_SINGLE));
                    add.push((a + 1, PARA_SKIP));
                    add.push((a + 2, PARA_SINGLE));
                    add.push((r2, PARA_DOUBLE));
                }
            }
        }
    }
    for (i, f) in add {
        res[i].flags |= f;
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 4: assignment
//──────────────────────────────────────────────────────────────────────────────

fn assign_helices(res: &mut [R]) {
    let n = res.len();
    // Clean internal helical residues (H-bonds on r-1,r,r+1, geometry permits).
    for a in PAD..n - PAD {
        if res[a].real
            && res[a - 1].flags & HELIX_HB != 0
            && res[a].flags & HELIX_HB != 0
            && res[a + 1].flags & HELIX_HB != 0
            && res[a].flags & PHIPSI_NOT_HELIX == 0
        {
            res[a].ss = b'H';
        }
    }
    // Fill a single H-bond-less residue inside an otherwise-helical span ('h').
    for a in PAD..n - PAD {
        if res[a].real
            && res[a - 2].flags & HELIX_HB != 0
            && res[a - 1].flags & HELIX_HB != 0
            && res[a - 1].flags & PHIPSI_HELIX != 0
            && res[a].flags & PHIPSI_HELIX != 0
            && res[a + 1].flags & HELIX_HB != 0
            && res[a + 1].flags & PHIPSI_HELIX != 0
            && res[a + 2].flags & HELIX_HB != 0
        {
            res[a].ss = b'h';
        }
    }
    for a in PAD..n - PAD {
        if res[a].real && res[a].ss == b'h' {
            res[a].flags |= HELIX_HB;
            res[a].ss = b'H';
        }
    }
    // Extend helix ends using geometry (forward then backward, in one pass).
    for a in PAD..n - PAD {
        if !res[a].real {
            continue;
        }
        let (f, fp1, fp2, fm1, fm2) = (
            res[a].flags,
            res[a + 1].flags,
            res[a + 2].flags,
            res[a - 1].flags,
            res[a - 2].flags,
        );
        let ss_next = res[a + 1].ss;
        let ss_prev = res[a - 1].ss;
        let geo = |x: u32| x & HELIX_HB != 0 && x & PHIPSI_HELIX != 0;
        if geo(f) && geo(fp1) && geo(fp2) && ss_next == b'H' {
            res[a].ss = b'H';
        }
        if geo(f) && geo(fm1) && geo(fm2) && ss_prev == b'H' {
            res[a].ss = b'H';
        }
    }
}

fn assign_strands(res: &mut [R]) {
    let n = res.len();
    for a in PAD..n - PAD {
        if !res[a].real {
            continue;
        }
        let fm1 = res[a - 1].flags;
        let f = res[a].flags;
        let fp1 = res[a + 1].flags;

        // Antiparallel.
        if f & ANTI_DOUBLE != 0 && f & PHIPSI_NOT_STRAND == 0 {
            res[a].ss = b'S';
        }
        if f & ANTI_BULGE != 0 && fp1 & ANTI_BULGE != 0 {
            res[a].ss = b'S';
            res[a + 1].ss = b'S';
        }
        if fm1 & ANTI_DOUBLE != 0
            && f & ANTI_SKIP != 0
            && f & PHIPSI_NOT_STRAND == 0
            && fp1 & (ANTI_SINGLE | ANTI_DOUBLE) != 0
        {
            res[a].ss = b'S';
        }
        if fm1 & (ANTI_SINGLE | ANTI_DOUBLE) != 0
            && f & ANTI_SKIP != 0
            && f & PHIPSI_NOT_STRAND == 0
            && fp1 & ANTI_DOUBLE != 0
        {
            res[a].ss = b'S';
        }
        // Open antiparallel ladder, geometry-supported.
        if fm1 & (ANTI_SINGLE | ANTI_DOUBLE) != 0
            && fm1 & PHIPSI_STRAND != 0
            && f & PHIPSI_STRAND != 0
            && fp1 & (ANTI_SINGLE | ANTI_DOUBLE) != 0
            && fp1 & PHIPSI_STRAND != 0
        {
            res[a - 1].ss = b'S';
            res[a].ss = b'S';
            res[a + 1].ss = b'S';
        }

        // Parallel.
        if f & PARA_DOUBLE != 0 && f & PHIPSI_NOT_STRAND == 0 {
            res[a].ss = b'S';
        }
        if fm1 & PARA_DOUBLE != 0
            && f & PARA_SKIP != 0
            && f & PHIPSI_NOT_STRAND == 0
            && fp1 & (PARA_SINGLE | PARA_DOUBLE) != 0
        {
            res[a].ss = b'S';
        }
        if fm1 & (PARA_SINGLE | PARA_DOUBLE) != 0
            && f & PARA_SKIP != 0
            && f & PHIPSI_NOT_STRAND == 0
            && fp1 & PARA_DOUBLE != 0
        {
            res[a].ss = b'S';
        }
        if fm1 & (PARA_SINGLE | PARA_DOUBLE) != 0
            && fm1 & PHIPSI_STRAND != 0
            && f & PARA_SKIP != 0
            && f & PHIPSI_STRAND != 0
            && fp1 & (PARA_SINGLE | PARA_DOUBLE) != 0
            && fp1 & PHIPSI_STRAND != 0
        {
            res[a - 1].ss = b'S';
            res[a].ss = b'S';
            res[a + 1].ss = b'S';
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 5: cleanup (min length 3; terminal strand residues must be paired)
//──────────────────────────────────────────────────────────────────────────────

fn cleanup(res: &mut [R]) {
    let n = res.len();
    let mut repeat = true;
    while repeat {
        repeat = false;
        for a in PAD..n - PAD {
            if !res[a].real {
                continue;
            }
            let ss = res[a].ss;
            let ssm1 = res[a - 1].ss;
            let ssp1 = res[a + 1].ss;
            let ssp2 = res[a + 2].ss;

            // No 2-residue segments.
            if (ss == b'S' && ssp1 == b'S' && ssm1 != b'S' && ssp2 != b'S')
                || (ss == b'H' && ssp1 == b'H' && ssm1 != b'H' && ssp2 != b'H')
            {
                res[a].ss = b'L';
                res[a + 1].ss = b'L';
                repeat = true;
                continue;
            }
            // No 1-residue segments.
            if (ss == b'S' && ssm1 != b'S' && ssp1 != b'S')
                || (ss == b'H' && ssm1 != b'H' && ssp1 != b'H')
            {
                res[a].ss = b'L';
                repeat = true;
                continue;
            }

            // Terminal strand residues must have a real strand-paired partner.
            if ss == b'S' && (ssm1 != b'S' || ssp1 != b'S') {
                let mut found = res[a].acc.iter().any(|&p| res[p].ss == b'S')
                    || res[a].don.iter().any(|&p| res[p].ss == b'S');
                // A "skip" residue may borrow a neighbor's partner.
                if !found && res[a].flags & (ANTI_SKIP | PARA_SKIP) != 0 {
                    if ssp1 == b'S' {
                        found = res[a + 1].acc.iter().any(|&p| res[p].ss == b'S');
                    }
                    if !found && ssm1 == b'S' {
                        found = res[a - 1].don.iter().any(|&p| res[p].ss == b'S');
                    }
                }
                if !found {
                    res[a].ss = b'L';
                    repeat = true;
                }
            }
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Geometry helpers
//──────────────────────────────────────────────────────────────────────────────

/// Dihedral A-B-C-D in degrees (same convention as molar's DSSP `dihedral_gmx`,
/// which is validated against `gmx dssp` — i.e. standard IUPAC φ/ψ signs).
fn dihedral(a: &Pos, b: &Pos, c: &Pos, d: &Pos) -> Float {
    let ba = a - b;
    let cd = d - c;
    let cb = b - c;
    let cbxba = cb.cross(&ba);
    let cbxcd = cb.cross(&cd);
    let cbxcbxcd = cb.cross(&cbxcd);
    let vdot1 = cbxcd.dot(&cbxcd);
    let vdot2 = cbxcbxcd.dot(&cbxcbxcd);
    if vdot1 > 0.0 && vdot2 > 0.0 {
        let x = cbxba.dot(&cbxcd) / vdot1.sqrt();
        let y = cbxba.dot(&cbxcbxcd) / vdot2.sqrt();
        y.atan2(x).to_degrees()
    } else {
        360.0
    }
}

/// Smallest absolute circular difference between two angles (degrees, ≤180).
fn circ_delta(a: Float, target: Float) -> Float {
    let d = (a - target).abs();
    if d > 180.0 {
        360.0 - d
    } else {
        d
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Compare against PyMOL's `dss` output for 2lao. PyMOL (and the PDB author
    /// records) put the strand at 178–180/181, where molar's DSSP over-extends
    /// to 178–185. Reference PyMOL `dss` string (178–188) = "EEE~~~~~~~H".
    #[test]
    fn dss_2lao() {
        // Reference: PyMOL 3.1 `dss` on 2lao (`cmd.dss`; H/S/L → H/E/~).
        const PYMOL: &str = "~~~~EEEEEE~~~~~~~~EEE~~EEEE~HHHHHHHHHHHHH~~EEEEEE~~HHHHHHHHHH~~~~EEE~~~~~~HHHHHH~~~~~~~~~~~EEEEEE~~~~~~~~HHHH~~~EEEEE~~~HHHHHHHHHHHHH~~EEEEE~~HHHHHHHHHH~~~~EEEEEHHHHHHH~~~~HHHH~EEE~~~~~~~HHHH~~~~~~~~~~~~HHHHHHHHHHHHHHHHH~HHHHHHHHH~~~~~~~~";

        let sys = System::from_file("tests/2lao.pdb").expect("load 2lao.pdb");
        let sel = sys.select_bound("protein").expect("select protein");
        let s = Dss::new(&sel).ss_string();
        assert_eq!(s.len(), PYMOL.len(), "residue count");

        let agree = s.bytes().zip(PYMOL.bytes()).filter(|(a, b)| a == b).count();
        let pct = 100.0 * agree as f32 / PYMOL.len() as f32;
        println!("PyMOL-port: {s}");
        println!("agreement with PyMOL dss: {agree}/{} = {pct:.1}%", PYMOL.len());

        // We match PyMOL essentially exactly; require very high agreement.
        assert!(pct >= 98.0, "only {pct:.1}% agreement with PyMOL dss");
        // The motivating case: 178-180 is a short strand, 182-188 is NOT strand
        // (DSSP over-extends 178-185).
        assert_eq!(&s[177..188], "EEE~~~~~~~H");
    }
}
