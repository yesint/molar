//! DSSP secondary structure assignment.
//!
//! Implements the Kabsch & Sander (1983) algorithm for assigning secondary structure
//! to protein residues based on backbone hydrogen bond patterns.
//!
//! # References
//! - W. Kabsch, C. Sander, Biopolymers 22:2577 (1983)
//! - Gromacs `gmx dssp` implementation

use crate::prelude::*;
use std::collections::{BTreeMap, HashSet};

//──────────────────────────────────────────────────────────────────────────────
// Public types
//──────────────────────────────────────────────────────────────────────────────

/// Secondary structure code assigned to a single residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SS {
    /// Alpha helix (4-turn helix, NH(i) → C=O(i-4)). Output: `H`
    AlphaHelix,
    /// 3₁₀ helix (3-turn helix, NH(i) → C=O(i-3)). Output: `G`
    Helix310,
    /// π helix (5-turn helix, NH(i) → C=O(i-5)). Output: `I`
    PiHelix,
    /// Poly-proline II helix (torsion-angle based). Output: `P`
    PolyProline,
    /// Extended beta strand in beta sheet (ladder). Output: `E`
    BetaSheet,
    /// Isolated beta bridge residue. Output: `B`
    BetaBridge,
    /// Hydrogen-bonded turn. Output: `T`
    Turn,
    /// Bend (Cα angle ≥ 70°). Output: `S`
    Bend,
    /// Loop / coil (no special designation). Output: `~`
    #[default]
    Coil,
    /// Break (missing backbone atoms). Output: `=`
    Break,
}

impl SS {
    /// Single-character code used by DSSP / Gromacs output.
    pub fn to_char(self) -> char {
        match self {
            SS::AlphaHelix  => 'H',
            SS::Helix310    => 'G',
            SS::PiHelix     => 'I',
            SS::PolyProline => 'P',
            SS::BetaSheet   => 'E',
            SS::BetaBridge  => 'B',
            SS::Turn        => 'T',
            SS::Bend        => 'S',
            SS::Coil        => '~',
            SS::Break       => '=',
        }
    }

    /// Parse a single DSSP character back to an `SS` variant.
    pub fn from_char(c: char) -> Option<Self> {
        Some(match c {
            'H' => SS::AlphaHelix,
            'G' => SS::Helix310,
            'I' => SS::PiHelix,
            'P' => SS::PolyProline,
            'E' => SS::BetaSheet,
            'B' => SS::BetaBridge,
            'T' => SS::Turn,
            'S' => SS::Bend,
            '~' | 'C' => SS::Coil,
            '=' => SS::Break,
            _ => return None,
        })
    }

    /// Assignment priority: higher wins; `Break` is never overwritten.
    fn priority(self) -> u8 {
        match self {
            SS::Break       => 255,
            SS::AlphaHelix  => 8,
            SS::BetaBridge  => 7,
            SS::BetaSheet   => 6,
            SS::Helix310    => 5,
            SS::PiHelix     => 4,
            SS::Turn        => 3,
            SS::Bend        => 2,
            SS::PolyProline => 1,
            SS::Coil        => 0,
        }
    }

    fn try_assign(&mut self, new: SS) {
        if new.priority() > self.priority() {
            *self = new;
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Internal types
//──────────────────────────────────────────────────────────────────────────────

/// Backbone atom indices and reconstructed/explicit H position for one residue.
struct BackboneResidue {
    n:  usize,       // atom index of amide N
    ca: usize,       // atom index of Cα
    c:  usize,       // atom index of carbonyl C
    o:  usize,       // atom index of carbonyl O
    h:  Option<Pos>, // amide H position (explicit or reconstructed)
}

//──────────────────────────────────────────────────────────────────────────────
// Algorithm entry point
//──────────────────────────────────────────────────────────────────────────────

pub(crate) fn compute_dssp(sel: &impl AtomPosAnalysis) -> Vec<SS> {
    let mut backbone = extract_backbone(sel);
    reconstruct_h(&mut backbone, sel);

    let n = backbone.len();
    let mut ss = vec![SS::Coil; n];

    // Mark breaks immediately
    for (i, b) in backbone.iter().enumerate() {
        if b.is_none() {
            ss[i] = SS::Break;
        }
    }

    let hbond = compute_hbonds(&backbone, sel);

    detect_helices(&backbone, &hbond, &mut ss);
    detect_beta(&backbone, &hbond, &mut ss);
    detect_bends(&backbone, sel, &mut ss);
    detect_polyproline(&backbone, sel, &mut ss);

    ss
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 1: Backbone extraction
//──────────────────────────────────────────────────────────────────────────────

fn extract_backbone(sel: &impl AtomPosAnalysis) -> Vec<Option<BackboneResidue>> {
    struct ResEntry {
        n:  Option<usize>,
        ca: Option<usize>,
        c:  Option<usize>,
        o:  Option<usize>,
        h:  Option<Pos>,
    }

    let mut by_res: BTreeMap<usize, ResEntry> = BTreeMap::new();

    for p in sel.iter_particle() {
        let entry = by_res.entry(p.atom.resindex).or_insert(ResEntry {
            n: None, ca: None, c: None, o: None, h: None,
        });
        match p.atom.name.as_str() {
            "N"  => entry.n  = Some(p.id),
            "CA" => entry.ca = Some(p.id),
            "C"  => entry.c  = Some(p.id),
            // Accept standard and C-terminal oxygen naming
            "O" | "OT1" | "OXT" => {
                if entry.o.is_none() { entry.o = Some(p.id); }
            }
            // Accept common amide H naming conventions
            "H" | "HN" | "1H" | "H1" => entry.h = Some(*p.pos),
            _ => {}
        }
    }

    by_res.into_values().map(|e| {
        match (e.n, e.ca, e.c, e.o) {
            (Some(n), Some(ca), Some(c), Some(o)) =>
                Some(BackboneResidue { n, ca, c, o, h: e.h }),
            _ => None,
        }
    }).collect()
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 2: Hydrogen reconstruction
//──────────────────────────────────────────────────────────────────────────────

/// Reconstruct amide H positions for all non-proline residues.
///
/// Formula matches Gromacs `gmx dssp -hmode dssp`:
///   H_i = N_i + normalize(C_{i-1} − O_{i-1}) × 0.1 nm
///
/// The mixed-unit scaling (C-O vector in nm, distance in Å) in Gromacs effectively
/// places H at N + unit_vector × 0.1 nm (= 1.0 Å), where the direction is
/// from O toward C of the preceding residue — pointing along the trans-peptide amide.
///
/// This function overwrites any previously stored H position so that all bonds are
/// computed from reconstructed H, consistent with the Gromacs -hmode dssp reference.
fn reconstruct_h(backbone: &mut Vec<Option<BackboneResidue>>, sel: &impl AtomPosAnalysis) {
    for i in 1..backbone.len() {
        if backbone[i].is_none() { continue; }

        let (prev_c, prev_o) = match backbone[i - 1].as_ref() {
            Some(p) => (p.c, p.o),
            None => continue, // preceding residue is a Break — no reconstruction possible
        };

        let n_pos  = get_pos(sel, backbone[i].as_ref().unwrap().n);
        let c_prev = get_pos(sel, prev_c);
        let o_prev = get_pos(sel, prev_o);

        // Gromacs -hmode dssp: H = N + normalize(C_prev − O_prev) × 0.1 nm
        let v = c_prev - o_prev; // vector from O toward C of the previous residue
        let norm = v.norm();
        if norm > 1e-6 {
            let h_pos = n_pos + v / norm * 0.1; // 0.1 nm = 1.0 Å
            backbone[i].as_mut().unwrap().h = Some(h_pos);
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 3: Hydrogen bond detection
//──────────────────────────────────────────────────────────────────────────────

/// Electrostatic H-bond factor: 0.084 × 33.2 (converts to nm, kcal/mol).
/// Original formula uses Å: E = 0.084 × 332 × (1/r_ON + 1/r_CH − 1/r_OH − 1/r_CN).
/// For nm: multiply by 0.1, giving factor = 0.084 × 33.2 = 2.7888.
const HBOND_FACTOR: f32 = 0.084 * 33.2; // kcal/mol·nm
const HBOND_THRESHOLD: f32 = -0.5;      // kcal/mol

/// Returns a set of `(donor_local_idx, acceptor_local_idx)` pairs where
/// the amide N–H of `donor` residue forms an H-bond with the C=O of `acceptor`.
fn compute_hbonds(backbone: &[Option<BackboneResidue>], sel: &impl AtomPosAnalysis)
    -> HashSet<(usize, usize)>
{
    let n = backbone.len();
    let mut hbond = HashSet::new();

    for donor in 0..n {
        let (dn, dh) = match backbone[donor].as_ref() {
            Some(b) => (get_pos(sel, b.n), match b.h {
                Some(h) => h,
                None => continue, // no H — cannot donate
            }),
            None => continue,
        };

        for acceptor in 0..n {
            if donor == acceptor { continue; }
            if donor.abs_diff(acceptor) < 2 { continue; } // skip adjacent residues

            let (ac, ao) = match backbone[acceptor].as_ref() {
                Some(b) => (get_pos(sel, b.c), get_pos(sel, b.o)),
                None => continue,
            };

            let e = hbond_energy(dn, &dh, ac, ao);
            if e < HBOND_THRESHOLD {
                hbond.insert((donor, acceptor));
            }
        }
    }

    hbond
}

#[inline]
fn hbond_energy(donor_n: &Pos, donor_h: &Pos, acceptor_c: &Pos, acceptor_o: &Pos) -> f32 {
    let r_on = (acceptor_o - donor_n).norm();
    let r_ch = (acceptor_c - donor_h).norm();
    let r_oh = (acceptor_o - donor_h).norm();
    let r_cn = (acceptor_c - donor_n).norm();

    // Guard against degenerate geometry
    if r_oh < 1e-4 || r_on < 1e-4 || r_ch < 1e-4 || r_cn < 1e-4 {
        return 0.0;
    }

    HBOND_FACTOR * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn)
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 4: Helix detection
//──────────────────────────────────────────────────────────────────────────────

fn detect_helices(
    backbone: &[Option<BackboneResidue>],
    hbond: &HashSet<(usize, usize)>,
    ss: &mut Vec<SS>,
) {
    let n = backbone.len();

    // Pre-compute n-turn start positions for all 3 helix types
    let mut n_turn_at = [vec![false; n], vec![false; n], vec![false; n]];
    for (t, n_turn) in [(0usize, 3usize), (1, 4), (2, 5)] {
        for i in 0..n {
            if i + n_turn < n
                && backbone[i].is_some()
                && backbone[i + n_turn].is_some()
                && hbond.contains(&(i + n_turn, i))
            {
                n_turn_at[t][i] = true;
            }
        }
    }

    // Process helices in Gromacs order: alpha (4) first, then 3-10 (3), then pi (5).
    // Alpha always assigns; 3-10 and pi have an "empty check" matching Gromacs behavior:
    //   3-10 (G, priority 5): skip if any residue in range already has E, B, or H (priority ≥ 6)
    //   pi   (I, priority 4): skip if any residue in range already has G, E, B, or H (priority ≥ 5)
    for &(t, n_turn, helix_code) in &[
        (1usize, 4usize, SS::AlphaHelix),
        (0usize, 3usize, SS::Helix310),
        (2usize, 5usize, SS::PiHelix),
    ] {
        let min_priority_block: u8 = match helix_code {
            SS::AlphaHelix => u8::MAX, // no block: always assign
            SS::Helix310   => SS::BetaSheet.priority(), // block if any residue has E/B/H (≥ 6)
            SS::PiHelix    => SS::Helix310.priority(),  // block if any residue has G/E/B/H (≥ 5)
            _              => u8::MAX,
        };

        let turns = &n_turn_at[t];

        // A proper helix requires two consecutive n-turns:
        // turns[i] AND turns[i+1] → residues i+1 … i+n_turn get the helix code.
        for i in 0..n {
            if turns[i] && i + 1 < n && turns[i + 1] {
                let lo = i + 1;
                let hi = (i + n_turn).min(n - 1);
                // Empty check: skip if any residue in range is already ≥ min_priority_block
                if (lo..=hi).any(|k| ss[k].priority() >= min_priority_block) {
                    continue;
                }
                for k in lo..=hi {
                    ss[k].try_assign(helix_code);
                }
            }
        }

        // Mark interior turn residues as T (NOT the acceptor i nor the donor i+n_turn).
        // Only residues not already assigned a helix get T (via try_assign priority).
        for i in 0..n {
            if turns[i] {
                for k in (i + 1)..(i + n_turn).min(n) {
                    ss[k].try_assign(SS::Turn);
                }
            }
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 5: Beta bridge / sheet detection
// Algorithm matches Gromacs 2025 dssp.cpp analyzeBridgesAndStrandsPatterns()
//──────────────────────────────────────────────────────────────────────────────

fn detect_beta(
    backbone: &[Option<BackboneResidue>],
    hbond: &HashSet<(usize, usize)>,
    ss: &mut Vec<SS>,
) {
    let n = backbone.len();
    if n < 5 { return; }

    // Per-residue bridge partner lists (matches Gromacs per-residue storage)
    let mut ap_partners:  Vec<Vec<usize>> = vec![vec![]; n];
    let mut par_partners: Vec<Vec<usize>> = vec![vec![]; n];

    // Bridge detection: i in 1..n-4, j in i+3..n-1
    // Requires backbone[i-1], backbone[i], backbone[i+1],
    //          backbone[j-1], backbone[j], backbone[j+1] all present (no chain break)
    for i in 1..n.saturating_sub(4) {
        if backbone[i - 1].is_none() || backbone[i].is_none() || backbone[i + 1].is_none() {
            continue;
        }
        for j in (i + 3)..n.saturating_sub(1) {
            if backbone[j - 1].is_none() || backbone[j].is_none() || backbone[j + 1].is_none() {
                continue;
            }

            // Anti-parallel bridge (i, j) — Gromacs / Kabsch & Sander:
            //   Pattern B: hbond(i+1, j-1) AND hbond(j+1, i-1)
            //   Pattern A: hbond(j, i)     AND hbond(i, j)      [mutual]
            let ap = (hbond.contains(&(i + 1, j - 1)) && hbond.contains(&(j + 1, i - 1)))
                  || (hbond.contains(&(j, i))          && hbond.contains(&(i, j)));

            if ap {
                ap_partners[i].push(j);
                ap_partners[j].push(i);
            }

            // Parallel bridge (i, j) — Gromacs / Kabsch & Sander:
            //   Pattern A: hbond(i+1, j) AND hbond(j, i-1)
            //   Pattern B: hbond(j+1, i) AND hbond(i, j-1)
            let par = (hbond.contains(&(i + 1, j)) && hbond.contains(&(j, i - 1)))
                   || (hbond.contains(&(j + 1, i)) && hbond.contains(&(i, j - 1)));

            if par {
                par_partners[i].push(j);
                par_partners[j].push(i);
            }
        }
    }

    // Strand detection — matches Gromacs second loop in analyzeBridgesAndStrandsPatterns()
    // For each i, check i+gap for gap in {1, 2}:
    //   if both residues have bridges of the same type, and their partners are within 6,
    //   mark all residues between partners as E, and all residues i..=i+gap as E.
    let has_break = |k: usize| -> bool {
        k == 0 || k + 1 >= n || backbone[k - 1].is_none() || backbone[k + 1].is_none()
    };

    for i in 1..n.saturating_sub(1) {
        for gap in 1..=2_usize {
            let ij = i + gap;
            if ij >= n { continue; }
            if has_break(i) || has_break(ij) { continue; }

            for (pi, pij) in [(&ap_partners[i], &ap_partners[ij]),
                              (&par_partners[i], &par_partners[ij])] {
                if pi.is_empty() || pij.is_empty() { continue; }
                for &ip in pi.iter() {
                    for &jp in pij.iter() {
                        let delta = ip.abs_diff(jp);
                        if delta < 6 {
                            // Mark E: all residues between the two partners (second strand)
                            let lo = ip.min(jp);
                            let hi = ip.max(jp);
                            for k in lo..=hi {
                                ss[k].try_assign(SS::BetaSheet);
                            }
                            // Mark E: residues i through i+gap (first strand)
                            for k in i..=ij {
                                ss[k].try_assign(SS::BetaSheet);
                            }
                        }
                    }
                }
            }
        }
    }

    // Isolated bridges: has bridge partners but was not assigned E → assign B
    for i in 1..n.saturating_sub(1) {
        if backbone[i].is_none() { continue; }
        if ss[i] != SS::BetaSheet && (!ap_partners[i].is_empty() || !par_partners[i].is_empty()) {
            ss[i].try_assign(SS::BetaBridge);
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 6: Bend detection
//──────────────────────────────────────────────────────────────────────────────

fn detect_bends(
    backbone: &[Option<BackboneResidue>],
    sel: &impl AtomPosAnalysis,
    ss: &mut Vec<SS>,
) {
    let n = backbone.len();

    for i in 2..n.saturating_sub(2) {
        if backbone[i - 2].is_none() || backbone[i].is_none() || backbone[i + 2].is_none() {
            continue;
        }

        let ca_m2 = get_pos(sel, backbone[i - 2].as_ref().unwrap().ca);
        let ca_0  = get_pos(sel, backbone[i].as_ref().unwrap().ca);
        let ca_p2 = get_pos(sel, backbone[i + 2].as_ref().unwrap().ca);

        let v1 = ca_0  - ca_m2;
        let v2 = ca_p2 - ca_0;

        let n1 = v1.norm();
        let n2 = v2.norm();
        if n1 < 1e-6 || n2 < 1e-6 { continue; }

        let cos_angle = (v1.dot(&v2) / (n1 * n2)).clamp(-1.0, 1.0);
        let angle_deg = cos_angle.acos().to_degrees();

        if angle_deg >= 70.0 {
            ss[i].try_assign(SS::Bend);
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Stage 7: PolyProline II helix detection (torsion angle-based)
// Matches Gromacs `gmx dssp -polypro` (default behavior, PPStretches::Default)
//──────────────────────────────────────────────────────────────────────────────

/// Detect PolyProline II helices from backbone φ/ψ torsion angles.
///
/// Assigns `SS::PolyProline` (P) to any Coil residue that belongs to a run of
/// three consecutive residues with φ ∈ [-104°, -46°] and ψ ∈ [116°, 174°].
/// This matches the Gromacs default (`-ppstretch default`, `-polypro true`).
fn detect_polyproline(
    backbone: &[Option<BackboneResidue>],
    sel: &impl AtomPosAnalysis,
    ss: &mut Vec<SS>,
) {
    let n = backbone.len();

    // Phi and psi angles in degrees; 360.0 marks invalid / not computed.
    let mut phi = vec![360.0f32; n];
    let mut psi = vec![360.0f32; n];

    // phi[i] = dihedral(C[i-1], N[i], CA[i], C[i]) — needs backbone[i-1] and backbone[i]
    // psi[i] = dihedral(N[i], CA[i], C[i], N[i+1]) — needs backbone[i] and backbone[i+1]
    for i in 1..n.saturating_sub(1) {
        if backbone[i - 1].is_none() || backbone[i].is_none() { continue; }
        let prev = backbone[i - 1].as_ref().unwrap();
        let curr = backbone[i].as_ref().unwrap();
        let c_prev  = get_pos(sel, prev.c);
        let n_curr  = get_pos(sel, curr.n);
        let ca_curr = get_pos(sel, curr.ca);
        let c_curr  = get_pos(sel, curr.c);
        phi[i] = dihedral_gmx(c_prev, n_curr, ca_curr, c_curr);
        if backbone[i + 1].is_some() {
            let n_next = get_pos(sel, backbone[i + 1].as_ref().unwrap().n);
            psi[i] = dihedral_gmx(n_curr, ca_curr, c_curr, n_next);
        }
    }

    // PP ranges (Gromacs / libcifpp values): phi ∈ [-75±29], psi ∈ [145±29]
    let phi_min = -75.0_f32 - 29.0; // -104°
    let phi_max = -75.0_f32 + 29.0; // -46°
    let psi_min = 145.0_f32 - 29.0; // 116°
    let psi_max = 145.0_f32 + 29.0; // 174°

    // Default stretch = 3 consecutive residues.
    // Gromacs loop: for i in 1..n while i+3 < n, checks i, i+1, i+2.
    for i in 1..n {
        if i + 3 >= n { break; }
        let in_phi = |k: usize| phi[k] >= phi_min && phi[k] <= phi_max;
        let in_psi = |k: usize| psi[k] >= psi_min && psi[k] <= psi_max;
        if !in_phi(i) || !in_phi(i + 1) || !in_phi(i + 2) { continue; }
        if !in_psi(i) || !in_psi(i + 1) || !in_psi(i + 2) { continue; }
        // Assign PP only to residues currently Coil (priority 0 < PP priority 1)
        ss[i].try_assign(SS::PolyProline);
        ss[i + 1].try_assign(SS::PolyProline);
        ss[i + 2].try_assign(SS::PolyProline);
    }
}

/// Dihedral angle A-B-C-D using the Gromacs formula (degrees).
/// Returns 360.0 for degenerate geometry.
#[inline]
fn dihedral_gmx(a: &Pos, b: &Pos, c: &Pos, d: &Pos) -> f32 {
    let vec_ba = a - b; // B→A
    let vec_cd = d - c; // C→D
    let vec_cb = b - c; // C→B
    let cbxba      = vec_cb.cross(&vec_ba);
    let cbxcd      = vec_cb.cross(&vec_cd);
    let cbxcbxcd   = vec_cb.cross(&cbxcd);
    let vdot1 = cbxcd.dot(&cbxcd);       // |CB×CD|²
    let vdot2 = cbxcbxcd.dot(&cbxcbxcd); // |CB×(CB×CD)|²
    if vdot1 > 0.0 && vdot2 > 0.0 {
        let x = cbxba.dot(&cbxcd)    / vdot1.sqrt();
        let y = cbxba.dot(&cbxcbxcd) / vdot2.sqrt();
        y.atan2(x).to_degrees()
    } else {
        360.0
    }
}

//──────────────────────────────────────────────────────────────────────────────
// Helpers
//──────────────────────────────────────────────────────────────────────────────

#[inline]
fn get_pos(sel: &impl RandomPosProvider, idx: usize) -> &Pos {
    unsafe {sel.get_pos_unchecked(idx)}
}

//──────────────────────────────────────────────────────────────────────────────
// Tests
//──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Known geometry: ideal alpha-helix H-bond geometry
    // N at origin; H points toward O (−z direction); O and C below.
    #[test]
    fn hbond_energy_alpha_helix() {
        // N at 0; H at 0.101 nm below N (pointing toward O); O at 0.29 nm below N
        let donor_n = Pos::new(0.0, 0.0, 0.0);
        let donor_h = Pos::new(0.0, 0.0, -0.101); // H points toward O
        let acceptor_o = Pos::new(0.0, 0.0, -0.290);
        let acceptor_c = Pos::new(0.0, 0.0, -0.390);

        let e = hbond_energy(&donor_n, &donor_h, &acceptor_c, &acceptor_o);
        // Should be a good H-bond (energy < -0.5 kcal/mol)
        assert!(e < -0.5, "Expected H-bond energy < -0.5 kcal/mol, got {:.3}", e);
    }

    #[test]
    fn hbond_energy_bad_geometry() {
        // Very long distance: no H-bond expected
        let donor_n = Pos::new(0.0, 0.0, 0.0);
        let donor_h = Pos::new(0.0, 0.0, 0.101);
        let acceptor_o = Pos::new(0.0, 0.0, -1.0); // 1 nm away
        let acceptor_c = Pos::new(0.0, 0.0, -1.1);

        let e = hbond_energy(&donor_n, &donor_h, &acceptor_c, &acceptor_o);
        assert!(e > HBOND_THRESHOLD, "No H-bond expected at long distance, got {:.3}", e);
    }

    /// Load PDB, run DSSP, compare with reference `.dat` file.
    ///
    /// `strip_breaks` removes `=` characters from the reference before comparing.
    /// Gromacs marks chain-boundary residues as `=` via its C–N distance check;
    /// our backbone extraction simply skips those residues entirely, so the counts
    /// differ when `strip_breaks` is false.
    fn check_dssp(pdb: &str, dat: &str, threshold: f32, strip_breaks: bool) -> anyhow::Result<()> {
        let sys = System::from_file(pdb)?;
        let sel = sys.select_bound("protein")?;
        let ss_string = sel.dssp_string();

        let raw = std::fs::read_to_string(dat)?;
        let expected: String = if strip_breaks {
            raw.trim().chars().filter(|&c| c != '=').collect()
        } else {
            raw.trim().to_string()
        };

        assert_eq!(ss_string.len(), expected.len(),
            "Residue count mismatch: got {}, expected {}", ss_string.len(), expected.len());

        let matches = ss_string.chars().zip(expected.chars())
            .filter(|(a, b)| a == b)
            .count();

        let accuracy = matches as f32 / ss_string.len() as f32;
        println!("DSSP accuracy: {:.1}% ({}/{} residues)", accuracy * 100.0, matches, ss_string.len());
        println!("Got:      {}", ss_string);
        println!("Expected: {}", expected);

        assert!(accuracy >= threshold,
            "DSSP accuracy {:.1}% below {:.0}% threshold", accuracy * 100.0, threshold * 100.0);
        Ok(())
    }

    #[test]
    fn dssp_protein_pdb() -> anyhow::Result<()> {
        check_dssp("tests/protein.pdb", "tests/protein_dssp.dat", 0.98, false)
    }

    #[test]
    fn dssp_2lao() -> anyhow::Result<()> {
        check_dssp("tests/2lao.pdb", "tests/2lao_dssp.dat", 0.95, false)
    }

    #[test]
    fn dssp_7pbd() -> anyhow::Result<()> {
        // Gromacs marks chain-boundary residues as `=`; we skip them → strip before comparing.
        check_dssp("tests/7pbd.pdb", "tests/7pbd_dssp.dat", 0.95, true)
    }

}
