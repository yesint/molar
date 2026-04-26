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

/// Backbone atom **local** indices (for `RandomPosProvider`) and H position.
enum BackboneResidue {
    Valid { n: usize, ca: usize, c: usize, o: usize, h: Option<Pos> },
    Break,
}

impl BackboneResidue {
    #[inline]
    fn is_valid(&self) -> bool { matches!(self, Self::Valid { .. }) }
    #[inline]
    fn is_break(&self) -> bool { matches!(self, Self::Break) }
}

//──────────────────────────────────────────────────────────────────────────────
// Public struct
//──────────────────────────────────────────────────────────────────────────────

/// DSSP secondary structure assignment result for a selection.
///
/// Constructed via [`MeasureAtomPos::dssp`].
pub struct Dssp {
    backbone: Vec<BackboneResidue>,  // local selection indices
    hbond:    HashSet<(usize, usize)>,
    ss:       Vec<SS>,
}

impl Dssp {
    /// Run the full DSSP algorithm on `sel` and return the result.
    pub fn new(sel: &(impl ParticleIterProvider + PosProvider)) -> Self {
        let backbone = Self::extract_backbone(sel);
        let ss = backbone.iter().map(|b| match b {
            BackboneResidue::Break => SS::Break,
            BackboneResidue::Valid { .. } => SS::Coil,
        }).collect();

        let mut this = Self { backbone, hbond: HashSet::new(), ss };
        this.reconstruct_h(sel);
        this.hbond = this.compute_hbonds(sel);
        this.detect_helices();
        this.detect_beta();
        this.detect_bends(sel);
        this.detect_polyproline(sel);
        this
    }

    //──────────────────────────────────────────────────────────────────────────
    // Public accessors
    //──────────────────────────────────────────────────────────────────────────

    /// Per-residue secondary structure codes, ordered by residue index.
    pub fn ss(&self) -> &[SS] { &self.ss }

    /// Compact string of DSSP codes, one character per residue.
    pub fn ss_string(&self) -> String { self.ss.iter().map(|s| s.to_char()).collect() }

    /// Number of residues.
    pub fn len(&self) -> usize { self.ss.len() }

    /// Returns `true` if the result contains no residues.
    pub fn is_empty(&self) -> bool { self.ss.is_empty() }

    //──────────────────────────────────────────────────────────────────────────
    // Stage 1: Backbone extraction
    //──────────────────────────────────────────────────────────────────────────

    fn extract_backbone(sel: &(impl ParticleIterProvider + PosProvider)) -> Vec<BackboneResidue> {
        struct ResEntry {
            n:  Option<usize>,
            ca: Option<usize>,
            c:  Option<usize>,
            o:  Option<usize>,
            h:  Option<Pos>,
        }

        let mut by_res: BTreeMap<usize, ResEntry> = BTreeMap::new();

        for (local_idx, p) in sel.iter_particle().enumerate() {
            let entry = by_res.entry(p.atom.resindex).or_insert(ResEntry {
                n: None, ca: None, c: None, o: None, h: None,
            });
            match p.atom.name.as_str() {
                "N"  => entry.n  = Some(local_idx),  // LOCAL — correct for get_pos_unchecked
                "CA" => entry.ca = Some(local_idx),
                "C"  => entry.c  = Some(local_idx),
                // Accept standard and C-terminal oxygen naming
                "O" | "OT1" | "OXT" => {
                    if entry.o.is_none() { entry.o = Some(local_idx); }
                }
                // Accept common amide H naming conventions; store pos directly
                "H" | "HN" | "1H" | "H1" => entry.h = Some(*p.pos),
                _ => {}
            }
        }

        by_res.into_values().map(|e| {
            match (e.n, e.ca, e.c, e.o) {
                (Some(n), Some(ca), Some(c), Some(o)) =>
                    BackboneResidue::Valid { n, ca, c, o, h: e.h },
                _ => BackboneResidue::Break,
            }
        }).collect()
    }

    //──────────────────────────────────────────────────────────────────────────
    // Stage 2: Hydrogen reconstruction
    //──────────────────────────────────────────────────────────────────────────

    /// Reconstruct amide H positions for all non-proline residues.
    ///
    /// Formula matches Gromacs `gmx dssp -hmode dssp`:
    ///   H_i = N_i + normalize(C_{i-1} − O_{i-1}) × 0.1 nm
    ///
    /// Overwrites any previously stored H so all bonds are computed from
    /// reconstructed H, consistent with the Gromacs -hmode dssp reference.
    fn reconstruct_h(&mut self, sel: &impl PosProvider) {
        for i in 1..self.backbone.len() {
            // Extract local indices (immutable borrows end here via Copy)
            let (prev_c, prev_o, n_idx) = match (&self.backbone[i - 1], &self.backbone[i]) {
                (BackboneResidue::Valid { c, o, .. }, BackboneResidue::Valid { n, .. })
                    => (*c, *o, *n),
                _ => continue,
            };

            let n_pos  = *get_pos(sel, n_idx);
            let c_prev = *get_pos(sel, prev_c);
            let o_prev = *get_pos(sel, prev_o);

            // Gromacs -hmode dssp: H = N + normalize(C_prev − O_prev) × 0.1 nm
            let v = c_prev - o_prev;
            let norm = v.norm();
            if norm > 1e-6 {
                let h_pos = n_pos + v / norm * 0.1;
                if let BackboneResidue::Valid { h, .. } = &mut self.backbone[i] {
                    *h = Some(h_pos);
                }
            }
        }
    }

    //──────────────────────────────────────────────────────────────────────────
    // Stage 3: Hydrogen bond detection
    //──────────────────────────────────────────────────────────────────────────

    fn compute_hbonds(&self, sel: &impl PosProvider) -> HashSet<(usize, usize)> {
        let n = self.backbone.len();
        let mut hbond = HashSet::new();

        for donor in 0..n {
            let (dn_idx, dh) = match &self.backbone[donor] {
                BackboneResidue::Valid { n, h: Some(h), .. } => (*n, *h),
                _ => continue, // no H or Break — cannot donate
            };
            let dn = get_pos(sel, dn_idx);

            for acceptor in 0..n {
                if donor == acceptor { continue; }
                if donor.abs_diff(acceptor) < 2 { continue; } // skip adjacent

                let (ac_idx, ao_idx) = match &self.backbone[acceptor] {
                    BackboneResidue::Valid { c, o, .. } => (*c, *o),
                    BackboneResidue::Break => continue,
                };

                let ac = get_pos(sel, ac_idx);
                let ao = get_pos(sel, ao_idx);

                let e = hbond_energy(dn, &dh, ac, ao);
                if e < HBOND_THRESHOLD {
                    hbond.insert((donor, acceptor));
                }
            }
        }

        hbond
    }

    //──────────────────────────────────────────────────────────────────────────
    // Stage 4: Helix detection
    //──────────────────────────────────────────────────────────────────────────

    fn detect_helices(&mut self) {
        let n = self.backbone.len();

        // Pre-compute n-turn start positions for all 3 helix types
        let mut n_turn_at = [vec![false; n], vec![false; n], vec![false; n]];
        for (t, n_turn) in [(0usize, 3usize), (1, 4), (2, 5)] {
            for i in 0..n {
                if i + n_turn < n
                    && self.backbone[i].is_valid()
                    && self.backbone[i + n_turn].is_valid()
                    && self.hbond.contains(&(i + n_turn, i))
                {
                    n_turn_at[t][i] = true;
                }
            }
        }

        // Process helices in Gromacs order: alpha (4) first, then 3-10 (3), then pi (5).
        for &(t, n_turn, helix_code) in &[
            (1usize, 4usize, SS::AlphaHelix),
            (0usize, 3usize, SS::Helix310),
            (2usize, 5usize, SS::PiHelix),
        ] {
            let min_priority_block: u8 = match helix_code {
                SS::AlphaHelix => u8::MAX,
                SS::Helix310   => SS::BetaSheet.priority(),
                SS::PiHelix    => SS::Helix310.priority(),
                _              => u8::MAX,
            };

            let turns = &n_turn_at[t];

            for i in 0..n {
                if turns[i] && i + 1 < n && turns[i + 1] {
                    let lo = i + 1;
                    let hi = (i + n_turn).min(n - 1);
                    if (lo..=hi).any(|k| self.ss[k].priority() >= min_priority_block) {
                        continue;
                    }
                    for k in lo..=hi {
                        self.ss[k].try_assign(helix_code);
                    }
                }
            }

            // Mark interior turn residues as T
            for i in 0..n {
                if turns[i] {
                    for k in (i + 1)..(i + n_turn).min(n) {
                        self.ss[k].try_assign(SS::Turn);
                    }
                }
            }
        }
    }

    //──────────────────────────────────────────────────────────────────────────
    // Stage 5: Beta bridge / sheet detection
    // Algorithm matches Gromacs 2025 dssp.cpp analyzeBridgesAndStrandsPatterns()
    //──────────────────────────────────────────────────────────────────────────

    fn detect_beta(&mut self) {
        let n = self.backbone.len();
        if n < 5 { return; }

        let mut ap_partners:  Vec<Vec<usize>> = vec![vec![]; n];
        let mut par_partners: Vec<Vec<usize>> = vec![vec![]; n];

        for i in 1..n.saturating_sub(4) {
            if self.backbone[i - 1].is_break() || self.backbone[i].is_break() || self.backbone[i + 1].is_break() {
                continue;
            }
            for j in (i + 3)..n.saturating_sub(1) {
                if self.backbone[j - 1].is_break() || self.backbone[j].is_break() || self.backbone[j + 1].is_break() {
                    continue;
                }

                let ap = (self.hbond.contains(&(i + 1, j - 1)) && self.hbond.contains(&(j + 1, i - 1)))
                      || (self.hbond.contains(&(j, i))          && self.hbond.contains(&(i, j)));
                if ap {
                    ap_partners[i].push(j);
                    ap_partners[j].push(i);
                }

                let par = (self.hbond.contains(&(i + 1, j)) && self.hbond.contains(&(j, i - 1)))
                       || (self.hbond.contains(&(j + 1, i)) && self.hbond.contains(&(i, j - 1)));
                if par {
                    par_partners[i].push(j);
                    par_partners[j].push(i);
                }
            }
        }

        for i in 1..n.saturating_sub(1) {
            for gap in 1..=2_usize {
                let ij = i + gap;
                if ij >= n { continue; }

                // has_break(i): i==0 || i+1>=n || backbone[i-1] or backbone[i+1] is Break
                let i_break  = i  == 0 || i  + 1 >= n
                    || self.backbone[i  - 1].is_break() || self.backbone[i  + 1].is_break();
                let ij_break = ij == 0 || ij + 1 >= n
                    || self.backbone[ij - 1].is_break() || self.backbone[ij + 1].is_break();
                if i_break || ij_break { continue; }

                for (pi, pij) in [(&ap_partners[i], &ap_partners[ij]),
                                  (&par_partners[i], &par_partners[ij])] {
                    if pi.is_empty() || pij.is_empty() { continue; }
                    for &ip in pi.iter() {
                        for &jp in pij.iter() {
                            if ip.abs_diff(jp) < 6 {
                                let lo = ip.min(jp);
                                let hi = ip.max(jp);
                                for k in lo..=hi {
                                    self.ss[k].try_assign(SS::BetaSheet);
                                }
                                for k in i..=ij {
                                    self.ss[k].try_assign(SS::BetaSheet);
                                }
                            }
                        }
                    }
                }
            }
        }

        for i in 1..n.saturating_sub(1) {
            if self.backbone[i].is_break() { continue; }
            if self.ss[i] != SS::BetaSheet
                && (!ap_partners[i].is_empty() || !par_partners[i].is_empty())
            {
                self.ss[i].try_assign(SS::BetaBridge);
            }
        }
    }

    //──────────────────────────────────────────────────────────────────────────
    // Stage 6: Bend detection
    //──────────────────────────────────────────────────────────────────────────

    fn detect_bends(&mut self, sel: &impl PosProvider) {
        let n = self.backbone.len();

        for i in 2..n.saturating_sub(2) {
            let (ca_m2_idx, ca_0_idx, ca_p2_idx) = match (
                &self.backbone[i - 2],
                &self.backbone[i],
                &self.backbone[i + 2],
            ) {
                (BackboneResidue::Valid { ca: m2, .. },
                 BackboneResidue::Valid { ca: z,  .. },
                 BackboneResidue::Valid { ca: p2, .. }) => (*m2, *z, *p2),
                _ => continue,
            };

            let ca_m2 = get_pos(sel, ca_m2_idx);
            let ca_0  = get_pos(sel, ca_0_idx);
            let ca_p2 = get_pos(sel, ca_p2_idx);

            let v1 = ca_0  - ca_m2;
            let v2 = ca_p2 - ca_0;

            let n1 = v1.norm();
            let n2 = v2.norm();
            if n1 < 1e-6 || n2 < 1e-6 { continue; }

            let cos_angle = (v1.dot(&v2) / (n1 * n2)).clamp(-1.0, 1.0);
            let angle_deg = cos_angle.acos().to_degrees();

            if angle_deg >= 70.0 {
                self.ss[i].try_assign(SS::Bend);
            }
        }
    }

    //──────────────────────────────────────────────────────────────────────────
    // Stage 7: PolyProline II helix detection (torsion angle-based)
    // Matches Gromacs `gmx dssp -polypro` (default behavior, PPStretches::Default)
    //──────────────────────────────────────────────────────────────────────────

    fn detect_polyproline(&mut self, sel: &impl PosProvider) {
        let n = self.backbone.len();

        let mut phi = vec![360.0; n];
        let mut psi = vec![360.0; n];

        for i in 1..n.saturating_sub(1) {
            let (prev_c_idx, n_idx, ca_idx, c_idx) = match (&self.backbone[i - 1], &self.backbone[i]) {
                (BackboneResidue::Valid { c: pc, .. },
                 BackboneResidue::Valid { n, ca, c, .. }) => (*pc, *n, *ca, *c),
                _ => continue,
            };

            let c_prev  = get_pos(sel, prev_c_idx);
            let n_curr  = get_pos(sel, n_idx);
            let ca_curr = get_pos(sel, ca_idx);
            let c_curr  = get_pos(sel, c_idx);
            phi[i] = dihedral_gmx(c_prev, n_curr, ca_curr, c_curr);

            if let BackboneResidue::Valid { n: nn, .. } = &self.backbone[i + 1] {
                let n_next = get_pos(sel, *nn);
                psi[i] = dihedral_gmx(n_curr, ca_curr, c_curr, n_next);
            }
        }

        let phi_min = -75.0 - 29.0; // -104°
        let phi_max = -75.0 + 29.0; // -46°
        let psi_min = 145.0 - 29.0; // 116°
        let psi_max = 145.0 + 29.0; // 174°

        let in_phi = |k: usize| phi[k] >= phi_min && phi[k] <= phi_max;
        let in_psi = |k: usize| psi[k] >= psi_min && psi[k] <= psi_max;

        for i in 1..n {
            if i + 3 >= n { break; }
            if !in_phi(i) || !in_phi(i + 1) || !in_phi(i + 2) { continue; }
            if !in_psi(i) || !in_psi(i + 1) || !in_psi(i + 2) { continue; }
            self.ss[i].try_assign(SS::PolyProline);
            self.ss[i + 1].try_assign(SS::PolyProline);
            self.ss[i + 2].try_assign(SS::PolyProline);
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
// H-bond energy constants and helper
//──────────────────────────────────────────────────────────────────────────────

/// Electrostatic H-bond factor: 0.084 × 33.2 (nm units, kcal/mol).
const HBOND_FACTOR: Float = 0.084 * 33.2;
const HBOND_THRESHOLD: Float = -0.5; // kcal/mol

#[inline]
fn hbond_energy(donor_n: &Pos, donor_h: &Pos, acceptor_c: &Pos, acceptor_o: &Pos) -> Float {
    let r_on = (acceptor_o - donor_n).norm();
    let r_ch = (acceptor_c - donor_h).norm();
    let r_oh = (acceptor_o - donor_h).norm();
    let r_cn = (acceptor_c - donor_n).norm();

    if r_oh < 1e-4 || r_on < 1e-4 || r_ch < 1e-4 || r_cn < 1e-4 {
        return 0.0;
    }

    HBOND_FACTOR * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn)
}

/// Dihedral angle A-B-C-D using the Gromacs formula (degrees).
/// Returns 360.0 for degenerate geometry.
#[inline]
fn dihedral_gmx(a: &Pos, b: &Pos, c: &Pos, d: &Pos) -> Float {
    let vec_ba = a - b;
    let vec_cd = d - c;
    let vec_cb = b - c;
    let cbxba      = vec_cb.cross(&vec_ba);
    let cbxcd      = vec_cb.cross(&vec_cd);
    let cbxcbxcd   = vec_cb.cross(&cbxcd);
    let vdot1 = cbxcd.dot(&cbxcd);
    let vdot2 = cbxcbxcd.dot(&cbxcbxcd);
    if vdot1 > 0.0 && vdot2 > 0.0 {
        let x = cbxba.dot(&cbxcd)    / vdot1.sqrt();
        let y = cbxba.dot(&cbxcbxcd) / vdot2.sqrt();
        y.atan2(x).to_degrees()
    } else {
        360.0
    }
}

#[inline]
fn get_pos(sel: &impl PosProvider, local_idx: usize) -> &Pos {
    unsafe { sel.get_pos_unchecked(local_idx) }
}

//──────────────────────────────────────────────────────────────────────────────
// Tests
//──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Known geometry: ideal alpha-helix H-bond geometry
    #[test]
    fn hbond_energy_alpha_helix() {
        let donor_n = Pos::new(0.0, 0.0, 0.0);
        let donor_h = Pos::new(0.0, 0.0, -0.101);
        let acceptor_o = Pos::new(0.0, 0.0, -0.290);
        let acceptor_c = Pos::new(0.0, 0.0, -0.390);

        let e = hbond_energy(&donor_n, &donor_h, &acceptor_c, &acceptor_o);
        assert!(e < -0.5, "Expected H-bond energy < -0.5 kcal/mol, got {:.3}", e);
    }

    #[test]
    fn hbond_energy_bad_geometry() {
        let donor_n = Pos::new(0.0, 0.0, 0.0);
        let donor_h = Pos::new(0.0, 0.0, 0.101);
        let acceptor_o = Pos::new(0.0, 0.0, -1.0);
        let acceptor_c = Pos::new(0.0, 0.0, -1.1);

        let e = hbond_energy(&donor_n, &donor_h, &acceptor_c, &acceptor_o);
        assert!(e > HBOND_THRESHOLD, "No H-bond expected at long distance, got {:.3}", e);
    }

    /// Load PDB, run DSSP, compare with reference `.dat` file.
    fn check_dssp(pdb: &str, dat: &str, threshold: Float, strip_breaks: bool) -> anyhow::Result<()> {
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

        let accuracy = matches as Float / ss_string.len() as Float;
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
        check_dssp("tests/7pbd.pdb", "tests/7pbd_dssp.dat", 0.95, true)
    }
}
