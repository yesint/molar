//! Secondary-structure types shared by all assignment algorithms.
//!
//! Holds the per-residue [`SS`] code, the [`SsAlgorithm`] selector, and the
//! unified [`SsResult`] returned by [`crate::measure::Measure::ss_compute`].
//! The concrete algorithms live in [`crate::dssp`] (Kabsch–Sander / GROMACS) and
//! [`crate::dss`] (PyMOL `dss`).

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
            SS::AlphaHelix => 'H',
            SS::Helix310 => 'G',
            SS::PiHelix => 'I',
            SS::PolyProline => 'P',
            SS::BetaSheet => 'E',
            SS::BetaBridge => 'B',
            SS::Turn => 'T',
            SS::Bend => 'S',
            SS::Coil => '~',
            SS::Break => '=',
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
    pub(crate) fn priority(self) -> u8 {
        match self {
            SS::Break => 255,
            SS::AlphaHelix => 8,
            SS::BetaBridge => 7,
            SS::BetaSheet => 6,
            SS::Helix310 => 5,
            SS::PiHelix => 4,
            SS::Turn => 3,
            SS::Bend => 2,
            SS::PolyProline => 1,
            SS::Coil => 0,
        }
    }

    /// Overwrite with `new` only if it has strictly higher priority.
    pub(crate) fn try_assign(&mut self, new: SS) {
        if new.priority() > self.priority() {
            *self = new;
        }
    }
}

/// Which secondary-structure assignment algorithm to run.
///
/// See [`crate::measure::Measure::ss_compute`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SsAlgorithm {
    /// Canonical Kabsch–Sander DSSP (proper β ladders + bounded β-bulge).
    #[default]
    Dssp,
    /// GROMACS-flavored DSSP — `analyzeBridgesAndStrandsPatterns`, which fills
    /// the residue span between nearby bridge partners and tends to over-extend
    /// β-strands. Kept for fidelity with `gmx dssp`.
    DsspGmx,
    /// PyMOL `dss`: fast H-bond-pattern + φ/ψ method, 3-state (H/E/coil).
    Dss,
}

impl SsAlgorithm {
    /// Canonical lower-case name.
    pub fn as_str(self) -> &'static str {
        match self {
            SsAlgorithm::Dssp => "dssp",
            SsAlgorithm::DsspGmx => "dssp_gmx",
            SsAlgorithm::Dss => "dss",
        }
    }
}

impl std::fmt::Display for SsAlgorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for SsAlgorithm {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "dssp" | "vanilla" => Ok(SsAlgorithm::Dssp),
            "dssp_gmx" | "dssp-gmx" | "gmx" => Ok(SsAlgorithm::DsspGmx),
            "dss" | "pymol" => Ok(SsAlgorithm::Dss),
            other => Err(format!(
                "unknown SS algorithm '{other}' (expected dssp|dssp_gmx|dss)"
            )),
        }
    }
}

/// Per-residue secondary-structure assignment, ordered by ascending residue
/// index (one entry per distinct `resindex` in the selection). The uniform
/// result of [`crate::measure::Measure::ss_compute`] across all algorithms.
pub struct SsResult {
    ss: Vec<SS>,
}

impl SsResult {
    pub(crate) fn new(ss: Vec<SS>) -> Self {
        Self { ss }
    }

    /// Per-residue codes.
    pub fn ss(&self) -> &[SS] {
        &self.ss
    }

    /// Compact string, one character per residue.
    pub fn ss_string(&self) -> String {
        self.ss.iter().map(|s| s.to_char()).collect()
    }

    pub fn len(&self) -> usize {
        self.ss.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ss.is_empty()
    }

    /// Consume into the underlying vector.
    pub fn into_vec(self) -> Vec<SS> {
        self.ss
    }
}
