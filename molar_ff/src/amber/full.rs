use super::*;
use crate::FFType;
use molar::prelude::*;

#[derive(Debug, thiserror::Error)]
pub enum AmberError {
    #[error(transparent)]
    Bonded(#[from] BondedError),
    #[error(transparent)]
    Nonbonded(#[from] NonbondedError),
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AmberEnergy {
    pub bonded: BondedEnergy,
    pub nonbonded: NonbondedEnergy,
}
impl AmberEnergy {
    pub fn total(self) -> f64 {
        self.bonded.total() + self.nonbonded.total()
    }
}
#[derive(Clone, Debug)]
pub struct AmberEvaluation {
    pub energy: AmberEnergy,
    pub forces: Vec<[f64; 3]>,
}

/// Bonded GAFF/GAFF2 plus cutoff LJ/Coulomb. Partial charges are read from the
/// source; callers choose the charge model (e.g. `ApplyCharges::apply_charges`).
/// Explicit hydrogens and complete bond orders are required. No implicit solvent.
#[derive(Clone, Debug)]
pub struct AmberModel {
    bonded: BondedModel,
    nonbonded: NonbondedModel,
}
impl AmberModel {
    pub fn gaff<S: AtomProvider + BondProvider>(
        source: &S,
        ff: FFType,
        options: NonbondedOptions,
    ) -> Result<Self, AmberError> {
        let bonded = BondedModel::gaff(source, ff)?;
        Self::assemble(source, bonded, AmberParameters::builtin(ff), options)
    }
    pub fn from_typed<S: AtomProvider + BondProvider>(
        source: &S,
        parameters: &AmberParameters,
        policy: ParameterPolicy,
        options: NonbondedOptions,
    ) -> Result<Self, AmberError> {
        let bonded = BondedModel::from_typed(source, parameters, policy)?;
        Self::assemble(source, bonded, parameters, options)
    }
    fn assemble<S: AtomProvider>(
        source: &S,
        bonded: BondedModel,
        parameters: &AmberParameters,
        options: NonbondedOptions,
    ) -> Result<Self, AmberError> {
        let charges: Vec<_> = source.iter_atoms().map(|a| a.get_charge() as f64).collect();
        let nonbonded = NonbondedModel::from_bonded(&bonded, parameters, &charges, options)?;
        Ok(Self { bonded, nonbonded })
    }
    pub fn bonded(&self) -> &BondedModel {
        &self.bonded
    }
    pub fn nonbonded(&self) -> &NonbondedModel {
        &self.nonbonded
    }
    pub fn atom_indices(&self) -> &[usize] {
        self.bonded.atom_indices()
    }
    pub fn evaluate(&self, positions: &[[f64; 3]]) -> Result<AmberEvaluation, EvaluationError> {
        self.evaluate_selected(positions, None)
    }
    /// Omit every interaction entirely within the fixed group. The returned
    /// energy differs from the full energy by the fixed-group constant; movable
    /// atom forces are identical to a full evaluation. Fixed atoms receive only
    /// reaction forces from terms involving movable atoms.
    pub fn evaluate_movable(
        &self,
        positions: &[[f64; 3]],
        movable: &[bool],
    ) -> Result<AmberEvaluation, EvaluationError> {
        self.evaluate_selected(positions, Some(movable))
    }
    fn evaluate_selected(
        &self,
        positions: &[[f64; 3]],
        movable: Option<&[bool]>,
    ) -> Result<AmberEvaluation, EvaluationError> {
        let bonded = match movable {
            Some(m) => self.bonded.evaluate_movable(positions, m)?,
            None => self.bonded.evaluate(positions)?,
        };
        let nonbonded = match movable {
            Some(m) => self.nonbonded.evaluate_movable(positions, m)?,
            None => self.nonbonded.evaluate(positions)?,
        };
        let forces: Vec<[f64; 3]> = bonded
            .forces
            .iter()
            .zip(nonbonded.forces)
            .map(|(b, n)| std::array::from_fn(|i| b[i] + n[i]))
            .collect();
        let energy = AmberEnergy {
            bonded: bonded.energy,
            nonbonded: nonbonded.energy,
        };
        if !energy.total().is_finite() || forces.iter().flatten().any(|f| !f.is_finite()) {
            return Err(EvaluationError::Overflow);
        }
        Ok(AmberEvaluation { energy, forces })
    }
}

/// Reusable energy-only evaluator for poses with the same fixed environment.
pub struct AmberScorer<'a> {
    bonded: super::energy::ActiveBonded<'a>,
    pairs: super::nonbonded::ActivePairs<'a, false>,
    active: Vec<usize>,
    fixed: Vec<(usize, [f64; 3])>,
}

impl AmberModel {
    /// Prepare active bonded terms and index the fixed atoms once for a series of poses.
    pub fn movable_scorer(
        &self,
        positions: &[[f64; 3]],
        movable: &[bool],
    ) -> Result<AmberScorer<'_>, EvaluationError> {
        let n = self.atom_indices().len();
        if positions.len() != n {
            return Err(EvaluationError::Size {
                expected: n,
                actual: positions.len(),
            });
        }
        if movable.len() != n {
            return Err(EvaluationError::MaskSize {
                expected: n,
                actual: movable.len(),
            });
        }
        Ok(AmberScorer {
            bonded: super::energy::ActiveBonded::new(&self.bonded, movable),
            pairs: super::nonbonded::ActivePairs::new(&self.nonbonded, positions, movable)?,
            active: (0..n).filter(|&i| movable[i]).collect(),
            fixed: (0..n)
                .filter(|&i| !movable[i])
                .map(|i| (i, positions[i]))
                .collect(),
        })
    }
}
impl AmberScorer<'_> {
    /// Evaluate a trial in model-local order. Movable atoms may change arbitrarily
    /// between calls. Fixed coordinates must match those supplied at construction.
    pub fn energy(&mut self, positions: &[[f64; 3]]) -> Result<AmberEnergy, EvaluationError> {
        let n = self.active.len() + self.fixed.len();
        if positions.len() != n {
            return Err(EvaluationError::Size {
                expected: n,
                actual: positions.len(),
            });
        }
        if let Some(&(i, _)) = self.fixed.iter().find(|&&(i, p)| positions[i] != p) {
            return Err(EvaluationError::FixedAtomMoved(i));
        }
        let mut result = NonbondedEvaluation {
            energy: Default::default(),
            forces: Vec::new(),
        };
        self.pairs.evaluate(positions, &self.active, &mut result)?;
        let energy = AmberEnergy {
            bonded: self.bonded.energy(positions)?,
            nonbonded: result.energy,
        };
        if !energy.total().is_finite() {
            return Err(EvaluationError::Overflow);
        }
        Ok(energy)
    }
}
