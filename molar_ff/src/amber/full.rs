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
        let bonded = self.bonded.evaluate(positions)?;
        let nonbonded = self.nonbonded.evaluate(positions)?;
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
