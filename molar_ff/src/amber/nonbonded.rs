//! Cutoff Lennard-Jones and Coulomb interactions in nonperiodic Cartesian space.
use super::{
    AmberParameters, BondedModel, EvaluationError, LennardJonesParameter, ParameterOrigin,
};
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

/// Coulomb constant in kJ mol^-1 nm e^-2.
pub const COULOMB_CONSTANT: f64 = 138.935_457_644_382;

/// LJ uses a quintic potential switch; Coulomb uses a shifted force.
/// Amber 1-4 exceptions retain their full, scaled interaction at all distances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonbondedOptions {
    /// Outer cutoff in nm. No periodic boundaries or long-range correction.
    pub cutoff: f64,
    /// LJ switch starts here in nm; must be positive and less than cutoff.
    pub switch_distance: f64,
}
impl Default for NonbondedOptions {
    fn default() -> Self {
        Self {
            cutoff: 1.2,
            switch_distance: 1.0,
        }
    }
}
impl NonbondedOptions {
    fn validate(self) -> Result<(), NonbondedError> {
        if !self.cutoff.is_finite()
            || !self.switch_distance.is_finite()
            || self.switch_distance <= 0.0
            || self.switch_distance >= self.cutoff
        {
            return Err(NonbondedError::InvalidCutoff);
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum NonbondedError {
    #[error("require finite 0 < switch distance < cutoff (nm)")]
    InvalidCutoff,
    #[error("expected {expected} partial charges, got {actual}")]
    ChargeCount { expected: usize, actual: usize },
    #[error("non-finite partial charge at source atom {0}")]
    NonFiniteCharge(usize),
    #[error("missing Lennard-Jones parameters (source atom, type): {0:?}")]
    MissingParameters(Vec<(usize, String)>),
}
#[derive(Clone, Debug)]
pub struct NonbondedAtom {
    pub parameter: LennardJonesParameter,
    /// Partial charge in elementary-charge units.
    pub charge: f64,
    pub origin: ParameterOrigin,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonbondedException {
    /// Model-local atom indices, in ascending order.
    pub atoms: [usize; 2],
    pub lj_scale: f64,
    pub coulomb_scale: f64,
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NonbondedEnergy {
    pub lennard_jones: f64,
    pub coulomb: f64,
}
impl NonbondedEnergy {
    pub fn total(self) -> f64 {
        self.lennard_jones + self.coulomb
    }
}
#[derive(Clone, Debug)]
pub struct NonbondedEvaluation {
    pub energy: NonbondedEnergy,
    pub forces: Vec<[f64; 3]>,
}

/// Immutable parameters and graph exclusions. Evaluation rebuilds a sparse cell
/// list so moving atoms can never outrun a stale neighbor list.
#[derive(Clone, Debug)]
pub struct NonbondedModel {
    atoms: Vec<NonbondedAtom>,
    exceptions: Vec<NonbondedException>,
    excluded: HashSet<[usize; 2]>,
    options: NonbondedOptions,
}
impl NonbondedModel {
    /// Charges and atoms must be in the bonded model's local order. Distances of
    /// one/two bonds are excluded; shortest distance three gets GAFF's 1/2 LJ
    /// and 1/1.2 Coulomb scales. Multiple torsions never duplicate a 1-4 pair.
    pub fn from_bonded(
        bonded: &BondedModel,
        parameters: &AmberParameters,
        charges: &[f64],
        options: NonbondedOptions,
    ) -> Result<Self, NonbondedError> {
        options.validate()?;
        let n = bonded.atom_indices().len();
        if charges.len() != n {
            return Err(NonbondedError::ChargeCount {
                expected: n,
                actual: charges.len(),
            });
        }
        let mut atoms = Vec::with_capacity(n);
        let mut missing = vec![];
        for (i, (t, &charge)) in bonded.atom_types().iter().zip(charges).enumerate() {
            if !charge.is_finite() {
                return Err(NonbondedError::NonFiniteCharge(bonded.atom_indices()[i]));
            }
            if let Some((parameter, origin)) = parameters.lennard_jones(t) {
                atoms.push(NonbondedAtom {
                    parameter,
                    charge,
                    origin: origin.clone(),
                });
            } else {
                missing.push((bonded.atom_indices()[i], t.clone()));
            }
        }
        if !missing.is_empty() {
            return Err(NonbondedError::MissingParameters(missing));
        }
        let mut adjacent = vec![vec![]; n];
        for bond in bonded.bonds() {
            let [a, b] = bond.atoms;
            adjacent[a].push(b);
            adjacent[b].push(a);
        }
        let mut exceptions = vec![];
        for start in 0..n {
            let mut distances = BTreeMap::from([(start, 0)]);
            let mut queue = VecDeque::from([start]);
            while let Some(atom) = queue.pop_front() {
                let depth = distances[&atom];
                if depth == 3 {
                    continue;
                }
                for &other in &adjacent[atom] {
                    if let std::collections::btree_map::Entry::Vacant(entry) =
                        distances.entry(other)
                    {
                        entry.insert(depth + 1);
                        queue.push_back(other);
                    }
                }
            }
            for (other, depth) in distances {
                if other <= start {
                    continue;
                }
                exceptions.push(NonbondedException {
                    atoms: [start, other],
                    lj_scale: if depth == 3 { 0.5 } else { 0.0 },
                    coulomb_scale: if depth == 3 { 1.0 / 1.2 } else { 0.0 },
                });
            }
        }
        let excluded = exceptions.iter().map(|e| e.atoms).collect();
        Ok(Self {
            atoms,
            exceptions,
            excluded,
            options,
        })
    }
    pub fn atoms(&self) -> &[NonbondedAtom] {
        &self.atoms
    }
    pub fn exceptions(&self) -> &[NonbondedException] {
        &self.exceptions
    }
    pub fn options(&self) -> NonbondedOptions {
        self.options
    }

    /// Energies in kJ/mol and forces in kJ/mol/nm. Coordinates are in nm, local
    /// order, and must already be unwrapped. Singular active pairs return errors.
    pub fn evaluate(&self, positions: &[[f64; 3]]) -> Result<NonbondedEvaluation, EvaluationError> {
        let n = self.atoms.len();
        if positions.len() != n {
            return Err(EvaluationError::Size {
                expected: n,
                actual: positions.len(),
            });
        }
        let mut cells: HashMap<[i32; 3], Vec<usize>> = HashMap::new();
        let mut keys = Vec::with_capacity(n);
        for (i, position) in positions.iter().enumerate() {
            if position.iter().any(|v| !v.is_finite()) {
                return Err(EvaluationError::NonFinite(i));
            }
            let cell = position.map(|v| (v / self.options.cutoff).floor());
            // Leave room for adjacent cell offsets; do not saturate huge coordinates.
            if cell
                .iter()
                .any(|v| *v <= i32::MIN as f64 || *v >= i32::MAX as f64)
            {
                return Err(EvaluationError::Overflow);
            }
            let key = cell.map(|v| v as i32);
            keys.push(key);
            cells.entry(key).or_default().push(i);
        }
        let mut result = NonbondedEvaluation {
            energy: NonbondedEnergy::default(),
            forces: vec![[0.0; 3]; n],
        };
        for (i, &[x, y, z]) in keys.iter().enumerate() {
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(neighbors) = cells.get(&[x + dx, y + dy, z + dz]) {
                            for &j in neighbors {
                                if j <= i || self.excluded.contains(&[i, j]) {
                                    continue;
                                }
                                self.add_pair([i, j], positions, 1.0, 1.0, true, &mut result)?;
                            }
                        }
                    }
                }
            }
        }
        // 1-4 exceptions are evaluated separately and are NOT subject to the cutoff.
        for e in &self.exceptions {
            if e.lj_scale != 0.0 || e.coulomb_scale != 0.0 {
                self.add_pair(
                    e.atoms,
                    positions,
                    e.lj_scale,
                    e.coulomb_scale,
                    false,
                    &mut result,
                )?;
            }
        }
        if !result.energy.total().is_finite()
            || result.forces.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(EvaluationError::Overflow);
        }
        Ok(result)
    }
    fn add_pair(
        &self,
        [i, j]: [usize; 2],
        positions: &[[f64; 3]],
        lj_scale: f64,
        coulomb_scale: f64,
        switched: bool,
        result: &mut NonbondedEvaluation,
    ) -> Result<(), EvaluationError> {
        let a = &self.atoms[i];
        let b = &self.atoms[j];
        let d: [f64; 3] = std::array::from_fn(|axis| positions[i][axis] - positions[j][axis]);
        let r2 = d.iter().map(|v| v * v).sum::<f64>();
        if switched && r2 >= self.options.cutoff.powi(2) {
            return Ok(());
        }
        let epsilon = a.parameter.epsilon.sqrt() * b.parameter.epsilon.sqrt() * lj_scale;
        let qq = COULOMB_CONSTANT * a.charge * b.charge * coulomb_scale;
        if epsilon == 0.0 && qq == 0.0 {
            return Ok(());
        }
        if r2 < 1e-24 {
            return Err(EvaluationError::CoincidentNonbonded { atoms: [i, j] });
        }
        let r = r2.sqrt();
        let ratio6 = ((a.parameter.radius + b.parameter.radius) / r).powi(6);
        let lj = epsilon * (ratio6 * ratio6 - 2.0 * ratio6);
        let lj_force = 12.0 * epsilon * (ratio6 * ratio6 - ratio6) / r;
        let (s, ds) = if switched {
            switch(r, self.options)
        } else {
            (1.0, 0.0)
        };
        let (coulomb, coulomb_force) = if switched {
            let rc = self.options.cutoff;
            (
                qq * (1.0 / r - 1.0 / rc + (r - rc) / (rc * rc)),
                qq * (1.0 / r2 - 1.0 / (rc * rc)),
            )
        } else {
            (qq / r, qq / r2)
        };
        result.energy.lennard_jones += s * lj;
        result.energy.coulomb += coulomb;
        let coefficient = (s * lj_force - ds * lj + coulomb_force) / r;
        for (axis, delta) in d.into_iter().enumerate() {
            let f = coefficient * delta;
            result.forces[i][axis] += f;
            result.forces[j][axis] -= f;
        }
        Ok(())
    }
}
fn switch(r: f64, options: NonbondedOptions) -> (f64, f64) {
    if r <= options.switch_distance {
        return (1.0, 0.0);
    }
    if r >= options.cutoff {
        return (0.0, 0.0);
    }
    let width = options.cutoff - options.switch_distance;
    let x = (r - options.switch_distance) / width;
    (
        1.0 - 10.0 * x.powi(3) + 15.0 * x.powi(4) - 6.0 * x.powi(5),
        -30.0 * x.powi(2) * (1.0 - x).powi(2) / width,
    )
}
