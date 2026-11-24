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
    kernel_atoms: Vec<PairAtom>,
    exceptions: Vec<NonbondedException>,
    excluded: HashSet<[usize; 2]>,
    options: NonbondedOptions,
}
// Pack the three scalars used by the inner pair loop next to one another.
#[derive(Clone, Copy, Debug)]
struct PairAtom {
    radius: f64,
    sqrt_epsilon: f64,
    charge: f64,
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
        let kernel_atoms = atoms
            .iter()
            .map(|a| PairAtom {
                radius: a.parameter.radius,
                sqrt_epsilon: a.parameter.epsilon.sqrt(),
                charge: a.charge,
            })
            .collect();
        Ok(Self {
            kernel_atoms,
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
        self.evaluate_selected(positions, None)
    }
    /// Evaluate moving–moving and moving–fixed interactions only. Neighbor searches
    /// start from movable atoms, so fixed–fixed pairs are never enumerated.
    /// Fixed-only 1–4 exceptions are also skipped before geometry evaluation.
    pub fn evaluate_movable(
        &self,
        positions: &[[f64; 3]],
        movable: &[bool],
    ) -> Result<NonbondedEvaluation, EvaluationError> {
        if movable.len() != self.atoms.len() {
            return Err(EvaluationError::MaskSize {
                expected: self.atoms.len(),
                actual: movable.len(),
            });
        }
        self.evaluate_selected(positions, Some(movable))
    }
    fn evaluate_selected(
        &self,
        positions: &[[f64; 3]],
        movable: Option<&[bool]>,
    ) -> Result<NonbondedEvaluation, EvaluationError> {
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
            if movable.is_some_and(|mask| !mask[i]) {
                continue;
            }
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(neighbors) = cells.get(&[x + dx, y + dy, z + dz]) {
                            for &j in neighbors {
                                let pair = [i.min(j), i.max(j)];
                                if j == i
                                    || (j < i && movable.is_none_or(|mask| mask[j]))
                                    || self.excluded.contains(&pair)
                                {
                                    continue;
                                }
                                self.add_pair::<true>(
                                    pair,
                                    positions,
                                    1.0,
                                    1.0,
                                    true,
                                    &mut result,
                                )?;
                            }
                        }
                    }
                }
            }
        }
        // 1-4 exceptions are evaluated separately and are NOT subject to the cutoff.
        for e in &self.exceptions {
            if movable.is_some_and(|mask| e.atoms.iter().all(|&i| !mask[i])) {
                continue;
            }
            if e.lj_scale != 0.0 || e.coulomb_scale != 0.0 {
                self.add_pair::<true>(
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
    fn add_pair<const FORCES: bool>(
        &self,
        [i, j]: [usize; 2],
        positions: &[[f64; 3]],
        lj_scale: f64,
        coulomb_scale: f64,
        switched: bool,
        result: &mut NonbondedEvaluation,
    ) -> Result<(), EvaluationError> {
        let a = &self.kernel_atoms[i];
        let b = &self.kernel_atoms[j];
        let d: [f64; 3] = std::array::from_fn(|axis| positions[i][axis] - positions[j][axis]);
        let r2 = d.iter().map(|v| v * v).sum::<f64>();
        if switched && r2 >= self.options.cutoff.powi(2) {
            return Ok(());
        }
        let epsilon = a.sqrt_epsilon * b.sqrt_epsilon * lj_scale;
        let qq = COULOMB_CONSTANT * a.charge * b.charge * coulomb_scale;
        if epsilon == 0.0 && qq == 0.0 {
            return Ok(());
        }
        if r2 < 1e-24 {
            return Err(EvaluationError::CoincidentNonbonded { atoms: [i, j] });
        }
        let r = r2.sqrt();
        let ratio6 = ((a.radius + b.radius) / r).powi(6);
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
        if !FORCES {
            return Ok(());
        }
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
// Only the movable coordinates may change during this workspace's lifetime.
// Build at cutoff + skin; two atoms each travelling half the skin cannot
// bring an unlisted pair inside the physical cutoff.
pub(super) struct ActivePairs<'a, const FORCES: bool = true> {
    model: &'a NonbondedModel,
    fixed: FixedGrid,
    ordinary: Vec<[usize; 2]>,
    scaled: Vec<NonbondedException>,
    reference: Vec<[f64; 3]>,
    fixed_rows: Vec<FixedNeighbors>,
    moving_heads: HashMap<[i32; 3], usize>,
    next_in_cell: Vec<usize>,
    pub(super) rebuilds: usize,
}
// Keep cell buffers when coordinates or the movable mask change between runs.
#[derive(Default)]
struct FixedGrid {
    cells: HashMap<[i32; 3], usize>,
    bins: Vec<Vec<usize>>,
}

impl FixedGrid {
    fn clear(&mut self) {
        self.cells.clear();
    }

    fn insert(&mut self, cell: [i32; 3], atom: usize) {
        let next = self.cells.len();
        let index = *self.cells.entry(cell).or_insert_with(|| {
            if next == self.bins.len() {
                self.bins.push(Vec::new());
            }
            self.bins[next].clear();
            next
        });
        self.bins[index].push(atom);
    }

    fn get(&self, cell: &[i32; 3]) -> Option<&[usize]> {
        self.cells.get(cell).map(|&i| self.bins[i].as_slice())
    }
}

// Only one atom in each of these pairs can move, so the full skin is safe.
#[derive(Default)]
struct FixedNeighbors {
    center: Option<[f64; 3]>,
    neighbors: Vec<usize>,
    cached: Option<([f64; 3], NonbondedEnergy)>,
}

impl FixedNeighbors {
    fn energy(
        &mut self,
        i: usize,
        positions: &[[f64; 3]],
        model: &NonbondedModel,
        grid: &FixedGrid,
    ) -> Result<NonbondedEnergy, EvaluationError> {
        let point = positions[i];
        if let Some((previous, energy)) = self.cached
            && previous == point
        {
            return Ok(energy);
        }
        let width = model.options.cutoff + SKIN;
        if self
            .center
            .is_none_or(|c| squared_distance(point, c) >= SKIN * SKIN)
        {
            let [x, y, z] = key(point, width, i)?;
            self.neighbors.clear();
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(atoms) = grid.get(&[x + dx, y + dy, z + dz]) {
                            for &j in atoms {
                                if squared_distance(point, positions[j]) <= width * width
                                    && !model.excluded.contains(&[i.min(j), i.max(j)])
                                {
                                    self.neighbors.push(j);
                                }
                            }
                        }
                    }
                }
            }
            self.center = Some(point);
        }
        let mut sum = NonbondedEvaluation {
            energy: Default::default(),
            forces: Vec::new(),
        };
        for &j in &self.neighbors {
            model.add_pair::<false>([i.min(j), i.max(j)], positions, 1., 1., true, &mut sum)?;
        }
        self.cached = Some((point, sum.energy));
        Ok(sum.energy)
    }
}

fn squared_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).powi(2)).sum()
}

const SKIN: f64 = 0.1;

fn key(p: [f64; 3], width: f64, index: usize) -> Result<[i32; 3], EvaluationError> {
    if p.iter().any(|x| !x.is_finite()) {
        return Err(EvaluationError::NonFinite(index));
    }
    let cell = p.map(|x| (x / width).floor());
    if cell
        .iter()
        .any(|x| *x <= i32::MIN as f64 || *x >= i32::MAX as f64)
    {
        return Err(EvaluationError::Overflow);
    }
    Ok(cell.map(|x| x as i32))
}
impl<'a, const FORCES: bool> ActivePairs<'a, FORCES> {
    pub(super) fn new(
        model: &'a NonbondedModel,
        positions: &[[f64; 3]],
        mask: &[bool],
    ) -> Result<Self, EvaluationError> {
        let mut pairs = Self::empty(model);
        pairs.reset(positions, mask)?;
        Ok(pairs)
    }
    pub(super) fn empty(model: &'a NonbondedModel) -> Self {
        Self {
            model,
            fixed: FixedGrid::default(),
            ordinary: Vec::new(),
            scaled: Vec::new(),
            reference: Vec::new(),
            fixed_rows: Vec::new(),
            moving_heads: HashMap::new(),
            next_in_cell: Vec::new(),
            rebuilds: 0,
        }
    }
    pub(super) fn reset(
        &mut self,
        positions: &[[f64; 3]],
        mask: &[bool],
    ) -> Result<(), EvaluationError> {
        self.reference.clear();
        self.ordinary.clear();
        self.scaled.clear();
        self.fixed_rows.clear();
        self.fixed.clear();
        self.rebuilds = 0;
        for (i, &p) in positions.iter().enumerate() {
            key(p, self.model.options.cutoff, i)?;
            if !mask[i] {
                self.fixed
                    .insert(key(p, self.model.options.cutoff + SKIN, i)?, i);
            }
        }
        self.scaled
            .extend(self.model.exceptions.iter().copied().filter(|e| {
                (e.lj_scale != 0. || e.coulomb_scale != 0.) && e.atoms.iter().any(|&i| mask[i])
            }));
        Ok(())
    }
    pub(super) fn evaluate(
        &mut self,
        p: &[[f64; 3]],
        active: &[usize],
        result: &mut NonbondedEvaluation,
    ) -> Result<(), EvaluationError> {
        let width = self.model.options.cutoff + SKIN;
        let stale = self.reference.len() != active.len()
            || active.iter().zip(&self.reference).any(|(&i, old)| {
                let d2 = squared_distance(p[i], *old);
                !d2.is_finite() || d2 >= (0.5 * SKIN).powi(2)
            });
        if stale {
            self.reference.clear();
            self.moving_heads.clear();
            self.next_in_cell.resize(active.len(), usize::MAX);
            // Insert in reverse to visit moving atoms in their original order.
            for (slot, &i) in active.iter().enumerate().rev() {
                key(p[i], self.model.options.cutoff, i)?;
                let cell = key(p[i], width, i)?;
                self.next_in_cell[slot] =
                    self.moving_heads.insert(cell, slot).unwrap_or(usize::MAX);
            }
            self.ordinary.clear();
            for &i in active {
                let [x, y, z] = key(p[i], width, i)?;
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for dz in -1..=1 {
                            let cell = [x + dx, y + dy, z + dz];
                            let fixed = if FORCES { self.fixed.get(&cell) } else { None };
                            let fixed = fixed.into_iter().flatten();
                            let moving = std::iter::successors(
                                self.moving_heads.get(&cell).copied(),
                                |&slot| {
                                    let next = self.next_in_cell[slot];
                                    (next != usize::MAX).then_some(next)
                                },
                            )
                            .map(|slot| active[slot])
                            .filter(|&j| j > i);
                            for j in fixed.copied().chain(moving) {
                                let pair = [i.min(j), i.max(j)];
                                let d2 = squared_distance(p[i], p[j]);
                                if d2 <= width * width && !self.model.excluded.contains(&pair) {
                                    self.ordinary.push(pair);
                                }
                            }
                        }
                    }
                }
            }
            self.reference.extend(active.iter().map(|&i| p[i]));
            self.rebuilds += 1;
        }
        // Unchanged coordinates retain the same energy against a fixed environment.
        if !FORCES {
            self.fixed_rows
                .resize_with(active.len(), FixedNeighbors::default);
            for (row, &i) in self.fixed_rows.iter_mut().zip(active) {
                let energy = row.energy(i, p, self.model, &self.fixed)?;
                result.energy.lennard_jones += energy.lennard_jones;
                result.energy.coulomb += energy.coulomb;
            }
        }
        for &pair in &self.ordinary {
            self.model
                .add_pair::<FORCES>(pair, p, 1., 1., true, result)?;
        }
        for e in &self.scaled {
            self.model.add_pair::<FORCES>(
                e.atoms,
                p,
                e.lj_scale,
                e.coulomb_scale,
                false,
                result,
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod cached_pairs_tests {
    use super::*;
    use crate::{FFType, amber::AmberModel};
    use molar::prelude::*;

    #[test]
    fn rejected_trials_and_masked_motion_match_the_reference() {
        let sys = System::from_file("tests/data/gaff_ref/sdf/Ethanol.sdf").unwrap();
        let model = AmberModel::gaff(&sys, FFType::Gaff2, NonbondedOptions::default()).unwrap();
        let p: Vec<_> = sys
            .iter_pos()
            .map(|q| [q.x as f64, q.y as f64, q.z as f64])
            .collect();
        let active = [0, 3];
        let mask: Vec<_> = (0..p.len()).map(|i| active.contains(&i)).collect();
        let nb = model.nonbonded();
        let mut cached = ActivePairs::<true>::new(nb, &p, &mask).unwrap();
        for shift in [0., 0.001, 1.25, 0.002, -1.30, 0.0] {
            let mut trial = p.clone();
            trial[0][0] += shift;
            trial[3][0] -= shift;
            let mut actual = NonbondedEvaluation {
                energy: Default::default(),
                forces: vec![[0.; 3]; p.len()],
            };
            cached.evaluate(&trial, &active, &mut actual).unwrap();
            let reference = nb.evaluate_movable(&trial, &mask).unwrap();
            let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * (1. + b.abs());
            assert!(close(actual.energy.total(), reference.energy.total()));
            for (a, b) in actual
                .forces
                .iter()
                .flatten()
                .zip(reference.forces.iter().flatten())
            {
                assert!(close(*a, *b));
            }
        }
        assert_eq!(cached.rebuilds, 5);
    }
}
