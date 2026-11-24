//! Preconditioned L-BFGS with a displacement-limited Wolfe line search.
use super::{AmberEnergy, AmberModel, EvaluationError, NonbondedEvaluation};
use std::collections::VecDeque;
type V = [f64; 3];

#[derive(Clone, Copy, Debug)]
pub struct MinimizeOptions {
    pub max_steps: usize,
    /// Maximum force on any movable atom, kJ/mol/nm.
    pub force_tolerance: f64,
    /// Maximum atom displacement for a trial step, nm.
    pub max_displacement: f64,
    /// Number of L-BFGS correction pairs (1..=64).
    pub history: usize,
}
impl Default for MinimizeOptions {
    fn default() -> Self {
        Self {
            max_steps: 1000,
            force_tolerance: 10.0,
            max_displacement: 0.01,
            history: 8,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StopReason {
    Converged,
    StepLimit,
    LineSearchStalled,
}
#[derive(Clone, Debug)]
pub struct MinimizeResult {
    pub steps: usize,
    /// Includes the initial evaluation and rejected line-search trials.
    pub evaluations: usize,
    pub neighbor_rebuilds: usize,
    pub skipped_updates: usize,
    pub initial_energy: AmberEnergy,
    pub final_energy: AmberEnergy,
    pub max_force: f64,
    pub reason: StopReason,
}
impl MinimizeResult {
    pub fn converged(&self) -> bool {
        self.reason == StopReason::Converged
    }
}
#[derive(Debug, thiserror::Error)]
pub enum MinimizeError {
    #[error("invalid minimizer options")]
    InvalidOptions,
    #[error("expected {expected} movable flags, got {actual}")]
    MaskSize { expected: usize, actual: usize },
    #[error(transparent)]
    Evaluation(#[from] EvaluationError),
}
#[derive(Default)]
struct Correction {
    s: Vec<V>,
    y: Vec<V>,
    rho: f64,
    scale: f64,
}
fn dot(a: &[V], b: &[V]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (0..3).map(|i| a[i] * b[i]).sum::<f64>())
        .sum()
}
fn norm(v: &V) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn max_norm(v: &[V]) -> f64 {
    v.iter().map(norm).fold(0.0, f64::max)
}
fn axpy(target: &mut [V], x: &[V], a: f64) {
    for (target, x) in target.iter_mut().zip(x) {
        for i in 0..3 {
            target[i] += a * x[i];
        }
    }
}
/// Reusable storage for independent Amber minimizations of one model. Each call
/// accepts a new coordinate set and movable mask and starts a fresh L-BFGS history.
pub struct AmberMinimizer<'a> {
    model: &'a AmberModel,
    pairs: super::nonbonded::ActivePairs<'a>,
    bonded: super::energy::ActiveBonded<'a>,
    scratch: NonbondedEvaluation,
    buffers: MinimizeBuffers,
}

#[derive(Default)]
struct MinimizeBuffers {
    active: Vec<usize>,
    trial: Vec<V>,
    gradient: Vec<V>,
    next_gradient: Vec<V>,
    direction: Vec<V>,
    s: Vec<V>,
    y: Vec<V>,
    stiffness: Vec<f64>,
    diagonal: Vec<f64>,
    alphas: Vec<f64>,
    history: VecDeque<Correction>,
    spare: Vec<Correction>,
}

impl AmberModel {
    /// Allocate reusable storage for a sequence of minimizations.
    pub fn minimizer(&self) -> AmberMinimizer<'_> {
        AmberMinimizer {
            model: self,
            pairs: super::nonbonded::ActivePairs::empty(self.nonbonded()),
            bonded: Default::default(),
            scratch: NonbondedEvaluation {
                energy: Default::default(),
                forces: Vec::new(),
            },
            buffers: Default::default(),
        }
    }
    /// Minimize coordinates in place. Only marked atoms move. Interactions entirely within the fixed
    /// group are excluded from evaluation and the reported energies. Every accepted step passes the Armijo energy
    /// decrease test. Failed trial evaluations are rejected, never sanitized to
    /// zero force. The returned force and energy describe the final coordinates.
    pub fn minimize(
        &self,
        positions: &mut [V],
        movable: &[bool],
        options: MinimizeOptions,
    ) -> Result<MinimizeResult, MinimizeError> {
        self.minimize_observed(positions, movable, options, |_, _, _| {})
    }

    /// Observe accepted iterates without restarting the minimizer. The callback
    /// receives the accepted step number, full coordinates and masked energy.
    pub fn minimize_observed(
        &self,
        positions: &mut [V],
        movable: &[bool],
        options: MinimizeOptions,
        accepted_step: impl FnMut(usize, &[V], AmberEnergy),
    ) -> Result<MinimizeResult, MinimizeError> {
        self.minimizer()
            .minimize_observed(positions, movable, options, accepted_step)
    }
}

impl AmberMinimizer<'_> {
    /// Run an independent minimization and observe each accepted step.
    pub fn minimize_observed(
        &mut self,
        positions: &mut [V],
        movable: &[bool],
        options: MinimizeOptions,
        mut accepted_step: impl FnMut(usize, &[V], AmberEnergy),
    ) -> Result<MinimizeResult, MinimizeError> {
        if !options.force_tolerance.is_finite()
            || options.force_tolerance <= 0.0
            || !options.max_displacement.is_finite()
            || options.max_displacement <= 0.0
            || options.history == 0
            || options.history > 64
        {
            return Err(MinimizeError::InvalidOptions);
        }
        if movable.len() != positions.len() {
            return Err(MinimizeError::MaskSize {
                expected: positions.len(),
                actual: movable.len(),
            });
        }
        if self.model.atom_indices().len() != positions.len() {
            return Err(EvaluationError::Size {
                expected: self.model.atom_indices().len(),
                actual: positions.len(),
            }
            .into());
        }
        let MinimizeBuffers {
            active,
            trial,
            gradient: g,
            next_gradient: next_g,
            direction,
            s,
            y,
            stiffness,
            diagonal,
            alphas,
            history,
            spare,
        } = &mut self.buffers;
        active.clear();
        active.extend((0..positions.len()).filter(|&i| movable[i]));
        trial.clear();
        trial.extend_from_slice(positions);
        for buffer in [&mut *g, &mut *next_g, &mut *direction, &mut *s, &mut *y] {
            buffer.resize(active.len(), [0.; 3]);
        }
        alphas.resize(options.history, 0.);
        spare.extend(history.drain(..));
        let bonded = &mut self.bonded;
        let pairs = &mut self.pairs;
        let scratch = &mut self.scratch;
        bonded.reset(self.model.bonded(), movable);
        pairs.reset(positions, movable)?;
        scratch.forces.resize(positions.len(), [0.; 3]);
        let mut evaluations = 0;
        let mut evaluate = |p: &[V], gradient: &mut [V]| -> Result<AmberEnergy, EvaluationError> {
            evaluations += 1;
            scratch.forces.fill([0.; 3]);
            scratch.energy = Default::default();
            let bonded_energy = bonded.evaluate(p, &mut scratch.forces)?;
            pairs.evaluate(p, active, scratch)?;
            let energy = AmberEnergy {
                bonded: bonded_energy,
                nonbonded: scratch.energy,
            };
            if !energy.total().is_finite()
                || scratch.forces.iter().flatten().any(|v| !v.is_finite())
            {
                return Err(EvaluationError::Overflow);
            }
            for (out, &i) in gradient.iter_mut().zip(active.iter()) {
                *out = scratch.forces[i].map(|f| -f);
            }
            Ok(energy)
        };
        let initial_energy = evaluate(positions, g)?;
        let mut energy = initial_energy;
        // Isotropic per-atom approximation to the bond-stretch Hessian.
        // Normalize by the mean stiffness, leaving an unbonded system at unit scale.
        stiffness.clear();
        stiffness.resize(positions.len(), 0.);
        for bond in self.model.bonded().bonds() {
            for i in bond.atoms {
                if movable[i] {
                    stiffness[i] += 2. * bond.parameter.k / 3.;
                }
            }
        }
        let mean = active.iter().map(|&i| stiffness[i]).sum::<f64>() / active.len().max(1) as f64;
        diagonal.clear();
        diagonal.extend(active.iter().map(|&i| {
            if mean > 0. && stiffness[i] > 0. {
                (mean / stiffness[i]).clamp(0.1, 10.)
            } else {
                1.
            }
        }));
        let mut skipped_updates = 0;
        let mut steps = 0;
        let reason = loop {
            if max_norm(g) <= options.force_tolerance {
                break StopReason::Converged;
            }
            if steps >= options.max_steps {
                break StopReason::StepLimit;
            }
            // Standard L-BFGS two-loop recursion: H_k * gradient.
            direction.copy_from_slice(g);
            for (i, pair) in history.iter().enumerate().rev() {
                let alpha = pair.rho * dot(&pair.s, direction);
                alphas[i] = alpha;
                axpy(direction, &pair.y, -alpha);
            }
            let gamma = history.back().map_or(1., |p| p.scale);
            for (v, weight) in direction.iter_mut().zip(diagonal.iter()) {
                for x in v {
                    *x *= gamma * weight;
                }
            }
            for (i, pair) in history.iter().enumerate() {
                let beta = pair.rho * dot(&pair.y, direction);
                axpy(direction, &pair.s, alphas[i] - beta);
            }
            for v in direction.iter_mut() {
                for x in v {
                    *x = -*x;
                }
            }
            let mut slope = dot(g, direction);
            if !slope.is_finite()
                || slope >= 0.0
                || direction.iter().flatten().any(|v| !v.is_finite())
            {
                spare.extend(history.drain(..));
                for (d, v) in direction.iter_mut().zip(g.iter()) {
                    *d = v.map(|x| -x);
                }
                slope = -dot(g, g);
            }
            let mut search = |direction: &[V], slope: f64| {
                let cap = options.max_displacement / max_norm(direction);
                let mut alpha = cap.min(1.);
                let mut lower = 0.;
                let mut upper = f64::INFINITY;
                for _ in 0..40 {
                    let mut changed = false;
                    for (&i, d) in active.iter().zip(direction.iter()) {
                        trial[i] = std::array::from_fn(|k| positions[i][k] + alpha * d[k]);
                        changed |= trial[i] != positions[i];
                    }
                    if !changed {
                        break;
                    }
                    if let Ok(value) = evaluate(trial, next_g) {
                        if value.total() < energy.total()
                            && value.total() <= energy.total() + 1e-4 * alpha * slope
                        {
                            // Accept a capped Armijo step when the curvature condition
                            // would require exceeding the displacement limit.
                            if alpha == cap || dot(next_g, direction) >= 0.9 * slope {
                                return Some(value);
                            }
                            lower = alpha;
                            alpha = if upper.is_finite() {
                                0.5 * (lower + upper)
                            } else {
                                (2. * alpha).min(cap)
                            };
                            continue;
                        }
                    }
                    upper = alpha;
                    alpha = 0.5 * (lower + upper);
                }
                None
            };
            let accepted = search(direction, slope).or_else(|| {
                // A poor inverse Hessian must not prevent a valid steepest-descent step.
                if direction
                    .iter()
                    .zip(g.iter())
                    .all(|(d, g)| (0..3).all(|k| d[k] == -g[k]))
                {
                    return None;
                }
                spare.extend(history.drain(..));
                for (d, v) in direction.iter_mut().zip(g.iter()) {
                    *d = v.map(|x| -x);
                }
                search(direction, -dot(g, g))
            });
            let Some(next) = accepted else {
                break StopReason::LineSearchStalled;
            };

            for (slot, &i) in active.iter().enumerate() {
                s[slot] = std::array::from_fn(|k| trial[i][k] - positions[i][k]);
                y[slot] = std::array::from_fn(|k| next_g[slot][k] - g[slot][k]);
            }
            let sy = dot(s, y);
            let y_metric: f64 = y
                .iter()
                .zip(diagonal.iter())
                .map(|(v, weight)| weight * v.iter().map(|x| x * x).sum::<f64>())
                .sum();
            let scale = sy / y_metric;
            if scale.is_finite()
                && scale > 0.
                && (1. / sy).is_finite()
                && sy.is_finite()
                && sy > 1e-10 * dot(s, s).sqrt() * dot(y, y).sqrt()
                && sy > 0.0
            {
                let mut entry = if history.len() == options.history {
                    history.pop_front().unwrap()
                } else {
                    spare.pop().unwrap_or_default()
                };
                entry.s.resize(active.len(), [0.; 3]);
                entry.y.resize(active.len(), [0.; 3]);
                std::mem::swap(&mut entry.s, s);
                std::mem::swap(&mut entry.y, y);
                entry.rho = 1. / sy;
                entry.scale = scale;
                history.push_back(entry);
            } else {
                skipped_updates += 1;
            }
            for &i in active.iter() {
                positions[i] = trial[i];
            }
            energy = next;
            std::mem::swap(g, next_g);
            steps += 1;
            accepted_step(steps, positions, energy);
        };
        Ok(MinimizeResult {
            steps,
            initial_energy,
            final_energy: energy,
            evaluations,
            neighbor_rebuilds: pairs.rebuilds,
            skipped_updates,
            max_force: max_norm(g),
            reason,
        })
    }
}
