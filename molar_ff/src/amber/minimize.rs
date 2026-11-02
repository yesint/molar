//! Limited-memory BFGS with an energy-decreasing backtracking line search.
use super::{AmberEnergy, AmberEvaluation, AmberModel, EvaluationError};
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
struct Correction {
    s: Vec<V>,
    y: Vec<V>,
    rho: f64,
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
fn gradient(evaluation: &AmberEvaluation, movable: &[bool]) -> Vec<V> {
    evaluation
        .forces
        .iter()
        .zip(movable)
        .map(|(f, &m)| if m { f.map(|v| -v) } else { [0.0; 3] })
        .collect()
}
fn axpy(target: &mut [V], x: &[V], a: f64) {
    for (target, x) in target.iter_mut().zip(x) {
        for i in 0..3 {
            target[i] += a * x[i];
        }
    }
}
impl AmberModel {
    /// Minimize coordinates in place. Only marked atoms move; all atoms still
    /// contribute to the energy. Every accepted step passes the Armijo energy
    /// decrease test. Failed trial evaluations are rejected, never sanitized to
    /// zero force. The returned force and energy describe the final coordinates.
    pub fn minimize(
        &self,
        positions: &mut [V],
        movable: &[bool],
        options: MinimizeOptions,
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
        let mut evaluation = self.evaluate(positions)?;
        let initial_energy = evaluation.energy;
        let mut g = gradient(&evaluation, movable);
        let mut history: VecDeque<Correction> = VecDeque::new();
        let mut steps = 0;
        let reason = loop {
            if max_norm(&g) <= options.force_tolerance {
                break StopReason::Converged;
            }
            if steps >= options.max_steps {
                break StopReason::StepLimit;
            }
            // Standard L-BFGS two-loop recursion: H_k * gradient.
            let mut direction = g.clone();
            let mut alphas = vec![0.0; history.len()];
            for (i, pair) in history.iter().enumerate().rev() {
                let alpha = pair.rho * dot(&pair.s, &direction);
                alphas[i] = alpha;
                axpy(&mut direction, &pair.y, -alpha);
            }
            let gamma = history
                .back()
                .map_or(1.0, |p| 1.0 / (p.rho * dot(&p.y, &p.y)));
            for v in &mut direction {
                for x in v {
                    *x *= gamma;
                }
            }
            for (i, pair) in history.iter().enumerate() {
                let beta = pair.rho * dot(&pair.y, &direction);
                axpy(&mut direction, &pair.s, alphas[i] - beta);
            }
            for v in &mut direction {
                for x in v {
                    *x = -*x;
                }
            }
            let mut slope = dot(&g, &direction);
            if !slope.is_finite()
                || slope >= 0.0
                || direction.iter().flatten().any(|v| !v.is_finite())
            {
                history.clear();
                direction = g.iter().map(|v| v.map(|x| -x)).collect();
                slope = -dot(&g, &g);
            }
            let search = |direction: &[V], slope: f64| {
                let mut alpha = (options.max_displacement / max_norm(direction)).min(1.0);
                for _ in 0..40 {
                    let trial: Vec<V> = positions
                        .iter()
                        .zip(direction)
                        .zip(movable)
                        .map(|((p, d), &m)| {
                            if m {
                                std::array::from_fn(|i| p[i] + alpha * d[i])
                            } else {
                                *p
                            }
                        })
                        .collect();
                    if let Ok(value) = self.evaluate(&trial) {
                        if value.energy.total() < evaluation.energy.total()
                            && value.energy.total()
                                <= evaluation.energy.total() + 1e-4 * alpha * slope
                        {
                            return Some((trial, value));
                        }
                    }
                    alpha *= 0.5;
                }
                None
            };
            let accepted = search(&direction, slope).or_else(|| {
                // A poor inverse Hessian must not prevent a valid steepest-descent step.
                history.clear();
                let steepest: Vec<_> = g.iter().map(|v| v.map(|x| -x)).collect();
                search(&steepest, -dot(&g, &g))
            });
            let Some((trial, next)) = accepted else {
                break StopReason::LineSearchStalled;
            };
            let next_g = gradient(&next, movable);
            let s: Vec<V> = trial
                .iter()
                .zip(positions.iter())
                .map(|(a, b)| std::array::from_fn(|i| a[i] - b[i]))
                .collect();
            let y: Vec<V> = next_g
                .iter()
                .zip(&g)
                .map(|(a, b)| std::array::from_fn(|i| a[i] - b[i]))
                .collect();
            let sy = dot(&s, &y);
            if sy.is_finite() && sy > 1e-10 * dot(&s, &s).sqrt() * dot(&y, &y).sqrt() && sy > 0.0 {
                if history.len() == options.history {
                    history.pop_front();
                }
                history.push_back(Correction {
                    s,
                    y,
                    rho: 1.0 / sy,
                });
            }
            positions.copy_from_slice(&trial);
            evaluation = next;
            g = next_g;
            steps += 1;
        };
        Ok(MinimizeResult {
            steps,
            initial_energy,
            final_energy: evaluation.energy,
            max_force: max_norm(&g),
            reason,
        })
    }
}
