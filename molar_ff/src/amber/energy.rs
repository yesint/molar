use super::{BondedModel, TermKind};
type V = [f64; 3];
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn mul(a: V, s: f64) -> V {
    a.map(|v| v * s)
}
fn dot(a: V, b: V) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BondedEnergy {
    pub bonds: f64,
    pub angles: f64,
    pub propers: f64,
    pub impropers: f64,
}
impl BondedEnergy {
    pub fn total(&self) -> f64 {
        self.bonds + self.angles + self.propers + self.impropers
    }
}
#[derive(Clone, Debug)]
pub struct BondedEvaluation {
    pub energy: BondedEnergy,
    pub forces: Vec<V>,
}
#[derive(Clone, Debug, thiserror::Error, PartialEq)]
pub enum EvaluationError {
    #[error("expected {expected} positions, got {actual}")]
    Size { expected: usize, actual: usize },
    #[error("expected {expected} movable flags, got {actual}")]
    MaskSize { expected: usize, actual: usize },
    #[error("non-finite coordinate at local atom {0}")]
    NonFinite(usize),
    #[error("undefined {kind:?} geometry at local atoms {atoms:?}")]
    Degenerate { kind: TermKind, atoms: Vec<usize> },
    #[error("coincident nonbonded atoms {atoms:?}")]
    CoincidentNonbonded { atoms: [usize; 2] },
    #[error("fixed atom {0} moved after the scoring workspace was created")]
    FixedAtomMoved(usize),
    #[error("force-field evaluation overflowed")]
    Overflow,
}
impl BondedModel {
    /// Evaluate nonperiodic Cartesian bonded energy and analytic forces. Coordinates
    /// must be in model-local order (see `atom_indices`) and already unwrapped.
    /// Invalid/undefined geometry returns an error, not a partial force array.
    pub fn evaluate(&self, positions: &[V]) -> Result<BondedEvaluation, EvaluationError> {
        self.evaluate_selected(positions, None)
    }
    /// Evaluate only terms containing at least one movable atom. Fixed-only terms
    /// are skipped before any geometry calculation, including degenerate terms.
    /// Forces include reactions on fixed atoms from the retained interactions.
    pub fn evaluate_movable(
        &self,
        positions: &[V],
        movable: &[bool],
    ) -> Result<BondedEvaluation, EvaluationError> {
        if movable.len() != self.atom_indices.len() {
            return Err(EvaluationError::MaskSize {
                expected: self.atom_indices.len(),
                actual: movable.len(),
            });
        }
        self.evaluate_selected(positions, Some(movable))
    }
    fn evaluate_selected(
        &self,
        positions: &[V],
        movable: Option<&[bool]>,
    ) -> Result<BondedEvaluation, EvaluationError> {
        if positions.len() != self.atom_indices.len() {
            return Err(EvaluationError::Size {
                expected: self.atom_indices.len(),
                actual: positions.len(),
            });
        }
        for (i, p) in positions.iter().enumerate() {
            if !p.iter().all(|v| v.is_finite()) {
                return Err(EvaluationError::NonFinite(i));
            }
        }
        let mut result = BondedEvaluation {
            energy: BondedEnergy::default(),
            forces: vec![[0.0; 3]; positions.len()],
        };
        let active = |ids: &[usize]| movable.is_none_or(|m| ids.iter().any(|&i| m[i]));
        result.energy = evaluate_terms::<true>(
            positions,
            &mut result.forces,
            self.bonds.iter().filter(|t| active(&t.atoms)),
            self.angles.iter().filter(|t| active(&t.atoms)),
            self.propers
                .iter()
                .map(|t| (TermKind::Proper, t))
                .chain(self.impropers.iter().map(|t| (TermKind::Improper, t)))
                .filter(|(_, t)| active(&t.atoms)),
        )?;
        if !result.energy.total().is_finite()
            || result.forces.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(EvaluationError::Overflow);
        }
        Ok(result)
    }
}
use super::{AngleParameter, BondParameter, BondedTerm, FourierTerm};
type BondTerm = BondedTerm<2, BondParameter>;
type AngleTerm = BondedTerm<3, AngleParameter>;
type TorsionTerm = BondedTerm<4, Vec<FourierTerm>>;
#[derive(Default)]
pub(super) struct ActiveBonded<'a> {
    bonds: Vec<&'a BondTerm>,
    angles: Vec<&'a AngleTerm>,
    torsions: Vec<(TermKind, &'a TorsionTerm)>,
}
impl<'a> ActiveBonded<'a> {
    pub(super) fn new(model: &'a BondedModel, movable: &[bool]) -> Self {
        let mut terms = Self::default();
        terms.reset(model, movable);
        terms
    }
    pub(super) fn reset(&mut self, model: &'a BondedModel, movable: &[bool]) {
        let active = |ids: &[usize]| ids.iter().any(|&i| movable[i]);
        self.bonds.clear();
        self.bonds
            .extend(model.bonds.iter().filter(|t| active(&t.atoms)));
        self.angles.clear();
        self.angles
            .extend(model.angles.iter().filter(|t| active(&t.atoms)));
        self.torsions.clear();
        self.torsions.extend(
            model
                .propers
                .iter()
                .map(|t| (TermKind::Proper, t))
                .chain(model.impropers.iter().map(|t| (TermKind::Improper, t)))
                .filter(|(_, t)| active(&t.atoms)),
        );
    }
    pub(super) fn evaluate(
        &self,
        positions: &[V],
        force: &mut [V],
    ) -> Result<BondedEnergy, EvaluationError> {
        self.evaluate_mode::<true>(positions, force)
    }
    pub(super) fn energy(&self, positions: &[V]) -> Result<BondedEnergy, EvaluationError> {
        self.evaluate_mode::<false>(positions, &mut [])
    }
    fn evaluate_mode<const FORCES: bool>(
        &self,
        positions: &[V],
        force: &mut [V],
    ) -> Result<BondedEnergy, EvaluationError> {
        evaluate_terms::<FORCES>(
            positions,
            force,
            self.bonds.iter().copied(),
            self.angles.iter().copied(),
            self.torsions.iter().copied(),
        )
    }
}
fn evaluate_terms<'a, const FORCES: bool>(
    positions: &[V],
    force: &mut [V],
    bonds: impl Iterator<Item = &'a BondTerm>,
    angles: impl Iterator<Item = &'a AngleTerm>,
    torsions: impl Iterator<Item = (TermKind, &'a TorsionTerm)>,
) -> Result<BondedEnergy, EvaluationError> {
    let mut energy = BondedEnergy::default();
    for term in bonds {
        let [a, b] = term.atoms;
        let d = sub(positions[a], positions[b]);
        let r = dot(d, d).sqrt();
        if r < 1e-12 {
            return Err(degenerate(TermKind::Bond, &term.atoms));
        }
        let dr = r - term.parameter.length;
        energy.bonds += term.parameter.k * dr * dr;
        if !FORCES {
            continue;
        }
        let f = mul(d, -2.0 * term.parameter.k * dr / r);
        force[a] = add(force[a], f);
        force[b] = sub(force[b], f);
    }
    for term in angles {
        let [a, b, c] = term.atoms;
        let u = sub(positions[a], positions[b]);
        let v = sub(positions[c], positions[b]);
        let (uu, vv) = (dot(u, u), dot(v, v));
        if uu < 1e-24 || vv < 1e-24 {
            return Err(degenerate(TermKind::Angle, &term.atoms));
        }
        let uv = (uu * vv).sqrt();
        let cosine = (dot(u, v) / uv).clamp(-1.0, 1.0);
        let sine = dot(cross(u, v), cross(u, v)).sqrt() / uv;
        let theta = sine.atan2(cosine);
        let d = theta - term.parameter.angle;
        energy.angles += term.parameter.k * d * d;
        if sine < 1e-12 {
            if d.abs() < 1e-10 {
                continue;
            }
            return Err(degenerate(TermKind::Angle, &term.atoms));
        }
        if !FORCES {
            continue;
        }
        let factor = -2.0 * term.parameter.k * d / sine;
        let fa = mul(sub(mul(u, cosine / uu), mul(v, 1.0 / uv)), factor);
        let fc = mul(sub(mul(v, cosine / vv), mul(u, 1.0 / uv)), factor);
        force[a] = add(force[a], fa);
        force[c] = add(force[c], fc);
        force[b] = sub(force[b], add(fa, fc));
    }
    for (kind, term) in torsions {
        if term.parameter.iter().all(|p| p.k == 0.0) {
            continue;
        }
        let [a, b, c, d] = term.atoms;
        let b1 = sub(positions[b], positions[a]);
        let b2 = sub(positions[c], positions[b]);
        let b3 = sub(positions[d], positions[c]);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let (n1s, n2s, b2s) = (dot(n1, n1), dot(n2, n2), dot(b2, b2));
        if n1s < 1e-28 || n2s < 1e-28 || b2s < 1e-24 {
            return Err(degenerate(kind, &term.atoms));
        }
        let length = b2s.sqrt();
        let phi = (length * dot(b1, n2)).atan2(dot(n1, n2));
        let mut e = 0.0;
        let mut derivative = 0.0;
        for p in &term.parameter {
            let arg = p.periodicity as f64 * phi - p.phase;
            e += p.k * (1.0 + arg.cos());
            if FORCES {
                derivative -= p.k * p.periodicity as f64 * arg.sin();
            }
        }
        if kind == TermKind::Proper {
            energy.propers += e;
        } else {
            energy.impropers += e;
        }
        if !FORCES {
            continue;
        }
        let g0 = mul(n1, -length / n1s);
        let g3 = mul(n2, length / n2s);
        let t1 = dot(b1, b2) / b2s;
        let t3 = dot(b3, b2) / b2s;
        let g1 = add(mul(g0, -1.0 - t1), mul(g3, t3));
        let g2 = add(mul(g0, t1), mul(g3, -1.0 - t3));
        for (i, g) in [(a, g0), (b, g1), (c, g2), (d, g3)] {
            force[i] = sub(force[i], mul(g, derivative));
        }
    }
    Ok(energy)
}
fn degenerate(kind: TermKind, atoms: &[usize]) -> EvaluationError {
    EvaluationError::Degenerate {
        kind,
        atoms: atoms.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amber::{BondedTerm, FourierTerm, Resolution, ResolutionMethod};
    #[test]
    fn signed_phase_matches_openmm_and_invalid_geometry_is_explicit() {
        let mut model = BondedModel {
            atom_indices: vec![0, 1, 2, 3],
            atom_types: vec!["c3".into(); 4],
            bonds: vec![],
            angles: vec![],
            impropers: vec![],
            propers: vec![BondedTerm {
                atoms: [0, 1, 2, 3],
                parameter: vec![FourierTerm {
                    k: 2.5,
                    periodicity: 1,
                    phase: 0.37,
                }],
                resolution: Resolution {
                    method: ResolutionMethod::Exact,
                    matched_types: vec![],
                    penalty: 0.0,
                    origins: vec![],
                },
            }],
        };
        let p = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        ];
        let r = model.evaluate(&p).unwrap();
        assert!(
            (r.energy.propers - 2.5 * (1.0 + (-std::f64::consts::FRAC_PI_2 - 0.37).cos())).abs()
                < 1e-12
        );
        assert!((r.forces[0][2] + 2.5 * 0.37_f64.cos()).abs() < 1e-12);
        assert!(matches!(
            model.evaluate(&p[..3]),
            Err(EvaluationError::Size { .. })
        ));
        assert!(matches!(
            model.evaluate(&[[0.0; 3]; 4]),
            Err(EvaluationError::Degenerate { .. })
        ));
        let mut bad = p;
        bad[1][0] = f64::NAN;
        assert_eq!(
            model.evaluate(&bad).unwrap_err(),
            EvaluationError::NonFinite(1)
        );
        model.propers[0].parameter[0].k = 0.0;
        assert_eq!(model.evaluate(&[[0.0; 3]; 4]).unwrap().energy.total(), 0.0);
    }
}
