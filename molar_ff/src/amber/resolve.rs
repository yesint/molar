//! Data-driven parameter completion using Amber's PARMCHK correspondence scores.
//! This is an independent implementation; no Amber executable is used at runtime.
use super::parameters::*;
use crate::FFType;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default)]
pub(crate) struct AtomRule {
    pub improper: bool,
    group: i32,
    equivalent_group: i32,
    pub element: u8,
    // bl, blf, central angle, central angle k, terminal angle, terminal angle k,
    // central torsion, terminal torsion, general/improper.
    corr: BTreeMap<String, ([f64; 9], bool)>,
}
#[derive(Clone, Debug, Default)]
pub(crate) struct Rules {
    pub atoms: BTreeMap<String, AtomRule>,
    weights: BTreeMap<String, f64>,
    angle_elements: BTreeMap<u8, (f64, f64)>,
    bond_elements: BTreeMap<(u8, u8), (f64, f64)>,
    bond_exponent: f64,
}
impl Rules {
    pub fn builtin(ff: FFType) -> Self {
        let mut result = Self::default();
        let mut current = String::new();
        let text = include_str!("../../assets/amber/PARMCHK.DAT");
        for line in text.lines() {
            let f: Vec<_> = line.split_whitespace().collect();
            match f.first().copied() {
                Some("PARM") if f.len() >= 7 => {
                    current = f[1].into();
                    let mut rule = AtomRule {
                        improper: f[2] == "1",
                        group: f[3].parse().unwrap(),
                        equivalent_group: f[5].parse().unwrap(),
                        element: f[6].parse().unwrap(),
                        ..Default::default()
                    };
                    rule.corr.insert(current.clone(), ([0.0; 9], true));
                    result.atoms.insert(current.clone(), rule);
                }
                Some("EQUA") => {
                    result
                        .atoms
                        .get_mut(&current)
                        .unwrap()
                        .corr
                        .insert(f[1].into(), ([0.0; 9], true));
                }
                Some("CORR") if f.len() >= 11 => {
                    let v = std::array::from_fn(|i| f[i + 2].parse().unwrap());
                    result
                        .atoms
                        .get_mut(&current)
                        .unwrap()
                        .corr
                        .insert(f[1].into(), (v, false));
                }
                Some(key)
                    if key.starts_with("WEIGHT_")
                        || key.starts_with("DEFAULT_")
                        || key == "THRESHOLD_BA" =>
                {
                    result.weights.insert(key.into(), f[1].parse().unwrap());
                }
                _ => (),
            }
        }
        let defaults = [
            "BL", "BLF", "BA_CTR", "BAF_CTR", "BA", "BAF", "TOR_CTR", "TOR",
        ];
        for rule in result.atoms.values_mut() {
            for (v, _) in rule.corr.values_mut() {
                for (i, name) in defaults.iter().enumerate() {
                    if v[i] < 0.0 {
                        v[i] = result.weights[&format!("DEFAULT_{name}")];
                    }
                }
                v[6] = v[6] * result.weights["DEFAULT_FRACT1"]
                    + v[8] * result.weights["DEFAULT_FRACT2"];
            }
        }
        let empirical = match ff {
            FFType::Gaff => include_str!("../../assets/amber/PARM_BLBA_GAFF.DAT"),
            FFType::Gaff2 => include_str!("../../assets/amber/PARM_BLBA_GAFF2.DAT"),
        };
        for line in empirical.lines() {
            let f: Vec<_> = line.split_whitespace().collect();
            if f.get(1) == Some(&"PC") {
                result.bond_exponent = f[2].parse().unwrap();
            }
            if f.get(1) == Some(&"BL") && f.len() >= 8 {
                let a: u8 = f[3].parse().unwrap();
                let b: u8 = f[5].parse().unwrap();
                result.bond_elements.insert(
                    (a.min(b), a.max(b)),
                    (f[6].parse().unwrap(), f[7].parse().unwrap()),
                );
            }
            if f.get(1) == Some(&"BA") && f.len() >= 6 {
                result.angle_elements.insert(
                    f[3].parse().unwrap(),
                    (f[4].parse().unwrap(), f[5].parse().unwrap()),
                );
            }
        }
        result
    }
    fn weight(&self, name: &str) -> f64 {
        self.weights[name]
    }
    fn replacement(
        &self,
        from: &str,
        to: &str,
        kind: TermKind,
        position: usize,
    ) -> Option<(f64, bool)> {
        if from == to {
            return Some((0.0, true));
        }
        let (v, equivalent) = self.atoms.get(from)?.corr.get(to)?;
        let score = match kind {
            TermKind::Bond => v[0] * self.weight("WEIGHT_BL") + v[1] * self.weight("WEIGHT_BLF"),
            TermKind::Angle if position == 1 => {
                (v[2] * self.weight("WEIGHT_BA") + v[3] * self.weight("WEIGHT_BAF"))
                    * self.weight("WEIGHT_BA_CTR")
            }
            TermKind::Angle => v[4] * self.weight("WEIGHT_BA") + v[5] * self.weight("WEIGHT_BAF"),
            TermKind::Proper if position == 1 || position == 2 => {
                v[6] * self.weight("WEIGHT_TOR_CTR")
            }
            TermKind::Proper => v[7],
            TermKind::Improper => {
                v[8] + if position == 2 {
                    self.weight("WEIGHT_IMPROPER")
                } else {
                    0.0
                }
            }
        };
        Some((score, *equivalent))
    }
    fn group_penalty(&self, types: &[String]) -> f64 {
        let groups: Vec<_> = types
            .iter()
            .filter_map(|t| self.atoms.get(t).map(|r| r.group))
            .collect();
        if groups.windows(2).any(|w| w[0] != w[1]) {
            self.weight("WEIGHT_GROUP")
        } else {
            0.0
        }
    }
    fn conjugation_penalty(&self, from: &[String], to: &[String], a: usize, b: usize) -> f64 {
        let group = |t: &str| self.atoms.get(t).map_or(0, |r| r.equivalent_group);
        let (x, y, u, v) = (
            group(&from[a]),
            group(&from[b]),
            group(&to[a]),
            group(&to[b]),
        );
        if x == 0 && y == 0 {
            return 0.0;
        }
        if (x.abs() + y.abs() == 3) != (u.abs() + v.abs() == 3) {
            self.weight("WEIGHT_EQUTYPE")
        } else if x + y == 0 && u < 0 && v < 0 {
            0.5 * self.weight("WEIGHT_EQUTYPE")
        } else {
            0.0
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Resolved {
    pub value: Value,
    pub resolution: Resolution,
    pub order: Vec<usize>,
}

impl AmberParameters {
    pub(crate) fn direct(&self, kind: TermKind, query: &[String]) -> Option<&Entry<Value>> {
        let index = self.index.get(&(kind, canonical(kind, query)))?;
        self.entries.get(&kind)?.get(*index)
    }
    pub(crate) fn resolve(
        &self,
        kind: TermKind,
        query: &[String],
        automatic: bool,
    ) -> Option<Resolved> {
        let permutations = permutations(kind, query.len());
        if let Some(entry) = self.direct(kind, query) {
            let usable = !matches!(&entry.value, Value::Bond(p) if p.k == 0.0)
                && !matches!(&entry.value, Value::Angle(p) if p.k == 0.0);
            if usable {
                if let Some(order) = permutations.iter().find(|order| {
                    order
                        .iter()
                        .enumerate()
                        .all(|(i, &j)| entry.types[i] == query[j])
                }) {
                    return Some(Resolved {
                        value: entry.value.clone(),
                        order: order.clone(),
                        resolution: Resolution {
                            method: ResolutionMethod::Exact,
                            matched_types: entry.types.clone(),
                            penalty: 0.0,
                            origins: entry.origins.clone(),
                        },
                    });
                }
            }
        }
        // Stage order matches parmchk2: exact, equivalent-specific, generic,
        // equivalent-generic, corresponding-specific, corresponding-generic.
        let mut best: Option<((u8, f64, usize, usize), Resolved)> = None;
        for (entry_index, entry) in self.entries.get(&kind).into_iter().flatten().enumerate() {
            if matches!(&entry.value, Value::Bond(p) if p.k == 0.0)
                || matches!(&entry.value, Value::Angle(p) if p.k == 0.0)
            {
                continue;
            }
            let wildcards = entry.types.iter().filter(|s| s.as_str() == "X").count();
            for order in &permutations {
                let actual: Vec<_> = order.iter().map(|&i| query[i].clone()).collect();
                let mut score = 0.0;
                let mut equivalents = true;
                let mut changed = false;
                let mut changes = 0;
                let mut valid = true;
                for (pos, (a, b)) in actual.iter().zip(&entry.types).enumerate() {
                    if b == "X" {
                        continue;
                    }
                    if a != b {
                        changed = true;
                        changes += 1;
                    }
                    match self.rules.replacement(a, b, kind, pos) {
                        Some((s, e)) if automatic || a == b => {
                            score += s;
                            equivalents &= e;
                        }
                        _ => {
                            valid = false;
                            break;
                        }
                    }
                }
                if !valid {
                    continue;
                }
                let stage = if !changed {
                    if wildcards == 0 { 0 } else { 2 }
                } else if equivalents {
                    if wildcards == 0 { 1 } else { 3 }
                } else if wildcards == 0 {
                    4
                } else {
                    5
                };
                if changed && !equivalents {
                    score += self.rules.group_penalty(&entry.types);
                    match kind {
                        TermKind::Bond => {
                            score += self.rules.conjugation_penalty(&actual, &entry.types, 0, 1)
                        }
                        TermKind::Proper => {
                            score += self.rules.conjugation_penalty(&actual, &entry.types, 1, 2)
                        }
                        TermKind::Improper => {
                            score += wildcards as f64 * self.rules.weight("WEIGHT_X")
                        }
                        _ => (),
                    }
                }
                if kind == TermKind::Angle
                    && stage >= 4
                    && score > self.rules.weight("THRESHOLD_BA")
                {
                    continue;
                }
                let rank = (stage, score + wildcards as f64 * 1e-6, changes, entry_index);
                if best.as_ref().is_none_or(|(b, _)| rank < *b) {
                    let method = if !changed {
                        if wildcards == 0 {
                            ResolutionMethod::Exact
                        } else {
                            ResolutionMethod::Wildcard
                        }
                    } else if equivalents {
                        ResolutionMethod::Equivalent
                    } else {
                        ResolutionMethod::Analogy
                    };
                    best = Some((
                        rank,
                        Resolved {
                            value: entry.value.clone(),
                            order: order.clone(),
                            resolution: Resolution {
                                method,
                                penalty: score,
                                matched_types: entry.types.clone(),
                                origins: entry.origins.clone(),
                            },
                        },
                    ));
                }
            }
        }
        if let Some((_, result)) = best {
            return Some(result);
        }
        if !automatic {
            return None;
        }
        if kind == TermKind::Bond {
            return self.empirical_bond(query);
        }
        if kind == TermKind::Angle {
            return self.empirical_angle(query);
        }
        if kind == TermKind::Improper
            && self.rules.atoms.get(&query[2]).is_some_and(|r| r.improper)
            && query.iter().all(|t| self.rules.atoms.contains_key(t))
        {
            return Some(Resolved {
                value: Value::Torsion(vec![FourierTerm {
                    k: 1.1 * 4.184,
                    phase: std::f64::consts::PI,
                    periodicity: 2,
                }]),
                order: (0..4).collect(),
                resolution: Resolution {
                    method: ResolutionMethod::DefaultImproper,
                    matched_types: query.to_vec(),
                    penalty: 0.0,
                    origins: vec![ParameterOrigin {
                        file: "parmchk2 default improper".into(),
                        line: 0,
                    }],
                },
            });
        }
        None
    }

    fn empirical_bond(&self, query: &[String]) -> Option<Resolved> {
        let entry = self.direct(TermKind::Bond, query)?;
        let Value::Bond(p) = entry.value else {
            return None;
        };
        let a = self.rules.atoms.get(&query[0])?.element;
        let b = self.rules.atoms.get(&query[1])?.element;
        let (_, constant) = self.rules.bond_elements.get(&(a.min(b), a.max(b)))?;
        let k = (constant - self.rules.bond_exponent * (p.length * 10.0).ln()).exp() * 418.4;
        if !k.is_finite() || k <= 0.0 {
            return None;
        }
        let mut origins = entry.origins.clone();
        origins.push(ParameterOrigin {
            file: match self.ff {
                FFType::Gaff => "PARM_BLBA_GAFF.DAT",
                FFType::Gaff2 => "PARM_BLBA_GAFF2.DAT",
            }
            .into(),
            line: 0,
        });
        Some(Resolved {
            value: Value::Bond(BondParameter {
                k,
                length: p.length,
            }),
            order: vec![0, 1],
            resolution: Resolution {
                method: ResolutionMethod::EmpiricalBond,
                matched_types: query.to_vec(),
                penalty: 0.0,
                origins,
            },
        })
    }

    fn empirical_angle(&self, query: &[String]) -> Option<Resolved> {
        let options: Vec<Vec<String>> = query
            .iter()
            .map(|q| {
                self.rules
                    .atoms
                    .get(q)
                    .map(|r| r.corr.keys().cloned().collect())
                    .unwrap_or_else(|| vec![q.clone()])
            })
            .collect();
        let mut best: Option<(f64, Resolved)> = None;
        for a in &options[0] {
            for b in &options[1] {
                for c in &options[2] {
                    let types = vec![a.clone(), b.clone(), c.clone()];
                    let mut score = self.rules.group_penalty(&types);
                    for i in 0..3 {
                        score += self
                            .rules
                            .replacement(&query[i], &types[i], TermKind::Angle, i)?
                            .0;
                    }
                    if best.as_ref().is_some_and(|(s, _)| score >= *s) {
                        continue;
                    }
                    let Some(left) =
                        self.direct(TermKind::Angle, &[a.clone(), b.clone(), a.clone()])
                    else {
                        continue;
                    };
                    let Some(right) =
                        self.direct(TermKind::Angle, &[c.clone(), b.clone(), c.clone()])
                    else {
                        continue;
                    };
                    let (Value::Angle(l), Value::Angle(r)) = (&left.value, &right.value) else {
                        unreachable!()
                    };
                    let Some(bl) = self.resolve(TermKind::Bond, &[a.clone(), b.clone()], true)
                    else {
                        continue;
                    };
                    let Some(br) = self.resolve(TermKind::Bond, &[b.clone(), c.clone()], true)
                    else {
                        continue;
                    };
                    let (Value::Bond(blp), Value::Bond(brp)) = (&bl.value, &br.value) else {
                        unreachable!()
                    };
                    let elements: Option<Vec<_>> = types
                        .iter()
                        .map(|t| {
                            self.rules
                                .atoms
                                .get(t)
                                .and_then(|r| self.rules.angle_elements.get(&r.element))
                        })
                        .collect();
                    let Some(e) = elements else { continue };
                    let theta = 0.5 * (l.angle + r.angle);
                    let (r1, r2) = (blp.length * 10.0, brp.length * 10.0);
                    let d = ((r1 - r2) / (r1 + r2)).powi(2);
                    let k = 143.9 * e[0].1 * e[1].0 * e[2].1 * (-2.0 * d).exp()
                        / ((r1 + r2) * theta.sqrt())
                        * 4.184;
                    if !k.is_finite() || k <= 0.0 {
                        continue;
                    }
                    let mut origins = left.origins.clone();
                    origins.extend(right.origins.clone());
                    origins.extend(bl.resolution.origins);
                    origins.extend(br.resolution.origins);
                    origins.push(ParameterOrigin {
                        file: match self.ff {
                            FFType::Gaff => "PARM_BLBA_GAFF.DAT",
                            FFType::Gaff2 => "PARM_BLBA_GAFF2.DAT",
                        }
                        .into(),
                        line: 0,
                    });
                    best = Some((
                        score,
                        Resolved {
                            value: Value::Angle(AngleParameter { k, angle: theta }),
                            order: (0..3).collect(),
                            resolution: Resolution {
                                method: ResolutionMethod::EmpiricalAngle,
                                matched_types: types,
                                penalty: score,
                                origins,
                            },
                        },
                    ));
                }
            }
        }
        best.map(|(_, r)| r)
    }
}

fn permutations(kind: TermKind, n: usize) -> Vec<Vec<usize>> {
    if kind == TermKind::Improper {
        vec![
            vec![0, 1, 2, 3],
            vec![0, 3, 2, 1],
            vec![1, 0, 2, 3],
            vec![1, 3, 2, 0],
            vec![3, 0, 2, 1],
            vec![3, 1, 2, 0],
        ]
    } else {
        vec![(0..n).collect(), (0..n).rev().collect()]
    }
}
