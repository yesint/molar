use super::resolve::Rules;
use crate::FFType;
use std::{collections::BTreeMap, sync::OnceLock};

/// Functional form of a bonded interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum TermKind {
    Bond,
    Angle,
    Proper,
    Improper,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondParameter {
    /// kJ/mol/nm²; E = k (r - length)².
    pub k: f64,
    /// Equilibrium length in nm.
    pub length: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngleParameter {
    /// kJ/mol/radian²; E = k (theta - angle)².
    pub k: f64,
    /// Radians.
    pub angle: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FourierTerm {
    /// kJ/mol, including Amber's IDIVF division; E = k [1 + cos(n phi - phase)].
    pub k: f64,
    pub periodicity: u32,
    /// Radians.
    pub phase: f64,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParameterOrigin {
    pub file: String,
    pub line: usize,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ResolutionMethod {
    Exact,
    Wildcard,
    Equivalent,
    Analogy,
    EmpiricalBond,
    EmpiricalAngle,
    DefaultImproper,
}
/// Audit trail for every interaction, including the actual lookup types and score.
#[derive(Clone, Debug, PartialEq)]
pub struct Resolution {
    pub method: ResolutionMethod,
    pub matched_types: Vec<String>,
    pub penalty: f64,
    pub origins: Vec<ParameterOrigin>,
}
impl Resolution {
    pub fn is_estimated(&self) -> bool {
        !matches!(
            self.method,
            ResolutionMethod::Exact | ResolutionMethod::Wildcard
        )
    }
}
#[derive(Clone, Debug)]
pub(crate) struct Entry<T> {
    pub types: Vec<String>,
    pub value: T,
    pub origins: Vec<ParameterOrigin>,
}
#[derive(Clone, Debug)]
pub(crate) enum Value {
    Bond(BondParameter),
    Angle(AngleParameter),
    Torsion(Vec<FourierTerm>),
}

/// Amber atom Lennard-Jones data. Pair Rmin is the sum of radii and pair epsilon
/// is the geometric mean. Zero epsilon/radius is valid (e.g. GAFF hydroxyl H).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LennardJonesParameter {
    /// Rmin / 2 in nm, not sigma.
    pub radius: f64,
    /// Well depth in kJ/mol.
    pub epsilon: f64,
}

/// A read-only parameter database. Built-ins are parsed once and shared between models.
/// Files and correspondence rules are pinned together in `assets/amber/source.json`.
#[derive(Clone, Debug)]
pub struct AmberParameters {
    pub(crate) nonbonded: BTreeMap<String, (LennardJonesParameter, ParameterOrigin)>,
    pub(crate) entries: BTreeMap<TermKind, Vec<Entry<Value>>>,
    pub(crate) rules: Rules,
    pub(crate) index: BTreeMap<(TermKind, Vec<String>), usize>,
    pub(crate) ff: FFType,
}
#[derive(Clone, Debug, thiserror::Error, PartialEq)]
#[error("{file}:{line}: {message}")]
pub struct ParameterError {
    pub file: String,
    pub line: usize,
    pub message: String,
}

impl AmberParameters {
    pub fn builtin(ff: FFType) -> &'static Self {
        static GAFF: OnceLock<AmberParameters> = OnceLock::new();
        static GAFF2: OnceLock<AmberParameters> = OnceLock::new();
        let (cell, text, name) = match ff {
            FFType::Gaff => (
                &GAFF,
                include_str!("../../assets/amber/gaff.dat"),
                "gaff.dat",
            ),
            FFType::Gaff2 => (
                &GAFF2,
                include_str!("../../assets/amber/gaff2.dat"),
                "gaff2.dat",
            ),
        };
        cell.get_or_init(|| Self::from_dat(text, name, ff).expect("validated bundled Amber data"))
    }

    /// Read Amber bonded and Rmin/epsilon nonbonded sections. Mass and legacy
    /// 10-12 hydrogen-bond sections are not evaluated.
    pub fn from_dat(text: &str, name: &str, ff: FFType) -> Result<Self, ParameterError> {
        let mut result = Self {
            entries: BTreeMap::new(),
            nonbonded: BTreeMap::new(),
            index: BTreeMap::new(),
            rules: Rules::builtin(ff),
            ff,
        };
        let mut blocks: Vec<Vec<(usize, &str)>> = vec![vec![]];
        for (i, line) in text.lines().enumerate() {
            if line.trim().is_empty() {
                if !blocks.last().unwrap().is_empty() {
                    blocks.push(vec![]);
                }
            } else {
                blocks.last_mut().unwrap().push((i + 1, line));
            }
        }
        if blocks.len() < 5 {
            return Err(error(name, 1, "missing Amber bonded sections"));
        }
        for (block, kind) in [
            (1, TermKind::Bond),
            (2, TermKind::Angle),
            (3, TermKind::Proper),
            (4, TermKind::Improper),
        ] {
            // The hydrophilic atom-type line directly precedes BOND in a parm.dat.
            let lines = if block == 1 && split_types(blocks[block][0].1, 2).is_none() {
                &blocks[block][1..]
            } else {
                &blocks[block][..]
            };
            result.read_section(lines, name, kind)?;
        }
        if let Some(start) = text
            .lines()
            .position(|l| l.trim_start().starts_with("MOD4"))
        {
            let header = text.lines().nth(start).unwrap();
            if header.split_whitespace().nth(1) != Some("RE") {
                return Err(error(
                    name,
                    start + 1,
                    "only Amber RE nonbonded parameters are supported",
                ));
            }
            for (i, line) in text.lines().enumerate().skip(start + 1) {
                if line.trim().is_empty() || line.trim() == "END" {
                    break;
                }
                result.read_nonbonded(line, name, i + 1)?;
            }
        }
        Ok(result)
    }

    /// Direct Lennard-Jones lookup; GAFF tables explicitly list all supported types.
    /// Unknown types are not assigned guessed van der Waals parameters.
    pub fn lennard_jones(
        &self,
        atom_type: &str,
    ) -> Option<(LennardJonesParameter, &ParameterOrigin)> {
        self.nonbonded
            .get(atom_type)
            .map(|(p, origin)| (*p, origin))
    }

    fn read_nonbonded(
        &mut self,
        line: &str,
        name: &str,
        line_no: usize,
    ) -> Result<(), ParameterError> {
        let line = line.trim();
        if line.is_empty() || line.starts_with(['#', '!']) {
            return Ok(());
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        let fail = || {
            error(
                name,
                line_no,
                "invalid or unresolved Lennard-Jones parameter",
            )
        };
        if fields.len() < 3 || line.to_ascii_lowercase().contains("attn") {
            return Err(fail());
        }
        let value = |i: usize| {
            fields[i]
                .replace(['D', 'd'], "E")
                .parse::<f64>()
                .ok()
                .filter(|v| v.is_finite() && *v >= 0.0)
                .ok_or_else(fail)
        };
        let parameter = LennardJonesParameter {
            radius: value(1)? * 0.1,
            epsilon: value(2)? * 4.184,
        };
        if !parameter.radius.is_finite() || !parameter.epsilon.is_finite() {
            return Err(fail());
        }
        self.nonbonded.insert(
            fields[0].into(),
            (
                parameter,
                ParameterOrigin {
                    file: name.into(),
                    line: line_no,
                },
            ),
        );
        Ok(())
    }

    /// Optional import for external Amber data. Not needed for automatic completion.
    /// Parsing is transactional; marked unresolved placeholders are rejected.
    pub fn with_frcmod(mut self, text: &str, name: &str) -> Result<Self, ParameterError> {
        let mut sections: BTreeMap<TermKind, Vec<(usize, &str)>> = BTreeMap::new();
        let mut current = None;
        let mut nonbonded = false;
        for (i, line) in text.lines().enumerate().skip(1) {
            let heading = line.trim();
            if matches!(
                heading,
                "BOND"
                    | "ANGLE"
                    | "DIHE"
                    | "DIHEDRAL"
                    | "IMPROPER"
                    | "IMPR"
                    | "MASS"
                    | "HBOND"
                    | "LJEDIT"
            ) {
                nonbonded = false;
            }
            if heading == "LJEDIT" {
                return Err(error(
                    name,
                    i + 1,
                    "pair-specific LJEDIT parameters are not supported",
                ));
            }
            match heading {
                "BOND" => current = Some(TermKind::Bond),
                "ANGLE" => current = Some(TermKind::Angle),
                "DIHE" | "DIHEDRAL" => current = Some(TermKind::Proper),
                "IMPROPER" | "IMPR" => current = Some(TermKind::Improper),
                "NONBON" | "NONB" => {
                    current = None;
                    nonbonded = true;
                }
                "MASS" | "HBOND" => current = None,
                "" => (),
                _ if !heading.starts_with(['#', '!']) => {
                    if nonbonded {
                        self.read_nonbonded(line, name, i + 1)?;
                    }
                    if let Some(kind) = current {
                        sections.entry(kind).or_default().push((i + 1, line));
                    }
                }
                _ => (),
            }
        }
        for (kind, lines) in sections {
            self.read_section(&lines, name, kind)?;
        }
        Ok(self)
    }

    fn read_section(
        &mut self,
        lines: &[(usize, &str)],
        name: &str,
        kind: TermKind,
    ) -> Result<(), ParameterError> {
        let mut pending: Option<(Vec<String>, Vec<FourierTerm>, Vec<ParameterOrigin>)> = None;
        for &(line_no, line) in lines {
            if line.trim_start().starts_with(['#', '!']) {
                continue;
            }
            let fail = |message| error(name, line_no, message);
            if line.to_ascii_lowercase().contains("attn") {
                return Err(fail("unresolved Amber placeholder"));
            }
            let n = match kind {
                TermKind::Bond => 2,
                TermKind::Angle => 3,
                _ => 4,
            };
            let (types, tail) =
                split_types(line, n).ok_or_else(|| fail("invalid atom-type tuple"))?;
            let count = match kind {
                TermKind::Bond | TermKind::Angle => 2,
                TermKind::Proper => 4,
                TermKind::Improper => 3,
            };
            let fields: Vec<_> = tail.split_whitespace().take(count).collect();
            if fields.len() != count {
                return Err(fail("missing numeric parameter"));
            }
            let values: Vec<f64> = fields
                .iter()
                .map(|f| {
                    f.replace(['D', 'd'], "E")
                        .parse::<f64>()
                        .ok()
                        .filter(|v| v.is_finite())
                        .ok_or_else(|| fail("non-finite or invalid number"))
                })
                .collect::<Result<_, _>>()?;
            let origin = ParameterOrigin {
                file: name.into(),
                line: line_no,
            };
            let value = match kind {
                TermKind::Bond => {
                    if values[0] < 0.0 || values[1] <= 0.0 {
                        return Err(fail("bond requires nonnegative k and positive length"));
                    }
                    Value::Bond(BondParameter {
                        k: values[0] * 418.4,
                        length: values[1] * 0.1,
                    })
                }
                TermKind::Angle => {
                    if values[0] < 0.0 || !(0.0..=180.0).contains(&values[1]) || values[1] == 0.0 {
                        return Err(fail("invalid angle constant or equilibrium angle"));
                    }
                    Value::Angle(AngleParameter {
                        k: values[0] * 4.184,
                        angle: values[1].to_radians(),
                    })
                }
                _ => {
                    let (div, k, phase, pn) = if kind == TermKind::Proper {
                        (values[0], values[1], values[2], values[3])
                    } else {
                        (1.0, values[0], values[1], values[2])
                    };
                    if div <= 0.0
                        || k < 0.0
                        || pn == 0.0
                        || pn.abs() > 100.0
                        || pn.abs().fract() != 0.0
                    {
                        return Err(fail("invalid Fourier divisor, amplitude or periodicity"));
                    }
                    let term = FourierTerm {
                        k: k * 4.184 / div,
                        periodicity: pn.abs() as u32,
                        phase: phase.to_radians(),
                    };
                    let (key, mut terms, mut origins) =
                        pending.take().unwrap_or((types.clone(), vec![], vec![]));
                    if key != types {
                        return Err(fail("torsion continuation changes atom types"));
                    }
                    terms.push(term);
                    origins.push(origin);
                    if pn < 0.0 {
                        pending = Some((key, terms, origins));
                        continue;
                    }
                    self.insert(
                        kind,
                        Entry {
                            types,
                            value: Value::Torsion(terms),
                            origins,
                        },
                    );
                    continue;
                }
            };
            self.insert(
                kind,
                Entry {
                    types,
                    value,
                    origins: vec![origin],
                },
            );
        }
        if pending.is_some() {
            return Err(error(
                name,
                lines.last().map_or(1, |v| v.0),
                "unterminated torsion continuation",
            ));
        }
        Ok(())
    }

    fn insert(&mut self, kind: TermKind, entry: Entry<Value>) {
        let entries = self.entries.entry(kind).or_default();
        // An override replaces a complete parameter/series with the same tuple.
        let key = canonical(kind, &entry.types);
        if let Some(&i) = self.index.get(&(kind, key.clone())) {
            entries[i] = entry;
        } else {
            self.index.insert((kind, key), entries.len());
            entries.push(entry);
        }
    }
}

pub(crate) fn canonical(kind: TermKind, types: &[String]) -> Vec<String> {
    let mut key = types.to_vec();
    if kind == TermKind::Improper {
        let mut outside = [key[0].clone(), key[1].clone(), key[3].clone()];
        outside.sort();
        key[0] = outside[0].clone();
        key[1] = outside[1].clone();
        key[3] = outside[2].clone();
    } else {
        let rev: Vec<_> = key.iter().rev().cloned().collect();
        if rev < key {
            key = rev;
        }
    }
    key
}
fn error(file: &str, line: usize, message: &str) -> ParameterError {
    ParameterError {
        file: file.into(),
        line,
        message: message.into(),
    }
}
fn split_types(line: &str, n: usize) -> Option<(Vec<String>, &str)> {
    let mut rest = line.trim_start();
    let mut types = vec![];
    for i in 0..n {
        if i + 1 < n {
            let (a, tail) = rest.split_once('-')?;
            let a = a.trim();
            if a.is_empty() || a.contains(char::is_whitespace) {
                return None;
            }
            types.push(a.into());
            rest = tail.trim_start();
        } else {
            let end = rest.find(char::is_whitespace)?;
            types.push(rest[..end].into());
            rest = &rest[end..];
        }
    }
    Some((types, rest))
}

#[cfg(test)]
mod tests {
    use super::*;
    const TEST_DAT: &str = "test\nc3 12.01\n\nhn\nc3-c3 300.0 1.50\n\nc3-c3-c3 50.0 109.5\n\nX -c3-c3-X  2 1.5 30.0 -3\nX -c3-c3-X  1 0.4 180.0 2\n\nX -X -c -o 10.5 180.0 2\n\nEND\n";
    #[test]
    fn units_continuations_reversal_and_overrides() {
        let parameters = AmberParameters::from_dat(TEST_DAT, "test.dat", FFType::Gaff2).unwrap();
        let query = |s: &str| s.split('-').map(str::to_owned).collect::<Vec<_>>();
        let bond = parameters
            .resolve(TermKind::Bond, &query("c3-c3"), false)
            .unwrap();
        let Value::Bond(b) = bond.value else { panic!() };
        assert_eq!(b.k, 300.0 * 418.4);
        assert!((b.length - 0.15).abs() < 1e-12);
        let torsion = parameters
            .resolve(TermKind::Proper, &query("hc-c3-c3-h1"), false)
            .unwrap();
        let Value::Torsion(t) = torsion.value else {
            panic!()
        };
        assert_eq!(t.len(), 2);
        assert_eq!(t[0].periodicity, 3);
        assert!((t[0].k - 3.138).abs() < 1e-12);
        assert_eq!(torsion.resolution.method, ResolutionMethod::Wildcard);
        let updated = parameters
            .clone()
            .with_frcmod(
                "overrides\nBOND\nc3-c3 400 1.6\nDIHE\nh1-c3-c3-hc 1 2.0 37 1\n",
                "extra",
            )
            .unwrap();
        let r = updated
            .resolve(TermKind::Proper, &query("hc-c3-c3-h1"), false)
            .unwrap();
        let Value::Torsion(t) = r.value else { panic!() };
        assert_eq!(t.len(), 1);
        assert_eq!(r.resolution.method, ResolutionMethod::Exact);
        assert_eq!(r.resolution.origins[0].file, "extra");
        assert!(
            parameters
                .clone()
                .with_frcmod("bad\nBOND\nc3-c3 0 0 ATTN, need revision\n", "bad")
                .is_err()
        );
        assert!(
            parameters
                .clone()
                .with_frcmod("bad\nDIHE\nX-c3-c3-X 0 1 0 2\n", "bad")
                .is_err()
        );
        assert!(
            AmberParameters::from_dat(
                &TEST_DAT.replace("1 0.4 180.0 2", "1 0.4 180.0 -2"),
                "bad",
                FFType::Gaff
            )
            .is_err()
        );
        assert!(
            parameters
                .with_frcmod("bad\nANGLE\nc3-c3-c3 NaN 109\n", "bad")
                .is_err()
        );
    }
    #[test]
    fn automatic_estimates_have_provenance_and_unknowns_stay_unresolved() {
        let parameters = AmberParameters::builtin(FFType::Gaff2);
        let query = vec!["not-a-type".into(), "c3".into()];
        assert!(parameters.resolve(TermKind::Bond, &query, true).is_none());
        let mut small = AmberParameters::from_dat(TEST_DAT, "test", FFType::Gaff2).unwrap();
        small=small.with_frcmod("empirical angle\nBOND\nn3-c3 300 1.45\nc3-os 320 1.43\nANGLE\nn3-c3-n3 65 110\nos-c3-os 70 108\n","test").unwrap();
        let q = vec!["n3".into(), "c3".into(), "os".into()];
        assert!(small.resolve(TermKind::Angle, &q, false).is_none());
        let r = small.resolve(TermKind::Angle, &q, true).unwrap();
        assert_eq!(r.resolution.method, ResolutionMethod::EmpiricalAngle);
        let Value::Angle(p) = r.value else { panic!() };
        assert!((p.angle - 109.0_f64.to_radians()).abs() < 1e-12);
        let expected =
            143.9 * 1.38 * 1.06 * 1.59 * (-2.0 * ((1.45_f64 - 1.43) / (1.45 + 1.43)).powi(2)).exp()
                / ((1.45 + 1.43) * 109.0_f64.to_radians().sqrt())
                * 4.184;
        assert!((p.k - expected).abs() < 1e-10);
        let q = vec!["c5".into(), "c6".into()];
        let r = small.resolve(TermKind::Bond, &q, true).unwrap();
        assert_eq!(r.resolution.method, ResolutionMethod::Equivalent);
        assert_eq!(r.resolution.matched_types, vec!["c3", "c3"]);
    }
    #[test]
    fn zero_stiffness_requires_completion_and_preserves_equilibrium_length() {
        // A blank hydrophilic line must not swallow the first bond entry.
        let dat = TEST_DAT
            .replace("hn\n", "")
            .replace("c3-c3 300.0 1.50", "br-br 0.0 2.337");
        let parameters = AmberParameters::from_dat(&dat, "test", FFType::Gaff2).unwrap();
        let query = vec!["br".into(), "br".into()];
        assert!(parameters.resolve(TermKind::Bond, &query, false).is_none());
        let resolved = parameters.resolve(TermKind::Bond, &query, true).unwrap();
        assert_eq!(resolved.resolution.method, ResolutionMethod::EmpiricalBond);
        let Value::Bond(p) = resolved.value else {
            panic!()
        };
        assert!((p.length - 0.2337).abs() < 1e-12);
        let expected = (8.61 - 4.68 * 2.337_f64.ln()).exp() * 418.4;
        assert!((p.k - expected).abs() < 1e-9);
        assert!(
            resolved
                .resolution
                .origins
                .iter()
                .any(|p| p.file == "PARM_BLBA_GAFF2.DAT")
        );
    }
}
