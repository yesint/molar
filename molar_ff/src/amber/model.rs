use super::{parameters::*, resolve::Resolved};
use crate::{ApplyFF, FFError, FFType};
use molar::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Automatic completion is the default. Strict mode is useful for validation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ParameterPolicy {
    #[default]
    Automatic,
    Strict,
}
#[derive(Clone, Debug)]
pub struct BondedTerm<const N: usize, T> {
    /// Local coordinate indices. Use `BondedModel::atom_indices` for source indices.
    pub atoms: [usize; N],
    pub parameter: T,
    pub resolution: Resolution,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MissingParameter {
    pub kind: TermKind,
    /// Source/global atom indices, suitable for a UI diagnostic.
    pub atoms: Vec<usize>,
    pub types: Vec<String>,
}
#[derive(Debug, thiserror::Error)]
pub enum BondedError {
    #[error(transparent)]
    Typing(#[from] FFError),
    #[error("atom {0} has no force-field type")]
    UntypedAtom(usize),
    #[error("invalid bond {0}-{1}: self-bond or duplicate bond")]
    InvalidBond(usize, usize),
    #[error("selection is not bond-complete: {0}-{1}")]
    OpenSelection(usize, usize),
    #[error("could not resolve bonded parameters: {0:?}")]
    Unresolved(Vec<MissingParameter>),
}
/// Immutable bonded model. Construction is atomic: unresolved interactions cannot
/// produce a partially usable model. Parameters are resolved once, outside minimization.
#[derive(Clone, Debug)]
pub struct BondedModel {
    pub(crate) atom_indices: Vec<usize>,
    pub(crate) atom_types: Vec<String>,
    pub(crate) bonds: Vec<BondedTerm<2, BondParameter>>,
    pub(crate) angles: Vec<BondedTerm<3, AngleParameter>>,
    pub(crate) propers: Vec<BondedTerm<4, Vec<FourierTerm>>>,
    pub(crate) impropers: Vec<BondedTerm<4, Vec<FourierTerm>>>,
}
impl BondedModel {
    /// Type and parameterize an entire molecule or bond-complete bound selection.
    /// Neither input atom types nor coordinates are changed. Explicit hydrogens and
    /// known bond orders are required, just as for `ApplyFF`.
    pub fn gaff<S: AtomProvider + BondProvider>(
        source: &S,
        ff: FFType,
    ) -> Result<Self, BondedError> {
        let (indices, pairs) = graph(source)?;
        let map: HashMap<_, _> = indices.iter().enumerate().map(|(i, &g)| (g, i)).collect();
        let mut top = Topology::default();
        for atom in source.iter_atoms() {
            top.atoms.push(&Atom::from(&atom));
        }
        for b in source.iter_bonds() {
            if let (Some(&a), Some(&c)) = (map.get(&b.pair()[0]), map.get(&b.pair()[1])) {
                top.bonds.push(&Bond::with_order(a, c, b.order()));
            }
        }
        top.apply_ff(ff)?;
        let types = top
            .iter_atoms()
            .map(|a| a.get_type_name().unwrap().to_owned())
            .collect();
        Self::build(
            indices,
            types,
            pairs,
            AmberParameters::builtin(ff),
            ParameterPolicy::Automatic,
        )
    }
    /// Parameterize types already stored on the source. Automatic completion does
    /// not require files, an AmberTools installation, or caller-supplied parameters.
    pub fn from_typed<S: AtomProvider + BondProvider>(
        source: &S,
        parameters: &AmberParameters,
        policy: ParameterPolicy,
    ) -> Result<Self, BondedError> {
        let (indices, pairs) = graph(source)?;
        let types = source
            .iter_atoms()
            .zip(&indices)
            .map(|(a, &i)| {
                a.get_type_name()
                    .map(str::to_owned)
                    .ok_or(BondedError::UntypedAtom(i))
            })
            .collect::<Result<_, _>>()?;
        Self::build(indices, types, pairs, parameters, policy)
    }
    fn build(
        atom_indices: Vec<usize>,
        atom_types: Vec<String>,
        pairs: Vec<[usize; 2]>,
        parameters: &AmberParameters,
        policy: ParameterPolicy,
    ) -> Result<Self, BondedError> {
        let mut model = Self {
            atom_indices,
            atom_types,
            bonds: vec![],
            angles: vec![],
            propers: vec![],
            impropers: vec![],
        };
        let mut adj = vec![vec![]; model.atom_indices.len()];
        for &[a, b] in &pairs {
            adj[a].push(b);
            adj[b].push(a);
        }
        for neighbors in &mut adj {
            neighbors.sort_unstable();
        }
        let mut requested: Vec<(TermKind, Vec<usize>)> =
            pairs.iter().map(|p| (TermKind::Bond, p.to_vec())).collect();
        for (center, neighbors) in adj.iter().enumerate() {
            for (i, &a) in neighbors.iter().enumerate() {
                for &b in &neighbors[i + 1..] {
                    requested.push((TermKind::Angle, vec![a, center, b]));
                }
            }
            if neighbors.len() == 3
                && parameters
                    .rules
                    .atoms
                    .get(&model.atom_types[center])
                    .is_some_and(|r| r.improper)
            {
                let mut outer = neighbors.clone();
                outer.sort_by_key(|&i| (&model.atom_types[i], model.atom_indices[i]));
                requested.push((
                    TermKind::Improper,
                    vec![outer[0], outer[1], center, outer[2]],
                ));
            }
        }
        for &[b, c] in &pairs {
            for &a in &adj[b] {
                for &d in &adj[c] {
                    if a != c && d != b && a != d {
                        requested.push((TermKind::Proper, vec![a, b, c, d]));
                    }
                }
            }
        }
        let mut cache: BTreeMap<(TermKind, Vec<String>), Option<Resolved>> = BTreeMap::new();
        let mut missing = vec![];
        for (kind, atoms) in requested {
            let types: Vec<_> = atoms.iter().map(|&i| model.atom_types[i].clone()).collect();
            let result = cache.entry((kind, types.clone())).or_insert_with(|| {
                parameters.resolve(kind, &types, policy == ParameterPolicy::Automatic)
            });
            let Some(r) = result else {
                missing.push(MissingParameter {
                    kind,
                    atoms: atoms.iter().map(|&i| model.atom_indices[i]).collect(),
                    types,
                });
                continue;
            };
            let ordered: Vec<_> = r.order.iter().map(|&i| atoms[i]).collect();
            let resolution = r.resolution.clone();
            match &r.value {
                Value::Bond(p) => model.bonds.push(BondedTerm {
                    atoms: ordered.try_into().unwrap(),
                    parameter: *p,
                    resolution,
                }),
                Value::Angle(p) => model.angles.push(BondedTerm {
                    atoms: ordered.try_into().unwrap(),
                    parameter: *p,
                    resolution,
                }),
                Value::Torsion(p) => {
                    let term = BondedTerm {
                        atoms: ordered.try_into().unwrap(),
                        parameter: p.clone(),
                        resolution,
                    };
                    if kind == TermKind::Proper {
                        model.propers.push(term);
                    } else {
                        model.impropers.push(term);
                    }
                }
            }
        }
        if !missing.is_empty() {
            return Err(BondedError::Unresolved(missing));
        }
        Ok(model)
    }
    pub fn atom_indices(&self) -> &[usize] {
        &self.atom_indices
    }
    pub fn atom_types(&self) -> &[String] {
        &self.atom_types
    }
    pub fn bonds(&self) -> &[BondedTerm<2, BondParameter>] {
        &self.bonds
    }
    pub fn angles(&self) -> &[BondedTerm<3, AngleParameter>] {
        &self.angles
    }
    pub fn propers(&self) -> &[BondedTerm<4, Vec<FourierTerm>>] {
        &self.propers
    }
    pub fn impropers(&self) -> &[BondedTerm<4, Vec<FourierTerm>>] {
        &self.impropers
    }
    pub fn estimated_terms(&self) -> impl Iterator<Item = (TermKind, &[usize], &Resolution)> {
        self.bonds
            .iter()
            .map(|t| (TermKind::Bond, t.atoms.as_slice(), &t.resolution))
            .chain(
                self.angles
                    .iter()
                    .map(|t| (TermKind::Angle, t.atoms.as_slice(), &t.resolution)),
            )
            .chain(
                self.propers
                    .iter()
                    .map(|t| (TermKind::Proper, t.atoms.as_slice(), &t.resolution)),
            )
            .chain(
                self.impropers
                    .iter()
                    .map(|t| (TermKind::Improper, t.atoms.as_slice(), &t.resolution)),
            )
            .filter(|(_, _, r)| r.is_estimated())
    }
}
fn graph<S: AtomProvider + BondProvider>(
    source: &S,
) -> Result<(Vec<usize>, Vec<[usize; 2]>), BondedError> {
    let indices: Vec<_> = source.iter_index().collect();
    let map: HashMap<_, _> = indices.iter().enumerate().map(|(i, &g)| (g, i)).collect();
    let mut pairs = BTreeSet::new();
    for b in source.iter_bonds() {
        let [a, c] = b.pair();
        match (map.get(&a), map.get(&c)) {
            (Some(&i), Some(&j)) => {
                let pair = [i.min(j), i.max(j)];
                if i == j || !pairs.insert(pair) {
                    return Err(BondedError::InvalidBond(a, c));
                }
            }
            (None, None) => (),
            _ => return Err(BondedError::OpenSelection(a, c)),
        }
    }
    Ok((indices, pairs.into_iter().collect()))
}
