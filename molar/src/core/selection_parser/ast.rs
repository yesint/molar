use num_traits::Bounded;
use regex::bytes::Regex;
use std::collections::HashSet;
use thiserror::Error;

use crate::prelude::*;

#[derive(Error, Debug)]
pub enum SelectionParserError {
    #[error("synatx error: {0}")]
    SyntaxError(String),

    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizes),

    #[error("no periodic box for pbc unwrapping")]
    PbcUnwrap,

    #[error("division by zero in math node eval")]
    DivisionByZero,

    #[error("sqrt of negative number in math node eval")]
    NegativeSqrt,

    #[error(transparent)]
    Measure(#[from] MeasureError),

    #[error("asked for atom {0} while inner expression selects {1} atoms")]
    OutOfBounds(usize, usize),
}

//##############################
//#  AST node types
//##############################

#[derive(Debug)]
pub(super) enum IntKeywordValue {
    Int(i32),
    IntRange(i32, i32),
}

#[derive(Debug)]
pub(super) enum StrKeywordValue {
    Str(String),
    Regex(Regex),
}

#[derive(Debug)]
pub(super) enum MathNode {
    Float(f32),
    Function(MathFunctionName, Box<Self>),
    X,
    Y,
    Z,
    Bfactor,
    Occupancy,
    Vdw,
    Mass,
    Charge,
    Plus(Box<Self>, Box<Self>),
    Minus(Box<Self>, Box<Self>),
    Mul(Box<Self>, Box<Self>),
    Div(Box<Self>, Box<Self>),
    Pow(Box<Self>, Box<Self>),
    Neg(Box<Self>),
    Dist(DistanceNode),
}

// Computes a vector value in various ways
#[derive(Debug)]
pub(super) enum VectorNode {
    Const(Pos),
    UnitConst(Pos),
    Com(Box<LogicalNode>, PbcDims),
    Cog(Box<LogicalNode>, PbcDims),
    NthAtomOf(Box<LogicalNode>, usize),
}

#[derive(Debug)]
pub(super) enum DistanceNode {
    Point(VectorNode, PbcDims),
    Line(VectorNode, VectorNode, PbcDims),
    LineDir(VectorNode, VectorNode, PbcDims),
    Plane(VectorNode, VectorNode, VectorNode, PbcDims),
    PlaneNormal(VectorNode, VectorNode, PbcDims),
}

pub(super) enum ComparisonOp {
    Eq,
    Neq,
    Leq,
    Lt,
    Geq,
    Gt,
}

#[derive(Debug)]
pub(super) enum ComparisonNode {
    // Simple
    Eq(MathNode, MathNode),
    Neq(MathNode, MathNode),
    Gt(MathNode, MathNode),
    Geq(MathNode, MathNode),
    Lt(MathNode, MathNode),
    Leq(MathNode, MathNode),
    // Chained left
    LtLt(MathNode, MathNode, MathNode),
    LeqLt(MathNode, MathNode, MathNode),
    LtLeq(MathNode, MathNode, MathNode),
    LeqLeq(MathNode, MathNode, MathNode),
    // Chained right
    GtGt(MathNode, MathNode, MathNode),
    GeqGt(MathNode, MathNode, MathNode),
    GtGeq(MathNode, MathNode, MathNode),
    GeqGeq(MathNode, MathNode, MathNode),
}

#[derive(Debug)]
pub(super) enum KeywordNode {
    Name(Vec<StrKeywordValue>),
    Resname(Vec<StrKeywordValue>),
    Chain(Vec<char>),

    Resid(Vec<IntKeywordValue>),
    Resindex(Vec<IntKeywordValue>),
    Index(Vec<IntKeywordValue>),
}

#[derive(Debug, PartialEq)]
pub(super) enum SameAttr {
    Residue,
    Chain,
}

#[derive(Debug, PartialEq)]
pub(super) struct WithinParams {
    pub(super) cutoff: f32,
    pub(super) pbc: PbcDims,
    pub(super) include_inner: bool,
}

#[derive(Debug)]
pub(super) enum LogicalNode {
    Not(Box<Self>),
    Or(Box<Self>, Box<Self>),
    And(Box<Self>, Box<Self>),
    Keyword(KeywordNode),
    Comparison(ComparisonNode),
    Same(SameAttr, Box<Self>),
    Within(WithinParams, Box<Self>),
    WithinPoint(WithinParams, VectorNode),
    All,
    Compound(CompoundNode),
}

pub(super) enum Keyword {
    Name,
    Resname,
    Resid,
    Resindex,
    Index,
    Occupancy,
    Bfactor,
    Chain,
    Residue,
}

#[derive(Debug, PartialEq, Clone)]
pub(super) enum MathFunctionName {
    Abs,
    Sqrt,
    Sin,
    Cos,
}

#[derive(Debug, PartialEq)]
pub(super) enum CompoundNode {
    Protein,
    Backbone,
    Sidechain,
    Water,
    NotWater,
    Hydrogen,
    NotHydrogen,
}

//##############################
//#  AST application stuff
//##############################

// Intermediate index type for applying AST
//type SubsetType = rustc_hash::FxHashSet<usize>;
pub(super) type SubsetType = Vec<usize>;

enum NodePayload {}

struct AstNode {
    payload: NodePayload,
    nodes: Vec<Box<AstNode>>,
}

#[derive(Clone)]
pub(super) struct EvaluationContext<'a> {
    topology: &'a Topology,
    state: &'a State,
    // This is a subset passed from outside
    // selection is completely within it
    // Parser is not changing subset
    global_subset: &'a SubsetType,
    // This is a context-dependent subset
    // (created for example by AND operation)
    // Selection can extend out of it
    // (like with WITHIN of BY)
    local_subset: Option<&'a SubsetType>,
}

impl<'a> EvaluationContext<'a> {
    pub(super) fn new(
        topology: &'a Topology,
        state: &'a State,
        global_subset: &'a SubsetType,
    ) -> Result<Self, SelectionParserError> {
        check_topology_state_sizes(topology, state)?;

        Ok(Self {
            topology,
            state,
            global_subset,
            local_subset: None,
        })
    }

    fn active_subset(&self) -> ActiveSubset {
        ActiveSubset {
            topology: &self.topology,
            state: &self.state,
            subset: self.local_subset.unwrap_or(self.global_subset),
        }
    }

    fn global_subset(&self) -> ActiveSubset {
        ActiveSubset {
            topology: &self.topology,
            state: &self.state,
            subset: self.global_subset,
        }
    }

    fn custom_subset(&self, custom: &'a SubsetType) -> ActiveSubset<'_> {
        ActiveSubset {
            topology: &self.topology,
            state: &self.state,
            subset: custom,
        }
    }

    fn clone_with_local_subset(&self, local_subset: &'a SubsetType) -> Self {
        Self {
            local_subset: Some(local_subset),
            ..*self
        }
    }
}

// Auxiliary struct representing a current active subset
// Created by ApplyData and contains eithr global or context subset
struct ActiveSubset<'a> {
    topology: &'a Topology,
    state: &'a State,
    subset: &'a SubsetType,
}

impl ActiveSubset<'_> {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>> {
        self.subset.iter().cloned().map(|i| unsafe {
            Particle {
                id: i,
                atom: self.topology.nth_atom_unchecked(i),
                pos: self.state.nth_pos_unchecked(i),
            }
        })
    }

    fn iter_ind_atom(&self) -> impl Iterator<Item = (usize, &Atom)> {
        self.subset
            .iter()
            .cloned()
            .map(|i| unsafe { (i, self.topology.nth_atom_unchecked(i)) })
    }
}

impl PosProvider for ActiveSubset<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.subset
            .iter()
            .map(|i| unsafe { self.state.nth_pos_unchecked(*i) })
    }
}

impl LenProvider for ActiveSubset<'_> {
    fn len(&self) -> usize {
        self.subset.len()
    }
}

impl RandomPosProvider for ActiveSubset<'_> {
    fn num_coords(&self) -> usize {
        self.subset.len()
    }
    
    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = *self.subset.get(i).unwrap();
        self.state.nth_pos_unchecked(ind)
    }
}

impl AtomProvider for ActiveSubset<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.subset
            .iter()
            .map(|i| unsafe { self.topology.nth_atom_unchecked(*i) })
    }
}

impl RandomAtomProvider for ActiveSubset<'_> {
    fn num_atoms(&self) -> usize {
        self.subset.len()
    }
    
    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = *self.subset.get(i).unwrap();
        self.topology.nth_atom_unchecked(ind)
    }
}

impl MassesProvider for ActiveSubset<'_> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.subset
            .iter()
            .map(|i| unsafe { self.topology.nth_atom_unchecked(*i).mass })
    }
}

impl BoxProvider for ActiveSubset<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl MeasurePos for ActiveSubset<'_> {}
impl MeasureMasses for ActiveSubset<'_> {}
impl MeasurePeriodic for ActiveSubset<'_> {}

//###################################
//#  AST nodes logic implementation
//###################################
impl VectorNode {
    fn get_vec(&self) -> &Pos {
        match self {
            Self::Const(v) => v,
            _ => unreachable!(),
        }
    }

    fn get_unit_vec(&self) -> &Pos {
        match self {
            Self::UnitConst(v) => v,
            _ => unreachable!(),
        }
    }

    fn precompute(&mut self, data: &EvaluationContext) -> Result<(), SelectionParserError> {
        match self {
            Self::Com(inner, dims) => {
                let res = inner.apply(data)?;
                let v = data.custom_subset(&res).center_of_mass_pbc_dims(*dims)?;
                *self = Self::Const(v);
            }
            Self::Cog(inner, dims) => {
                let res = inner.apply(data)?;
                let v = data
                    .custom_subset(&res)
                    .center_of_geometry_pbc_dims(*dims)?;
                *self = Self::Const(v);
            }
            Self::NthAtomOf(inner, i) => {
                let res = inner.apply(data)?;
                let v = data
                    .state
                    .nth_pos(*i)
                    .ok_or_else(|| SelectionParserError::OutOfBounds(*i, res.len()))?;
                *self = Self::Const(*v);
            }
            _ => (),
        }
        Ok(())
    }

    fn precompute_unit(&mut self, data: &EvaluationContext) -> Result<(), SelectionParserError> {
        match self {
            Self::Const(v) => {
                *self = Self::Const(Pos::from(v.coords.normalize()));
            }
            Self::UnitConst(_) => (),
            _ => {
                self.precompute(data)?;
                self.precompute_unit(data)?;
            }
        }
        Ok(())
    }
}

impl DistanceNode {
    pub fn precompute(&mut self, data: &EvaluationContext) -> Result<(), SelectionParserError> {
        match self {
            Self::Point(v1, _) => v1.precompute(data)?,
            Self::Line(v1, v2, _) => {
                v1.precompute(data)?;
                v2.precompute(data)?;
            }
            Self::Plane(v1, v2, v3, _) => {
                v1.precompute(data)?;
                v2.precompute(data)?;
                v3.precompute(data)?;
            }
            Self::LineDir(v, dir, _) => {
                v.precompute(data)?;
                dir.precompute_unit(data)?;
            }
            Self::PlaneNormal(v, n, _) => {
                v.precompute(data)?;
                n.precompute_unit(data)?;
            }
        }
        Ok(())
    }
}

impl LogicalNode {
    fn map_same_prop<T>(
        &self,
        data: &EvaluationContext,
        inner: &SubsetType,
        prop_fn: fn(&Atom) -> &T,
    ) -> SubsetType
    where
        T: Eq + std::hash::Hash + Copy,
    {
        // Collect all properties from the inner
        let mut properties = HashSet::<T>::new();
        let sub: ActiveSubset<'_> = data.custom_subset(inner);
        for at in sub.iter_atoms() {
            properties.insert(*prop_fn(at));
        }

        let mut res = SubsetType::default();
        // Now loop over *global* subset and add all atoms with the same property
        let sub = data.global_subset();
        for (i, at) in sub.iter_ind_atom() {
            for prop in properties.iter() {
                if prop_fn(at) == prop {
                    res.push(i);
                    break;
                }
            }
        }
        res
    }

    pub fn apply(&self, data: &EvaluationContext) -> Result<SubsetType, SelectionParserError> {
        use rustc_hash::FxHashSet;
        match self {
            Self::Not(node) => {
                // Here we always use global subset!
                let set1 = FxHashSet::<usize>::from_iter(data.global_subset.into_iter().cloned());
                let set2 = FxHashSet::<usize>::from_iter(node.apply(data)?.into_iter());
                Ok(set1.difference(&set2).cloned().collect())
            }

            Self::Or(a, b) => {
                let set1 = FxHashSet::<usize>::from_iter(a.apply(data)?.into_iter());
                let set2 = FxHashSet::<usize>::from_iter(b.apply(data)?.into_iter());
                Ok(set1.union(&set2).cloned().collect())
            }

            Self::And(a, b) => {
                let a_res = a.apply(data)?;
                // Create new instance of data and set a context subset to
                // the result of a
                let b_data = data.clone_with_local_subset(&a_res);

                let set1 = FxHashSet::<usize>::from_iter(a_res.iter().cloned());
                let set2 = FxHashSet::<usize>::from_iter(b.apply(&b_data)?.into_iter());
                Ok(set1.intersection(&set2).cloned().collect())
            }

            Self::Keyword(node) => node.apply(data),

            Self::Comparison(node) => node.apply(data),

            Self::Same(prop, node) => {
                let inner = node.apply(data)?;
                let res = match prop {
                    // Here we use the global subset!
                    SameAttr::Residue => self.map_same_prop(data, &inner, |at| &at.resindex),
                    SameAttr::Chain => self.map_same_prop(data, &inner, |at| &at.chain),
                };
                Ok(res)
            }

            Self::Within(prop, node) => {
                let inner = node.apply(data)?;

                let sub1 = data.global_subset();
                let sub2 = data.custom_subset(&inner);
                // Perform distance search
                let mut res: SubsetType = if prop.pbc == PBC_NONE {
                    // Non-periodic variant
                    // Find extents
                    let (mut lower, mut upper) = get_min_max(data.state, inner.iter().cloned());
                    lower.add_scalar_mut(-prop.cutoff - f32::EPSILON);
                    upper.add_scalar_mut(prop.cutoff + f32::EPSILON);

                    distance_search_within(
                        prop.cutoff, 
                        &sub1,
                        &sub2,
                        sub1.subset.iter().cloned(),
                        sub2.subset.iter().cloned(),
                        &lower,
                        &upper
                    )
                } else {
                    // Periodic variant
                    distance_search_within_pbc(
                        prop.cutoff,
                        &sub1,
                        &sub2,
                        sub1.subset.iter().cloned(),
                        sub2.subset.iter().cloned(),
                        data.state.get_box().unwrap(),
                        prop.pbc,
                    )
                };

                // Add inner if asked
                if prop.include_inner {
                    res.extend(inner);
                }
                Ok(res)
            }

            Self::WithinPoint(prop, point) => {
                let sub1 = data.global_subset();
                let pvec = point.get_vec();
                // Perform distance search
                let res: SubsetType = if prop.pbc == PBC_NONE {
                    // Non-periodic variant
                    // Find extents
                    let lower = pvec.coords.add_scalar(-prop.cutoff - f32::EPSILON);
                    let upper = pvec.coords.add_scalar(prop.cutoff + f32::EPSILON);

                    distance_search_within(
                        prop.cutoff, 
                        &sub1, 
                        pvec, 
                        sub1.subset.iter().cloned(),
                        0..1,
                        &lower, 
                        &upper
                    )
                } else {
                    // Periodic variant
                    distance_search_within_pbc(
                        prop.cutoff,
                        &sub1,
                        pvec,
                        sub1.subset.iter().cloned(),
                        0..1,
                        data.state.get_box().unwrap(),
                        prop.pbc,
                    )
                };
                Ok(res)
            }

            // All always works for global subset
            Self::All => Ok(data.global_subset.iter().cloned().collect()),

            Self::Compound(comp) => comp.apply(data),
        }
    }

    pub fn precompute(&mut self, data: &EvaluationContext) -> Result<(), SelectionParserError> {
        match self {
            Self::WithinPoint(_, p) => p.precompute(data)?,
            Self::Comparison(c) => c.precompute(data)?,
            _ => (),
        }
        Ok(())
    }

    pub fn is_coord_dependent(&self) -> bool {
        match self {
            Self::Not(node)
            | Self::Same(_, node) => node.is_coord_dependent(),
            
            Self::Comparison(node) => node.is_coord_dependent(),

            Self::Or(a, b)
            | Self::And(a, b) => a.is_coord_dependent() || b.is_coord_dependent(),

            Self::Within(_, _) 
            | Self::WithinPoint(_, _) => true,
                
            // All always works for global subset
            _ => false,
        }
    }
}

// impl VectorNode {
//     pub fn eval(&self, data: &ApplyData) -> Result<Pos, SelectionParserError> {
//         match self {
//             Self::Const(p) => Ok(*p),

//             Self::Com(pbc,node) => {
//                 let inner = node.apply(data)?;
//                 let c = if *pbc {
//                     data.custom_subset(&inner).center_of_mass().unwrap()
//                 } else {
//                     data.custom_subset(&inner).center_of_mass_pbc().unwrap()
//                 };
//                 Ok(c)
//             },

//             Self::Cog(pbc,node) => {
//                 let inner = node.apply(data)?;
//                 let c = if *pbc {
//                     data.custom_subset(&inner).center_of_geometry()
//                 } else {
//                     data.custom_subset(&inner).center_of_geometry_pbc()
//                 };
//                 Ok(c)
//             }
//         }
//     }
// }

impl CompoundNode {
    fn is_protein(atom: &Atom) -> bool {
        match atom.resname.as_str() {
            "GLY" | "ALA" | "VAL" | "PHE" | "PRO" | "MET" | "ILE" | "LEU" | "ASP" | "GLU"
            | "LYS" | "ARG" | "SER" | "THR" | "TYR" | "HIS" | "CYS" | "ASN" | "GLN" | "TRP"
            | "HSE" | "HSD" | "HSP" | "CYX" => true,
            _ => false,
        }
    }

    fn is_backbone(atom: &Atom) -> bool {
        if !Self::is_protein(atom) {
            return false;
        }
        match atom.name.as_str() {
            "C" | "N" | "O" | "CA" => true,
            _ => false,
        }
    }

    fn is_sidechain(atom: &Atom) -> bool {
        if !Self::is_protein(atom) {
            return false;
        }
        !Self::is_backbone(atom)
    }

    fn is_water(atom: &Atom) -> bool {
        match atom.resname.as_str() {
            "SOL" | "HOH" | "TIP3" | "TIP4" | "TIP5" => true,
            _ => false,
        }
    }

    fn is_hydrogen(atom: &Atom) -> bool {
        // Find first letter in file name
        if let Some(c) = atom.name.chars().find(char::is_ascii_alphabetic) {
            c == 'H'
        } else {
            false
        }
    }

    pub fn apply(&self, data: &EvaluationContext) -> Result<SubsetType, SelectionParserError> {
        let sub = data.active_subset();
        match self {
            Self::Protein => {
                let mut res = SubsetType::new();
                for (i, at) in sub.iter_ind_atom() {
                    if Self::is_protein(at) {
                        res.push(i)
                    }
                }
                Ok(res)
            }

            Self::Backbone => {
                let mut res = SubsetType::new();
                for (i, at) in sub.iter_ind_atom() {
                    if Self::is_backbone(at) {
                        res.push(i)
                    }
                }
                Ok(res)
            }

            Self::Sidechain => {
                let mut res = SubsetType::new();
                for (i, at) in sub.iter_ind_atom() {
                    if Self::is_sidechain(at) {
                        res.push(i)
                    }
                }
                Ok(res)
            }

            Self::Water => {
                let mut res = SubsetType::new();
                for (i, at) in sub.iter_ind_atom() {
                    if Self::is_water(at) {
                        res.push(i)
                    }
                }
                Ok(res)
            }

            Self::NotWater => {
                let mut res = SubsetType::new();
                for (i, at) in sub.iter_ind_atom() {
                    if !Self::is_water(at) {
                        res.push(i)
                    }
                }
                Ok(res)
            }

            Self::Hydrogen => {
                let mut res = SubsetType::new();
                for (i, at) in sub.iter_ind_atom() {
                    if Self::is_hydrogen(at) {
                        res.push(i)
                    }
                }
                Ok(res)
            }

            Self::NotHydrogen => {
                let mut res = SubsetType::new();
                for (i, at) in sub.iter_ind_atom() {
                    if !Self::is_hydrogen(at) {
                        res.push(i)
                    }
                }
                Ok(res)
            }
        }
    }
}

fn get_min_max(state: &State, iter: impl IndexIterator) -> (Vector3f, Vector3f) {
    let mut lower = Vector3f::max_value();
    let mut upper = Vector3f::min_value();
    for i in iter {
        let crd = unsafe { state.nth_pos_unchecked_mut(i) };
        for d in 0..3 {
            if crd[d] < lower[d] {
                lower[d] = crd[d]
            }
            if crd[d] > upper[d] {
                upper[d] = crd[d]
            }
        }
    }
    (lower, upper)
}

impl KeywordNode {
    fn map_str_values(
        &self,
        data: &EvaluationContext,
        values: &Vec<StrKeywordValue>,
        f: fn(&Atom) -> &String,
    ) -> SubsetType {
        let mut res = SubsetType::default();
        let sub = data.active_subset();
        for (ind, a) in sub.iter_ind_atom() {
            for val in values {
                match val {
                    StrKeywordValue::Str(s) => {
                        if s == f(a) {
                            res.push(ind);
                            break;
                        }
                    }
                    StrKeywordValue::Regex(r) => {
                        if r.is_match(f(a).as_bytes()) {
                            res.push(ind);
                            break;
                        }
                    }
                }
            }
        }
        res
    }

    fn map_int_values(
        &self,
        data: &EvaluationContext,
        values: &Vec<IntKeywordValue>,
        f: fn(&Atom, usize) -> i32,
    ) -> SubsetType {
        let mut res = SubsetType::default();
        let sub = data.active_subset();
        for (ind, a) in sub.iter_ind_atom() {
            for val in values {
                match *val {
                    IntKeywordValue::Int(v) => {
                        if v == f(a, ind) {
                            res.push(ind);
                            break;
                        }
                    }
                    IntKeywordValue::IntRange(b, e) => {
                        let val = f(a, ind);
                        if b <= val && val <= e {
                            res.push(ind);
                            break;
                        }
                    }
                }
            }
        }
        res
    }

    fn apply(&self, data: &EvaluationContext) -> Result<SubsetType, SelectionParserError> {
        match self {
            Self::Name(values) => Ok(self.map_str_values(data, values, |a| &a.name)),
            Self::Resname(values) => Ok(self.map_str_values(data, values, |a| &a.resname)),
            Self::Resid(values) => Ok(self.map_int_values(data, values, |a, _i| a.resid)),
            Self::Resindex(values) => {
                Ok(self.map_int_values(data, values, |a, _i| a.resindex as i32))
            }
            Self::Index(values) => Ok(self.map_int_values(data, values, |_a, i| i as i32)),
            Self::Chain(values) => {
                let mut res = SubsetType::default();
                let sub = data.active_subset();
                for (i, a) in sub.iter_ind_atom() {
                    for c in values {
                        if c == &a.chain {
                            res.push(i);
                        }
                    }
                }
                Ok(res)
            }
        }
    }
}

impl MathNode {
    fn closest_image(
        &self,
        point: &Pos,
        pbox: Option<&PeriodicBox>,
    ) -> Result<Option<Pos>, SelectionParserError> {
        match self {
            Self::Dist(d) => match d {
                DistanceNode::Point(target, dims)
                | DistanceNode::Line(target, _, dims)
                | DistanceNode::LineDir(target, _, dims)
                | DistanceNode::Plane(target, _, _, dims)
                | DistanceNode::PlaneNormal(target, _, dims) => {
                    if dims.any() {
                        Ok(Some(
                            pbox.ok_or_else(|| SelectionParserError::PbcUnwrap)?
                                .closest_image_dims(point, target.get_vec(), *dims),
                        ))
                    } else {
                        Ok(None)
                    }
                }
            },
            _ => Ok(None),
        }
    }

    fn eval(&self, atom: &Atom, pos: &Pos) -> Result<f32, SelectionParserError> {
        match self {
            Self::Float(v) => Ok(*v),
            Self::X => Ok(pos[0]),
            Self::Y => Ok(pos[1]),
            Self::Z => Ok(pos[2]),
            Self::Bfactor => Ok(atom.bfactor),
            Self::Occupancy => Ok(atom.occupancy),
            Self::Vdw => Ok(atom.vdw()),
            Self::Mass => Ok(atom.mass),
            Self::Charge => Ok(atom.charge),
            Self::Plus(a, b) => Ok(a.eval(atom, pos)? + b.eval(atom, pos)?),
            Self::Minus(a, b) => Ok(a.eval(atom, pos)? - b.eval(atom, pos)?),
            Self::Mul(a, b) => Ok(a.eval(atom, pos)? * b.eval(atom, pos)?),
            Self::Div(a, b) => {
                let b_val = b.eval(atom, pos)?;
                if b_val == 0.0 {
                    return Err(SelectionParserError::DivisionByZero);
                }
                Ok(a.eval(atom, pos)? / b_val)
            }
            Self::Pow(a, b) => Ok(a.eval(atom, pos)?.powf(b.eval(atom, pos)?)),
            Self::Neg(v) => Ok(-v.eval(atom, pos)?),
            Self::Function(func, v) => {
                let val = v.eval(atom, pos)?;
                match func {
                    MathFunctionName::Abs => Ok(val.abs()),
                    MathFunctionName::Sqrt => {
                        if val < 0.0 {
                            return Err(SelectionParserError::NegativeSqrt);
                        }
                        Ok(val.sqrt())
                    }
                    MathFunctionName::Sin => Ok(val.sin()),
                    MathFunctionName::Cos => Ok(val.cos()),
                }
            }
            Self::Dist(d) => {
                // Points are considered correctly unwrapped!
                match d {
                    DistanceNode::Point(p, _) => Ok((pos - p.get_vec()).norm()),
                    DistanceNode::Line(p1, p2, _) => {
                        let p1 = p1.get_vec();
                        let p2 = p2.get_vec();
                        let v = p2 - p1;
                        let w = pos - p1;
                        Ok((w - v * (w.dot(&v) / v.norm_squared())).norm())
                    }
                    DistanceNode::LineDir(p, dir, _) => {
                        let w = pos - p.get_vec();
                        let dir = dir.get_unit_vec();
                        Ok((w - dir.coords * w.dot(&dir.coords)).norm())
                    }
                    DistanceNode::Plane(p1, p2, p3, _) => {
                        let p1 = p1.get_vec();
                        let p2 = p2.get_vec();
                        let p3 = p3.get_vec();
                        // Plane normal
                        let n = (p2 - p1).cross(&(p3 - p1));
                        let w = pos - p1;
                        Ok((n * (w.dot(&n) / n.norm_squared())).norm())
                    }
                    DistanceNode::PlaneNormal(p, n, _) => {
                        let w = pos - p.get_vec();
                        let n = n.get_unit_vec();
                        Ok((n.coords * w.dot(&n.coords)).norm())
                    }
                }
            }
        }
    }

    fn is_coord_dependent(&self) -> bool {
        match self {
            Self::X | Self::Y | Self::Z | Self::Dist(_) => true,
            Self::Plus(a, b)
            | Self::Minus(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::Pow(a, b) => {
                a.is_coord_dependent() || b.is_coord_dependent()
            },
            Self::Neg(v) | Self::Function(_, v) => {
                v.is_coord_dependent()
            },
            _ => false,
        }
    }

    fn precompute(&mut self, data: &EvaluationContext) -> Result<(), SelectionParserError> {
        match self {
            Self::Dist(p) => p.precompute(data)?,
            _ => (),
        }
        Ok(())
    }
}

impl ComparisonNode {
    fn eval_op(
        data: &EvaluationContext,
        v1: &MathNode,
        v2: &MathNode,
        op: fn(f32, f32) -> bool,
    ) -> Result<SubsetType, SelectionParserError> {
        let mut res = SubsetType::default();
        let b = data.state.get_box();
        let sub = data.active_subset();

        for p in sub.iter_particle() {
            let pt1 = &v1.closest_image(p.pos, b)?.unwrap_or(*p.pos);
            let pt2 = &v2.closest_image(p.pos, b)?.unwrap_or(*p.pos);

            if op(v1.eval(p.atom, pt1)?, v2.eval(p.atom, pt2)?) {
                res.push(p.id);
            }
        }
        Ok(res)
    }

    fn eval_op_chained(
        data: &EvaluationContext,
        v1: &MathNode,
        v2: &MathNode,
        v3: &MathNode,
        op1: fn(f32, f32) -> bool,
        op2: fn(f32, f32) -> bool,
    ) -> Result<SubsetType, SelectionParserError> {
        let mut res = SubsetType::default();
        let b = data.state.get_box();
        let sub = data.active_subset();

        for p in sub.iter_particle() {
            let pt1 = &v1.closest_image(p.pos, b)?.unwrap_or(*p.pos);
            let pt2 = &v2.closest_image(p.pos, b)?.unwrap_or(*p.pos);
            let pt3 = &v3.closest_image(p.pos, b)?.unwrap_or(*p.pos);

            let mid = v2.eval(p.atom, pt2)?;
            if op1(v1.eval(p.atom, pt1)?, mid) && op2(mid, v3.eval(p.atom, pt3)?) {
                res.push(p.id);
            }
        }
        Ok(res)
    }

    fn apply(&self, data: &EvaluationContext) -> Result<SubsetType, SelectionParserError> {
        match self {
            // Simple
            Self::Eq(v1, v2) => Self::eval_op(data, v1, v2, |a, b| a == b),
            Self::Neq(v1, v2) => Self::eval_op(data, v1, v2, |a, b| a != b),
            Self::Gt(v1, v2) => Self::eval_op(data, v1, v2, |a, b| a > b),
            Self::Geq(v1, v2) => Self::eval_op(data, v1, v2, |a, b| a >= b),
            Self::Lt(v1, v2) => Self::eval_op(data, v1, v2, |a, b| a < b),
            Self::Leq(v1, v2) => Self::eval_op(data, v1, v2, |a, b| a <= b),
            // Chained left
            Self::LtLt(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a < b, |a, b| a < b)
            }
            Self::LtLeq(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a < b, |a, b| a <= b)
            }
            Self::LeqLt(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a <= b, |a, b| a < b)
            }
            Self::LeqLeq(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a <= b, |a, b| a <= b)
            }
            // Chained right
            Self::GtGt(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a > b, |a, b| a > b)
            }
            Self::GtGeq(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a > b, |a, b| a >= b)
            }
            Self::GeqGt(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a >= b, |a, b| a > b)
            }
            Self::GeqGeq(v1, v2, v3) => {
                Self::eval_op_chained(data, v1, v2, v3, |a, b| a >= b, |a, b| a >= b)
            }
        }
    }

    fn precompute(&mut self, data: &EvaluationContext) -> Result<(), SelectionParserError> {
        match self {
            // Simple
            Self::Eq(v1, v2)
            | Self::Neq(v1, v2)
            | Self::Gt(v1, v2)
            | Self::Geq(v1, v2)
            | Self::Lt(v1, v2)
            | Self::Leq(v1, v2) => {
                v1.precompute(data)?;
                v2.precompute(data)?;
            }
            // Chained
            Self::LtLt(v1, v2, v3)
            | Self::LtLeq(v1, v2, v3)
            | Self::LeqLt(v1, v2, v3)
            | Self::LeqLeq(v1, v2, v3)
            | Self::GtGt(v1, v2, v3)
            | Self::GtGeq(v1, v2, v3)
            | Self::GeqGt(v1, v2, v3)
            | Self::GeqGeq(v1, v2, v3) => {
                v1.precompute(data)?;
                v2.precompute(data)?;
                v3.precompute(data)?;
            }
        }
        Ok(())
    }
    
    fn is_coord_dependent(&self) -> bool {
        match self {
            // Simple
            Self::Eq(v1, v2)
            | Self::Neq(v1, v2)
            | Self::Gt(v1, v2)
            | Self::Geq(v1, v2)
            | Self::Lt(v1, v2)
            | Self::Leq(v1, v2) => v1.is_coord_dependent() || v2.is_coord_dependent(),
            // Chained
            Self::LtLt(v1, v2, v3)
            | Self::LtLeq(v1, v2, v3)
            | Self::LeqLt(v1, v2, v3)
            | Self::LeqLeq(v1, v2, v3)
            | Self::GtGt(v1, v2, v3)
            | Self::GtGeq(v1, v2, v3)
            | Self::GeqGt(v1, v2, v3)
            | Self::GeqGeq(v1, v2, v3) => v1.is_coord_dependent() || v2.is_coord_dependent() || v3.is_coord_dependent(),
        }
    }
}
