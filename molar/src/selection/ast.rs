use crate::{prelude::*, selection::sel};
use log::Log;
use num_traits::Bounded;
use regex::bytes::Regex;
use std::{borrow::Cow, collections::HashSet};
use thiserror::Error;

use super::utils::check_topology_state_sizes;

//##############################
//#  AST node types
//##############################

#[derive(Debug,Clone)]
pub(super) enum IntKeywordArg {
    Int(isize),
    IntRange(isize, isize),
}

#[derive(Debug,Clone)]
pub(super) enum StrKeywordArg {
    Str(String),
    Regex(Regex),
}

#[derive(Debug,Clone)]
pub(super) enum BinaryOperator {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
}

#[derive(Debug,Clone)]
pub(super) enum MathNode {
    Float(f32),
    Function(MathFunctionName, Box<Self>),
    X,
    Y,
    Z,
    Xof(VectorNode),
    Yof(VectorNode),
    Zof(VectorNode),
    Bfactor,
    Occupancy,
    Vdw,
    Mass,
    Charge,
    BinaryOp(Box<Self>, BinaryOperator, Box<Self>),
    UnaryMinus(Box<Self>),
    Dist(DistanceNode),
}

// Computes a vector value in various ways
#[derive(Debug,Clone)]
pub(super) enum VectorNode {
    Const(Pos),
    UnitConst(Pos),
    Com(Box<LogicalNode>, PbcDims),
    Cog(Box<LogicalNode>, PbcDims),
    NthAtomOf(Box<LogicalNode>, usize),
}

#[derive(Debug,Clone)]
pub(super) enum DistanceNode {
    Point(VectorNode, PbcDims),
    Line(VectorNode, VectorNode, PbcDims),
    LineDir(VectorNode, VectorNode, PbcDims),
    Plane(VectorNode, VectorNode, VectorNode, PbcDims),
    PlaneNormal(VectorNode, VectorNode, PbcDims),
}

#[derive(Debug,Clone)]
pub(super) enum ComparisonOp {
    Eq,
    Neq,
    Leq,
    Lt,
    Geq,
    Gt,
}

#[derive(Debug,Clone)]
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

#[derive(Debug,Clone)]
pub(super) enum KeywordNode {
    // String keywords
    Name(Vec<StrKeywordArg>),
    Resname(Vec<StrKeywordArg>),
    Chain(Vec<char>),
    // Int keywords
    Resid(Vec<IntKeywordArg>),
    Resindex(Vec<IntKeywordArg>),
    Index(Vec<IntKeywordArg>),
}

#[derive(Debug, PartialEq, Clone)]
pub(super) enum SameAttr {
    Residue,
    Chain,
}

#[derive(Debug, PartialEq, Clone)]
pub(super) struct WithinParams {
    pub(super) cutoff: f32,
    pub(super) pbc: PbcDims,
    pub(super) include_inner: bool,
}

// Top level
#[derive(Clone)]
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
    Chemical(ChemicalNode),
    Precomputed(Vec<usize>),
}

#[derive(Debug,Clone)]
pub(super) enum AtomAttr {
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

#[derive(Debug, PartialEq, Clone)]
pub(super) enum ChemicalNode {
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

#[derive(Clone)]
enum EvalSubset<'a> {
    Whole(std::ops::Range<usize>),
    Part(&'a [usize]),
}

enum EvalSubsetIter<'a> {
    Whole(std::ops::Range<usize>),
    Part(std::iter::Copied<std::slice::Iter<'a, usize>>),
}

impl<'a> Iterator for EvalSubsetIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            EvalSubsetIter::Whole(r) => r.next(),
            EvalSubsetIter::Part(it) => it.next(),
        }
    }
}

impl EvalSubset<'_> {
    fn iter(&self) -> EvalSubsetIter<'_> {
        match self {
            Self::Whole(r) => EvalSubsetIter::Whole(r.clone()),
            Self::Part(ind) => EvalSubsetIter::Part(ind.iter().copied()),
        }
    }

    fn len(&self) -> usize {
        match self {
            EvalSubset::Whole(r) => r.len(),
            EvalSubset::Part(ind) => ind.len(),
        }
    }

    unsafe fn get_unchecked(&self, i: usize) -> usize {
        match self {
            EvalSubset::Whole(_) => i,
            EvalSubset::Part(ind) => *ind.get_unchecked(i),
        }
    }
}

#[derive(Clone)]
pub(super) struct EvaluationContext<'a> {
    topology: &'a Topology,
    state: &'a State,
    // This is a subset passed from outside.
    // Selection is completely within it.
    // Parser is never changing it.
    global_subset: EvalSubset<'a>,
    cur_subset: EvalSubset<'a>,
}

impl<'a> EvaluationContext<'a> {
    pub(super) fn new_whole(
        topology: &'a Topology,
        state: &'a State,
    ) -> Result<Self, SelectionParserError> {
        check_topology_state_sizes(topology, state)?;

        Ok(Self {
            topology,
            state,
            global_subset: EvalSubset::Whole(0..state.len()),
            cur_subset: EvalSubset::Whole(0..state.len()),
        })
    }

    pub(super) fn new_part(
        topology: &'a Topology,
        state: &'a State,
        part: &'a [usize],
    ) -> Result<Self, SelectionParserError> {
        check_topology_state_sizes(topology, state)?;

        Ok(Self {
            topology,
            state,
            global_subset: EvalSubset::Part(part),
            cur_subset: EvalSubset::Part(part),
        })
    }

    fn with_custom_subset(&'a self, custom_subset: &'a [usize]) -> Self {
        Self {
            cur_subset: EvalSubset::Part(custom_subset),
            ..self.clone()
        }
    }

    fn with_global_subset(&'a self) -> Self {
        Self {
            cur_subset: self.global_subset.clone(),
            ..self.clone()
        }
    }

    fn iter(&self) -> EvalSubsetIter<'_> {
        self.cur_subset.iter()
    }
}


impl EvaluationContext<'_> {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>> {
        self.iter().map(|i| unsafe {
            Particle {
                id: i,
                atom: self.topology.get_atom_unchecked(i),
                pos: self.state.get_pos_unchecked(i),
            }
        })
    }

    fn iter_ind_atom(&self) -> impl Iterator<Item = (usize, &Atom)> {
        self.iter()
            .map(|i| unsafe { (i, self.topology.get_atom_unchecked(i)) })
    }
}

impl PosIterProvider for EvaluationContext<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.iter()
            .map(|i| unsafe { self.state.get_pos_unchecked(i) })
    }
}

impl LenProvider for EvaluationContext<'_> {
    fn len(&self) -> usize {
        self.cur_subset.len()
    }
}

impl RandomPosProvider for EvaluationContext<'_> {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.cur_subset.get_unchecked(i);
        self.state.get_pos_unchecked(ind)
    }
}

impl AtomIterProvider for EvaluationContext<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self
            .iter()
            .map(|i| unsafe { self.topology.get_atom_unchecked(i) })
    }
}

impl RandomAtomProvider for EvaluationContext<'_> {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.cur_subset.get_unchecked(i);
        self.topology.get_atom_unchecked(ind)
    }
}

impl BoxProvider for EvaluationContext<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl MeasurePos for EvaluationContext<'_> {}
impl MeasureMasses for EvaluationContext<'_> {}
impl MeasurePeriodic for EvaluationContext<'_> {}

//###################################
//#  AST nodes logic implementation
//###################################

pub(super) trait Evaluate {
    fn is_state_dependent(&self) -> bool;
    fn apply<'a>(&'a mut self, data: &'a EvaluationContext)
        -> Result<Cow<'a, [usize]>, SelectionParserError>;
}

impl DistanceNode {
    fn closest_image(&mut self, point: &mut Pos, data: &EvaluationContext) {
        if let Some(pbox) = data.state.get_box() {
            match self {
                Self::Point(target, dims)
                | Self::Line(target, _, dims)
                | Self::LineDir(target, _, dims)
                | Self::Plane(target, _, _, dims)
                | Self::PlaneNormal(target, _, dims) => {
                    if dims.any() {
                        *point =
                            pbox.closest_image_dims(point, target.get_vec(data).unwrap(), *dims)
                    }
                }
            }
        }
    }
}

impl VectorNode {
    fn get_vec(&mut self, data: &EvaluationContext) -> Result<&Pos, SelectionParserError> {
        match self {
            Self::Const(v) => Ok(v),
            Self::UnitConst(v) => Ok(v),
            Self::Com(inner, dims) => {
                let res = inner.apply(data)?;
                let v = if *dims == PBC_NONE {
                    data.with_custom_subset(&res).center_of_mass()?
                } else {
                    data.with_custom_subset(&res).center_of_mass_pbc_dims(*dims)?
                };
                *self = Self::Const(v);
                self.get_vec(data)
            }
            Self::Cog(inner, dims) => {
                let res = inner.apply(data)?;
                let v = if *dims == PBC_NONE {
                    data.with_custom_subset(&res).center_of_geometry()
                } else {
                    data.with_custom_subset(&res)
                        .center_of_geometry_pbc_dims(*dims)?
                };
                *self = Self::Const(v);
                self.get_vec(data)
            }
            Self::NthAtomOf(inner, i) => {
                let res = inner.apply(data)?;
                let v = data
                    .state
                    .get_pos(*i)
                    .ok_or_else(|| SelectionParserError::OutOfBounds(*i, res.len()))?;
                *self = Self::Const(*v);
                self.get_vec(data)
            }
        }
    }

    fn get_unit_vec(&mut self, data: &EvaluationContext) -> Result<&Pos, SelectionParserError> {
        match self {
            Self::UnitConst(v) => Ok(v),
            _ => {
                *self = Self::UnitConst(Pos::from(self.get_vec(data)?.coords.normalize()));
                self.get_unit_vec(data)
            }
        }
    }

    fn is_state_dependent(&self) -> bool {
        match self {
            Self::Cog(node, _) | Self::Com(node, _) | Self::NthAtomOf(node,_) => node.is_state_dependent(),
            Self::Const(_) | Self::UnitConst(_) => false,
        }
    }
}

impl LogicalNode {
    fn map_same_attr<T>(
        data: &EvaluationContext,
        inner: &[usize],
        prop_fn: fn(&Atom) -> &T,
    ) -> Vec<usize>
    where
        T: Eq + std::hash::Hash + Copy,
    {
        // Collect all properties from the inner
        let mut properties = HashSet::<T>::new();
        let sub = data.with_custom_subset(inner);
        for at in sub.iter_atoms() {
            properties.insert(*prop_fn(at));
        }

        let mut res = vec![];
        // Now loop over *global* subset and add all atoms with the same property
        let sub = data.with_global_subset();
        for (i, at) in sub.iter_ind_atom() {
            let cur_prop = prop_fn(at);
            if properties.contains(cur_prop) {
                res.push(i);
            }
        }
        res
    }
}

use std::fmt;

impl fmt::Debug for LogicalNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LogicalNode::Precomputed(v) => {
                if v.len() > 10 {
                    write!(f, "Precomputed {:?}...", &v[..10])
                } else {
                    write!(f, "Precomputed {:?}", v)
                }
            }

            LogicalNode::Not(a) => f.debug_tuple("Not").field(a).finish(),
            LogicalNode::Or(a, b) => f.debug_tuple("Or").field(a).field(b).finish(),
            LogicalNode::And(a, b) => f.debug_tuple("And").field(a).field(b).finish(),
            LogicalNode::Keyword(k) => f.debug_tuple("Keyword").field(k).finish(),
            LogicalNode::Comparison(c) => f.debug_tuple("Comparison").field(c).finish(),
            LogicalNode::Same(attr, n) => f.debug_tuple("Same").field(attr).field(n).finish(),
            LogicalNode::Within(p, n) => f.debug_tuple("Within").field(p).field(n).finish(),
            LogicalNode::WithinPoint(p, v) => f.debug_tuple("WithinPoint").field(p).field(v).finish(),
            LogicalNode::All => f.write_str("All"),
            LogicalNode::Chemical(c) => f.debug_tuple("Chemical").field(c).finish(),
        }
    }
}

impl Evaluate for LogicalNode {
    fn is_state_dependent(&self) -> bool {
        match self {
            Self::All | Self::Chemical(_) | Self::Keyword(_) | Self::Precomputed(_) => false,
            Self::Within(_, _) => true,
            Self::WithinPoint(_, _) => true,
            Self::Not(v) | Self::Same(_, v) => v.is_state_dependent(),
            Self::And(a, b) | Self::Or(a, b) => a.is_state_dependent() && b.is_state_dependent(),
            Self::Comparison(v) => v.is_state_dependent(),
        }
    }

    fn apply<'a>(
        &'a mut self,
        data: &'a EvaluationContext,
    ) -> Result<Cow<'a, [usize]>, SelectionParserError> {
        use rustc_hash::FxHashSet;
        match self {
            Self::Precomputed(v) => Ok(Cow::from(&*v)),

            Self::Not(node) => {
                // Here we always use global subset!
                let set1 = FxHashSet::from_iter(data.iter());
                let set2 = FxHashSet::from_iter(node.apply(data)?.into_iter().cloned());
                let res = set1.difference(&set2).cloned().collect();
                
                if node.is_state_dependent() {
                    Ok(Cow::from(res)) // Owned res
                } else {
                    *self = Self::Precomputed(res);
                    self.apply(data) // will give borrowed res
                }
            }

            Self::Or(a, b) => {
                let set1 = FxHashSet::from_iter(a.apply(data)?.into_iter().cloned());
                let set2 = FxHashSet::from_iter(b.apply(data)?.into_iter().cloned());
                let res = set1.union(&set2).cloned().collect();
                
                if a.is_state_dependent() || b.is_state_dependent() {
                    Ok(Cow::from(res)) // Owned res
                } else {
                    *self = Self::Precomputed(res);
                    self.apply(data) // will give borrowed res
                }
            }

            Self::And(a, b) => {
                let a_res = a.apply(data)?;
                // Create new instance of data and set a context subset to
                // the result of a.
                let b_data = data.with_custom_subset(&a_res);

                let set1 = FxHashSet::from_iter(a_res.iter().cloned());
                let set2 = FxHashSet::from_iter(b.apply(&b_data)?.into_iter().cloned());
                let res = set1.intersection(&set2).cloned().collect();
                
                if a.is_state_dependent() || b.is_state_dependent() {
                    Ok(Cow::from(res)) // Owned res
                } else {
                    *self = Self::Precomputed(res);
                    self.apply(data) // will give borrowed res
                }
            }

            Self::Keyword(node) => { 
                *self = Self::Precomputed(node.apply(data)?.into_owned());
                self.apply(data) // will give borrowed res
            }

            Self::Comparison(node) => {
                let res = node.apply(data)?.into_owned();
                if node.is_state_dependent() {
                    Ok(Cow::from(res)) // Owned res
                } else {
                    *self = Self::Precomputed(res);
                    self.apply(data) // will give borrowed res
                }
            }

            Self::Same(attr, node) => {
                let inner_res = node.apply(data)?;
                let res = match attr {
                    // Here we use the global subset!
                    SameAttr::Residue => Self::map_same_attr(data, &inner_res, |at| &at.resindex),
                    SameAttr::Chain => Self::map_same_attr(data, &inner_res, |at| &at.chain),
                };
                
                if node.is_state_dependent() {
                    Ok(Cow::from(res)) // Owned res
                } else {
                    *self = Self::Precomputed(res);
                    self.apply(data) // will give borrowed res
                }
            }

            Self::Within(params, node) => {
                // Inner expr have to be evaluated in global context!
                let glob_data = data.with_global_subset();
                let inner_res = node.apply(&glob_data)?;

                let sub1 = data;
                let sub2 = data.with_custom_subset(&inner_res);
                // Perform distance search
                let mut res: Vec<usize> = if params.pbc == PBC_NONE {
                    // Non-periodic variant
                    // Find extents
                    let (mut lower, mut upper) = get_min_max(data.state, inner_res.iter().cloned());
                    lower.add_scalar_mut(-params.cutoff - f32::EPSILON);
                    upper.add_scalar_mut(params.cutoff + f32::EPSILON);

                    distance_search_within(
                        params.cutoff,
                        sub1,
                        &sub2,
                        sub1.iter(),
                        sub2.iter(),
                        &lower,
                        &upper,
                    )
                } else {
                    // Periodic variant
                    distance_search_within_pbc(
                        params.cutoff,
                        sub1,
                        &sub2,
                        sub1.iter(),
                        sub2.iter(),
                        data.state.require_box()?,
                        params.pbc,
                    )
                };

                // Add inner if asked
                if params.include_inner {
                    res.extend(inner_res.iter());
                }
                Ok(Cow::from(res))
            }

            Self::WithinPoint(prop, point) => {
                let sub1 = data.with_global_subset();
                let pvec = point.get_vec(data)?;
                // Perform distance search
                let res: Vec<usize> = if prop.pbc == PBC_NONE {
                    // Non-periodic variant
                    // Find extents
                    let lower = pvec.coords.add_scalar(-prop.cutoff - f32::EPSILON);
                    let upper = pvec.coords.add_scalar(prop.cutoff + f32::EPSILON);

                    distance_search_within(
                        prop.cutoff,
                        &sub1,
                        pvec,
                        sub1.iter(),
                        0..1,
                        &lower,
                        &upper,
                    )
                } else {
                    // Periodic variant
                    distance_search_within_pbc(
                        prop.cutoff,
                        &sub1,
                        pvec,
                        sub1.iter(),
                        0..1,
                        data.state.require_box()?,
                        prop.pbc,
                    )
                };
                Ok(Cow::from(res))
            }

            // All always works in global subset
            Self::All => Ok(Cow::from_iter(data.with_global_subset().iter())),

            Self::Chemical(c) => { 
                *self = Self::Precomputed(c.apply(data)?.into_owned());
                self.apply(data) // will give borrowed res
            },
        }
    }

    // pub fn is_coord_dependent(&self) -> bool {
    //     match self {
    //         Self::Not(node) | Self::Same(_, node) => node.is_coord_dependent(),

    //         Self::Comparison(node) => node.is_coord_dependent(),

    //         Self::Or(a, b) | Self::And(a, b) => a.is_coord_dependent() || b.is_coord_dependent(),

    //         Self::Within(_, _) | Self::WithinPoint(_, _) => true,

    //         // All always works for global subset
    //         _ => false,
    //     }
    // }
}

impl ChemicalNode {
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
            "SOL" | "HOH" | "TIP3" | "TIP4" | "TIP5" | "OPC" => true,
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
}

impl Evaluate for ChemicalNode {
    fn is_state_dependent(&self) -> bool {
        false
    }

    fn apply(&mut self, data: &EvaluationContext)
            -> Result<Cow<'_, [usize]>, SelectionParserError> {
        
        let mut res = vec![];
        match self {
            Self::Protein => {
                for (i, at) in data.iter_ind_atom() {
                    if Self::is_protein(at) {
                        res.push(i)
                    }
                }
            }

            Self::Backbone => {
                for (i, at) in data.iter_ind_atom() {
                    if Self::is_backbone(at) {
                        res.push(i)
                    }
                }
            }

            Self::Sidechain => {
                for (i, at) in data.iter_ind_atom() {
                    if Self::is_sidechain(at) {
                        res.push(i)
                    }
                }
            }

            Self::Water => {
                for (i, at) in data.iter_ind_atom() {
                    if Self::is_water(at) {
                        res.push(i)
                    }
                }
            }

            Self::NotWater => {
                for (i, at) in data.iter_ind_atom() {
                    if !Self::is_water(at) {
                        res.push(i)
                    }
                }
            }

            Self::Hydrogen => {
                for (i, at) in data.iter_ind_atom() {
                    if Self::is_hydrogen(at) {
                        res.push(i)
                    }
                }
            }

            Self::NotHydrogen => {
                for (i, at) in data.iter_ind_atom() {
                    if !Self::is_hydrogen(at) {
                        res.push(i)
                    }
                }
            }
        }

        Ok(Cow::from(res))
    }
}

fn get_min_max(state: &State, iter: impl IndexIterator) -> (Vector3f, Vector3f) {
    let mut lower = Vector3f::max_value();
    let mut upper = Vector3f::min_value();
    for i in iter {
        let crd = unsafe { state.get_pos_unchecked(i) };
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
    fn map_str_args(
        &self,
        data: &EvaluationContext,
        args: &[StrKeywordArg],
        f: fn(&Atom) -> &String,
    ) -> Vec<usize> {
        let mut res = vec![];
        
        for (ind, a) in data.iter_ind_atom() {
            for arg in args {
                match arg {
                    StrKeywordArg::Str(s) => {
                        if s == f(a) {
                            res.push(ind);
                            break;
                        }
                    }
                    StrKeywordArg::Regex(r) => {
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

    fn map_int_args(
        &self,
        data: &EvaluationContext,
        args: &Vec<IntKeywordArg>,
        f: fn(&Atom, usize) -> isize,
    ) -> Vec<usize> {
        let mut res = vec![];
        
        for (ind, a) in data.iter_ind_atom() {
            for arg in args {
                match *arg {
                    IntKeywordArg::Int(v) => {
                        if v == f(a, ind) {
                            res.push(ind);
                            break;
                        }
                    }
                    IntKeywordArg::IntRange(b, e) => {
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
}

impl Evaluate for KeywordNode {
    fn is_state_dependent(&self) -> bool {
        false
    }

    fn apply(&mut self, data: &EvaluationContext)
            -> Result<Cow<'_, [usize]>, SelectionParserError> {
        let res = match &*self {
            Self::Name(args) => self.map_str_args(data, args, |a| &a.name),
            Self::Resname(args) => self.map_str_args(data, args, |a| &a.resname),
            Self::Resid(args) => self.map_int_args(data, args, |a, _i| a.resid as isize),
            Self::Resindex(args) => self.map_int_args(data, args, |a, _i| a.resindex as isize),
            Self::Index(args) => self.map_int_args(data, args, |_a, i| i as isize),
            Self::Chain(args) => {
                let mut res = vec![];
                
                for (i, a) in data.iter_ind_atom() {
                    for c in args {
                        if *c == a.chain {
                            res.push(i);
                        }
                    }
                }
                res
            }
        };

        Ok(Cow::from(res))
    }
}

impl MathNode {
    fn is_state_dependent(&self) -> bool {
        match self {
            Self::Float(_) => false,
            Self::X | Self::Y | Self::Z => true,
            Self::Xof(vec) |
            Self::Yof(vec) |
            Self::Zof(vec) => vec.is_state_dependent(),
            Self::Bfactor | Self::Occupancy | Self::Vdw | Self::Mass | Self::Charge => false,
            Self::BinaryOp(a, _, b) => a.is_state_dependent() || b.is_state_dependent(),
            Self::UnaryMinus(v) |
            Self::Function(_, v) => v.is_state_dependent(),
            Self::Dist(_) => true,
        }
    }

    fn eval(&mut self, at: usize, data: &EvaluationContext) -> Result<f32, SelectionParserError> {
        let atom = unsafe { data.topology.get_atom_unchecked(at) };
        match self {
            Self::Float(v) => Ok(*v),
            Self::X => Ok(unsafe { data.state.get_pos_unchecked(at) }[0]),
            Self::Y => Ok(unsafe { data.state.get_pos_unchecked(at) }[1]),
            Self::Z => Ok(unsafe { data.state.get_pos_unchecked(at) }[2]),
            Self::Xof(vec) => Ok(vec.get_vec(data)?.x),
            Self::Yof(vec) => Ok(vec.get_vec(data)?.y),
            Self::Zof(vec) => Ok(vec.get_vec(data)?.z),
            Self::Bfactor => Ok(atom.bfactor),
            Self::Occupancy => Ok(atom.occupancy),
            Self::Vdw => Ok(atom.vdw()),
            Self::Mass => Ok(atom.mass),
            Self::Charge => Ok(atom.charge),
            Self::BinaryOp(a, op, b) => {
                let ret = match op {
                    BinaryOperator::Add => a.eval(at, data)? + b.eval(at, data)?,
                    BinaryOperator::Sub => a.eval(at, data)? - b.eval(at, data)?,
                    BinaryOperator::Mul => a.eval(at, data)? * b.eval(at, data)?,
                    BinaryOperator::Pow => a.eval(at, data)?.powf(b.eval(at, data)?),
                    BinaryOperator::Div => {
                        let b_val = b.eval(at, data)?;
                        if b_val == 0.0 {
                            return Err(SelectionParserError::DivisionByZero);
                        }
                        a.eval(at, data)? / b_val
                    }
                };
                // Precomute if possible
                if !a.is_state_dependent() && !b.is_state_dependent() {
                    *self = Self::Float(ret);
                }
                Ok(ret)
            }

            Self::UnaryMinus(v) => {
                let ret = -v.eval(at, data)?;
                // Precomute if possible
                if !v.is_state_dependent() {
                    *self = Self::Float(ret);
                }
                Ok(ret)
            }

            Self::Function(func, v) => {
                use MathFunctionName as M;
                let val = v.eval(at, data)?;
                let ret = match func {
                    M::Abs => val.abs(),
                    M::Sqrt => {
                        if val < 0.0 {
                            return Err(SelectionParserError::NegativeSqrt);
                        }
                        val.sqrt()
                    }
                    M::Sin => val.sin(),
                    M::Cos => val.cos(),
                };
                // Precomute if possible
                if !v.is_state_dependent() {
                    *self = Self::Float(ret);
                }
                Ok(ret)
            }

            Self::Dist(node) => {
                let mut pos = unsafe { data.state.get_pos_unchecked(at) }.clone();
                // Point should be unwrapped first!
                node.closest_image(&mut pos, data);

                match node {
                    DistanceNode::Point(p, _) => Ok((pos - p.get_vec(data)?).norm()),
                    DistanceNode::Line(p1, p2, _) => {
                        let p1 = p1.get_vec(data)?;
                        let p2 = p2.get_vec(data)?;
                        let v = p2 - p1;
                        let w = pos - p1;
                        Ok((w - v * (w.dot(&v) / v.norm_squared())).norm())
                    }
                    DistanceNode::LineDir(p, dir, _) => {
                        let w = pos - p.get_vec(data)?;
                        let dir = dir.get_unit_vec(data)?.coords;
                        Ok((w - dir * w.dot(&dir)).norm())
                    }
                    DistanceNode::Plane(p1, p2, p3, _) => {
                        let p1 = p1.get_vec(data)?;
                        let p2 = p2.get_vec(data)?;
                        let p3 = p3.get_vec(data)?;
                        // Plane normal
                        let n = (p2 - p1).cross(&(p3 - p1));
                        let w = pos - p1;
                        Ok((n * (w.dot(&n) / n.norm_squared())).norm())
                    }
                    DistanceNode::PlaneNormal(p, n, _) => {
                        let w = pos - p.get_vec(data)?;
                        let n = n.get_unit_vec(data)?.coords;
                        Ok((n * w.dot(&n)).norm())
                    }
                }
            }
        }
    }

    // fn is_coord_dependent(&self) -> bool {
    //     match self {
    //         Self::X | Self::Y | Self::Z | Self::Dist(_) => true,
    //         Self::Add(a, b)
    //         | Self::Sub(a, b)
    //         | Self::Mul(a, b)
    //         | Self::Div(a, b)
    //         | Self::Pow(a, b) => a.is_coord_dependent() || b.is_coord_dependent(),
    //         Self::Neg(v) | Self::Function(_, v) => v.is_coord_dependent(),
    //         _ => false,
    //     }
    // }
}

impl ComparisonNode {
    fn eval_op(
        data: &EvaluationContext,
        v1: &mut MathNode,
        v2: &mut MathNode,
        op: fn(f32, f32) -> bool,
    ) -> Result<Vec<usize>, SelectionParserError> {
        let mut res = vec![];
        
        for p in data.iter_particle() {
            if op(v1.eval(p.id, data)?, v2.eval(p.id, data)?) {
                res.push(p.id);
            }
        }
        Ok(res)
    }

    fn eval_op_chained(
        data: &EvaluationContext,
        v1: &mut MathNode,
        v2: &mut MathNode,
        v3: &mut MathNode,
        op1: fn(f32, f32) -> bool,
        op2: fn(f32, f32) -> bool,
    ) -> Result<Vec<usize>, SelectionParserError> {
        let mut res = vec![];
        
        for p in data.iter_particle() {
            let mid = v2.eval(p.id, data)?;
            if op1(v1.eval(p.id, data)?, mid) && op2(mid, v3.eval(p.id, data)?) {
                res.push(p.id);
            }
        }
        Ok(res)
    }
}

impl Evaluate for ComparisonNode {
    fn is_state_dependent(&self) -> bool {
        match self {
            // Simple
            Self::Eq(v1, v2) |
            Self::Neq(v1, v2) |
            Self::Gt(v1, v2) |
            Self::Geq(v1, v2) |
            Self::Lt(v1, v2)|
            Self::Leq(v1, v2) => v1.is_state_dependent() || v2.is_state_dependent(),
            // Chained left
            Self::LtLt(v1, v2, v3) |
            Self::LtLeq(v1, v2, v3) |
            Self::LeqLt(v1, v2, v3) | 
            Self::LeqLeq(v1, v2, v3) | 
            // Chained right
            Self::GtGt(v1, v2, v3) |
            Self::GtGeq(v1, v2, v3) | 
            Self::GeqGt(v1, v2, v3) | 
            Self::GeqGeq(v1, v2, v3) => v1.is_state_dependent() || v2.is_state_dependent() || v3.is_state_dependent(),
        }
    }

    fn apply(&mut self, data: &EvaluationContext)
            -> Result<Cow<'_, [usize]>, SelectionParserError> {
        let res = match self {
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
        }?;

        Ok(Cow::from(res))
    }

    // fn is_coord_dependent(&self) -> bool {
    //     match self {
    //         // Simple
    //         Self::Eq(v1, v2)
    //         | Self::Neq(v1, v2)
    //         | Self::Gt(v1, v2)
    //         | Self::Geq(v1, v2)
    //         | Self::Lt(v1, v2)
    //         | Self::Leq(v1, v2) => v1.is_coord_dependent() || v2.is_coord_dependent(),
    //         // Chained
    //         Self::LtLt(v1, v2, v3)
    //         | Self::LtLeq(v1, v2, v3)
    //         | Self::LeqLt(v1, v2, v3)
    //         | Self::LeqLeq(v1, v2, v3)
    //         | Self::GtGt(v1, v2, v3)
    //         | Self::GtGeq(v1, v2, v3)
    //         | Self::GeqGt(v1, v2, v3)
    //         | Self::GeqGeq(v1, v2, v3) => {
    //             v1.is_coord_dependent() || v2.is_coord_dependent() || v3.is_coord_dependent()
    //         }
    //     }
    // }
}

#[derive(Error, Debug)]
pub enum SelectionParserError {
    #[error("syntax error: {0}")]
    SyntaxError(String),

    #[error("selection has incompatible topology and state: {0}")]
    DifferentSizes(#[from] TopologyStateSizesError),

    #[error("periodic selection for non-periodic system: {0}")]
    PbcUnwrap(#[from] PeriodicBoxError),

    #[error("division by zero")]
    DivisionByZero,

    #[error("sqrt of negative number")]
    NegativeSqrt,

    #[error("selection measure operation failed: {0}")]
    Measure(#[from] MeasureError),

    #[error("asked for atom {0} while selection inner expression has {1} atoms")]
    OutOfBounds(usize, usize),
}
