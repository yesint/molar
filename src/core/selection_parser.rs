use nalgebra::Unit;
use num_traits::Bounded;
use regex::bytes::Regex;
use sorted_vec::SortedSet;
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
}

//##############################
//#  AST node types
//##############################

#[derive(Debug, PartialEq)]
pub(crate) enum IntKeywordValue {
    Int(i32),
    IntRange(i32, i32),
}

#[derive(Debug)]
pub(crate) enum StrKeywordValue {
    Str(String),
    Regex(Regex),
}

#[derive(Debug, PartialEq)]
pub(crate) enum MathNode {
    Float(f32),
    Function(MathFunctionName, Box<Self>),
    X,
    Y,
    Z,
    Bfactor,
    Occupancy,
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
pub(crate) enum VectorNode {
    Const(Vector3f),
    NthOf(LogicalNode),
    Com(LogicalNode),
    Cog(LogicalNode),
}

#[derive(Debug, PartialEq)]
pub(crate) enum DistanceNode {
    Point(Pos, PbcDims),
    Line(Pos, Pos, PbcDims),
    LineDir(Pos, Unit<Vector3f>, PbcDims),
    Plane(Pos, Pos, Pos, PbcDims),
    PlaneNormal(Pos, Unit<Vector3f>, PbcDims),
}

enum ComparisonOp {
    Eq,
    Neq,
    Leq,
    Lt,
    Geq,
    Gt,
}

#[derive(Debug, PartialEq)]
pub(crate) enum ComparisonNode {
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
pub(crate) enum KeywordNode {
    Name(Vec<StrKeywordValue>),
    Resname(Vec<StrKeywordValue>),
    Chain(Vec<char>),

    Resid(Vec<IntKeywordValue>),
    Resindex(Vec<IntKeywordValue>),
    Index(Vec<IntKeywordValue>),
}

#[derive(Debug, PartialEq)]
pub(crate) enum SameProp {
    Residue,
    Chain,
}

#[derive(Debug, PartialEq)]
pub struct WithinProp {
    cutoff: f32,
    pbc: PbcDims,
    include_inner: bool,
}

#[derive(Debug)]
pub(crate) enum LogicalNode {
    Not(Box<Self>),
    Or(Box<Self>, Box<Self>),
    And(Box<Self>, Box<Self>),
    Keyword(KeywordNode),
    Comparison(ComparisonNode),
    Same(SameProp, Box<Self>),
    Within(WithinProp, Box<Self>),
    All,
}

enum Keyword {
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

#[derive(Debug, PartialEq)]
pub(crate) enum MathFunctionName {
    Abs,
    Sqrt,
    Sin,
    Cos,
}

//##############################
//#  AST application stuff
//##############################

// Intermediate index type for applying AST
//type SubsetType = rustc_hash::FxHashSet<usize>;
type SubsetType = Vec<usize>;

#[derive(Clone)]
pub(crate) struct ApplyData<'a> {
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
    context_subset: Option<&'a SubsetType>,
}

impl<'a> ApplyData<'a> {
    fn new(
        topology: &'a Topology,
        state: &'a State,
        global_subset: &'a SubsetType,
    ) -> Result<Self, SelectionParserError> {
        check_topology_state_sizes(topology, state)?;

        Ok(Self {
            topology,
            state,
            global_subset,
            context_subset: None,
        })
    }

    fn active_subset(&self) -> ActiveSubset {
        ActiveSubset {
            topology: &self.topology,
            state: &self.state,
            subset: self.context_subset.unwrap_or(self.global_subset),
        }
    }

    fn global_subset(&self) -> ActiveSubset {
        ActiveSubset {
            topology: &self.topology,
            state: &self.state,
            subset: self.global_subset,
        }
    }

    fn custom_subset(&'a self, custom: &'a SubsetType) -> ActiveSubset {
        ActiveSubset {
            topology: &self.topology,
            state: &self.state,
            subset: custom,
        }
    }

    fn clone_with_context(&self, context_subset: &'a SubsetType) -> Self {
        Self {
            context_subset: Some(context_subset),
            ..*self
        }
    }

    fn len(&self) -> usize {
        self.global_subset.len()
    }

    // Returns iterator over subset.
    // If context is set, iterates over context
    // otherwise iterates over global subset
    fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        match self.context_subset {
            Some(cont) => cont.iter().cloned(),
            None => self.global_subset.iter().cloned(),
        }
    }

    fn iter_ind_atom_pos(&self) -> impl Iterator<Item = (usize, &Atom, &Pos)> {
        self.iter().map(|i| unsafe {
            (
                i,
                self.topology.nth_atom_unchecked(i),
                self.state.nth_pos_unchecked(i),
            )
        })
    }

    fn iter_ind_atom(&self) -> impl Iterator<Item = (usize, &Atom)> {
        self.iter()
            .map(|i| unsafe { (i, self.topology.nth_atom_unchecked(i)) })
    }

    fn iter_atom_from_index(
        &self,
        index: impl Iterator<Item = usize>,
    ) -> impl Iterator<Item = &Atom> {
        index.map(|i| unsafe { self.topology.nth_atom_unchecked(i) })
    }

    fn iter_ind_pos_from_index<'b>(
        &'b self,
        index: impl ExactSizeIterator<Item = usize> + 'b,
    ) -> impl ExactSizeIterator<Item = (usize, Pos)> + 'b {
        index.map(|i| unsafe { (i, *self.state.nth_pos_unchecked(i)) })
    }

    fn iter_ind_atom_from_index(
        &self,
        index: impl ExactSizeIterator<Item = usize>,
    ) -> impl ExactSizeIterator<Item = (usize, &Atom)> {
        index.map(|i| unsafe { (i, self.topology.nth_atom_unchecked(i)) })
    }
}

// Auxiliary struct representing a current active subset
// Created by ApplyData and contains eithr global or context subset
struct ActiveSubset<'a> {
    topology: &'a Topology,
    state: &'a State,
    subset: &'a SubsetType,
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

impl RandomPos for ActiveSubset<'_> {
    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = *self.subset.get(i).unwrap();
        self.state.nth_pos_unchecked(ind)
    }
}

impl AtomsProvider for ActiveSubset<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.subset
            .iter()
            .map(|i| unsafe { self.topology.nth_atom_unchecked(*i) })
    }
}

impl RandomAtom for ActiveSubset<'_> {
    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = *self.subset.get(i).unwrap();
        self.topology.nth_atom_unchecked(ind)
    }
}

//###################################
//#  AST nodes logic implementation
//###################################

impl LogicalNode {
    /*
    pub fn is_coord_dependent(&self) -> bool {
        match self {
            Self::Not(node) => node.is_coord_dependent(),
            Self::Or(a,b) |
            Self::And(a,b) => {
                a.is_coord_dependent() || b.is_coord_dependent()
            },
            Self::Keyword(node) => node.is_coord_dependent(),
            Self::Comparison(node) =>  node.is_coord_dependent(),
            Self::All => false,
            Self::Same(_,node) => node.is_coord_dependent(),
            Self::Within(_,node) => node.is_coord_dependent(),
        }
    }
    */

    fn map_same_prop<'a, T>(
        &self,
        data: &'a ApplyData,
        inner: &SubsetType,
        prop_fn: fn(&'a Atom) -> &'a T,
    ) -> SubsetType
    where
        T: Eq + std::hash::Hash + Copy,
    {
        // Collect all properties from the inner
        let mut properties = HashSet::<T>::new();
        for at in data.iter_atom_from_index(inner.iter().cloned()) {
            properties.insert(*prop_fn(at));
        }

        let mut res = SubsetType::default();
        // Now loop over *global* subset and add all atoms with the same property
        for (i, at) in data.iter_ind_atom_from_index(data.global_subset.iter().cloned()) {
            for prop in properties.iter() {
                if prop_fn(at) == prop {
                    res.push(i);
                    break;
                }
            }
        }
        res
    }

    pub fn apply(&self, data: &ApplyData) -> Result<SubsetType, SelectionParserError> {
        match self {
            Self::Not(node) => {
                // Here we always use global subset!
                let set1 = rustc_hash::FxHashSet::<usize>::from_iter(data.global_subset.into_iter().cloned());
                let set2 = rustc_hash::FxHashSet::<usize>::from_iter(node.apply(data)?.into_iter());
                Ok(set1.difference(&set2).cloned().collect())
            }

            Self::Or(a, b) => {
                let set1 = rustc_hash::FxHashSet::<usize>::from_iter(a.apply(data)?.into_iter());
                let set2 = rustc_hash::FxHashSet::<usize>::from_iter(b.apply(data)?.into_iter());
                Ok(set1.union(&set2).cloned().collect())
            },

            Self::And(ref a, ref b) => {
                // Check which operand is coord-dependent and put it first
                // This could be optimized to precomputed index and will
                // make re-evaluation of AST faster
                //if !a.is_coord_dependent() && b.is_coord_dependent() {
                //    mem::swap(a,  &mut b);
                //}

                let a_res = a.apply(data)?;
                // Create new instance of data and set a context subset to
                // the result of a
                let b_data = data.clone_with_context(&a_res);

                let set1 = rustc_hash::FxHashSet::<usize>::from_iter(a_res.iter().cloned());
                let set2 = rustc_hash::FxHashSet::<usize>::from_iter(b.apply(&b_data)?.into_iter());
                Ok(set1.intersection(&set2).cloned().collect())
            }

            Self::Keyword(node) => node.apply(data),

            Self::Comparison(node) => node.apply(data),

            Self::Same(prop, node) => {
                let inner = node.apply(data)?;
                let res = match prop {
                    // Here we use the global subset!
                    SameProp::Residue => self.map_same_prop(data, &inner, |at| &at.resindex),
                    SameProp::Chain => self.map_same_prop(data, &inner, |at| &at.chain),
                };
                Ok(res)
            }

            Self::Within(prop, node) => {
                let inner = node.apply(data)?;
                let mut res: SubsetType = Default::default();
                // Perform distance search
                if prop.pbc == PBC_NONE {
                    // Non-periodic variant
                    // Find extents
                    let (mut lower, mut upper) = get_min_max(data.state, inner.iter().cloned());
                    lower.add_scalar_mut(-prop.cutoff - f32::EPSILON);
                    upper.add_scalar_mut(prop.cutoff + f32::EPSILON);

                    let sub1 = data.global_subset();
                    let sub2 = data.custom_subset(&inner);
                    //println!("sub: {} {}",sub1.subset.len(),sub2.subset.len());
                    //println!("inner: {:?}",sub2.subset);                    
                    res = distance_search_within(prop.cutoff, &sub1, &sub2, &lower, &upper);
                    
                    // let searcher = DistanceSearcherDouble::new(
                    //     prop.cutoff,
                    //     // Global subset is used!
                    //     data.iter_ind_pos_from_index(data.global_subset.iter().cloned()),
                    //     data.iter_ind_pos_from_index(inner.iter().cloned()),
                    //     &lower,
                    //     &upper,
                    // );
                    // res = searcher.search();
                } else {
                    // Periodic variant
                    let searcher = DistanceSearcherDouble::new_periodic(
                        prop.cutoff,
                        // Global subset is used!
                        data.iter_ind_pos_from_index(data.global_subset.iter().cloned()),
                        data.iter_ind_pos_from_index(inner.iter().cloned()),
                        data.state.get_box().unwrap(),
                        &prop.pbc,
                    );
                    res = searcher.search();
                };

                // Add inner if asked
                if prop.include_inner {
                    res.extend(inner);
                }
                Ok(res)
            }

            // All always works for global subset
            Self::All => Ok(data.global_subset.iter().cloned().collect()),
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
    pub fn is_coord_dependent(&self) -> bool {
        false
    }

    fn map_str_values(
        &self,
        data: &ApplyData,
        values: &Vec<StrKeywordValue>,
        f: fn(&Atom) -> &String,
    ) -> SubsetType {
        let mut res = SubsetType::default();
        for (ind, a) in data.iter_ind_atom() {
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
        data: &ApplyData,
        values: &Vec<IntKeywordValue>,
        f: fn(&Atom, usize) -> i32,
    ) -> SubsetType {
        let mut res = SubsetType::default();
        for (ind, a) in data.iter_ind_atom() {
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

    fn apply(&self, data: &ApplyData) -> Result<SubsetType, SelectionParserError> {
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
                for (i, a) in data.iter_ind_atom() {
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
    /*
    pub fn is_coord_dependent(&self) -> bool {
        match self {
            Self::Float(_) | Self::Bfactor | Self::Occupancy => false,
            Self::X | Self::Y | Self::Z => true,
            Self::Plus(a, b) |
            Self::Minus(a, b) |
            Self::Mul(a, b) |
            Self::Div(a, b) |
            Self::Pow(a, b) => {
                a.is_coord_dependent() || b.is_coord_dependent()
            },
            Self::Neg(v) => v.is_coord_dependent(),
            Self::Function(_,v) => v.is_coord_dependent(),
            Self::Dist(_) => true,
        }
    }
    */

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
                    if dims[0] || dims[1] || dims[2] {
                        Ok(Some(
                            pbox.ok_or_else(|| SelectionParserError::PbcUnwrap)?
                                .closest_image_dims(point, target, dims),
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
                    DistanceNode::Point(p, _) => Ok((pos - p).norm()),
                    DistanceNode::Line(p1, p2, _) => {
                        let v = p2 - p1;
                        let w = pos - p1;
                        Ok((w - v * (w.dot(&v) / v.norm_squared())).norm())
                    }
                    DistanceNode::LineDir(p, dir, _) => {
                        let w = pos - p;
                        // dir is already normalized
                        Ok((w - dir.into_inner() * w.dot(dir)).norm())
                    }
                    DistanceNode::Plane(p1, p2, p3, _) => {
                        // Plane normal
                        let n = (p2 - p1).cross(&(p3 - p1));
                        let w = pos - p1;
                        Ok((n * (w.dot(&n) / n.norm_squared())).norm())
                    }
                    DistanceNode::PlaneNormal(p, n, _) => {
                        // Plane normal is already normalized
                        let w = pos - p;
                        Ok((n.into_inner() * w.dot(&n)).norm())
                    }
                }
            }
        }
    }
}

impl ComparisonNode {
    /*
    pub fn is_coord_dependent(&self) -> bool {
        match self {
            Self::Eq(v1, v2) |
            Self::Neq(v1, v2) |
            Self::Gt(v1, v2) |
            Self::Geq(v1, v2) |
            Self::Lt(v1, v2) |
            Self::Leq(v1, v2) => {
                v1.is_coord_dependent() || v2.is_coord_dependent()
            },
            // Chained left
            Self::LtLt(v1, v2, v3) |
            Self::LtLeq(v1, v2, v3) |
            Self::LeqLt(v1, v2, v3) |
            Self::LeqLeq(v1, v2, v3) |
            // Chained right
            Self::GtGt(v1, v2, v3) |
            Self::GtGeq(v1, v2, v3) |
            Self::GeqGt(v1, v2, v3) |
            Self::GeqGeq(v1, v2, v3) => {
                v1.is_coord_dependent() || v2.is_coord_dependent() || v3.is_coord_dependent()
            }
        }
    }
    */

    fn eval_op(
        data: &ApplyData,
        v1: &MathNode,
        v2: &MathNode,
        op: fn(f32, f32) -> bool,
    ) -> Result<SubsetType, SelectionParserError> {
        let mut res = SubsetType::default();
        let b = data.state.get_box();

        for (i, atom, pos) in data.iter_ind_atom_pos() {
            let p1 = &v1.closest_image(pos, b)?.unwrap_or(*pos);
            let p2 = &v2.closest_image(pos, b)?.unwrap_or(*pos);

            if op(v1.eval(atom, p1)?, v2.eval(atom, p2)?) {
                res.push(i);
            }
        }
        Ok(res)
    }

    fn eval_op_chained(
        data: &ApplyData,
        v1: &MathNode,
        v2: &MathNode,
        v3: &MathNode,
        op1: fn(f32, f32) -> bool,
        op2: fn(f32, f32) -> bool,
    ) -> Result<SubsetType, SelectionParserError> {
        let mut res = SubsetType::default();
        let b = data.state.get_box();

        for (i, atom, pos) in data.iter_ind_atom_pos() {
            let p1 = &v1.closest_image(pos, b)?.unwrap_or(*pos);
            let p2 = &v2.closest_image(pos, b)?.unwrap_or(*pos);
            let p3 = &v3.closest_image(pos, b)?.unwrap_or(*pos);

            let mid = v2.eval(atom, p2)?;
            if op1(v1.eval(atom, p1)?, mid) && op2(mid, v3.eval(atom, p3)?) {
                res.push(i);
            }
        }
        Ok(res)
    }

    fn apply(&self, data: &ApplyData) -> Result<SubsetType, SelectionParserError> {
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
}

//##############################
//#  Grammar
//##############################

peg::parser! {
    grammar selection_parser() for str {
        // Optional whitespace
        rule _ = (" " / "\t")*
        // Mandatory whitespace
        rule __ = (" " / "\t")+
        // Mandatory whitespace unless followed by paren
        rule ___ = _ &"(" / __

        rule uint() -> u32
            = n:$(['0'..='9']+)
            { n.parse().unwrap() }

        rule int()  -> i32
            = n:$(("-"/"+")? uint())
            { n.parse().unwrap() }

        rule float_val() -> f32
            = n:$((int() ("." uint())? / ("-"/"+") "." uint()) (("e"/"E") int())?)
            { n.parse().unwrap() }

        rule float() -> MathNode
            = v:float_val() { MathNode::Float(v) }

        // Keywords
        rule keyword_name() -> Keyword = "name" {Keyword::Name}
        rule keyword_resid() -> Keyword = "resid" {Keyword::Resid}
        rule keyword_resindex() -> Keyword = "resindex" {Keyword::Resindex}
        rule keyword_resname() -> Keyword = "resname" {Keyword::Resname}
        rule keyword_index() -> Keyword = "index" {Keyword::Index}
        rule keyword_chain() -> Keyword = "chain" {Keyword::Chain}
        rule keyword_residue() -> Keyword = "residue" {Keyword::Residue}
        rule keyword_occupancy() -> Keyword = ("occupancy" / "occ") {Keyword::Occupancy}
        rule keyword_bfactor() -> Keyword = ("bfactor" / "beta") {Keyword::Bfactor}

        rule int_keyword_expr() -> KeywordNode
            = s:(keyword_resid() / keyword_resindex() / keyword_index()) __
              v:((int_range() / int_single()) ++ __)
            {
                match s {
                    Keyword::Resid => KeywordNode::Resid(v),
                    Keyword::Resindex => KeywordNode::Resindex(v),
                    Keyword::Index => KeywordNode::Index(v),
                    _ => unreachable!()
                }
            }

        rule int_range() -> IntKeywordValue
            = i1:int() _ ":" _ i2:int()
            { IntKeywordValue::IntRange(i1,i2) }

        rule int_single() -> IntKeywordValue
            = i:int()
            { IntKeywordValue::Int(i) }

        rule chain_keyword_expr() -> KeywordNode
        = s:keyword_chain() __ v:(['a'..='z' | 'A'..='Z' | '0'..='9'] ++ __)
        {
            KeywordNode::Chain(v)
        }

        rule str_keyword_expr() -> KeywordNode
        = s:(keyword_name() / keyword_resname()) __
            v:((str_value() / regex_value()) ++ __)
        {?
            match s {
                Keyword::Name => Ok(KeywordNode::Name(v)),
                Keyword::Resname => Ok(KeywordNode::Resname(v)),
                _ => Err(unreachable!())
            }
        }

        rule keyword_expr() -> KeywordNode
        = v:(int_keyword_expr() / str_keyword_expr() / chain_keyword_expr()) {v}

        rule regex_value() -> StrKeywordValue
        = "'" s:$((!"'" [_])+) "'"
        {?
            match Regex::new(&format!("^{s}$")) {
                Ok(r) => Ok(StrKeywordValue::Regex(r)),
                Err(_) => Err("Invalid regex value {}"),
            }
        }

        rule str_value() -> StrKeywordValue
        = !("and"/"or") s:$((![' '|'\''|'"'] [_])+)
        { StrKeywordValue::Str(s.to_owned()) }


        // 3-vector value
        rule vec3() -> Pos = vec3_spaces() / vec3_comas()

        rule vec3_spaces() -> Pos
        = x:float_val() __ y:float_val() __ z:float_val() {Pos::new(x, y, z)}

        rule vec3_comas() -> Pos
        = "[" _ x:float_val() _ "," _ y:float_val() _ "," _ z:float_val() _ "]" {
            Pos::new(x, y, z)
        }

        //rule nth_of() -> Pos
        //= logical_expr()

        // Distance
        rule distance() -> DistanceNode
        = distance_point() / distance_line_2points() / distance_line_point_dir()
        / distance_plane_3points() / distance_plane_point_normal()

        rule distance_point() -> DistanceNode
        = "dist" __ b:pbc_expr()? "point" __ p:vec3() {
            DistanceNode::Point(p,b.unwrap_or(PBC_NONE))
        }

        rule distance_line_2points() -> DistanceNode
        = "dist" __ b:pbc_expr()? "line" __ p1:vec3() __ p2:vec3() {
            DistanceNode::Line(p1,p2,b.unwrap_or(PBC_NONE))
        }

        rule distance_line_point_dir() -> DistanceNode
        = "dist" __ b:pbc_expr()? "line" __ p:vec3() __ "dir" __ dir:vec3() {
            DistanceNode::LineDir(p,Unit::new_normalize(dir.coords),b.unwrap_or(PBC_NONE))
        }

        rule distance_plane_3points() -> DistanceNode
        = "dist" __ b:pbc_expr()? "plane" __ p1:vec3() __ p2:vec3() __ p3:vec3() {
            DistanceNode::Plane(p1,p2,p3,b.unwrap_or(PBC_NONE))
        }

        rule distance_plane_point_normal() -> DistanceNode
        = "dist" __ b:pbc_expr()? "plane" __ p:vec3() __ "normal" __ n:vec3() {
            DistanceNode::PlaneNormal(p,Unit::new_normalize(n.coords),b.unwrap_or(PBC_NONE))
        }


        rule abs_function() -> MathFunctionName = "abs" {MathFunctionName::Abs}
        rule sqrt_function() -> MathFunctionName = "sqrt" {MathFunctionName::Sqrt}
        rule sin_function() -> MathFunctionName = "sin" {MathFunctionName::Sin}
        rule cos_function() -> MathFunctionName = "cos" {MathFunctionName::Cos}

        rule math_function_name() -> MathFunctionName
        = abs_function() / sqrt_function() / sin_function() / cos_function()

        // Math
        rule math_expr() -> MathNode
        = precedence!{
            x:(@) _ "+" _ y:@ { MathNode::Plus(Box::new(x),Box::new(y)) }
            x:(@) _ "-" _ y:@ { MathNode::Minus(Box::new(x),Box::new(y)) }
                    "-" _ v:@ { MathNode::Neg(Box::new(v)) }
            --
            x:(@) _ "*" _ y:@ { MathNode::Mul(Box::new(x),Box::new(y)) }
            x:(@) _ "/" _ y:@ { MathNode::Div(Box::new(x),Box::new(y)) }
            --
            x:@ _ "^" _ y:(@) { MathNode::Pow(Box::new(x),Box::new(y)) }
            --
            v:float() {v}
            ['x'|'X'] { MathNode::X }
            ['y'|'Y'] { MathNode::Y }
            ['z'|'Z'] { MathNode::Z }
            d:distance() {MathNode::Dist(d)}
            keyword_occupancy() { MathNode::Occupancy }
            keyword_bfactor() { MathNode::Bfactor }
            f:math_function_name() _ "(" _ e:math_expr() _ ")" {
                MathNode::Function(f,Box::new(e))
            }
            "(" _ e:math_expr() _ ")" { e }
        }

        // Comparisons
        rule comparison_op_eq() -> ComparisonOp = "==" {ComparisonOp::Eq}
        rule comparison_op_neq() -> ComparisonOp = "!=" {ComparisonOp::Neq}
        rule comparison_op_leq() -> ComparisonOp = "<=" {ComparisonOp::Leq}
        rule comparison_op_lt() -> ComparisonOp = "<" {ComparisonOp::Lt}
        rule comparison_op_geq() -> ComparisonOp = ">=" {ComparisonOp::Geq}
        rule comparison_op_gt() -> ComparisonOp = ">" {ComparisonOp::Gt}

        // Simple comparison
        rule comparison_expr() -> ComparisonNode =
            a:math_expr() _
            op:(comparison_op_eq()/comparison_op_neq()/
                comparison_op_leq()/comparison_op_lt()/
                comparison_op_geq()/comparison_op_gt()) _
            b:math_expr() {
                use ComparisonOp as C;
                match op {
                    C::Eq => { ComparisonNode::Eq(a,b) },
                    C::Neq => { ComparisonNode::Neq(a,b) },
                    C::Leq => { ComparisonNode::Leq(a,b) },
                    C::Lt => { ComparisonNode::Lt(a,b) },
                    C::Geq => { ComparisonNode::Geq(a,b) },
                    C::Gt => { ComparisonNode::Gt(a,b) },
                    _ => unreachable!(),
                }
            }

        // Chained comparison
        rule comparison_expr_chained() -> ComparisonNode
        = comparison_expr_chained_l() / comparison_expr_chained_r()

        rule comparison_expr_chained_l() -> ComparisonNode
        =   a:math_expr() _
            op1:(comparison_op_leq()/comparison_op_lt()) _
            b:math_expr() _
            op2:(comparison_op_leq()/comparison_op_lt()) _
            c:math_expr()
        {
            use ComparisonOp as C;
            match (op1,op2) {
                (C::Lt,C::Lt) => { ComparisonNode::LtLt(a,b,c) },
                (C::Lt,C::Leq) => { ComparisonNode::LtLeq(a,b,c) },
                (C::Leq,C::Lt) => { ComparisonNode::LeqLt(a,b,c) },
                (C::Leq,C::Leq) => { ComparisonNode::LeqLeq(a,b,c) },
                _ => unreachable!(),
            }
        }

        rule comparison_expr_chained_r() -> ComparisonNode
        =   a:math_expr() _
            op1:(comparison_op_geq()/comparison_op_gt()) _
            b:math_expr() _
            op2:(comparison_op_geq()/comparison_op_gt()) _
            c:math_expr()
        {
            use ComparisonOp as C;
            match (op1,op2) {
                (C::Gt,C::Gt) => { ComparisonNode::GtGt(a,b,c) },
                (C::Gt,C::Geq) => { ComparisonNode::GtGeq(a,b,c) },
                (C::Geq,C::Gt) => { ComparisonNode::GeqGt(a,b,c) },
                (C::Geq,C::Geq) => { ComparisonNode::GeqGeq(a,b,c) },
                _ => unreachable!(),
            }
        }

        // "Same" expressions
        rule same_expr() -> SameProp
        = "same" __ t:(keyword_residue() / keyword_chain()) __ "as" {
            match t {
                Keyword::Residue => SameProp::Residue,
                Keyword::Chain => SameProp::Chain,
                _ => unreachable!(),
            }
        }

        // Single PBC dimention
        rule pbc_dim() -> bool
        = v:$("1" / "0" / "y" / "n") {
            match v {
                "1" | "y" => true,
                "0" | "n" => false,
                _ => unreachable!()
            }
        }

        // PBC for within
        rule pbc_expr() -> [bool;3]
        = pbc_with_dims() / pbc_full_no_dims() / pbc_none_no_dims()

        rule pbc_full_no_dims() -> [bool;3]
        = "pbc" __
        {PBC_FULL}

        rule pbc_none_no_dims() -> [bool;3]
        = "nopbc" __
        {PBC_NONE}

        rule pbc_with_dims() -> [bool;3]
        = "pbc" __ p:(pbc_dim()*<3>) __ {
            [p[0],p[1],p[2]]
        }

        // Within
        rule within_expr() -> WithinProp
        = "within" __ d:float() __ p:pbc_expr()? s:$(("self" __)?) "of" {
            if let MathNode::Float(cutoff) = d {
                let pbc = match p {
                    Some(dims) => dims,
                    None => PBC_NONE,
                };
                let include_inner = !s.is_empty();
                WithinProp {cutoff, pbc, include_inner}
            } else {
                unreachable!()
            }
        }

        // Logic
        pub rule logical_expr() -> LogicalNode
        = precedence!{
            // Binary
            x:(@) _ "or" _ y:@ { LogicalNode::Or(Box::new(x),Box::new(y)) }
            x:(@) _ "and" _ y:@ { LogicalNode::And(Box::new(x),Box::new(y)) }
            // Unary prefixes
            "not" ___ v:@ { LogicalNode::Not(Box::new(v)) }
            t:same_expr() ___ v:@ { LogicalNode::Same(t,Box::new(v)) }
            p:within_expr() ___ v:@ {LogicalNode::Within(p,Box::new(v))}
            --
            v:keyword_expr() { LogicalNode::Keyword(v) }
            v:comparison_expr() { LogicalNode::Comparison(v) }
            v:comparison_expr_chained() { LogicalNode::Comparison(v) }
            "all" _ { LogicalNode::All }
            "(" _ e:logical_expr() _ ")" { e }
        }
    } // grammar
} // parser

//##############################
//#  Public interface
//##############################

// Alias for top-level rule
pub struct SelectionExpr {
    ast: LogicalNode,
    sel_str: String,
}

impl SelectionExpr {
    pub fn get_str(&self) -> &str {
        &self.sel_str
    }
}

impl TryFrom<&str> for SelectionExpr {
    type Error = SelectionParserError;
    fn try_from(value: &str) -> std::prelude::v1::Result<Self, Self::Error> {
        Ok(Self {
            ast: selection_parser::logical_expr(value).map_err(|e| {
                let s = format!(
                    "\n{}\n{}^\nExpected {}",
                    value,
                    "-".repeat(e.location.column - 1),
                    e.expected
                );
                SelectionParserError::SyntaxError(s)
            })?,
            sel_str: value.to_owned(),
        })
    }
}

impl SelectionExpr {
    pub fn new(s: &str) -> Result<Self, SelectionParserError> {
        Ok(s.try_into()?)
    }

    pub fn apply_whole(
        &self,
        topology: &Topology,
        state: &State,
    ) -> Result<SortedSet<usize>, SelectionParserError> {
        let subset = SubsetType::from_iter(0..topology.num_atoms());
        let data = ApplyData::new(topology, state, &subset)?;
        let index = Vec::<usize>::from_iter(self.ast.apply(&data)?.into_iter());
        Ok(index.into())
    }

    pub fn apply_subset(
        &self,
        topology: &Topology,
        state: &State,
        subset: impl Iterator<Item = usize>,
    ) -> Result<SortedSet<usize>, SelectionParserError> {
        let subset = SubsetType::from_iter(subset);
        let data = ApplyData::new(topology, state, &subset)?;
        let index = self.ast.apply(&data)?.into_iter().collect::<Vec<usize>>();
        Ok(index.into())
    }
}

//##############################
//#  Tests
//##############################

#[cfg(test)]
mod tests {
    use triomphe::UniqueArc;

    use super::{SelectionExpr, State, Topology};
    use crate::io::*;

    #[test]
    fn within_syntax_test() {
        let _ast: SelectionExpr = "within 0.5 pbc yyy of resid 555".try_into().unwrap();
    }

    fn read_test_pdb() -> (UniqueArc<Topology>, UniqueArc<State>) {
        let mut h = FileHandler::open("tests/albumin.pdb").unwrap();
        let structure = h.read_topology().unwrap();
        let state = h.read_state().unwrap().unwrap();
        (structure, state)
    }

    fn read_test_pdb2() -> (UniqueArc<Topology>, UniqueArc<State>) {
        let mut h = FileHandler::open("tests/albumin.pdb").unwrap();
        let structure = h.read_topology().unwrap();
        let state = h.read_state().unwrap().unwrap();
        (structure, state)
    }

    fn get_selection_index(sel_str: &str) -> Vec<usize> {
        let topst = read_test_pdb();
        let ast: SelectionExpr = sel_str.try_into().expect("Error generating AST");
        ast.apply_whole(&topst.0, &topst.1)
            .expect("Error applying AST")
            .to_vec()
    }

    fn get_selection_index2(sel_str: &str) -> Vec<usize> {
        let ast: SelectionExpr = sel_str.try_into().expect("Error generating AST");
        let topst = read_test_pdb2();
        ast.apply_whole(&topst.0, &topst.1)
            .expect("Error applying AST")
            .to_vec()
    }

    #[test]
    #[should_panic]
    fn test_invalid_syntax() {
        let _ast: SelectionExpr = "resname A B C D and resid a:6".try_into().unwrap();
    }

    #[test]
    fn test_sqrt() {
        let topst = read_test_pdb2();

        let ast: SelectionExpr = "sqrt (x^2)<5^2".try_into().expect("Error generating AST");
        let vec1 = ast
            .apply_whole(&topst.0, &topst.1)
            .expect("Error applying AST");

        let ast: SelectionExpr = "x<25".try_into().expect("Error generating AST");
        let vec2 = ast
            .apply_whole(&topst.0, &topst.1)
            .expect("Error applying AST");

        assert_eq!(vec1.len(), vec2.len());
    }

    #[test]
    fn test_dist_syntax() {
        let _ast: SelectionExpr = "dist point 1.9 2.9 3.8 > 0.4"
            .try_into()
            .expect("Error generating AST");
    }

    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/generated_selection_tests.in"
    ));

    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/generated_pteros_tests.in"
    ));
}
