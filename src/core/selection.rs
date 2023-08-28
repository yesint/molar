use std::default;
use std::path::Iter;

use anyhow::{bail, Result};
use ascii::{AsciiString,AsciiStr};
use nalgebra::{Point3, Vector3};
use regex::bytes::Regex;

#[derive(Debug, PartialEq)]
pub enum IntKeywordValue {
    Int(i32),
    IntRange(i32, i32),
}

#[derive(Debug)]
pub enum StrKeywordValue {
    Str(AsciiString),
    Regex(Regex),
}

#[derive(Debug, PartialEq)]
pub enum MathNode {
    Float(f32),
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
}

#[derive(Debug, PartialEq)]
pub enum ComparisonNode {    
    Eq(Box<MathNode>, Box<MathNode>),
    Neq(Box<MathNode>, Box<MathNode>),
    Gt(Box<MathNode>, Box<MathNode>),
    Geq(Box<MathNode>, Box<MathNode>),
    Lt(Box<MathNode>, Box<MathNode>),
    Leq(Box<MathNode>, Box<MathNode>),
    Lt3(Box<MathNode>, Box<MathNode>, Box<MathNode>),
    Gt3(Box<MathNode>, Box<MathNode>, Box<MathNode>),
}

#[derive(Debug)]
pub enum KeywordNode {
    Name(Vec<StrKeywordValue>),
    Resname(Vec<StrKeywordValue>),
    Chain(Vec<StrKeywordValue>),

    Resid(Vec<IntKeywordValue>),
    Resindex(Vec<IntKeywordValue>),
    Index(Vec<IntKeywordValue>),
}

#[derive(Debug)]
pub enum LogicalNode {
    Not(Box<Self>),
    Or(Box<Self>, Box<Self>),
    And(Box<Self>, Box<Self>),
    Keyword(KeywordNode),
    Comparison(ComparisonNode),
}

use super::structure::{Structure, self};
use super::state::State;
use super::atom::Atom;
use std::collections::HashSet;

type SubsetType = HashSet<usize>;

#[derive(Debug, Clone)]
struct ApplyData<'a> {
    structure: &'a Structure,
    state: &'a State,
    subset: SubsetType,    
}

impl<'a> ApplyData<'a> {
    fn new(structure: &'a Structure, state: &'a State) -> Result<Self> {
        if structure.atoms.len() != state.coords.len() {
            bail!("There are {} atoms but {} positions",structure.atoms.len(),state.coords.len());
        }

        Ok(Self {
            structure,
            state, 
            subset: SubsetType::from_iter(0..structure.atoms.len()),
        })
    }

    fn new_with_subset(structure: &'a Structure, state: &'a State, subset: &SubsetType) -> Result<Self> {
        if structure.atoms.len() != state.coords.len() {
            bail!("There are {} atoms but {} positions",structure.atoms.len(),state.coords.len());
        }

        Ok(Self {
            structure,
            state, 
            subset: subset.clone(),
        })
    }


    fn len(&self) -> usize {
        self.subset.len()
    }
}

/*
struct ApplyDataIterator<'a> {
    atom: &'a Atom,
    pos: &'a [f32;3],
}

impl<'a> Iterator for ApplyData<'a> {
    type Item = ApplyDataIterator<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        
            let ind = match self.subset {
                Some(s) => s.iter().next().unwrap(),
                None => &self.iter_pos,
            }

            ApplyDataIterator {
                atom: &self.structure.atoms[ind],
                pos: &self.state.coords[ind],
            }
        
    }
}
*/

impl LogicalNode {
    /*
    fn is_coord_dependent(&self) -> bool {
        match self {
            Self::Not(node) => node.is_coord_dependent(),
            Self::Or(a,b) => a.is_coord_dependent() || b.is_coord_dependent(),
            Self::And(a,b) => a.is_coord_dependent() || b.is_coord_dependent(),
            Self::Keyword(node) => node.is_coord_dependent(),
            Self::Comparison(node) =>  node.is_coord_dependent()
        }
    }
    */

    fn apply(&self,data: &ApplyData) -> Result<SubsetType> {
        match self {
            Self::Not(node) => node.apply(data),
            Self::Or(a,b) => {
                Ok(a.apply(data)?.union(&b.apply(data)?).cloned().collect())
            },
            Self::And(a,b) => {
                let a_res = a.apply(data)?;
                let b_data = ApplyData::new_with_subset(data.structure, data.state, &a_res)?;
                Ok(a_res.intersection(&b.apply(&b_data)?).cloned().collect())
            },
            Self::Keyword(node) => node.apply(data),
            Self::Comparison(node) => node.apply(data),
        }
    }
}


impl KeywordNode {

    fn map_str_values(&self, data: &ApplyData, values: &Vec<StrKeywordValue>, f: fn(&Atom)->&AsciiString) -> SubsetType {
        let mut res = SubsetType::new();
        for (i,a) in data.structure.atoms.iter().enumerate() {
            for val in values {
                match val {
                    StrKeywordValue::Str(s) => {
                        if s == f(a) { res.insert(i); }
                    }
                    StrKeywordValue::Regex(r) => {
                        if r.is_match(f(a).as_bytes()) { res.insert(i); }
                    }
                }
            }
        }
        res
    }

    fn map_int_values(&self, data: &ApplyData, values: &Vec<IntKeywordValue>, f: fn(&Atom,usize)->i32) -> SubsetType {
        let mut res = SubsetType::new();
        for (i,a) in data.structure.atoms.iter().enumerate() {
            for val in values {
                match val {
                    IntKeywordValue::Int(v) => {
                        if *v == f(a,i) { res.insert(i); }
                    }
                    IntKeywordValue::IntRange(b,e) => {
                        if *b<f(a,i) && f(a,i)<*e { res.insert(i); }
                    }
                }
            }
        }
        res
    }

    fn apply(&self,data: &ApplyData) -> Result<SubsetType> {
        
        match self {
            Self::Name(values) => {
                Ok(self.map_str_values(data, values, |a:&Atom| &a.name ))
            },
            Self::Resname(values) => {
                Ok(self.map_str_values(data, values, |a:&Atom| &a.resname ))
            },
            Self::Resid(values) => {
                Ok(self.map_int_values(data, values, |a:&Atom,_i:usize| a.resid as i32))
            },
            Self::Resindex(values) => {
                Ok(self.map_int_values(data, values, |a:&Atom,_i:usize| a.resindex as i32 ))
            },
            Self::Index(values) => {
                Ok(self.map_int_values(data, values, |_a:&Atom,i:usize| i as i32))
            },
            Self::Chain(values) => {
                let mut res = SubsetType::new();
                // Sanity check
                //for val in values {
                //    if val.len() !=1 { bail!("Chain has to be one character, given {val}") }
                //}
                for (i,a) in data.structure.atoms.iter().enumerate() {
                    for val in values {
                        match val {                            
                            StrKeywordValue::Str(s) => {
                                if s.chars().next().unwrap() == a.chain { res.insert(i); }
                            }
                            StrKeywordValue::Regex(r) => {
                                bail!("Can't use regex with chains, given {r}")
                            }
                        }
                    }
                }
                Ok(res)
            }
        }
    }
}


impl MathNode {
    fn eval(&self,data: &ApplyData, i: usize) -> Result<f32> {
        match self {
            Self::Float(v) => Ok(v),
            Self::X => Ok(data.state.coords[i][0]),
            Self::Y => Ok(data.state.coords[i][1]),
            Self::Z => Ok(data.state.coords[i][2]),
            Self::Bfactor => Ok(data.structure.atoms[i].bfactor),
            Self::Occupancy => Ok(data.structure.atoms[i].occupancy),
            Self::Plus(a,b) => Ok(a.eval(data,i)?+b.eval(data,i)?),
            Self::Minus(a,b) => Ok(a.eval(data,i)?-b.eval(data,i)?),
            Self::Mul(a,b) => Ok(a.eval(data,i)?*b.eval(data,i)?),
            Self::Div(a,b) => {
                let b_val = b.eval(data,i)?;
                if b_val==0.0 {bail!("Division by zero at atom {i}")}
                Ok(a.eval(data,i)?/b_val)
            },
            Self::Pow(a,b) => Ok(a.eval(data,i)?.powf(b.eval(data,i)?)),
            Self::Neg(v) => Ok(-v.eval(data,i)?),            
        }
    }
}


impl ComparisonNode {
    
    fn apply(&self,data: &ApplyData) -> Result<SubsetType> {
        let mut res = SubsetType::new();
        match self {            
            Self::Eq(a,b)=>{
                for (i,at) in data.structure.atoms.iter().enumerate() {                    
                    if a.eval(data,i)? == b.eval(data,i)? { res.insert(i); }
                }
                Ok(res)
            }
        }
    }
}
//fn apply_ast(ast: &LogicalNode) -> Result<Vec<u32>> {
    //use array_tool::vec::{Intersect,Union};
    
//}


peg::parser! {
    grammar selection_parser() for str {

        rule _ = (" " / "\t")* // Optional whitespace
        rule __ = (" " / "\t")+ //Mandatory whitespace

        rule uint() -> u32
            = n:$(['0'..='9']+)
            { n.parse().unwrap() }

        rule int()  -> i32
            = n:$(("-"/"+")? uint())
            { n.parse().unwrap() }

        rule float() -> MathNode
            = n:$((int() ("." int())? / ("-"/"+") "." int()) (("e"/"E") int())?)
            { MathNode::Float(n.parse().unwrap()) }

        pub rule int_keyword_expr() -> KeywordNode
            = s:$("resid" / "resindex" / "index") __
              v:((int_range() / int_single()) ++ __)
            {
                match s {
                    "resid" => KeywordNode::Resid(v),
                    "resindex" => KeywordNode::Resindex(v),
                    "index" => KeywordNode::Index(v),
                    _ => unreachable!()
                }
            }

        rule int_range() -> IntKeywordValue
            = i1:int() ":" i2:int()
            { IntKeywordValue::IntRange(i1,i2) }

        rule int_single() -> IntKeywordValue
            = i:int()
            { IntKeywordValue::Int(i) }

        pub rule str_keyword_expr() -> KeywordNode
        = s:$("name" / "resname" / "chain") __
            v:((str_value() / regex_value()) ++ __)
        {?
            match s {
                "name" => Ok(KeywordNode::Name(v)),
                "resname" => Ok(KeywordNode::Resname(v)),
                "chain" => {
                    for el in &v[..] {
                        match el {
                            StrKeywordValue::Str(s) => {
                                if s.len() != 1 {
                                    return Err("Chain has to be a single char")
                                }
                            }
                            StrKeywordValue::Regex(_) => {
                                return Err("Chain can't match a regex")
                            }
                        }
                    }
                    Ok(KeywordNode::Chain(v))
                },
                _ => Err(unreachable!())
            }
        }

        rule keyword_expr() -> KeywordNode
        = v:(int_keyword_expr() / str_keyword_expr()) {v}

        rule regex_value() -> StrKeywordValue
        = "'" s:$((!"'" [_])+) "'"
        {?
            match Regex::new(&format!("^{s}$")) {
                Ok(r) => Ok(StrKeywordValue::Regex(r)),
                Err(_) => Err("Invalid regex value"),
            }            
        }

        rule str_value() -> StrKeywordValue
        = !("and"/"or") s:$((![' '|'\''|'"'] [_])+)
        { StrKeywordValue::Str(AsciiString::from_ascii(s).unwrap()) }

        // Math
        pub rule math_expr() -> MathNode
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
            ("occupancy" / "occ") { MathNode::Occupancy }
            ("bfactor" / "beta") { MathNode::Bfactor }
            // TODO:
            // v:distance_expr() {v}
            "(" _ e:math_expr() _ ")" { e }
          }


        // Comparisons
        pub rule comparison_expr() -> ComparisonNode
        = v:(
            //comparison_gt3() / comparison_lt3() /
            comparison_eq() / comparison_neq() /
            comparison_gt() / comparison_geq() /
            comparison_lt() / comparison_leq() 
            ) {v}

        rule comparison_eq() -> ComparisonNode
        = a:math_expr() _ "==" _ b:math_expr()
          { ComparisonNode::Eq(
                Box::new(a),
                Box::new(b)
            ) }
        
        rule comparison_neq() -> ComparisonNode
        = a:math_expr() _ "!=" _ b:math_expr()
        { ComparisonNode::Neq(
                Box::new(a),
                Box::new(b)
            ) }
        
        rule comparison_gt() -> ComparisonNode
        = a:math_expr() _ ">" _ b:math_expr()
        { ComparisonNode::Gt(
                Box::new(a),
                Box::new(b)
            ) }
        
        rule comparison_geq() -> ComparisonNode
        = a:math_expr() _ ">=" _ b:math_expr()
        { ComparisonNode::Geq(
                Box::new(a),
                Box::new(b)
            ) }

        rule comparison_lt() -> ComparisonNode
        = a:math_expr() _ "<" _ b:math_expr()
        { ComparisonNode::Lt(
                Box::new(a),
                Box::new(b)
            ) }
        
        rule comparison_leq() -> ComparisonNode
        = a:math_expr() _ "<=" _ b:math_expr()
        { ComparisonNode::Leq(
                Box::new(a),
                Box::new(b)
            ) }

        rule comparison_gt3() -> ComparisonNode
        = a:math_expr() _ ">" _ b:math_expr() _ ">" _ c:math_expr()
        { ComparisonNode::Gt3(
                Box::new(a),
                Box::new(b),
                Box::new(c)
            ) }
        
        rule comparison_lt3() -> ComparisonNode
        = a:math_expr() _ "<" _ b:math_expr() _ "<" _ c:math_expr()
        { ComparisonNode::Lt3(
                Box::new(a),
                Box::new(b),
                Box::new(c)
            ) }
    
        // Logic

        pub rule logical_expr() -> LogicalNode
        = precedence!{
            x:(@) _ "or" _ y:@ { LogicalNode::Or(Box::new(x),Box::new(y)) }
            x:(@) _ "and" _ y:@ { LogicalNode::And(Box::new(x),Box::new(y)) }
            "not" _ v:@ { LogicalNode::Not(Box::new(v)) }
            //TODO:
            // v:by_expr() {LogicalNode::By(v)}
            // v:within_expr() {LogicalNode::Within(v)}
            --
            v:keyword_expr() { LogicalNode::Keyword(v) }
            v:comparison_expr() { LogicalNode::Comparison(v) }
            "(" _ e:logical_expr() _ ")" { e }
        }

    } // grammar
} // parser

#[cfg(test)]
mod tests {
    use super::selection_parser;

    #[test]
    pub fn test_int_keyword_expr() {
        let res = selection_parser::int_keyword_expr("index 1  2 3:4 5");
        println!("{:?}", res);
    }

    #[test]
    pub fn test_str_keyword_expr() {
        let res = selection_parser::str_keyword_expr("chain Cf A");
        println!("{:?}", res);
    }

    #[test]
    pub fn test_logical() {
        let res = selection_parser::logical_expr("name C B and not(name A Z or x>7)");
        println!("{:?}", res);
    }

    #[test]
    pub fn test_math() {
        let res = selection_parser::math_expr("- 1  + - ( -2  *5 )");
        println!("{:?}", res);
    }

    #[test]
    pub fn test_comparison() {
        let res = selection_parser::comparison_expr("(x +4 )< x");
        println!("{:#?}", res);

        
    }
}
