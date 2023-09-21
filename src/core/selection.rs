use anyhow::{bail, Result};
use ascii::{AsciiString,AsciiStr};
use nalgebra::{Point3, Vector3};
use regex::bytes::Regex;

use super::structure::Structure;
use super::state::State;
use super::atom::Atom;
use std::collections::HashSet;


//##############################
//#  AST node types 
//##############################


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
    Eq(MathNode, MathNode),
    Neq(MathNode, MathNode),
    Gt(MathNode, MathNode),
    Geq(MathNode, MathNode),
    Lt(MathNode, MathNode),
    Leq(MathNode, MathNode),
    //Lt3(Box<MathNode>, Box<MathNode>, Box<MathNode>),
    //Gt3(Box<MathNode>, Box<MathNode>, Box<MathNode>),
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
pub enum SameProp {
    Residue,
    Chain,
}

#[derive(Debug)]
pub enum LogicalNode {
    Not(Box<Self>),
    Or(Box<Self>, Box<Self>),
    And(Box<Self>, Box<Self>),
    Keyword(KeywordNode),
    Comparison(ComparisonNode),
    Same(SameProp,Box<Self>),
}


//##############################
//#  AST application stuff
//##############################

// Intermediate index type for applying AST
type SubsetType = HashSet<usize>;

#[derive(Debug, Clone)]
pub struct ApplyData<'a> {
    structure: &'a Structure,
    state: &'a State,
    subset: SubsetType,
      
}

impl<'a> ApplyData<'a> {    
    fn new(structure: &'a Structure, state: &'a State, subset: &SubsetType) -> Result<Self> {
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

//###################################
//#  AST nodes logic implementation
//###################################

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

    fn map_same_prop<'a,T>(&self, data: &'a ApplyData, inner: &SubsetType, prop_fn: fn(&'a Atom)->&'a T) -> SubsetType 
    where T: Eq+std::hash::Hash+Copy {
        // Collect all properties from the inner
        let mut properties = HashSet::<T>::new();
        for el in inner.iter().cloned() {
            properties.insert(*prop_fn(&data.structure.atoms[el]));
        }
        
        // Now loop over current cubset and add all atoms with the same property
        let mut res = SubsetType::new();
        for el in data.subset.iter().cloned() {
            for prop in properties.iter() {
                if prop_fn(&data.structure.atoms[el]) == prop {
                    res.insert(el);
                    break;
                }
            }
        }
        res
    }


    pub fn apply(&self,data: &ApplyData) -> Result<SubsetType> {
        match self {
            Self::Not(node) => {
                Ok(data.subset.difference(&node.apply(data)?).cloned().collect())
            },
            Self::Or(a,b) => {
                Ok(a.apply(data)?.union(&b.apply(data)?).cloned().collect())
            },
            Self::And(a,b) => {
                let a_res = a.apply(data)?;
                let b_data = ApplyData::new(data.structure, data.state, &a_res)?;
                Ok(a_res.intersection(&b.apply(&b_data)?).cloned().collect())
            },
            Self::Keyword(node) => node.apply(data),
            Self::Comparison(node) => node.apply(data),
            Self::Same(prop,node) => {
                let inner = node.apply(data)?;
                let res = match prop {
                    SameProp::Residue => self.map_same_prop(data, &inner, |at| &at.resindex),
                    SameProp::Chain => self.map_same_prop(data, &inner, |at| &at.chain),
                };
                Ok(res)
            },
        }
    }
}


impl KeywordNode {

    fn map_str_values(&self, data: &ApplyData, values: &Vec<StrKeywordValue>, f: fn(&Atom)->&AsciiString) -> SubsetType {
        let mut res = SubsetType::new();
        for ind in data.subset.iter().cloned() {
            let a = &data.structure.atoms[ind];
            for val in values {
                match val {
                    StrKeywordValue::Str(s) => {
                        if s == f(a) {
                            res.insert(ind);
                            break;
                        }
                    }
                    StrKeywordValue::Regex(r) => {
                        if r.is_match(f(a).as_bytes()) {
                            res.insert(ind);
                            break;
                        }
                    }
                }
            }
        }
        res
    }

    fn map_int_values(&self, data: &ApplyData, values: &Vec<IntKeywordValue>, f: fn(&Atom,usize)->i32) -> SubsetType {
        let mut res = SubsetType::new();
        for ind in data.subset.iter().cloned() {
            let a = &data.structure.atoms[ind];
            for val in values {
                match *val {
                    IntKeywordValue::Int(v) => {
                        if v == f(a,ind) {
                            res.insert(ind);
                            break;
                        }
                    }
                    IntKeywordValue::IntRange(b,e) => {
                        let val = f(a,ind);
                        if b<=val && val<=e {
                            res.insert(ind);
                            break;
                        }
                    }
                }
            }
        }
        res
    }

    fn apply(&self,data: &ApplyData) -> Result<SubsetType> {  
        match self {
            Self::Name(values) => {
                Ok(self.map_str_values(data, values, |a| &a.name ))
            },
            Self::Resname(values) => {
                Ok(self.map_str_values(data, values, |a| &a.resname ))
            },
            Self::Resid(values) => {
                Ok(self.map_int_values(data, values, |a,_i| a.resid))
            },
            Self::Resindex(values) => {
                Ok(self.map_int_values(data, values, |a,_i| a.resindex as i32 ))
            },
            Self::Index(values) => {
                Ok(self.map_int_values(data, values, |_a,i| i as i32))
            },
            Self::Chain(values) => {
                let mut res = SubsetType::new();
                for (i,a) in data.structure.atoms.iter().enumerate() {
                    for val in values {
                        match val {                            
                            StrKeywordValue::Str(s) => {
                                if s.chars().next().unwrap() == a.chain { res.insert(i); }
                            },
                            _ => unreachable!()
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
            Self::Float(v) => Ok(*v),
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
    
    fn eval_op(data: &ApplyData, v1: &MathNode, v2: &MathNode, op: fn(f32,f32)->bool) -> Result<SubsetType> {
        let mut res = SubsetType::new();
        for i in data.subset.iter().cloned() {
            if op(v1.eval(data,i)?, v2.eval(data,i)?) { res.insert(i); }
        }
        Ok(res)
    }

    fn apply(&self,data: &ApplyData) -> Result<SubsetType> {
        match self {            
            Self::Eq(v1,v2) => 
                Self::eval_op(data, v1, v2, |a,b| a==b),
            Self::Neq(v1,v2) => 
                Self::eval_op(data, v1, v2, |a,b| a!=b),
            Self::Gt(v1,v2) => 
                Self::eval_op(data, v1, v2, |a,b| a>b),
            Self::Geq(v1,v2) => 
                Self::eval_op(data, v1, v2, |a,b| a>=b),
            Self::Lt(v1,v2) => 
                Self::eval_op(data, v1, v2, |a,b| a<b),
            Self::Leq(v1,v2) => 
                Self::eval_op(data, v1, v2, |a,b| a<=b),
        }
    }
}

//##############################
//#  Grammar
//##############################

peg::parser! {
    grammar selection_parser() for str {

        rule _ = (" " / "\t")* // Optional whitespace
        rule __ = (" " / "\t")+ // Mandatory whitespace

        rule uint() -> u32
            = n:$(['0'..='9']+)
            { n.parse().unwrap() }

        rule int()  -> i32
            = n:$(("-"/"+")? uint())
            { n.parse().unwrap() }

        rule float() -> MathNode
            = n:$((int() ("." uint())? / ("-"/"+") "." uint()) (("e"/"E") int())?)
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
        = a:math_expr() _ op:$("=="/"!="/"<="/"<"/">="/">") _ b:math_expr() {
            match op {
                "==" => { ComparisonNode::Eq(a,b) },
                "!=" => { ComparisonNode::Neq(a,b) },
                "<=" => { ComparisonNode::Leq(a,b) },
                "<" => { ComparisonNode::Lt(a,b) },
                ">=" => { ComparisonNode::Geq(a,b) },
                ">" => { ComparisonNode::Gt(a,b) },
                _ => unreachable!(),
            }
        }

        // By expressions
        pub rule same_expr() -> SameProp
        = "same" __ t:$(("residue" / "chain")) __ "as" {
            match t {
                "residue" => SameProp::Residue,
                "chain" => SameProp::Chain,
                _ => unreachable!(),
            }
        }

        // Logic
        pub rule logical_expr() -> LogicalNode
        = precedence!{
            x:(@) _ "or" _ y:@ { LogicalNode::Or(Box::new(x),Box::new(y)) }
            x:(@) _ "and" _ y:@ { LogicalNode::And(Box::new(x),Box::new(y)) }
            "not" _ v:@ { LogicalNode::Not(Box::new(v)) }
            t:same_expr() _ v:@ { LogicalNode::Same(t,Box::new(v)) }
            //TODO:
            // v:within_expr() {LogicalNode::Within(v)}
            --
            v:keyword_expr() { LogicalNode::Keyword(v) }
            v:comparison_expr() { LogicalNode::Comparison(v) }
            "(" _ e:logical_expr() _ ")" { e }
        }

    } // grammar
} // parser

//##############################
//#  Public interface
//##############################

// Alias for top-level rule
pub type SelectionAst = LogicalNode;

pub fn generate_ast(sel_str: &str) -> Result<SelectionAst> {
    Ok(selection_parser::logical_expr(sel_str)?)
}

pub fn apply_ast_whole(ast: &SelectionAst, structure: &Structure, state: &State) -> Result<Vec<usize>> {
    let data = ApplyData {
        structure,
        state,
        subset: SubsetType::from_iter(0..structure.atoms.len())
    };
    let mut index = Vec::<usize>::from_iter(ast.apply(&data)?.into_iter());
    index.sort();
    Ok( index )
}

pub fn apply_ast_subset(ast: &SelectionAst, structure: &Structure, state: &State, subset: &Vec<usize>) -> Result<Vec<usize>> {
    let data = ApplyData {
        structure,
        state,
        subset: SubsetType::from_iter(subset.iter().cloned())
    };
    Ok( Vec::<usize>::from_iter(ast.apply(&data)?.into_iter()) )
}

//##############################
//#  Tests
//##############################

#[cfg(test)]
mod tests {
    use super::{selection_parser, generate_ast, apply_ast_whole};
    use crate::{io::*, core::Structure,core::State};
    use lazy_static::lazy_static;

    fn read_test_pdb() -> (Structure,State) {
        let mut h = FileHandler::new_reader("triclinic.pdb").unwrap();
        let structure = h.read_structure().unwrap();
        let state = h.read_next_state().unwrap().unwrap();
        (structure,state)
    }

    // Read the test PDB file once and provide the content for tests
    lazy_static! {
        static ref SS: (Structure,State) = read_test_pdb();
    }

    fn get_selection_index(sel_str: &str) -> Vec<usize> {
        let ast = generate_ast(sel_str).expect("Error generating AST");
        apply_ast_whole(&ast, &SS.0, &SS.1).expect("Error applying AST")
    }

    //include!()

    //----------------------------------------------------------------------
    
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
        let res = selection_parser::logical_expr("name C B and not(name A Z or x+1>7)");
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

    #[test]
    pub fn test_apply() {
        use crate::io::FileHandler;
        let mut h = FileHandler::new_reader("colored.pdb").unwrap();
        let structure = h.read_structure().unwrap();
        let state = h.read_next_state().unwrap().unwrap();

        let ast = generate_ast("name N and resid 1:5 and x<20").expect("Error generating AST");
        let mut index = apply_ast_whole(&ast, &structure, &state).expect("Error applying");
        index.sort();

        println!("index: {:?}",index);
    }

 
}
