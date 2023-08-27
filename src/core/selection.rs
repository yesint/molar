use std::default;

use anyhow::{bail, Result};
use ascii::AsciiString;
use nalgebra::{Point3, Vector3};
use regex::bytes::Regex;

#[derive(Debug, PartialEq)]
pub enum IntKeywordValue {
    Int(i32),
    IntRange(i32, i32),
}

#[derive(Debug, PartialEq)]
pub enum StrKeywordValue {
    Str(AsciiString),
    Regex(AsciiString),
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
    Math(MathNode),
    Eq(Box<Self>, Box<Self>),
    Neq(Box<Self>, Box<Self>),
    Gt(Box<Self>, Box<Self>),
    Geq(Box<Self>, Box<Self>),
    Lt(Box<Self>, Box<Self>),
    Leq(Box<Self>, Box<Self>),
    Lt3(Box<Self>, Box<Self>, Box<Self>),
    Gt3(Box<Self>, Box<Self>, Box<Self>),
}

#[derive(Debug, PartialEq)]
pub enum KeywordNode {
    Name(Vec<StrKeywordValue>),
    Resname(Vec<StrKeywordValue>),
    Chain(Vec<StrKeywordValue>),

    Resid(Vec<IntKeywordValue>),
    Resindex(Vec<IntKeywordValue>),
    Index(Vec<IntKeywordValue>),
}

#[derive(Debug, PartialEq)]
pub enum LogicalNode {
    Not(Box<Self>),
    Or(Box<Self>, Box<Self>),
    And(Box<Self>, Box<Self>),
    Keyword(KeywordNode),
    Comparison(ComparisonNode),
}

use super::atom::Atom;
use super::state::State;
use super::structure::{self, Structure};
use std::collections::HashSet;

type SubsetType = HashSet<usize>;

struct Selection<'a> {
    structure: &'a Structure,
    state: &'a State,
    subset: Option<SubsetType>,
    cur_iter_index: usize,
}

impl<'a> Selection<'a> {
    fn new(structure: &'a Structure, state: &'a State) -> Result<Self> {
        if structure.atoms.len() != state.coords.len() {
            bail!(
                "There are {} atoms but {} positions",
                structure.atoms.len(),
                state.coords.len()
            );
        }

        Ok(Self {
            structure,
            state,
            subset: None,
            cur_iter_index: 0,
        })
    }

    fn len(&self) -> usize {
        match &self.subset {
            Some(s) => s.len(),
            None => self.structure.atoms.len()
        }
    }
}

struct SelectionAtomHandler<'a> {
    atom: &'a Atom,
    pos: &'a [f32; 3],
}

impl<'a> Iterator for Selection<'a> {
    type Item = SelectionAtomHandler<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur_iter_index == self.len() {
            None
        } else {
            let ind = match &self.subset {
                Some(s) => self.subset.iter().next(),
                None => self.cur_iter_index,
            };
            let ret = Self::Item {
                atom: &self.structure.atoms[ind],
                pos: &self.state.coords[ind],
            };
            
            self.cur_iter_index += 1;
            Some(ret)
        }
    }
}

/*
impl LogicalNode {
    fn apply(&self,data: &AstApplyData) -> Result<HashSet<usize>> {
        match self {
            Self::Not(node) => node.apply(data),
            Self::Or(a,b) => {
                Ok(a.apply(data)?.union(&b.apply(data)?).cloned().collect())
            },
            Self::And(a,b) => {
                Ok(a.apply(data)?.intersection(&b.apply(data)?).cloned().collect())
            },
            Self::Keyword(node) => node.apply(data),
            Self::Comparison(node) => node.apply(data)
        }
    }
}

impl KeywordNode {
    fn apply(&self,data: &AstApplyData) -> Result<HashSet<usize>> {
        let mut res = HashSet::<usize>::new();
        match self {
            Self::Chain(node) => {
                for (i,a) in data.structure.atoms.iter().enumerate() {
                    for val in node {
                        match val {
                            StrKeywordValue::Str(s) => {
                                if s.chars().next().unwrap() == a.chain {
                                    res.insert(i);
                                }
                            }
                            StrKeywordValue::Regex(_) => unreachable!()
                        }

                    }
                }
                Ok(res)
            }
        }
    }
}
*/
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
        { StrKeywordValue::Regex(AsciiString::from_ascii(s).unwrap()) }

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
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b))
            ) }

        rule comparison_neq() -> ComparisonNode
        = a:math_expr() _ "!=" _ b:math_expr()
        { ComparisonNode::Neq(
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b))
            ) }

        rule comparison_gt() -> ComparisonNode
        = a:math_expr() _ ">" _ b:math_expr()
        { ComparisonNode::Gt(
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b))
            ) }

        rule comparison_geq() -> ComparisonNode
        = a:math_expr() _ ">=" _ b:math_expr()
        { ComparisonNode::Geq(
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b))
            ) }

        rule comparison_lt() -> ComparisonNode
        = a:math_expr() _ "<" _ b:math_expr()
        { ComparisonNode::Lt(
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b))
            ) }

        rule comparison_leq() -> ComparisonNode
        = a:math_expr() _ "<=" _ b:math_expr()
        { ComparisonNode::Leq(
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b))
            ) }

        rule comparison_gt3() -> ComparisonNode
        = a:math_expr() _ ">" _ b:math_expr() _ ">" _ c:math_expr()
        { ComparisonNode::Gt3(
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b)),
                Box::new(ComparisonNode::Math(c))
            ) }

        rule comparison_lt3() -> ComparisonNode
        = a:math_expr() _ "<" _ b:math_expr() _ "<" _ c:math_expr()
        { ComparisonNode::Lt3(
                Box::new(ComparisonNode::Math(a)),
                Box::new(ComparisonNode::Math(b)),
                Box::new(ComparisonNode::Math(c))
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
