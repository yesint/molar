use sorted_vec::SortedSet;
use crate::prelude::*;

mod grammar;
mod ast;

use ast::{EvaluationContext, LogicalNode, SubsetType};
pub use ast::SelectionParserError;

//##############################
//#  Public interface
//##############################

#[derive(Debug)]
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
            ast: grammar::selection_parser::logical_expr(value).map_err(|e| {
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
        &mut self,
        topology: &Topology,
        state: &State,
    ) -> Result<SortedSet<usize>, SelectionParserError> {
        let subset = SubsetType::from_iter(0..topology.num_atoms());
        let data = EvaluationContext::new(topology, state, &subset)?;
        self.ast.precompute(&data)?;
        let index = Vec::<usize>::from_iter(self.ast.apply(&data)?.into_iter());
        Ok(index.into())
    }

    pub fn apply_subset(
        &mut self,
        topology: &Topology,
        state: &State,
        subset: impl Iterator<Item = usize>,
    ) -> Result<SortedSet<usize>, SelectionParserError> {
        let subset = SubsetType::from_iter(subset);
        let data = EvaluationContext::new(topology, state, &subset)?;
        self.ast.precompute(&data)?;
        let index = self.ast.apply(&data)?.into_iter().collect::<Vec<usize>>();
        Ok(index.into())
    }
}

//##############################
//#  Tests
//##############################

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;
    use super::{SelectionExpr, State, Topology};
    use crate::io::*;
    
    static TOPST: LazyLock<(Topology,State)> = LazyLock::new(|| {
        let mut h = FileHandler::open("tests/albumin.pdb").unwrap();
        h.read().unwrap()
    });
    

    #[test]
    fn within_syntax_test() {
        let _ast: SelectionExpr = "within 0.5 pbc yyy of resid 555".try_into().unwrap();
    }

    fn get_selection_index(sel_str: &str) -> Vec<usize> {
        let mut ast: SelectionExpr = sel_str.try_into().expect("Error generating AST");
        ast.apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST")
            .to_vec()
    }

    fn get_selection_index2(sel_str: &str) -> Vec<usize> {
        let mut ast: SelectionExpr = sel_str.try_into().expect("Error generating AST");
        ast.apply_whole(&TOPST.0, &TOPST.1)
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
        let mut ast: SelectionExpr = "sqrt (x^2)<5^2".try_into().expect("Error generating AST");
        let vec1 = ast
            .apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST");

        let mut ast: SelectionExpr = "x<25".try_into().expect("Error generating AST");
        let vec2 = ast
            .apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST");

        assert_eq!(vec1.len(), vec2.len());
    }

    #[test]
    fn test_dist_syntax() {
        let _ast: SelectionExpr = "dist point 1.9 2.9 3.8 > 0.4"
            .try_into()
            .expect("Error generating AST");
    }

    #[test]
    fn within_from_point() {
        let _ast: SelectionExpr = "within 0.5 of com pbc 101 of protein"
            .try_into()
            .expect("Error generating AST");
    }

    #[test]
    fn debug_print() {
        let ast: SelectionExpr = "within 0.5 of com pbc 101 of protein"
            .try_into()
            .expect("Error generating AST");
        println!("{:?}",ast);
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
