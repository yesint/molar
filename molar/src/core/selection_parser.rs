use std::cell::RefCell;

use crate::prelude::*;
use sorted_vec::SortedSet;

mod ast;
mod grammar;

mod grammar3;

pub use ast::SelectionParserError;
use ast::{EvaluationContext, LogicalNode};

//##############################
//#  Public interface
//##############################

#[derive(Debug)]
pub struct SelectionExpr {
    ast: RefCell<LogicalNode>,
    sel_str: String,
}

impl SelectionExpr {
    pub fn get_str(&self) -> &str {
        &self.sel_str
    }
}

impl SelectionExpr {
    pub fn new(s: &str) -> Result<Self, SelectionParserError> {
        Ok(Self {
            ast: RefCell::new(grammar::selection_parser::logical_expr(s).map_err(|e| {
                let s = format!(
                    "\n{s}\n{}^\nExpected {}",
                    "-".repeat(e.location.column - 1),
                    e.expected
                );
                SelectionParserError::SyntaxError(s)
            })?),
            sel_str: s.to_owned(),
        })
    }

    pub fn apply_whole(
        &self,
        topology: &Topology,
        state: &State,
    ) -> Result<SortedSet<usize>, SelectionParserError> {
        let subset = (0..topology.num_atoms()).collect::<Vec<_>>();
        let data = EvaluationContext::new(topology, state, &subset)?;
        let mut ast = self.ast.borrow_mut();
        Ok(ast.apply(&data)?.into())
    }

    pub fn apply_subset(
        &self,
        topology: &Topology,
        state: &State,
        subset: &[usize],
    ) -> Result<SortedSet<usize>, SelectionParserError> {
        let data = EvaluationContext::new(topology, state, subset)?;
        let mut ast = self.ast.borrow_mut();
        Ok(ast.apply(&data)?.into())
    }
}

//##############################
//#  Tests
//##############################

#[cfg(test)]
mod tests {
    use super::{SelectionExpr, State, Topology};
    use crate::io::*;
    use std::sync::LazyLock;

    static TOPST: LazyLock<(Topology, State)> = LazyLock::new(|| {
        let mut h = FileHandler::open("tests/albumin.pdb").unwrap();
        h.read().unwrap()
    });

    #[test]
    fn within_syntax_test() {
        let _ast = SelectionExpr::new("within 0.5 pbc yyy of resid 555").unwrap();
    }

    fn get_selection_index(sel_str: &str) -> Vec<usize> {
        let ast = SelectionExpr::new(sel_str).expect("Error generating AST");
        ast.apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST")
            .to_vec()
    }

    fn get_selection_index2(sel_str: &str) -> Vec<usize> {
        let ast = SelectionExpr::new(sel_str).expect("Error generating AST");
        ast.apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST")
            .to_vec()
    }

    #[test]
    #[should_panic]
    fn test_invalid_syntax() {
        let _ast = SelectionExpr::new("resname A B C D and resid a:6").unwrap();
    }

    #[test]
    fn test_sqrt() {
        let ast = SelectionExpr::new("sqrt (x^2)<5^2").expect("Error generating AST");
        let vec1 = ast
            .apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST");

        let ast = SelectionExpr::new("x<25").expect("Error generating AST");
        let vec2 = ast
            .apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST");

        assert_eq!(vec1.len(), vec2.len());
    }

    #[test]
    fn test_dist_syntax() {
        let _ast = SelectionExpr::new("dist point 1.9 2.9 3.8 > 0.4")
            .expect("Error generating AST");
    }

    #[test]
    fn within_from_point() {
        let _ast = SelectionExpr::new("within 0.5 of com pbc 101 of protein")
            .expect("Error generating AST");
    }

    #[test]
    fn debug_print() {
        let ast = SelectionExpr::new("within 0.5 of com pbc 101 of protein")
            .expect("Error generating AST");
        println!("{:?}", ast);
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
