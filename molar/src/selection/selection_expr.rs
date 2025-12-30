use std::cell::RefCell;
use crate::{prelude::*, selection::ast::{Evaluate, SelectionParserError}};

//mod grammar3;

//##############################
//#  Public interface
//##############################

/// A compiled selection expression that can be used to select atoms based on various criteria
///
/// Selection expressions are parsed from strings and compiled into an abstract syntax tree (AST)
/// that can be evaluated against a molecular system to select atoms matching the criteria.
#[derive(Debug)]
pub struct SelectionExpr {
    // The internals of the parser AST are changed when parsed
    // but we don't want to expose it, so internal mutability. 
    ast: RefCell<super::ast::LogicalNode>,
    sel_str: String,
}

impl SelectionExpr {
    /// Returns the original selection expression string
    pub fn get_str(&self) -> &str {
        &self.sel_str
    }
}

impl SelectionExpr {
    /// Creates a new selection expression from a string
    ///
    /// # Arguments
    /// * `s` - The selection expression string to parse
    ///
    /// # Returns
    /// * `Ok(SelectionExpr)` if parsing succeeds
    /// * `Err(SelectionParserError)` if there's a syntax error in the expression
    ///
    /// # Examples
    /// ```
    /// # use molar::prelude::SelectionExpr;
    /// let expr = SelectionExpr::new("resname ALA").unwrap();
    /// ```
    pub fn new(s: impl AsRef<str>) -> Result<Self, SelectionParserError> {
        let s = s.as_ref().trim();
        Ok(Self {
            ast: RefCell::new(
                super::grammar::selection_parser::logical_expr(s).map_err(|e| {
                    let err_str = format!(
                        "\n{s}\n{}^\nExpected {}",
                        "-".repeat(e.location.column - 1),
                        e.expected
                    );
                    SelectionParserError::SyntaxError(err_str)
                })?,
            ),
            sel_str: s.to_owned(),
        })
    }

    /// Applies the selection expression to all atoms in the system
    pub(crate) fn apply_whole(
        &self,
        topology: &Topology,
        state: &State,
    ) -> Result<SVec, SelectionParserError> {
        let data = super::ast::EvaluationContext::new(topology, state, None)?;
        let mut ast = self.ast.borrow_mut();
        //println!("BEFORE:\n{ast:#?}\n");
        let ind = ast.apply(&data)?.into_owned();
        Ok(SVec::from_unsorted(ind))
    }

    /// Applies the selection expression to a subset of atoms
    pub(crate) fn apply_subset(
        &self,
        topology: &Topology,
        state: &State,
        subset: &[usize],
    ) -> Result<SVec, SelectionParserError> {
        let data = super::ast::EvaluationContext::new(topology, state, Some(subset))?;
        let mut ast = self.ast.borrow_mut();
        let ind = ast.apply(&data)?.into_owned();
        Ok(SVec::from_unsorted(ind))
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
        let ast = SelectionExpr::new("within 0.3 of x>2+2*3").expect("Error generating AST");
        let _vec1 = ast
            .apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST");

        let ast = SelectionExpr::new("x<25").expect("Error generating AST");
        let _vec2 = ast
            .apply_whole(&TOPST.0, &TOPST.1)
            .expect("Error applying AST");

        //assert_eq!(vec1.len(), vec2.len());
    }

    #[test]
    fn test_dist_syntax() {
        let _ast =
            SelectionExpr::new("dist point 1.9 2.9 3.8 > 0.4").expect("Error generating AST");
    }

    #[test]
    fn test_leading_space() {
        let _ast =
            SelectionExpr::new("  dist point 1.9 2.9 3.8 > 0.4").expect("Error generating AST");
    }

    #[test]
    fn test_trainling_space() {
        let _ast =
            SelectionExpr::new("dist point 1.9 2.9 3.8 > 0.4  ").expect("Error generating AST");
    }

    #[test]
    fn within_from_point() {
        let _ast = SelectionExpr::new("within 0.5 of com pbc 101 of protein")
            .expect("Error generating AST");
    }

    #[test]
    fn test_x_of() {
        let _ast = SelectionExpr::new("x < x of com of name CA")
            .expect("Error generating AST");
        // x < name(CA).com.x
        // pbc=yyy; resid(13 14 15).within(0.35).cog.dist < 3.0
    }

    #[test]
    fn test_nth_pos_of() {
        let _ast = SelectionExpr::new("x < x of pos 3 of name CA")
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
        "/tests/generated_vmd_tests.in"
    ));

    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/generated_pteros_tests.in"
    ));
}
