use std::cell::RefCell;
use crate::{prelude::*, selection::ast::{Evaluate, SelectionParserError, SyntaxError}};

/// Byte range of the token at/around `offset` (where the PEG parser failed), so the
/// whole offending word can be highlighted rather than a single character. If `offset`
/// lands on whitespace or past the end (an "expected more" failure), anchors on the
/// previous non-whitespace byte. Selection strings are single-line ASCII, so byte and
/// character positions coincide.
fn offending_token_span(s: &str, offset: usize) -> std::ops::Range<usize> {
    let b = s.as_bytes();
    let n = b.len();
    if n == 0 {
        return 0..0;
    }
    let is_ws = |c: u8| c == b' ' || c == b'\t';
    let mut i = offset.min(n - 1);
    while i > 0 && is_ws(b[i]) {
        i -= 1;
    }
    if is_ws(b[i]) {
        return offset..offset; // input is all whitespace up to the caret
    }
    let mut start = i;
    while start > 0 && !is_ws(b[start - 1]) {
        start -= 1;
    }
    let mut end = i + 1;
    while end < n && !is_ws(b[end]) {
        end += 1;
    }
    start..end
}

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
                    let offset = e.location.offset;
                    // Curate peg's raw "expected" set: strip the quotes it puts around
                    // literals ("and" -> and), drop empty/whitespace tokens that slip
                    // through, dedupe, and sort for a stable message.
                    let mut expected: Vec<String> = e
                        .expected
                        .tokens()
                        .map(|t| t.trim_matches('"').to_string())
                        .filter(|t| !t.trim().is_empty())
                        .collect();
                    expected.sort();
                    expected.dedup();
                    SelectionParserError::SyntaxError(SyntaxError {
                        input: s.to_owned(),
                        offset,
                        span: offending_token_span(s, offset),
                        expected,
                    })
                })?,
            ),
            sel_str: s.to_owned(),
        })
    }

    /// Applies selection expression to all atoms in the system
    pub(super) fn apply_whole (
        &self,
        sys: &(impl PosProvider + AtomProvider + BoxProvider + VelProvider + ForceProvider)
    ) -> Result<SVec, SelectionParserError> {
        let data = super::ast::EvalContext::new(sys, None)?;
        let mut ast = self.ast.borrow_mut();
        //println!("BEFORE:\n{ast:#?}\n");
        let ind = ast.apply(&data)?.into_owned();
        Ok(SVec::from_unsorted(ind))
    }

    /// Applies the selection expression to a subset of atoms
    pub(super) fn apply_subset (
        &self,
        sys: &(impl PosProvider + AtomProvider + BoxProvider + VelProvider + ForceProvider),
        subset: &[usize],
    ) -> Result<SVec, SelectionParserError> {
        let data = super::ast::EvalContext::new(sys, Some(subset))?;
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
    use super::SelectionExpr;
    use crate::System;
    use std::sync::LazyLock;

    #[test]
    fn syntax_error_has_word_span_and_curated_expected() {
        use crate::prelude::SelectionParserError;
        let Err(SelectionParserError::SyntaxError(e)) = SelectionExpr::new("chain A an resid 5")
        else {
            panic!("expected a syntax error");
        };
        // The span covers the whole offending word, not a single character.
        assert_eq!(&e.input[e.span.clone()], "an");
        // Whitespace tokens are curated out; only meaningful tokens remain.
        assert!(!e.expected.iter().any(|t| t.trim().is_empty()));
        assert!(e.expected.contains(&"and".to_string()));
        assert!(e.expected.contains(&"or".to_string()));
    }

    static SYS: LazyLock<System> = LazyLock::new(|| {
        System::from_file("tests/albumin.pdb").unwrap()
    });

    #[test]
    fn within_syntax_test() {
        let _ast = SelectionExpr::new("within 0.5 pbc yyy of resid 555").unwrap();
    }

    /// A bareword keyword/operator must be a whole word: a following identifier
    /// character means it is *not* that keyword. So "backboneand x<4" and
    /// "name CAand x<4" must NOT parse (the missing space is a syntax error),
    /// while their spaced forms must.
    #[test]
    fn keywords_require_word_boundary() {
        // Glued operator after a compound / a value: rejected.
        assert!(SelectionExpr::new("backboneand x<4").is_err());
        assert!(SelectionExpr::new("name CAand x<4").is_err());
        assert!(SelectionExpr::new("proteinor water").is_err());
        assert!(SelectionExpr::new("notprotein").is_err());

        // Properly spaced forms: accepted.
        assert!(SelectionExpr::new("backbone and x<4").is_ok());
        assert!(SelectionExpr::new("name CA and x<4").is_ok());
        assert!(SelectionExpr::new("protein or water").is_ok());
        assert!(SelectionExpr::new("not protein").is_ok());

        // Plain keywords and multi-value lists still parse.
        assert!(SelectionExpr::new("backbone").is_ok());
        assert!(SelectionExpr::new("name CA CB CG").is_ok());
        assert!(SelectionExpr::new("all").is_ok());

        // A name that merely *starts* with "and"/"or" is a valid value, not an
        // operator (the guard is whole-word) — previously these were rejected.
        assert!(SelectionExpr::new("name android orbital").is_ok());
    }

    fn get_selection_index(sel_str: &str) -> Vec<usize> {
        let ast = SelectionExpr::new(sel_str).expect("Error generating AST");
        ast.apply_whole(&*SYS)
            .expect("Error applying AST")
            .to_vec()
    }

    fn get_selection_index2(sel_str: &str) -> Vec<usize> {
        let ast = SelectionExpr::new(sel_str).expect("Error generating AST");
        ast.apply_whole(&*SYS)
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
            .apply_whole(&*SYS)
            .expect("Error applying AST");

        let ast = SelectionExpr::new("x<25").expect("Error generating AST");
        let _vec2 = ast
            .apply_whole(&*SYS)
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
