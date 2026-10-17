use molar::prelude::*;
use std::cell::Cell;

fn system(xs: &[Float], names: &[&str]) -> System {
    let mut top = Topology::default();
    for name in names {
        top.atoms.push(&Atom::new().with_name(name).with_mass(1.0));
    }
    let mut state = State::new_fake(xs.len());
    for (p, &x) in state.coords.iter_mut().zip(xs) {
        *p = Pos::new(x, 0.0, 0.0);
    }
    System::new(top, state).unwrap()
}

struct Resolver {
    id: SelectionResolverId,
    names: Vec<String>,
    source: System,
    second: System,
    domain: Vec<usize>,
    binds: Cell<usize>,
    missing: bool,
}
impl Resolver {
    fn new() -> Self {
        Self {
            id: SelectionResolverId::new(),
            names: vec![],
            source: system(&[5.0, 0.1], &["CB", "CA"]),
            second: system(&[2.1], &["CA"]),
            domain: vec![1],
            binds: Cell::new(0),
            missing: false,
        }
    }
}
impl SelectionResolver for Resolver {
    fn parse_atom(&self, input: &str, position: usize) -> Option<ForeignParse> {
        let rest = &input[position..];
        if !rest.starts_with('m') && !rest.starts_with("deleted") { return None; }
        let len = rest.char_indices().find(|(_, c)| c.is_whitespace() || matches!(c, ':' | '(' | ')' | '&' | '|'))
            .map_or(rest.len(), |(i, _)| i);
        Some(ForeignParse::matched(position + len, &rest[..len]))
    }
    fn identity(&self) -> SelectionResolverId {
        self.id
    }
    fn prepare(&mut self, variable: &str) -> Result<SelectionReferenceId, SelectionResolveError> {
        // This micro-grammar belongs to the caller, not to molar.
        if !["m1", "m1.r2[3]", "m2", "m1[0]", "m1.r_empty"].contains(&variable) {
            return Err(SelectionResolveError::new("unknown variable"));
        }
        self.names.push(variable.into());
        Ok(self.id.reference(self.names.len() - 1))
    }
    fn bind(
        &self,
        reference: SelectionReferenceId,
    ) -> Result<SelectionSourceView<'_>, SelectionResolveError> {
        self.binds.set(self.binds.get() + 1);
        if self.missing {
            return Err(SelectionResolveError::new("source was deleted"));
        }
        match self.names[reference.slot()].as_str() {
            "m1.r2[3]" | "m1[0]" => SelectionSourceView::from_parts(
                self.source.topology(),
                self.source.state(),
                Some(&self.domain),
            )
            .map_err(|e| SelectionResolveError::new(e.to_string())),
            "m1.r_empty" => SelectionSourceView::from_parts(
                self.source.topology(),
                self.source.state(),
                Some(&[]),
            )
            .map_err(|e| SelectionResolveError::new(e.to_string())),
            "m2" => Ok(self.second.selection_view()),
            _ => Ok(self.source.selection_view()),
        }
    }
}
fn indices(sel: Sel) -> Vec<usize> {
    sel.into_svec().to_vec()
}
fn default_system() -> System {
    system(&[0.0, 2.0, 5.0], &["CA", "CB", "CA"])
}

#[test]
fn foreign_indices_never_become_default_indices() {
    let default = default_system();
    let mut resolver = Resolver::new();
    let expr =
        SelectionExpr::new_with_resolver("within 0.5 of m1:(name CA)", &mut resolver).unwrap();
    assert_eq!(
        indices(default.select_with(&expr, &resolver).unwrap()),
        vec![0]
    );
}
#[test]
fn opaque_names_and_selected_domain() {
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver("within 0.5 of m1.r2[3]", &mut resolver).unwrap();
    assert_eq!(resolver.names, ["m1.r2[3]"]);
    assert_eq!(
        indices(default_system().select_with(&expr, &resolver).unwrap()),
        vec![0]
    );
    let expr = SelectionExpr::new_with_resolver("within 0.5 of m1[0]:(not name CB)", &mut resolver)
        .unwrap();
    assert_eq!(
        indices(default_system().select_with(&expr, &resolver).unwrap()),
        vec![0]
    );
}
#[test]
fn prepare_once_and_bind_once_per_evaluation() {
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver(
        "within 0.5 of (m1:(name CA) or m1:(name CB))",
        &mut resolver,
    )
    .unwrap();
    assert_eq!(resolver.names, ["m1"]);
    assert_eq!(resolver.binds.get(), 0);
    let default = default_system();
    assert_eq!(
        indices(default.select_with(&expr, &resolver).unwrap()),
        vec![0, 2]
    );
    assert_eq!(resolver.binds.get(), 1);
    default.select_with(&expr, &resolver).unwrap();
    assert_eq!(resolver.binds.get(), 2);
}
#[test]
fn mixed_source_reference_union() {
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver("within 0.5 of (m1.r2[3] or m2)", &mut resolver)
        .unwrap();
    assert_eq!(
        indices(default_system().select_with(&expr, &resolver).unwrap()),
        vec![0, 1]
    );
}
#[test]
fn rebinding_uses_current_coordinates_and_domain() {
    let default = default_system();
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver("within 0.5 of m1.r2[3]", &mut resolver).unwrap();
    assert_eq!(
        indices(default.select_with(&expr, &resolver).unwrap()),
        vec![0]
    );
    resolver.domain = vec![0];
    assert_eq!(
        indices(default.select_with(&expr, &resolver).unwrap()),
        vec![2]
    );
    let mut frame = resolver.source.state().clone();
    frame.coords[0] = Pos::new(2.1, 0.0, 0.0);
    resolver.source.set_state(frame).unwrap();
    assert_eq!(
        indices(default.select_with(&expr, &resolver).unwrap()),
        vec![1]
    );
}
#[test]
fn deletion_and_preparation_errors_have_variable_spans() {
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver("within 0.5 of m1", &mut resolver).unwrap();
    resolver.missing = true;
    match default_system().select_with(&expr, &resolver) {
        Err(SelectionError::Parser(SelectionParserError::Variable { name, span, .. })) => {
            assert_eq!(name, "m1");
            assert_eq!(&expr.get_str()[span], "m1");
        }
        _ => panic!("expected a positioned binding error"),
    }
    assert!(matches!(
        SelectionExpr::new_with_resolver("within 0.5 of deleted", &mut resolver),
        Err(SelectionParserError::Variable { .. })
    ));
}
#[test]
fn rejects_unrelated_resolver_and_missing_resolver() {
    let mut resolver = Resolver::new();
    let other = Resolver::new();
    let expr = SelectionExpr::new_with_resolver("within 0.5 of m1", &mut resolver).unwrap();
    assert!(matches!(
        default_system().select_with(&expr, &other),
        Err(SelectionError::Parser(
            SelectionParserError::ResolverMismatch
        ))
    ));
    assert!(matches!(
        SelectionExpr::new("m1"),
        Err(SelectionParserError::SyntaxError(_))
    ));
    assert!(matches!(
        default_system().select(&expr),
        Err(SelectionError::Parser(
            SelectionParserError::ResolverRequired { .. }
        ))
    ));
}
#[test]
fn public_foreign_result_and_same_residue_are_rejected() {
    let mut resolver = Resolver::new();
    for query in ["m1:(all)", "all or m1", "same residue as m1"] {
        let expr = SelectionExpr::new_with_resolver(query, &mut resolver).unwrap();
        assert!(
            matches!(
                default_system().select_with(&expr, &resolver),
                Err(SelectionError::Parser(SelectionParserError::SourceMismatch))
            ),
            "{query}"
        );
    }
}
#[test]
fn empty_references_and_foreign_self_include() {
    let mut resolver = Resolver::new();
    let expr =
        SelectionExpr::new_with_resolver("within 0.5 of m1.r_empty", &mut resolver).unwrap();
    assert!(matches!(
        default_system().select_with(&expr, &resolver),
        Err(SelectionError::EmptyExpr(_))
    ));
    let expr =
        SelectionExpr::new_with_resolver("within 0.5 self of m1.r2[3]", &mut resolver).unwrap();
    assert_eq!(
        indices(default_system().select_with(&expr, &resolver).unwrap()),
        vec![0]
    );
}
#[test]
fn compiled_subselection_and_mutable_binding() {
    let mut default = default_system();
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver("within 0.5 of m1", &mut resolver).unwrap();
    let domain = default.select([0usize, 1].as_slice()).unwrap();
    assert_eq!(
        indices(default.sub_select_with(&domain, &expr, &resolver).unwrap()),
        vec![0]
    );
    assert_eq!(
        default
            .bind(&domain)
            .select_bound_with(&expr, &resolver)
            .unwrap()
            .get_index_slice(),
        [0]
    );
    assert_eq!(
        default
            .select_bound_mut_with(&expr, &resolver)
            .unwrap()
            .get_index_slice(),
        [0, 2]
    );
}
#[test]
fn centers_and_nth_atom_use_reference_source() {
    let mut resolver = Resolver::new();
    for query in [
        "within 0.5 of com of m1.r2[3]",
        "within 0.5 of cog of m1.r2[3]",
        "within 0.5 of pos 0 of m1.r2[3]",
    ] {
        let expr = SelectionExpr::new_with_resolver(query, &mut resolver).unwrap();
        assert_eq!(
            indices(default_system().select_with(&expr, &resolver).unwrap()),
            vec![0],
            "{query}"
        );
    }
}
#[test]
fn legacy_expression_reuse_is_not_stale() {
    let expr = SelectionExpr::new("name CA").unwrap();
    assert_eq!(
        indices(system(&[0.0, 1.0], &["CA", "CB"]).select(&expr).unwrap()),
        [0]
    );
    assert_eq!(
        indices(system(&[0.0, 1.0], &["CB", "CA"]).select(&expr).unwrap()),
        [1]
    );
    let sys = default_system();
    let domain = sys.select([1usize, 2].as_slice()).unwrap();
    assert_eq!(indices(sys.sub_select(&domain, &expr).unwrap()), [2]);
}
#[test]
fn checked_custom_view_rejects_bad_arrays_and_indices() {
    let sys = default_system();
    assert!(
        SelectionSourceView::from_buffers(&sys.topology().atoms, &sys.state().coords[..1], None)
            .is_err()
    );
    for indices in [&[2usize, 1][..], &[1usize, 1][..], &[3usize][..]] {
        assert!(
            SelectionSourceView::from_parts(sys.topology(), sys.state(), Some(indices)).is_err()
        );
    }
    let view = SelectionSourceView::from_parts(sys.topology(), sys.state(), Some(&[0, 2])).unwrap();
    let expr = SelectionExpr::new("name CA").unwrap();
    assert_eq!(
        indices(expr.evaluate(&view, &Resolver::new()).unwrap()),
        [0, 2]
    );
}
#[test]
fn numeric_ranges_and_colons_in_names_stay_compatible() {
    let sys = system(&[0.0, 1.0, 2.0], &["a:b", "CB", "CA"]);
    assert_eq!(indices(sys.select("index 0:1").unwrap()), [0, 1]);
    assert_eq!(indices(sys.select("name a:b").unwrap()), [0]);
}

#[test]
fn periodic_within_uses_outer_box() {
    let mut default = system(&[0.1, 1.0], &["CA", "CB"]);
    let mut frame = default.state().clone();
    frame.pbox = Some(PeriodicBox::from_vectors_angles(2.0, 2.0, 2.0, 90.0, 90.0, 90.0).unwrap());
    default.set_state(frame).unwrap();
    let mut resolver = Resolver::new();
    resolver.source = system(&[1.9], &["CA"]);
    let expr = SelectionExpr::new_with_resolver("within 0.3 pbc of m1", &mut resolver).unwrap();
    assert_eq!(indices(default.select_with(&expr, &resolver).unwrap()), [0]);
    assert!(matches!(
        default_system().select_with(&expr, &resolver),
        Err(SelectionError::Parser(SelectionParserError::PbcUnwrap(_)))
    ));
}
#[test]
fn same_source_scope_and_external_binding_adapter() {
    struct Adapter<'a> {
        table: &'a Resolver,
        source: &'a System,
    }
    impl SelectionResolver for Adapter<'_> {
        fn identity(&self) -> SelectionResolverId {
            self.table.identity()
        }
        fn prepare(&mut self, _: &str) -> Result<SelectionReferenceId, SelectionResolveError> {
            unreachable!("binding adapter does not prepare")
        }
        fn bind(
            &self,
            _: SelectionReferenceId,
        ) -> Result<SelectionSourceView<'_>, SelectionResolveError> {
            Ok(self.source.selection_view())
        }
    }
    let sys = default_system();
    let mut table = Resolver::new();
    let expr = SelectionExpr::new_with_resolver("m1:(name CA)", &mut table).unwrap();
    let adapter = Adapter {
        table: &table,
        source: &sys,
    };
    assert_eq!(indices(sys.select_with(&expr, &adapter).unwrap()), [0, 2]);
    let domain = sys.select([1usize, 2].as_slice()).unwrap();
    assert_eq!(
        indices(sys.sub_select_with(&domain, &expr, &adapter).unwrap()),
        [2]
    );
}
#[test]
fn local_union_centers_and_sparse_nth_atom() {
    let sys = default_system();
    assert_eq!(
        indices(sys.select("within 0.1 of pos 0 of index 2").unwrap()),
        [2]
    );
    assert_eq!(
        indices(
            sys.select("within 0.6 of cog of (index 0 or index 2)")
                .unwrap()
        ),
        [1]
    );
}
#[test]
fn nested_scope_resets_domain_and_binds_nested_references() {
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver(
        "within 0.5 of m1:(name CA and within 0.5 of m2)",
        &mut resolver,
    )
    .unwrap();
    assert!(matches!(
        default_system().select_with(&expr, &resolver),
        Err(SelectionError::EmptyExpr(_))
    ));
    assert_eq!(resolver.binds.get(), 2);
}
#[test]
fn centers_rebind_between_evaluations() {
    let mut resolver = Resolver::new();
    let default = default_system();
    let expr =
        SelectionExpr::new_with_resolver("within 0.5 of com of m1.r2[3]", &mut resolver).unwrap();
    assert_eq!(indices(default.select_with(&expr, &resolver).unwrap()), [0]);
    let mut frame = resolver.source.state().clone();
    frame.coords[1] = Pos::new(2.1, 0.0, 0.0);
    resolver.source.set_state(frame).unwrap();
    assert_eq!(indices(default.select_with(&expr, &resolver).unwrap()), [1]);
}

#[test]
fn reference_intersections_keep_their_source() {
    let mut resolver = Resolver::new();
    let expr =
        SelectionExpr::new_with_resolver("within 0.5 of (m1 and m1:(name CA))", &mut resolver)
            .unwrap();
    assert_eq!(
        indices(default_system().select_with(&expr, &resolver).unwrap()),
        [0]
    );
    assert_eq!(resolver.binds.get(), 1);
}

#[test]
fn periodic_scalar_reference_error_is_returned_without_panicking() {
    let mut default = default_system();
    let mut frame = default.state().clone();
    frame.pbox =
        Some(PeriodicBox::from_vectors_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap());
    default.set_state(frame).unwrap();
    let mut resolver = Resolver::new();
    let expr =
        SelectionExpr::new_with_resolver("dist pbc point com of m1 < 0.5", &mut resolver).unwrap();
    resolver.missing = true;
    assert!(matches!(
        default.select_with(&expr, &resolver),
        Err(SelectionError::Parser(
            SelectionParserError::Variable { .. }
        ))
    ));
}

#[test]
fn unqualified_reference_predicates_keep_subset_semantics() {
    let sys = default_system();
    assert!(matches!(
        sys.select("within 0.1 of (index 0 and x < x of com of name CA)"),
        Err(SelectionError::EmptyExpr(_))
    ));
}

#[test]
fn bare_references_match_scopes_and_share_prepared_binding() {
    let mut resolver = Resolver::new();
    let expr = SelectionExpr::new_with_resolver(
        "within 0.5 of (m1.r2[3] or m1.r2[3])", &mut resolver,
    ).unwrap();
    assert_eq!(resolver.names, ["m1.r2[3]"]);
    assert_eq!(indices(default_system().select_with(&expr, &resolver).unwrap()), [0]);
    assert_eq!(resolver.binds.get(), 1);
}

#[test]
fn bare_references_support_groups_centers_and_boolean_operands() {
    for text in ["within 0.5 of m1[0]", "within 0.5 of com of m1.r2[3]",
        "within 0.5 of (m1 and m1:(name CA))"] {
        let mut resolver = Resolver::new();
        let expr = SelectionExpr::new_with_resolver(text, &mut resolver).unwrap();
        assert_eq!(indices(default_system().select_with(&expr, &resolver).unwrap()), [0]);
    }
}

#[test]
fn external_hooks_parse_nested_operands_and_evaluate_each_call_in_current_source() {
    struct Extension { id: SelectionResolverId, names: Vec<String>, peer: System, calls: Cell<usize>, last_size: Cell<usize>, invalid: bool }
    impl SelectionResolver for Extension {
        fn identity(&self) -> SelectionResolverId { self.id }
        fn parse_atom(&self, input: &str, position: usize) -> Option<ForeignParse> {
            let rest = &input[position..];
            if rest.starts_with("broken") {
                return Some(ForeignParse { end: position + 6, name: "broken".into(), error: Some("invalid external syntax".into()) });
            }
            if !rest.starts_with("peer") { return None; }
            Some(ForeignParse::matched(position + 4, "peer"))
        }
        fn parse_prefix(&self, input: &str, position: usize) -> Option<ForeignParse> {
            input[position..].starts_with("touches").then(|| ForeignParse::matched(position + 7, "touches"))
        }
        fn prepare(&mut self, text: &str) -> Result<SelectionReferenceId, SelectionResolveError> {
            self.names.push(text.into());
            Ok(self.id.reference(self.names.len() - 1))
        }
        fn bind(&self, reference: SelectionReferenceId) -> Result<SelectionSourceView<'_>, SelectionResolveError> {
            assert_eq!(self.names[reference.slot()], "peer");
            Ok(self.peer.selection_view())
        }
        fn evaluate_prefix(&self, reference: SelectionReferenceId, current: SelectionSourceView<'_>,
            operands: &[SelectionSourceView<'_>]) -> Result<Vec<usize>, SelectionResolveError> {
            assert_eq!(self.names[reference.slot()], "touches");
            if self.invalid { return Ok(vec![99]); }
            self.last_size.set(current.atom_storage().len());
            self.calls.set(self.calls.get() + 1);
            Ok(current.iter_particle().filter(|p| operands.iter().any(|view|
                view.iter_particle().any(|q| (p.pos - q.pos).norm() < 0.5)))
                .map(|p| p.id).collect())
        }
    }
    let mut resolver = Extension { id: SelectionResolverId::new(), names: vec![],
        peer: system(&[0.1], &["P"]), calls: Cell::new(0), last_size: Cell::new(0), invalid: false };
    for (text, expected) in [
        ("touches (peer or name CB)", vec![0, 1]),
        ("touches peer or name CA", vec![0, 2]),
        ("touches (touches (peer or name CB))", vec![0, 1]),
        ("index 1 and touches (peer or name CB)", vec![1]),
        ("touches name CA or touches name CB", vec![0, 1, 2]),
    ] {
        let expr = SelectionExpr::new_with_resolver(text, &mut resolver).unwrap();
        assert_eq!(indices(default_system().select_with(&expr, &resolver).unwrap()), expected, "{text}");
    }
    let expr = SelectionExpr::new_with_resolver("touches peer or touches name CB", &mut resolver).unwrap();
    resolver.calls.set(0);
    assert_eq!(indices(default_system().select_with(&expr, &resolver).unwrap()), [0, 1]);
    assert_eq!(resolver.calls.get(), 2, "a shared operation head must not cache different operands together");
    for text in ["peer", "peer or index 0", "index 0 or peer", "index 99 and peer"] {
        let expr = SelectionExpr::new_with_resolver(text, &mut resolver).unwrap();
        assert!(matches!(default_system().select_with(&expr, &resolver),
            Err(SelectionError::Parser(SelectionParserError::SourceMismatch))), "{text}");
    }
    assert!(matches!(SelectionExpr::new_with_resolver("broken", &mut resolver),
        Err(SelectionParserError::Variable { error, .. }) if error.message == "invalid external syntax"));
    let scoped = SelectionExpr::new_with_resolver("within 0.5 of peer:(touches peer)", &mut resolver).unwrap();
    assert_eq!(indices(default_system().select_with(&scoped, &resolver).unwrap()), [0]);
    assert_eq!(resolver.last_size.get(), 1, "a source filter must pass its own source to the prefix evaluator");
    let expr = SelectionExpr::new_with_resolver("touches peer", &mut resolver).unwrap();
    resolver.invalid = true;
    assert!(matches!(default_system().select_with(&expr, &resolver),
        Err(SelectionError::Parser(SelectionParserError::Variable { error, .. }))
            if error.message.contains("out-of-bounds")));
}
