# Selection sources and variable resolvers

A selection has one default source. Unqualified predicates use that source.
Variables supply other sources for reference operands. For example:

```text
protein and within 0.5 of m0.r1
within 0.5 of m1:(name CA)
within 0.5 of m1[0]
within 0.5 of m0.r2[3]:(name CA)
```

Molar treats `m0.r1`, `m1[0]` and `m0.r2[3]` as opaque strings. It does not parse
molecule, representation, group or frame identifiers. The caller owns that
micro-grammar. Molar recognizes the token before `:` and an optional parenthesized
selection after it. The token must not contain whitespace, `:`, parentheses, `&`
or `|`. A scope with no group, such as `m0.r1:`, means all atoms in that source's
allowed domain. The colon is optional for a bare variable, so `m0.r1` is
equivalent to `m0.r1:`. Bare tokens start with an ASCII letter or underscore and
exclude comparison operators and commas; selection keywords retain their existing
meaning. Use the colon form for an otherwise ambiguous variable name. Existing ranges such as `index 0:10` and colons in atom-name values
remain valid.

## Three separate roles

- `SelectionSource` is an adapter trait. A system, bound selection or custom object
  implements it to supply a view. It does not resolve variable names.
- `SelectionSourceView` is the concrete, borrowed data handle used by evaluation.
  It contains full atom storage and coordinate arrays, an optional subset, and
  optional bonds, box, velocities and forces. It owns no molecular data.
- `SelectionResolver` owns the prepared-reference table. Its `prepare` method
  interprets a variable name once during expression construction. Its `bind`
  method supplies a view of current data during each evaluation.

A detached `Sel` contains only indices. Bind it before using it as a source.
`System`, bound selection types, `SelBoundParts`, and views implement
`SelectionSource`. Custom sources can use checked `from_parts` or `from_buffers`
constructors. Subset indices must be sorted, unique and in bounds. Empty source
domains are allowed. Optional arrays must be absent or full length.

## Prepare once, bind on every evaluation

`SelectionExpr::new_with_resolver(text, &mut resolver)` calls `prepare` once per
distinct variable string in that expression. It stores the variable's prepared
ID and byte span in the AST. It does not call `bind` or retain borrowed source data.

`system.select_with(&expr, &resolver)` calls `bind` at most once per reached ID
in that evaluation. The view cache is discarded when evaluation returns. The
next evaluation binds again, so changes to frames, representation membership or
source existence take effect. The compiled AST is retained intact; evaluated
nodes are cached only within an evaluation.

Each prepared table has a `SelectionResolverId`. Create it once and preserve it
in evaluation adapters that borrow the table and the current scene. Use
`identity.reference(slot)` to create IDs for table entries. Molar rejects a
reference or evaluation resolver with a different table identity. The caller
must keep slots stable for the lifetime of compiled expressions.

This small example resolves one opaque name. A viewer's resolver can instead
prepare a parsed molecule/representation/frame description and bind it against
its current scene:

```rust
use molar::prelude::*;

struct References<'a> {
    identity: SelectionResolverId,
    other: &'a System,
}

impl SelectionResolver for References<'_> {
    fn identity(&self) -> SelectionResolverId { self.identity }

    fn prepare(&mut self, name: &str)
        -> Result<SelectionReferenceId, SelectionResolveError>
    {
        if name == "m0.r2[3]" {
            Ok(self.identity.reference(0))
        } else {
            Err(SelectionResolveError::new("unknown source reference"))
        }
    }

    fn bind(&self, reference: SelectionReferenceId)
        -> Result<SelectionSourceView<'_>, SelectionResolveError>
    {
        if reference != self.identity.reference(0) {
            return Err(SelectionResolveError::new("unknown prepared reference"));
        }
        Ok(self.other.selection_view())
    }
}

# fn example(default: &System, other: &System) -> Result<(), Box<dyn std::error::Error>> {
let mut resolver = References {
    identity: SelectionResolverId::new(),
    other,
};
let expr = SelectionExpr::new_with_resolver(
    "within 0.5 of m0.r2[3]:(name CA)",
    &mut resolver,
)?;
let selected = default.select_with(&expr, &resolver)?;
// Existing APIs also accept a resolved definition:
let bound = default.select_bound(expr.with_resolver(&resolver))?;
# let _ = (selected, bound);
# Ok(())
# }
```

`select_with` and related methods take a compiled expression. Construction needs
mutable access for preparation; evaluation needs only a shared binding resolver.
The methods are `select_with`, `select_bound_with`, `sub_select_with`, and
`select_bound_mut_with`. `expr.evaluate(&source, &resolver)` accepts a custom
default source. `expr.with_resolver(&resolver)` works with existing methods that
accept `SelectionDef`, without changing that trait.

`SelectionExpr::new` and existing selection methods remain available for
expressions without variables. Compiled expressions now also work in
subselections.

## Source and coordinate rules

- Public `Sel` results contain indices from the default source only. A foreign
  scope used as the final selection or ordinary Boolean filter is an error.
- `within` may consume reference scopes from other sources. Reference unions can
  contain several sources. Reference intersections require matching sources. Numeric atom indices from different sources are never
  combined as if they described the same atoms.
- A scoped group evaluates within that source's allowed domain. A scope resets
  the domain, including any candidate subset from an outer `and` expression.
- Source identity is shared atom storage plus coordinate-buffer identity. Aliases
  and subsets over the same data share source identity; distinct coordinate
  buffers are distinct sources.
- `same residue as` and `same chain as` require the current source. Equal residue
  or chain values do not establish identity across sources.
- `com`, `cog`, and `pos N of` use the reference source's coordinates and masses.
  These operations require a single source. Same-source reference unions are
  combined before measurement; mixed-source centers and positional references
  report an error. `pos N` addresses the sorted selected atoms.
- Distances and coordinates use nanometers. Sources must share a coordinate
  frame. Explicit periodic `within` uses the outer source's box, with no fallback
  to a reference source's box. Periodic centers use the measured source's box.
- `within ... self` includes reference atoms only when they belong to the current
  source and its allowed domain. It never adds foreign atom indices.
- Empty reference sets are valid. Final empty selections retain the existing
  `EmptyExpr` error. Missing velocities, forces or a required periodic box produce
  their existing errors. Sources without bonds behave like bondless systems.

Preparation and binding errors carry the variable string and its byte span in
the trimmed expression. Source deletion can therefore be reported at evaluation
without rebuilding a variable registry. The viewer is still responsible for
triggering evaluation when referenced data changes, and for detecting cycles if
binding requires recursive representation evaluation.

Mutable binding follows Rust's normal borrowing rules. If a resolver borrows the
system to be modified, first evaluate to a detached `Sel`, release the resolver's
borrows, then bind mutably.
