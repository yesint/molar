# pymolar documentation coverage audit

Audit scope: the Python interface in this checkout, reviewed on 2026-10-04.
The root README, both precision packages, Python type files, Rust binding docstrings, the
Sphinx generator, and the documentation workflow were checked. The audit used
Rust LSP navigation, source inspection, and a freshly built Python extension.
The external published site was not used as evidence for this checkout.

## Baseline

The native module exports 16 classes and 6 functions. All 22 objects and all
121 public class members have docstrings. The Python package adds
`AnalysisTask` and three PBC constants. Thus the main problem was missing
usage contracts and incorrect examples, not missing method summaries.

| Area | Gap before the change | Documentation added or corrected |
| --- | --- | --- |
| Discovery | README covered installation but had no capability map or reading order | Agent guide, `llms.txt`, README links and root example time fix |
| Selections | No complete query guide; local/global indices and different tuple endpoints were unclear | Query examples, endpoint rules, set operations, expansion, splitting |
| Data sharing | Frame assignment and frame swaps appeared interchangeable | Shared ownership, independent references, swap effects, selection lifetime |
| NumPy | Copy/view behavior, dtype, memory layout, and array lifetime were incomplete | Shapes, copies, views, Fortran storage for assignment, view lifetime |
| Analysis | Geometry methods had summaries but no complete workflows | Geometry, fitting, RMSD, contacts, PBC, secondary structure examples |
| IO | Seek-before-iteration, field selection, subset writes, optional formats, and retained context handlers were unclear | IO rules, current format table, NDX and output examples |
| Chemistry | Charge docs incorrectly used partial `charge` as formal-charge input; MOL2 was cited despite no handler | Separate formal charges, preparation stages, warnings, supported elements, charge limitations |
| Tasks | Hook lifecycle, current state source, skip/end rules, and multi-file timing were incomplete; f64 hooks had no docs | CLI task example and current behavior in both packages |
| API names | RMSD example called nonexistent `rmsd`; SASA docs called nonexistent methods | `rmsd_py` example and explicit SASA availability |
| Type files | Missing `Sasa` and `SysParticleIterator`; missing 12 methods; phantom `SasaResults` and `_pbc` methods; wrong return types and defaults | Both native type files updated; package type files added |
| Published reference | Generator supplied only API pages and omitted special methods and constants | Checked-in guides, operators, constants, downloadable types, hosted `llms.txt`, strict build option |

The 12 missing methods were `FileHandler.read_state_pick`, `skip_to_last`,
`write_state_pick`, and `Sel.apply_ff`, `apply_charges`, `dssp`, `dssp_string`,
`ss`, `ss_string`, `unwrap_simple`, `whole_chains`, and `whole_residues`.
`Particle.id`, selection indices, iterator item types, transform return types,
optional atom types, and `FileHandler.__next__` were also checked.

## Force-field API follow-up

The binding now exposes all public `molar_ff` operations on the object types
supported by Rust. `Topology` has `apply_ff` and `apply_charges`, matching
`System` and `Sel`. Preparation remains System-only, as in Rust.

The new `FFType` and `ChargeModel` choices work alongside existing strings.
`PrepareOptions` exposes all connectivity, bond-order, search, and hydrogen
settings. Nested option objects are mutable references. `FFError`,
`ChargeError`, and `BondPerceptionError` preserve error variants and fields
through `kind` and `details`, while remaining subclasses of `ValueError`.
The source charge documentation was corrected: partial charges preserve total
formal charge, including for ions.

Both type files and the dedicated force-field guide cover these additions.
The updated native API has 28 classes, 6 functions, and 154 public class
members, including inherited exception methods.

## Remaining implementation limits

These limits were recorded in the original documentation audit and remain
outside the force-field API additions.

- Python cannot start or update a SASA calculation. `Sasa` is only an exported
  result type. Membrane bindings are disabled. Several Rust-only capabilities
  are not Python methods.
- `Particle.pos = value` writes to the first coordinate slot. Use component
  setters or selection coordinates for an indexed edit.
- `Sel.set_same_*` changes the full backing topology. Use selected atom views
  for restricted edits.
- `FileHandler.write_state(sel)` and `write_topology(sel)` use backing data.
  Use `write(sel)` or `sel.save()` for a subset.
- The context manager does not close a retained handler. Release it to finish
  output. Iteration consumes the reader and disables other handler operations.
- Position views require careful lifetime control across data swaps and changes
  to atom count. Recreate particles and selections after structural changes.
- The coordinate setter requires the documented dtype and Fortran layout.
- System/state box assignment requires an existing box. Use an all-atom
  selection's box setter to add a box. Inspect `state.box` for absence.
- `AnalysisTask.state` can hold previous-frame data after swaps. Multi-file
  `--add-time` reads that object. Zero bounds are ignored, and end frame limits
  count processed frames after skipping.
- Python does not expose direct bond editing, the separate formal-charge field,
  or general force-field parameter export.

## Maintenance

Keep the checked-in guides, both native type files, and both package type
files synchronized with changes to bindings. Run the docs example check against
fresh wheels, then build Sphinx with `--strict`. The API reference is generated
from the installed extension, so an old wheel can hide new methods or retain
old docstrings. Prefer `--skip-install` after an explicit wheel installation.

## Validation results

- Built fresh single- and double-precision wheels.
- Executed all 10 checked-in Python examples with each wheel, including full
  trajectory loops and the CLI task example: 20 successful example checks.
- Checked both type files against 16 native classes, 6 functions, and all 121
  public class members; confirmed identical native type files.
- Built both HTML references with Sphinx warnings treated as errors.
- Confirmed coordinate copy/layout, shared state, frame assignment, indexed
  position assignment, and full-topology bulk-setter behavior at runtime.
- Checked Python syntax, documentation links, and Git whitespace errors.

### Force-field follow-up validation

- Built both precision wheels with the new bindings.
- Passed 71 Python API and regression tests, including both precision packages.
- Passed 12 documentation examples per package: 24 successful example checks.
- Checked all 28 classes, 6 functions, and 154 public members against both type files.
- Built both Sphinx references with warnings treated as errors.
- Added the force-field API tests to the documentation CI workflow.
