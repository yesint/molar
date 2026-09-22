# Bond perception and hydrogen addition plan

## Goal

MolAR must prepare a chemical graph when an input format gives only atoms and
coordinates, or when its bonds do not have chemical orders. The first target is
feature parity with the PDB preparation path used by Antechamber. Later work can
add more chemistry models without changing the public data flow.

The feature belongs in `molar`. Connectivity, bond orders, formal charges, and
hydrogen atoms are general molecular data. `molar_ff` consumes this data for
force-field atom typing and charge calculation. It must not own the only copy of
the preparation logic.

## Required behavior

The complete preparation path has five separate stages:

1. Perceive atom connectivity from coordinates.
2. Assign heavy-atom bond orders and formal charges from a connection table.
3. Calculate implicit-hydrogen counts and report chemical ambiguities.
4. Perceive rings and aromaticity from the assigned Kekule graph.
5. Optionally add explicit hydrogen atoms and their coordinates.

Each stage has a non-mutating calculation operation. A result is applied only
after validation. This makes failure atomic and lets applications inspect an
ambiguous result before they change a topology.

## Public API

The planned high-level API is:

```rust
let bonds = perceive_connectivity(&system, &ConnectivityOptions::default())?;
system.set_bonds(bonds)?;

let assignment = assign_bond_orders(system.topology(), &BondOrderOptions::default())?;
assignment.apply_to(system.topology_mut())?;

let addition = plan_hydrogen_addition(&system, &HydrogenOptions::default())?;
addition.apply_to(&mut system)?;
```

`System` convenience methods will provide the same operations when direct
mutable topology access is not public.

### Connectivity result

`perceive_connectivity` returns `BondStorage`. All created bonds have
`BondOrder::Unspecified`. It does not assign orders and it does not add atoms.

`ConnectivityOptions` controls:

- covalent-radius tolerance in nanometers;
- the minimum accepted atom distance;
- periodic dimensions;
- maximum-coordination cleanup;
- later, element-pair overrides and metal handling.

The implementation uses the existing cell-list distance search. It first makes
all radius-qualified candidate edges. It then removes the longest edges from
atoms that exceed a safe coordination limit. Candidate and final edge order is
deterministic.

### Bond assignment result

`BondAssignment` is parallel to the source topology:

```rust
pub struct BondAssignment {
    source_pairs: Vec<[usize; 2]>,
    bond_orders: Vec<Option<BondOrder>>,
    formal_charges: Vec<Option<i32>>,
    implicit_hydrogens: Vec<u8>,
    warnings: Vec<PerceptionWarning>,
}
```

`None` means that the existing topology value must stay unchanged. The source
pairs protect against applying a result after a structural edit. `apply_to`
checks all sizes and pairs before it writes any value. Bond-order writes do not
invalidate `BondAdjacency`.

Implicit-hydrogen counts are initially calculated result data. They are not a
new topology column. They can be calculated again from assigned Kekule bond
orders and formal charges. The assignment stores them so that callers can check
and use the exact solver result.

## Bond-order and formal-charge solver

The solver will follow the observable Antechamber `bondtype` behavior. We will
not copy source or tables until their exact license is verified. Reference tests
will compare output with Antechamber.

### Input modes

- `PreserveKnown`: keep each concrete input order and solve only unspecified
  bonds.
- `ReassignAll`: treat all non-aromatic input orders as variables.
- Aromatic input is Kekulized before the valence search, or rejected when it
  cannot be Kekulized.
- An optional total formal charge constrains the complete component or molecule.

### Valence states

Each atom gets an ordered set of candidate states. A state contains target
valence, formal charge, maximum coordination, implicit-hydrogen range, and a
penalty. Tables cover common organic elements first and expand to the same useful
range as Antechamber APS data.

Explicit hydrogen is part of the graph and consumes valence. Missing hydrogen is
a solver variable:

```text
target valence = assigned bond-order sum + implicit hydrogen count
```

Hydrogen bonds are constrained to single order. Hydrogen cannot have more than
one neighbor.

### Search

The solver operates per connected component:

1. Apply fixed bond orders and fixed formal charges.
2. Remove impossible valence states.
3. Propagate atoms or bonds that have only one possible value.
4. Select the most constrained remaining atom or bond.
5. Search its alternatives with branch-and-bound.
6. Minimize the sum of atom-state penalties and secondary deterministic costs.
7. Keep equally scored distinct solutions so ambiguity is visible.

The search has explicit limits for states, branches, and elapsed work. A limit
produces an error or warning, never a silent partial assignment.

### Special chemistry

General valence search alone cannot select all resonance, tautomer, and
protonation states. Ordered constraints will cover at least:

- carboxylate, nitro, phosphate, sulfate, and guanidinium groups;
- amides and other conjugated carbonyl groups;
- common aromatic five- and six-membered rings;
- nitrile, isocyanate, azide, and similar linear groups;
- standard amino-acid and nucleic-acid residue templates.

Templates take priority only when residue and atom identity is reliable. Unknown
residues use the general solver. Conflicting template and graph data produces a
diagnostic.

## Protein PDB files without hydrogen

PDB residue templates set standard polymer bonds, bond orders, and expected
formal charges. Geometry perception fills missing links and handles ligands.
Terminal state, histidine tautomer, acidic protonation, disulfides, and modified
residues can be ambiguous. Options will let the caller select a pH model or give
explicit residue-state overrides. The default must report unresolved choices.

PDB `CONECT` bonds and file bond orders are evidence. `PreserveKnown` keeps them.
The coordinate stage can merge perceived bonds with known bonds instead of
replacing them when requested.

## Adding hydrogens

Hydrogen addition is a separate structural transaction because it changes atom
indices, `State`, selections, molecule ranges, and adjacency.

`HydrogenAddition` contains new `Atom` rows, positions, parent atom indices,
new single bonds, and an old-to-new index map. The first version supports common
tetrahedral, trigonal, and linear geometry. Ring and conjugated geometry use the
local heavy-atom frame. A seeded deterministic orientation is used when rotation
around one bond is not fixed.

Application updates `Topology` and `State` together. It also updates molecule
ranges and invalidates adjacency. Existing velocities and forces either receive
zero entries for new atoms or cause an error, controlled by an option. Hydrogen
removal or replacement is a separate explicit mode.

Hydrogen addition does not select a protonation state. It materializes the
implicit-hydrogen counts selected by bond assignment.

## Diagnostics

Errors identify atom and bond indices and include the relevant element, current
order, valence, and component charge. Warnings include:

- two or more equal best assignments;
- a charge chosen without an explicit total-charge constraint;
- unsupported or weakly supported elements;
- coordination cleanup that removed a distance candidate;
- a residue template mismatch;
- missing coordinates for hydrogen placement.

No preparation warning is written only to a log. It is part of the returned
result.

## Test plan

### Unit tests

- Covalent radii, distance thresholds, minimum distance, periodic edges, and
  coordination cleanup.
- Assignment application, stale-result rejection, and no partial mutation on
  error.
- Each valence-state table row and special chemical constraint.
- Explicit and implicit hydrogen valence accounting.
- Hydrogen geometry for linear, trigonal, tetrahedral, aromatic, and terminal
  cases.

### Reference tests

- Small molecules with all major functional groups against Antechamber.
- The existing GAFF parity set after bonds are removed or orders are erased.
- SDF round trips with known formal charges and orders.
- Hydrogen-free PDB proteins, terminal variants, histidine states, waters, ions,
  disulfides, cofactors, and unknown ligands.
- mmCIF equivalents of the PDB cases.

Reference fixtures store the tool version, command, input charge, and expected
connection table. License-compatible generated output is checked into the test
data, not external program code.

### Invariants

- Perception is deterministic.
- A successful assignment leaves no unspecified order in its solved scope.
- The formal-charge sum matches the requested charge.
- Each explicit hydrogen has one single bond.
- Adding explicit hydrogen reduces the corresponding implicit count to zero.
- Applying a result after atom or bond structural change fails without mutation.
- Existing concrete orders survive `PreserveKnown` mode.

## Delivery sequence

1. Connectivity options, covalent radii, cell-list candidate search, cleanup,
   `System` integration, and tests.
2. Transactional `BondAssignment` storage, validation, application, and tests.
3. Common-organic valence states and the component search engine.
4. Antechamber reference harness and parity fixes.
5. Functional-group and aromatic constraints.
6. Protein and nucleic-acid residue templates.
7. Hydrogen-addition planning, geometry, and atomic `System` application.
8. `molar_ff` opt-in preparation wrapper; keep its low-level force-field apply
   operation strict.
9. Python bindings and end-user documentation.

The first two items create stable public boundaries. Later solver and hydrogen
work can improve without changing file handlers or force-field code.
