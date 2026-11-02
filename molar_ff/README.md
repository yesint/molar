# molar_ff

Force-field atom typing, partial charges, and bonded mechanics for [MolAR](https://github.com/yesint/molar).

Assigns [GAFF](https://ambermd.org/antechamber/gaff.html) (General Amber Force Field)
atom types to a molecule, reproducing AmberTools `antechamber`'s `atomtype` perception.

Typing is exposed through the [`ApplyFF`] trait, implemented for any MolAR object that
provides atoms and bonds (`System`, `Topology`, and bound selections):

```rust,ignore
use molar::prelude::*;
use molar_ff::{ApplyFF, FFType};

let mut sys = System::from_file("ligand.sdf")?;
sys.apply_ff(FFType::Gaff)?;               // types the whole system
// or type just a selection (treated as the molecule):
sys.select_bound_mut("resname LIG")?.apply_ff(FFType::Gaff)?;
```

The assigned type is written into each atom's `type_name`.

**Input requirement:** molecules must already carry bond orders (e.g. from an SDF/mol2
file). Inputs without bond orders (PDB/GRO) return an error — `molar_ff` does not perceive
bond orders.

Both GAFF and GAFF2 (`FFType::Gaff2`) are implemented and tested.

## Espaloma partial charges

`ApplyCharges` predicts Espaloma partial charges. This function is always present;
it does not require a Cargo feature. The old `espaloma` feature name remains as an
empty compatibility feature for downstream manifests.

The source network is in `assets/espaloma_charge.onnx`. Production code uses a
reviewed fixed network and generated raw weights. It does not link an ONNX parser or
a general graph executor. See `assets/README.md` for the weight format and the exact
regeneration procedure.

The implementation is suitable for large sparse molecular systems:

- Bond messages use `BondAdjacency` and need `O(E)` graph work. The implementation
  does not create an `n × n` adjacency matrix.
- Hidden-state memory is proportional to the atom count. Two `[n, 128]` buffers are
  live during a message layer.
- Native builds use SIMD matrix kernels. Large row sets also use Rayon tasks. Small
  molecules stay serial to prevent task overhead.
- Temporary neighbor storage is limited to one 512-atom block per active task.

The normal test suite compares raw network output and final charges with the Python
reference and with the full molecule corpus. An ignored 20,000-atom sparse test is
available for release-mode performance checks; the command is in `assets/README.md`.

## Automatic Amber bonded mechanics

`amber::BondedModel::gaff(&source, FFType::Gaff2)` types and parameterizes a molecule
without changing its atoms. It needs no AmberTools installation or parameter files.
Input must have explicit hydrogens and known bond orders. Use `PrepareForFF` first
if connectivity, bond orders, or hydrogens need to be prepared.

```rust,ignore
use molar::prelude::*;
use molar_ff::{FFType, amber::BondedModel};

let sys = System::from_file("ligand.sdf")?;
let model = BondedModel::gaff(&sys, FFType::Gaff2)?;
let positions: Vec<[f64; 3]> = sys.state().coords.iter()
    .map(|p| [p.x as f64, p.y as f64, p.z as f64]).collect();
let result = model.evaluate(&positions)?;
println!("Bonded energy: {} kJ/mol", result.energy.total());
// result.forces contains -dE/dx in kJ/mol/nm.
```

The model includes harmonic bonds and angles, multi-term proper torsions, and
periodic impropers. Coordinates are nm, angles are radians, and energies are
kJ/mol. Amber harmonic energies are `k * displacement²`, without a factor of 1/2.
Torsions use `k * (1 + cos(n * phi - phase))`; `k` includes Amber's IDIVF divisor.

Parameter completion is automatic:

1. Use specific terms, equivalent types, and wildcard terms from the bundled
   GAFF 1.81 or GAFF2 2.2.30 tables.
2. Use the weighted type correspondences in Amber's `PARMCHK.DAT` to find analogues.
3. Where possible, estimate angles from symmetric angles and bond lengths using
   Amber's empirical formula. A zero-stiffness bond with a known equilibrium length
   can use the empirical bond formula. Eligible planar centers can use Amber's
   default improper term.
4. If an interaction still has no parameter, return `BondedError::Unresolved` with
   **all** missing interactions, source atom indices, and types. No partial model or
   silent zero-force placeholder is returned.

Every term carries its resolution method, matched types, score, and source location.
`model.estimated_terms()` lists terms obtained through completion. Estimated values
remain approximations, not new fitted force-field parameters. See
[the parameter and validation notes](assets/amber/README.md) for exact behavior.

`BondedModel::from_typed` uses existing atom types. Its `ParameterPolicy::Automatic`
uses the same completion, while `Strict` permits only direct and wildcard matches.
`AmberParameters::with_frcmod` is an optional import API for existing Amber workflows;
it is not required by the automatic API. A *frcmod* file is just a parameter override
file; Amber's *parmchk2* program contains the parameter-completion logic.

Models also accept bond-complete bound selections. Term and coordinate indices are
local to the model; `atom_indices()` maps them to source indices. Selections that cut
bonds are rejected. Coordinate arrays must use model-local order and be unwrapped;
the evaluator does not apply periodic boundaries. Invalid or undefined geometry
returns an error.

This module supplies **bonded** energy and forces. It does not yet include Lennard-
Jones/electrostatic interactions, exclusions/1–4 scaling, a minimizer, or viewer
integration. It must not be treated as a complete GAFF force field. Improper topology
is inferred from the GAFF/GAFF2 planar-center flags in `PARMCHK.DAT`; arbitrary custom
force fields with different improper topology rules are outside this API's scope.

Tests compare all four energy components and forces with independent OpenMM
reference calculations, check finite-difference gradients, and cover missing-term
errors. A separate coverage audit exercises all 597 molecules in both force fields:

```sh
cargo test -p molar_ff
cargo test -p molar_ff --test amber_bonded automatic_corpus -- --ignored --nocapture
```
