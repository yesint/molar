# molar_ff

Force-field atom typing for [MolAR](https://github.com/yesint/molar).

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

GAFF is implemented and validated first; GAFF2 (`FFType::Gaff2`) is planned.

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
