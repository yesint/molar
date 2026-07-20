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
