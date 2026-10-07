# LAMMPS data files

MolAR reads and writes one polymer structure per `.data` file. It supports the
`molecular` atom style and the `bond` and `angle` styles with the same atom columns.
This support is intended for KG models.

## Reading

Supported data:

- Atom rows: `id molecule type x y z`, with optional `ix iy iz` integers on every row.
- Orthogonal and restricted triclinic boxes (`xy xz yz`).
- Numeric atom type IDs, molecule IDs, optional type masses, and bond connectivity.

The header must declare positive atom and atom type counts and all three box
bounds. Atom IDs can be sparse and rows can be out of order. MolAR sorts atoms by molecule
ID, then atom ID. Bond endpoints follow the new order. Positive molecule IDs become
`resid` values and contiguous molecule ranges. ID zero means an unassigned atom;
these atoms do not form a molecule range. Molecule IDs must fit a nonnegative `i32`.
Atoms have name `B`, residue name `MOL`, and no assigned chemical element.

MolAR stores unwrapped coordinates. It applies image flags with the box vectors:

```text
r_unwrapped = r + ix*A + iy*B + iz*C
```

It then subtracts the box origin. This translation preserves distances and box
shape. A six-column row retains its position relative to the box origin; MolAR
does not infer images from bond connectivity.

`Masses` can appear before or after `Atoms`. If present, it must supply a positive
mass for every declared type. If absent, each atom gets a mass of one input mass
unit. Bond type IDs are checked, then discarded. Chemical bond order remains
`BondOrder::Unspecified`.

The reader discards velocities, coefficients, angles, dihedrals, impropers, and
other recognized fixed-length sections. They are not restored during writing.
Unknown custom sections, type label sections, other atom styles, general
triclinic boxes, variable-length `Bodies`, compressed files, dump trajectories,
and restart files are not supported.

Use `read()` for topology and state together. You can also use `read_topology()`
and `read_state()` in either order. State iteration yields one frame, then EOF.

## Units

The default assigns **1 nm per input distance unit** and **1 atomic mass unit per
input mass unit**. This is a KG convention. It does not infer a physical scale
from LJ reduced units. The title comment, including any `units = ...` text, is
not used to select scales. Time is zero; velocities and forces are absent.

Set explicit positive, finite conversion scales when needed. For `real` or
`metal` distances, use a length scale of `0.1` to convert Angstroms to nanometers.
For an LJ model, supply the physical values of sigma and the reference mass.
Writing uses the inverse scales. Use the same options to read the output.

```rust
use molar::prelude::*;

let options = LammpsOptions { length_scale: 0.1, mass_scale: 1.0 };
let (top, state) = FileHandler::open_lammps("polymer.data", options)?.read()?;
let system = System::new(top, state)?;
FileHandler::create_lammps("output.data", options)?.write(&system)?;
```

`open_lammps` and `create_lammps` select the format directly, independent of the
extension. `from_lammps_reader` accepts a byte source with explicit options.
The usual `open`, `create`, and `from_reader("data", ...)` use default scales.

```python
import pymolar as mol

reader = mol.FileHandler("polymer.data", "r",
                         lammps_length_scale=0.1, lammps_mass_scale=1.0)
top, state = reader.read()
writer = mol.FileHandler("output.data", "w",
                         lammps_length_scale=0.1, lammps_mass_scale=1.0)
writer.write((top, state))
```

The Python scale keywords require a `.data` extension. `System("polymer.data")`
and `system.save("output.data")` use default scales.

## Writing

The writer requires a valid orthogonal or restricted triclinic periodic box,
positive masses, nonnegative molecule IDs, and equal atom and coordinate counts.
Atoms of one type must have the same mass. Missing atom types use type 1. Used
numeric atom type IDs are retained. Unused type IDs below the highest used ID
receive a mass of one input mass unit; other unused declarations are omitted.

The output has a zero box origin, consecutive atom and bond IDs, `Masses`,
`Atoms # molecular`, and, when present, `Bonds`. Coordinates are wrapped, with
calculated image flags that retain unwrapped positions. Image flags must fit
`i32`. All bonds use type 1. The output has no chemical bond orders or interaction
parameters. Set the required interaction coefficients in the LAMMPS input script.

Saving a selection includes only bonds whose two endpoints are selected. Their
endpoints refer to the output atom IDs. A second structure write is an error.
Separate writes require `write_topology()` before `write_state()`; the latter
writes the complete file. Output is flushed after a complete write.
