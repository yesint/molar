# Amber bonded data and automatic completion

The data files are copied unchanged from
[AmberClassic](https://github.com/Amber-MD/AmberClassic/tree/76b00b2e9d7dc41ab1ce7c7aa2fe52b3608ec002),
commit `76b00b2e9d7dc41ab1ce7c7aa2fe52b3608ec002`. `source.json` records the
upstream paths and SHA-256 checksums. Preserve these files and their notices when
updating them. `UPSTREAM-LICENSE` reproduces the upstream root license and its
component-specific scope; it does not replace the MolAR source license.

| File | Use |
| --- | --- |
| `gaff.dat` | GAFF 1.81, May 2017: bonded parameters |
| `gaff2.dat` | GAFF2 2.2.30, October 2025: bonded parameters |
| `PARMCHK.DAT` | Equivalent types, correspondence penalties, improper-center flags |
| `PARM_BLBA_GAFF.DAT` | GAFF empirical bond and angle constants |
| `PARM_BLBA_GAFF2.DAT` | GAFF2 empirical bond and angle constants |

The existing GAFF2 typing definition in MolAR is identical to the definition at
this upstream revision. The database reads bonded and Lennard-Jones sections. See
[nonbonded methods](NONBONDED.md) for cutoff mechanics and minimization.

## Resolution

The Rust resolver implements the automatic completion strategy in Amber's
[`parmchk2.c`](https://github.com/Amber-MD/AmberClassic/blob/76b00b2e9d7dc41ab1ce7c7aa2fe52b3608ec002/src/antechamber/parmchk2.c).
It runs in-process on native and WebAssembly targets. No external executable,
network service, or user-supplied `frcmod` is required.

Search order is exact, equivalent-specific, wildcard, equivalent-wildcard,
weighted analogue-specific, then weighted analogue-wildcard. Bond and proper
torsion scores include conjugation penalties. Angle-center replacements receive
the larger center penalty. Angle analogues above `THRESHOLD_BA` are rejected.
Ties favor fewer wildcards, fewer substitutions, then original table order.
This is deterministic, but is not a promise of byte-for-byte `parmchk2` output.

If direct transfer is unavailable, the resolver uses the empirical angle rule:

```text
theta = (theta(a-b-a) + theta(c-b-c)) / 2
d = ((r_ab - r_bc) / (r_ab + r_bc))^2
k = 143.9 Z_a C_b Z_c exp(-2d) / ((r_ab + r_bc) sqrt(theta))
```

The calculation above uses Angstrom bond lengths and radians; its kcal/mol
constant is converted to kJ/mol at the API boundary. The source symmetric angles
and bonds must exist or be resolvable. A zero-stiffness bond entry can use
`exp(constant - exponent * ln(r))` when it has a known equilibrium length and
supported element constants. The resolver does not invent missing bond lengths.
An eligible planar center with no resolved improper can use Amber's default
1.1 kcal/mol, phase 180 degrees, periodicity 2 improper.

Improper centers occupy the third position. All six peripheral permutations are
searched, and the model preserves the order used for the chosen parameter.
Negative periodicities continue a Fourier series; the final positive periodicity
ends it. Proper torsion amplitudes include division by IDIVF. A later explicit
entry for the same tuple replaces the whole previous series.

Every estimated term has a method, score, matched types, and source locations.
A source line of zero denotes a formula or rule rather than a parameter row.
Strict mode rejects estimated terms. Automatic mode fails atomically with all
unresolved interactions; neither mode silently omits a missing bonded term.

## Validation and regeneration

`tests/data/amber_bonded_reference.json` contains 18 reference calculations for
nine molecules in both force fields. It includes topology, perturbed coordinates,
four energy components, and all force components. The set covers exact, wildcard,
equivalent, analogue, default-improper, and empirical-angle resolution. OpenMM
8.6.1 on the Reference platform calculates the energies and forces; ParmEd 4.3.1
reads the independent Amber parameter files and `parmchk2` results.

The harness first compares Rust-selected parameters with `parmchk2` output.
Because `frcmod` rounds values, it then uses the original table precision for
energy calculations and independently evaluates the empirical-angle formula.
It does not use Rust energies or forces to generate expected values. Coordinates
and interaction ordering come from the export, so this tests parameter completion
and evaluation, not an independent reconstruction of LEaP topology.

Two upstream behaviors are deliberately not reproduced:

- At the pinned revision, `parmchk2` sets `improperparmnum2 = impropernum`, using
  the molecule's number of impropers as a parameter-table bound. The reference
  build applies `tools/parmchk2-reference.patch` to use `improperparmnum` instead.
  The production Rust resolver always searches the complete table.
- A positional peripheral comparison in `parmchk2` can miss an available direct
  improper wildcard. For example, a carbonyl center with c3/o/os neighbors has
  the published `X-X-c-o` term. The reference harness uses the base table for
  direct impropers across all peripheral permutations, before considering any
  `parmchk2` estimate.

To regenerate, build `parmchk2` from the pinned upstream source with the supplied
one-line patch. Its Amber home must contain the pinned `dat/antechamber/PARM*.DAT`
and `dat/leap/parm/gaff*.dat` files. In a separate Python environment install
`openmm==8.6.1`, `parmed==4.3.1`, and NumPy. From the workspace root:

```sh
MOLAR_AMBER_EXPORT=/tmp/amber-models.json cargo test -p molar_ff \
  --test amber_bonded export_reference_input -- --ignored
python molar_ff/tools/check_amber_bonded.py \
  --models /tmp/amber-models.json \
  --parmchk2 /path/to/patched/parmchk2 --amber-home /path/to/amber
cargo test -p molar_ff --test amber_bonded
cargo test -p molar_ff --test amber_bonded automatic_corpus -- --ignored --nocapture
```

The export uses the default coordinate precision; do not enable `f64` when
regenerating. The fixture itself tests both precision builds. Python and Amber
are only needed to regenerate references, not to run Rust tests or use the crate.

The full corpus audit resolves 595 of 597 molecules in each force field. Hydrogen
chloride and hydrogen fluoride remain explicit missing-bond errors. Passing this
audit means parameters are available; it does not establish that every estimated
parameter is accurate for every molecule.
