# Cutoff nonbonded mechanics

Lengths are nm; energies are kJ/mol; charges are elementary-charge units. For a
pair, `R = radius_i + radius_j`, `epsilon = sqrt(epsilon_i * epsilon_j)` and
`K = 138.935457644382 q_i q_j`:

```text
U_LJ(r) = epsilon [(R/r)^12 - 2 (R/r)^6]
U_C(r)  = K [1/r - 1/rc + (r-rc)/rc^2]       r < rc
F_C(r)  = K [1/r^2 - 1/rc^2]                 r < rc
```

Ordinary Coulomb energy and force vanish at `rc`. This is the undamped
shifted-force method (the zero-damping limit of damped shifted force), not a
quintic Coulomb potential switch and not PME. Position-independent electrostatic
self terms are not included; energies describe the stated pair potential.

Ordinary LJ is multiplied by `S=1-10x^3+15x^4-6x^5`, where
`x=(r-rs)/(rc-rs)` between the switch start `rs` and cutoff `rc`. Below `rs`,
`S=1`; at/above `rc`, energy and force are zero. Forces include the derivative
of S. Both first and second derivatives of the LJ potential are continuous.

Defaults are `rs=1.0 nm`, `rc=1.2 nm`. These are configurable cleanup-model choices,
**not Amber-recommended universal values**. A 0.8 nm real-space cutoff with PME
also includes reciprocal-space electrostatics and cannot be reproduced by cutting
off the total Coulomb interaction at 0.8 nm. Cutoff-only electrostatics can alter
conformational preferences, particularly for charged or large molecules.

Shortest graph distances of one/two bonds exclude both interactions. Distance
three gives the GAFF SCEE=1.2 and SCNB=2.0 scaling, once per pair. These special
pairs use the unmodified LJ and Coulomb forms at every distance. Rings use the
shortest path, so a longer three-bond path never overrides a shorter exclusion.
This API implements GAFF's fixed 1–4 factors; arbitrary per-torsion SCEE/SCNB
settings and pair-specific LJEDIT parameters are outside its scope. LJEDIT imports
are explicitly rejected.

## Sources and choice of cutoff

- [Amber energy documentation](https://amberhub.chpc.utah.edu/energy/) identifies
  8 Angstrom as the default direct-space Ewald cutoff.
- [Lee et al., JCTC 2016, doi:10.1021/acs.jctc.5b00935](https://pmc.ncbi.nlm.nih.gov/articles/PMC4712441/)
  compares CHARMM36 protocols, including 10–12 Angstrom force switching. This is
  not a universal GAFF switching prescription.
- [OpenMM standard forces](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#lennard-jones-interaction)
  defines the quintic LJ potential switch and Lorentz–Berthelot mixing.
- [Whitmore et al., JCTC 2025, doi:10.1021/acs.jctc.5c00646](https://pmc.ncbi.nlm.nih.gov/articles/PMC12515373/)
  distinguishes potential switching, force switching, and potential shifting and
  measures their different cutoff effects. Its free-energy results do not directly
  validate this vacuum cleanup model.
- [Fennell and Gezelter, JCP 2006, doi:10.1063/1.2206581](https://doi.org/10.1063/1.2206581)
  treats shifted-force electrostatics and pairwise alternatives to lattice sums.
- [Piana et al., PLOS ONE 2012, doi:10.1371/journal.pone.0039918](https://pmc.ncbi.nlm.nih.gov/articles/PMC3386949/)
  examines electrostatic treatment and 8–12 Angstrom cutoffs in protein folding;
  smooth cutoff forces do not by themselves guarantee Ewald-equivalent results.

## Independent validation

`tests/data/amber_nonbonded_reference.json` contains 36 cases generated with
OpenMM 8.6.1's Reference platform and ParmEd 4.3.1. OpenMM's `NonbondedForce`
computes LJ, mixing, switching, and graph exclusions/1–4 pairs. A separate
`CustomNonbondedForce` differentiates the stated shifted-force Coulomb energy;
`CustomBondForce` handles full, scaled 1–4 Coulomb terms. The generator imports
molecular topology/types/coordinates, not Rust energies, forces, or exclusions.

The cases include both GAFF versions, nine molecules, unlike-type pairs around
both cutoff boundaries, zero-LJ sites, and a 1–4 pair beyond the cutoff. Tests also
check finite differences, cell crossings, small rings, missing parameters,
invalid geometry, and final minimizer energy/force consistency.

From the workspace root (Python dependencies as in the bonded reference notes):

```sh
MOLAR_AMBER_EXPORT=/tmp/amber-models.json cargo test -p molar_ff \
  --test amber_bonded export_reference_input -- --ignored
python molar_ff/tools/check_amber_nonbonded.py --models /tmp/amber-models.json
cargo test -p molar_ff --test amber_nonbonded --test amber_minimize
```
