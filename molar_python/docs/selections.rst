Selections and query syntax
===========================

Index rules
-----------

.. code-block:: python

   import pymolar as mol

   sys = mol.System("molar/tests/protein.pdb")
   all_atoms = sys()
   first_three = sys((0, 3))       # global 0, 1, 2: end excluded
   sparse = sys([0, 2, 4, 6])     # global indices; sorted and unique
   two = sparse([1, 3])           # local 1, 3 -> global 2, 6
   middle = sparse((1, 2))        # local 1, 2 -> global 2, 4: end included
   last = sparse[-1]              # Particle with global id 6
   assert last.id == 6
   assert 6 in sparse             # membership uses global indices

``System`` tuple ranges exclude the end. ``Sel`` tuple ranges include the
end. Lists passed to ``Sel`` contain local indices. String ``index`` queries
always refer to global atom indices, including inside a sub-selection.
Do not use Python slice objects as selection arguments.

Queries
-------

Query keywords are case-sensitive. Values must match the input naming.
Use parentheses to make Boolean grouping explicit. The current parser gives
``and`` and ``or`` the same precedence. Names and residue names support lists
of exact values or regular expressions between slashes. Regular expressions
match the full value; shell wildcard syntax is not supported.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Query
     - Meaning
   * - ``all``
     - All atoms in the current search scope
   * - ``protein``, ``backbone``, ``sidechain``, ``water``
     - Built-in groups based on molecular names
   * - ``hydrogen``, ``noh``, ``polh``, ``apolh``, ``now``
     - Hydrogens, non-hydrogens, polar or apolar hydrogens, non-water
   * - ``name CA`` or ``name N CA C O``
     - One or more atom names
   * - ``name /C.*/``
     - Atom names that start with C
   * - ``resname ALA GLY``
     - Residue names
   * - ``resid 1:10`` or ``resid 1 4 8``
     - Residue identifiers; colon ranges include both endpoints
   * - ``resindex 0:9``
     - Internal residue indices
   * - ``index 0:9``
     - Global atom indices
   * - ``chain A B``
     - Chain identifiers
   * - ``protein and (name CA or name N)``
     - Boolean expression
   * - ``not water``
     - Complement within the query scope
   * - ``mass > 12`` or ``0 <= z < 2``
     - Numeric property or coordinate comparisons
   * - ``occupancy > 0.5`` or ``beta < 30``
     - Occupancy (alias ``occ``) or B-factor (alias ``bfactor``)
   * - ``within 0.35 of (resname LIG)``
     - Neighbors in nm; excludes the target atoms
   * - ``within 0.35 self of (resname LIG)``
     - Includes the target atoms
   * - ``within 0.35 pbc 110 of (resname LIG)``
     - Periodic neighbors in x and y; requires a box
   * - ``same residue as (name CA)``
     - Whole residues that contain selected atoms
   * - ``same chain as (name CA)``
     - Whole chains that contain selected atoms
   * - ``dist point [0, 0, 0] < 1``
     - Distance from a point in nm

Numeric expressions support ``x``, ``y``, ``z``, ``mass``, ``charge``,
``vdw``, occupancy and B-factor. ``vx``, ``vy``, ``vz`` and ``fx``, ``fy``,
``fz`` require velocity or force data in the state. Arithmetic supports
``+``, ``-``, ``*``, ``/``, ``^`` and ``abs()``, ``sqrt()``, ``sin()``,
``cos()``. Spatial expressions also support ``dist line`` and ``dist plane``,
``com of (query)``, and ``pos n of (query)``. Use the tested forms above for
routine tasks; do not assume that another library's selection syntax works.

Shared selections and set operations
------------------------------------

.. code-block:: python

   import pymolar as mol

   sys = mol.System("molar/tests/protein.pdb")
   protein = sys("protein")
   ca = protein("name CA")
   backbone = sys("backbone")
   union = protein | backbone
   common = protein & backbone
   sidechains = protein - backbone
   non_ca = ~ca
   residues = ca.whole_residues()
   chains = ca.whole_chains()
   residue_groups = protein.split_resindex()
   chain_groups = protein.split_chain()

Binary set operations require the same topology object. Separate systems
loaded from the same file do not meet that requirement. The result uses the
left selection's state. Complement, ``whole_residues()``, and
``whole_chains()`` can include atoms outside the original selection.

``split_molecule()`` uses topology connectivity. Check that bonds are present
before using it; coordinate proximity alone does not define a molecule.
String sub-selections search within their parent selection. Re-evaluate
spatial queries after each frame update. Recreate selections after atom
removal or addition.
