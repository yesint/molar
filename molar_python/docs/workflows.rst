Workflows
=========

The examples use ``pymolar``. For double precision, change the import to
``pymolar_f64``. Run examples with repository test files from the workspace
root. Use your own files in production.

Read, measure, edit, and save
-----------------------------

.. code-block:: python

   import numpy as np
   import pymolar as mol

   sys = mol.System("molar/tests/protein.pdb")
   protein = sys("protein")
   ca = protein("name CA")
   print(len(sys), len(ca), sys.time)
   print("center / nm:", ca.cog())
   print("center of mass / nm:", ca.com())
   print("radius of gyration / nm:", ca.gyration())
   moments, axes = ca.inertia()   # shapes (3,), (3, 3); moments in Da*nm^2
   lo, hi = ca.min_max()         # each shape (3,), in nm

   coords = ca.coords           # independent copy, shape (3, len(ca))
   coords[0, :] += 0.1
   ca.coords = np.asfortranarray(coords, dtype=ca.coords.dtype)
   ca.translate(np.array([0.0, 0.1, 0.0], dtype=ca.coords.dtype))
   for atom in ca.iter_atoms():
       atom.bfactor = 10.0       # changes only these selected atoms
   ca.save("ca.pdb")

Use known masses for mass-weighted measurements. Do not assume every format
supplies correct masses, atom types, charges, connectivity, or a periodic box.
``principal_transform(pbc=False)`` returns a transform for principal-axis
alignment; apply it with ``apply_transform()``.

Stream a trajectory and keep selections valid
---------------------------------------------

.. code-block:: python

   import pymolar as mol

   sys = mol.System("molar/tests/protein.pdb")
   ca = sys("protein and name CA")
   rows = []
   reader = mol.FileHandler("molar/tests/protein.xtc", "r")
   for frame in reader:
       sys.replace_state_deep(frame)
       rows.append((sys.time, ca.gyration()))
   print("processed frames:", len(rows))

The topology and trajectory must have the same atom count and atom order.
The frame swap also changes time and box. The input ``frame`` then holds the
old state; do not use its time as the new frame time. This pattern also keeps
other selections attached to the shared state current. Recreate spatial
selections inside the loop when their membership must follow moving atoms.

Use ``read_state_pick(coords=True, velocities=False, forces=False)`` before
iteration to select fields on formats that support it, such as TRR. This is
a single-frame read; iteration uses the ordinary frame reader. A state with
coordinates omitted is unsuitable for coordinate analysis.

Fit and compute RMSD
--------------------

.. code-block:: python

   import pymolar as mol

   reference = mol.System("molar/tests/protein.pdb")
   mobile = mol.System("molar/tests/protein.pdb")
   ref_ca = reference("protein and name CA")
   mobile_ca = mobile("protein and name CA")
   mobile_all = mobile()
   values = []
   for frame in mol.FileHandler("molar/tests/protein.xtc", "r"):
       mobile.replace_state_deep(frame)
       transform = mol.fit_transform(mobile_ca, ref_ca)
       mobile_all.apply_transform(transform)
       values.append((mobile.time, mol.rmsd_py(mobile_ca, ref_ca)))
   print("RMSD values / nm:", values[:3])

The reference uses independent state data. ``fit_transform`` returns an
opaque ``IsometryTransform`` object; the package does not export its class
or a constructor. Keep the object and pass it to ``apply_transform``.
``rmsd_mw`` is the mass-weighted measurement. Both RMSD functions use the
current coordinates and require equal-size, corresponding selections.

``fit_transform_matching`` first aligns sequences of atom names, then fits
matching atoms. It does not perform a general chemical graph match. Repeated
atom names can make correspondence ambiguous. Resolve atom correspondence
before fitting when sequence matching is unsuitable. Neither fitting nor
RMSD makes molecules whole across periodic boundaries by itself.

Periodic boxes and contacts
---------------------------

.. code-block:: python

   import numpy as np
   import pymolar as mol

   sys = mol.System()
   for name, x in [("C1", 0.05), ("C2", 1.95), ("C3", 1.0)]:
       atom = mol.Atom()
       atom.name = name
       atom.atomic_number = 6
       atom.mass = 12.011
       sys.append(atom, np.array([x, 0.0, 0.0], dtype=np.float32))
   atoms = sys()
   atoms.box = mol.PeriodicBox([2.0, 2.0, 2.0], [90.0, 90.0, 90.0])
   pairs, distances = mol.distance_search(0.2, atoms, dims=mol.PBC_FULL)
   assert pairs.shape == (1, 2)
   assert np.isclose(distances[0], 0.1, atol=1e-6)
   box = sys.state.box
   assert np.isclose(box.distance(sys[0].pos, sys[1].pos, mol.PBC_FULL), 0.1)

For double precision, use ``np.float64`` for the construction arrays above.
``PeriodicBox(matrix)`` expects box vectors as **columns**. A matrix from an
interface that stores box vectors as rows must be transposed. Its lengths
use nm and ``PeriodicBox(lengths, angles)`` uses degrees.
``wrap_point`` returns a point in the primary cell. ``to_box_coords`` and
``to_lab_coords`` convert between fractional and Cartesian coordinates.
``get_matrix()``, ``get_box_extents()``, ``get_lab_extents()``,
``to_vectors_angles()``, and ``is_triclinic()`` inspect box geometry.
``shortest_vector`` returns a minimum-image displacement. ``closest_image``
returns the image closest to a target. ``distance_squared`` returns nm squared.

Contact results have shapes ``(n_pairs, 2)`` and ``(n_pairs,)``. Pair indices
are global in their respective systems, not local selection indices. Do not
assume an order for pairs. No contacts gives shapes ``(0, 2)`` and ``(0,)``.
A numeric cutoff must be finite and positive. PBC uses the first selection's
box; ensure both inputs are in compatible coordinates and boxes.

``distance_search("vdw", sel1, sel2)`` uses sums of van der Waals radii.
It requires two selections and valid element assignments. Single-selection
``"vdw"`` search raises ``NotImplementedError``. Numeric single-selection
search returns atom pairs without self-pairs. A search between overlapping
selections can include identical atoms; use disjoint selections or filter
results if needed.

``unwrap_simple()`` changes coordinates by imaging each atom against the
first selected atom. Use it for a small molecule whose extent fits in the
box. It does not traverse bonds to unwrap a long chain.

Secondary structure
-------------------

.. code-block:: python

   import pymolar as mol

   sys = mol.System("molar/tests/protein.pdb")
   protein = sys("protein")
   codes = protein.dssp()
   assert "".join(codes) == protein.dssp_string()
   for residue, code in zip(protein.split_resindex(), codes):
       first = residue[0]
       print(first.chain, first.resid, code)
   print(protein.ss_string("dss"))

Use all protein backbone atoms, not a CA-only selection. The algorithms are
``"dssp"`` (default), ``"dssp_gmx"``, and ``"dss"``. DSSP codes are H
(alpha helix), G (3-10 helix), I (pi helix), P (poly-proline II helix), E
(beta strand), B (beta bridge), T (turn), S (bend), ``~`` (coil), and ``=``
(break or missing backbone atoms). Non-protein or incomplete residues can
produce break codes. The result has one code per residue, not per atom.

Prepare a molecule and assign atom types
----------------------------------------

.. code-block:: python

   import math
   from pathlib import Path
   import pymolar as mol

   # XYZ input coordinates are in Angstrom. Readers convert them to nm.
   lines = ["6", "benzene heavy atoms"]
   for k in range(6):
       angle = math.pi * k / 3.0
       lines.append(f"C {1.39 * math.cos(angle):.4f} {1.39 * math.sin(angle):.4f} 0")
   Path("benzene.xyz").write_text("\n".join(lines) + "\n")
   sys = mol.System("benzene.xyz")
   assert sys.perceive_connectivity() == 6
   warnings = sys.perceive_bond_orders(infer_hydrogens=True, total_charge=0)
   print("bond-order diagnostics:", warnings)
   assert sys.add_hydrogens() == 6
   sys.apply_ff("gaff2")
   sys.apply_charges("espaloma")
   for atom in sys.iter_atoms():
       print(atom.name, atom.type_name, atom.charge)
   sys.save("benzene.sdf")

``perceive_connectivity(tolerance=0.045, min_distance=0.04, pbc=False)``
replaces connectivity from elements and coordinates; bonds initially have
unspecified order. Distances use nm. ``perceive_bond_orders`` assigns orders
and formal charges and returns diagnostic warning strings. Its default
requires explicit hydrogens; set ``infer_hydrogens=True`` for a structure
without hydrogens. ``total_charge`` constrains a single-fragment assignment.
Review the returned warnings.

``prepare_for_ff(infer_hydrogens=True, add_hydrogens=True, total_charge=0)``
combines missing-connectivity perception, bond-order assignment, and hydrogen
addition. It returns ``None`` and logs bond-order warnings. Use the separate
steps when diagnostics must be stored. Hydrogen addition appends atoms;
create selections again afterwards.

``apply_ff`` accepts ``"gaff"`` or ``"gaff2"`` and writes ``type_name``.
It assigns atom types, not full bonded or nonbonded parameters.
``apply_charges`` accepts ``"espaloma"`` and writes partial ``charge``.
It reads the separate formal charges from topology, not the existing partial
``charge`` values. Its current charge equilibration makes the predicted
charges sum to zero over the input system or selection. Do not treat that
result as charge-conserving for an ion. Supported elements are H, C, N, O,
F, P, S, Cl, Br, and I. For multiple molecules, call charge prediction on each
complete molecular selection when that is the intended scope. Bonds with
unspecified order or bonds crossing a selection boundary cause ``ValueError``.

Formats, writing, and index files
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Extensions
     - Input
     - Use
   * - ``.pdb``, ``.ent``, ``.gro``, ``.xyz``, ``.cif``, ``.mmcif``
     - Structure; supported frame reads
     - ``System(path)``, structure output
   * - ``.sdf``, ``.sd``, ``.mol``
     - Structure with bond orders and formal charges
     - Molecule preparation and output
   * - ``.xtc``, ``.dcd``, ``.trr``
     - Trajectory states
     - Use a separate topology; trajectory output
   * - ``.itp``
     - Topology
     - ``read_topology()``; combine with a state
   * - ``.tpr``
     - Topology and state
     - Requires the GROMACS runtime plugin; read only
   * - ``.cpt``
     - Checkpoint state
     - Requires the GROMACS runtime plugin; read only
   * - ``.nc``, ``.ncdf``
     - AMBER trajectory states
     - Requires a wheel built with ``molar/netcdf``

The extension chooses the handler. Each format supports only its relevant
operations. Seek and selective-field methods are format-dependent; an
unsupported operation raises an IO error. There is no MOL2 handler in the
current dispatch table. See the README for optional format build instructions.

.. code-block:: python

   from pathlib import Path
   import pymolar as mol

   sys = mol.System("molar/tests/protein.pdb")
   ca = sys("protein and name CA")
   Path("ca.ndx").write_text(ca.to_gromacs_ndx("CA"))
   restored = mol.NdxFile("ca.ndx").get_group_as_sel("CA", sys)
   assert restored.index.tolist() == ca.index.tolist()
   with mol.FileHandler("subset.pdb", "w") as writer:
       writer.write(ca)
   del writer                    # release retained handler and finish output

``to_gromacs_ndx`` returns text, not a file. NDX uses one-based indices;
reading and writing convert them automatically. The index file must match
the system's atom order. ``write`` accepts ``System``, ``Sel``, or
``(Topology, State)``. ``write_topology`` also accepts ``Topology``;
``write_state`` also accepts ``State``. For a subset, use ``write(sel)``
or ``sel.save(path)``. With a trajectory writer, call ``write_state(sys)``
once per frame. ``write_state_pick`` has ``coords``, ``velocities``, and
``forces`` switches; support depends on format and source data.

``FileStats`` is a snapshot with ``elapsed_time`` (``datetime.timedelta``),
``frames_processed`` (integer), and ``cur_t`` (ps). Access ``stats`` before
converting a reader to an iterator.

Command-line analysis tasks
---------------------------

Save this example as ``radius.py``:

.. code-block:: python

   import pymolar as mol

   class Radius(mol.AnalysisTask):
       def register_args(self, parser):
           parser.add_argument("--selection", default="protein and name CA")

       def pre_process(self):
           self.atoms = self.src(self.args.selection)
           self.rows = []

       def process_frame(self):
           self.rows.append((self.src.time, self.atoms.gyration()))

       def post_process(self):
           for time_ps, radius_nm in getattr(self, "rows", []):
               print(time_ps, radius_nm)

   if __name__ == "__main__":
       Radius()

.. code-block:: sh

   python radius.py -f molar/tests/protein.pdb molar/tests/protein.xtc --skip 2

Construction runs the full task. ``register_args(parser)`` runs before
argument parsing. ``pre_process()`` runs once on the first processed frame.
``process_frame()`` runs on each selected frame. ``post_process()`` runs after
reading, even if there were no processed frames. During processing, use
``self.src`` for current coordinates and time. ``self.state`` is the input
swap object and can contain previous-frame data after a swap.

.. list-table::
   :header-rows: 1

   * - Argument
     - Current behavior
   * - ``-f``, ``--files``
     - Topology file first, followed by at least one trajectory
   * - ``-b``, ``--begin``
     - Frame index, or integer time with ``ps``, ``ns``, ``us`` suffix; applied to each trajectory
   * - ``-e``, ``--end``
     - Processed-frame count limit, or inclusive time limit; not an ordinary trajectory slice
   * - ``--skip``
     - Process every Nth eligible frame; default 1; use a positive integer
   * - ``--log``
     - Log every N processed frames; default 100; use a positive integer
   * - ``--add-time``
     - Adds offsets between trajectory files

A frame count measures processed frames after skipping. Zero-valued begin
and end limits are ignored by the current truth checks. Time suffixes require
integers. Prefer the explicit trajectory loop when exact multi-file time or
frame boundary control is required: the current ``--add-time`` implementation
reads offsets from the swap object. Do not infer the Rust task API from these
Python hooks.
