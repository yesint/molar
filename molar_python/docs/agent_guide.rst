Agent guide
===========

Scope and sources
-----------------

``pymolar`` is the Python interface to MolAR. Use it to read molecular
structures, process trajectories, select atoms, edit coordinates, measure
geometry, assign secondary structure, and prepare small molecules for atom
typing and charge prediction. This guide describes the Python interface in
this checkout. The Rust interface has additional capabilities.

Read this page first. Then read :doc:`selections` for query syntax and
:doc:`workflows` for complete examples. Use :doc:`api_reference` for method
signatures. The package also supplies ``molar.pyi`` for static inspection.
Use ``help(pymolar.Sel)`` to check an installed version. Do not infer Python
method names from the Rust API.

The HTML build includes :download:`API types <_static/molar.pyi>`,
:download:`the coverage audit <_static/coverage-audit.md>`, and
:download:`installation instructions <_static/README.md>`. The site root
also contains ``llms.txt`` for agent discovery.

Install and select precision
----------------------------

.. code-block:: sh

   python -m pip install pymolar
   # Optional separate package:
   python -m pip install pymolar-f64

Use ``import pymolar as mol`` for single precision. Use
``import pymolar_f64 as mol`` for double precision. The two packages can be
installed together. They use separate native types. Keep all objects in one
workflow in the same package. Obtain the NumPy dtype from ``sel.coords.dtype``.

Task to API map
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 45 25

   * - Task
     - API
     - Result or change
   * - Load a structure
     - ``System(path)``
     - Topology and first frame
   * - Load topology and coordinates separately
     - ``FileHandler.read_topology()``, ``read_state()``, ``System(top, st)``
     - System with equal atom counts
   * - Select atoms
     - ``sys(query)``, ``sys(indices)``, ``sys((start, end))``, ``sys()``
     - ``Sel`` with shared data
   * - Read or edit one atom
     - ``sys[i]``, ``sel[i]``, ``iter_atoms()``
     - Mutable particle or atom view
   * - Read or edit coordinates
     - ``sel.coords``, ``sel.translate(vector)``, ``sel.apply_transform(tr)``
     - Coordinate copy or explicit edit
   * - Process trajectory frames
     - Iterate ``FileHandler(path, "r")``; use ``sys.replace_state_deep(st)``
     - Existing selections see new frame
   * - Measure geometry
     - ``sel.com()``, ``cog()``, ``gyration()``, ``inertia()``, ``min_max()``
     - Centers, radius, moments, axes, bounds
   * - Fit structures and measure RMSD
     - ``fit_transform()``, ``fit_transform_matching()``, ``rmsd_py()``, ``rmsd_mw()``
     - Opaque transform or distance in nm
   * - Find contacts
     - ``distance_search(cutoff, sel1, sel2=None, dims=None)``
     - Global index pairs and distances
   * - Assign secondary structure
     - ``sel.dssp()``, ``dssp_string()``, ``ss(algo)``, ``ss_string(algo)``
     - One code per residue
   * - Work with fragments
     - ``split_resindex()``, ``split_chain()``, ``split_molecule()``
     - Lists of shared selections
   * - Expand selections
     - ``whole_residues()``, ``whole_chains()``
     - Atoms from the full system
   * - Build or change a system
     - ``System()``, ``append(atom, pos)``, ``append(sel)``, ``remove(arg)``
     - Changes atom count
   * - Prepare chemistry
     - ``perceive_connectivity()``, ``perceive_bond_orders()``, ``add_hydrogens()``, ``prepare_for_ff()``
     - Bonds, formal charges, hydrogens
   * - Assign atom types and partial charges
     - ``sys.apply_ff("gaff2")``, ``sel.apply_charges("espaloma")``
     - Changes atom properties
   * - Save structures or trajectories
     - ``sys.save(path)``, ``sel.save(path)``, ``FileHandler.write(data)``
     - Format chosen by extension
   * - Read or write GROMACS index groups
     - ``NdxFile(path).get_group_as_sel(name, sys)``, ``sel.to_gromacs_ndx(name)``
     - Converts between NDX and Python indices
   * - Run a command-line analysis
     - Subclass ``AnalysisTask``; construct it to run
     - Calls task hooks for frames

Data and units
--------------

* Coordinates and distances use **nm**. Time uses **ps**. Mass uses **Da**.
  Partial charge uses elementary charge units. Box angles use degrees.
  Readers convert file units, for example Angstrom coordinates in XYZ.
* ``Topology`` contains atoms and connectivity. ``State`` contains a frame.
  Obtain these objects from readers or system properties; they have no public
  Python constructors. ``System()`` creates an empty system.
* Python ``Sel`` holds indices and references to topology and state. It is
  already usable for analysis. There is no Python ``bind()`` step.
* ``Atom()`` creates a detached atom. ``particle.atom`` and ``iter_atoms()``
  return mutable ``AtomView`` objects. These views write to the topology.
  ``AtomView`` is a returned type, not a package-level constructor.
* ``resid`` is the residue identifier from the input. ``resindex`` identifies
  residues in internal order. Use ``split_resindex()`` to separate residues
  even when residue identifiers repeat across chains.
* ``type_name`` and ``type_id`` can be ``None``. The Python interface does not
  expose the separate integer formal-charge field or direct bond editing.

Rules for correct use
---------------------

1. Atom indices are zero-based. A selection has sorted, unique indices.
   ``sel[i]`` uses a local index; ``particle.id`` and ``sel.index`` are global.
   See :doc:`selections` for the different tuple range endpoints.
2. Query results contain fixed indices. Coordinate-dependent queries do not
   update by themselves. Evaluate the query again on each required frame.
   Avoid empty queries and ranges. Empty text-query results raise errors.
3. ``sel.coords`` is a copy with shape ``(3, len(sel))``. To write it back, use
   an array with the same dtype and **Fortran-contiguous** storage:
   ``sel.coords = np.asfortranarray(coords, dtype=sel.coords.dtype)``.
   The current setter reads columns from consecutive memory locations.
4. ``particle.pos`` and arrays from ``iter_pos()`` are writable views with
   shape ``(3,)``. Copy them for stored results. Release coordinate views
   before frame swaps, append, removal, or hydrogen addition. Recreate
   selections and particles after atom-count changes.
5. ``sys.state = frame`` changes only that system's state reference. Existing
   selections retain their old reference. ``replace_state_deep(frame)`` swaps
   frame data in the shared state, so existing selections see the new data.
   The input ``frame`` receives the old data. Read time from ``sys.time``
   after the swap. Atom order must match the topology.
6. PBC is off by default for centers, contacts, and geometry measurements.
   Pass ``dims=mol.PBC_FULL`` or ``pbc=True`` where supported. A box is required.
   Box matrix **columns** are box vectors. ``shortest_vector()`` and
   ``closest_image()`` default to all periodic dimensions.
7. Selection edits change shared data. ``sel.system`` gives the full backing
   system, not an independent system restricted to those atoms. Load a second
   structure for an independent reference. Never fit a trajectory to a
   reference that shares its state.
8. Use ``fit_transform(mobile, reference)`` with equal-length, corresponding
   atoms. It does not apply the transform. RMSD functions measure current
   coordinates without fitting. Apply the transform explicitly first.
9. Call seek methods before trajectory iteration. Iteration consumes the
   handler's reader. Read, write, seek, ``stats``, and ``file_name`` methods
   then raise ``TypeError``. Use a separate handler if needed.
10. Use complete molecules for atom typing and charge prediction. Inspect
    bond-order warnings before accepting chemistry assignments.

Current interface limits
------------------------

``rmsd_py`` is the exported RMSD name; ``rmsd`` is not exported.
``gyration(pbc=True)``, ``inertia(pbc=True)``, and
``principal_transform(pbc=True)`` are the Python calls. Separate methods with
``_pbc`` suffixes are not exported.

``Sasa`` is exported as a result type, but the current Python API supplies
neither a public constructor nor ``Sel.sasa()`` nor ``Sasa.update()``.
Do not generate a SASA workflow for this version. Membrane analysis bindings
are disabled. Rust analysis tasks, parallel selection methods, direct
velocity/force arrays, and general force-field parameter export are not
Python APIs.

Known behavior to account for: ``Particle.pos = value`` currently writes to
the first coordinate slot. Use ``particle.x``, ``y``, and ``z`` setters or
write a selection coordinate array instead. ``set_same_*`` methods currently
apply to the full backing topology; iterate selected atom views for restricted
property edits. ``FileHandler.write_state(sel)`` and ``write_topology(sel)``
use the backing state/topology. Use ``sel.save()`` or ``write(sel)`` for a subset.
``System.append(sel)`` copies atoms and coordinates but does not copy bonds.
Use a selection from the same backing system with ``remove(sel)``.
``FileHandler.__exit__`` does not close a retained handler; release it to
finish buffered output. These are implementation limits, not recommended
usage patterns.

For a missing box, inspect ``sys.state.box`` and catch ``AttributeError``.
Avoid ``sys.box`` or ``sel.box`` reads until a box exists. To add a new box,
set ``sys().box = PeriodicBox(...)``. The system and state box setters require
an existing box. Errors from file operations are usually ``OSError``;
selection creation uses ``TypeError`` on ``System`` and ``RuntimeError`` on
``Sel``. Selection algebra uses ``ValueError`` for incompatible topology or
empty intersection, difference, or complement. Treat error text as diagnostic
information, not as a stable API.
