Force-field API
===============

Python access to molar_ff
-------------------------

The Python interface exposes each public ``molar_ff`` operation on the object
types supported by Rust:

.. list-table::
   :header-rows: 1

   * - Operation
     - System
     - Topology
     - Sel
   * - ``apply_ff(ff)``
     - Yes
     - Yes
     - Yes
   * - ``apply_charges(model)``
     - Yes
     - Yes
     - Yes
   * - ``prepare_for_ff(..., options=None)``
     - Yes
     - Requires a System with coordinates
     - Prepare its backing System

Rust implements ``PrepareForFF`` on ``System`` only. Preparation can create
bonds and append hydrogens, so it needs a topology and coordinate state.
``Topology.apply_ff()`` and ``Topology.apply_charges()`` need no coordinates.
They edit the same topology seen by a system and its selections.

To prepare a ligand independently, load its structure as a separate System.
``sel.system`` gives the full backing system, so preparation on that object
acts on the full system. ``Sel.apply_ff`` and ``Sel.apply_charges`` require
complete molecules: a bond crossing the selection boundary is an error.
The low-level model and GAFF engine functions are private Rust implementation
details, not additional public operations.

Typed choices
-------------

``FFType.Gaff`` and ``FFType.Gaff2`` select GAFF and GAFF2 atom typing.
``ChargeModel.Espaloma`` selects the charge model. Existing case-insensitive
strings remain valid: ``"gaff"``, ``"gaff2"``, and ``"espaloma"``. Defaults
remain ``"gaff"`` and ``"espaloma"``. Use named enum values, not integers.
An unknown string raises ``ValueError``; an unsupported argument type raises
``TypeError``.

``apply_ff`` writes ``atom.type_name``. It does not assign a complete set of
force-field parameters. ``apply_charges`` writes the partial ``atom.charge``
and preserves the separate topology formal charges. The partial charges sum
to the total formal charge in scope, including for ions. Multiple fragments
are treated as one scope for charge equilibration; use complete per-molecule
selections when separate fragment charge totals are required.

Full preparation options
------------------------

``PrepareOptions`` exposes all Rust preparation fields through mutable Python
objects. Its fields are ``connectivity``, ``bond_orders``, and
``add_hydrogens``. Each default construction has independent nested objects.
Explicitly supplied nested objects are shared references. Thus edits such as
``options.bond_orders.limits.max_branches = 100000`` take effect.

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Field
     - Default
     - Meaning
   * - ``ConnectivityOptions.tolerance``
     - 0.045 nm
     - Added to the sum of covalent radii; finite and nonnegative
   * - ``ConnectivityOptions.minimum_distance``
     - 0.04 nm
     - Minimum bond distance; finite and nonnegative
   * - ``ConnectivityOptions.pbc``
     - ``[False, False, False]``
     - Periodic lattice directions; a box is required when any direction is enabled
   * - ``ConnectivityOptions.cleanup_overcoordination``
     - ``True``
     - Remove longest candidate edges from over-coordinated atoms
   * - ``BondOrderOptions.input_orders``
     - ``InputOrders.ReassignAll``
     - Standalone options reassign non-aromatic orders; PreserveKnown keeps known orders
   * - ``BondOrderOptions.hydrogens``
     - ``HydrogenPolicy.AllExplicit``
     - Require explicit H; InferFromGeometry allows missing H and requires coordinates
   * - ``BondOrderOptions.use_functional_groups``
     - ``True``
     - Use canonical orders for recognized functional groups
   * - ``BondOrderOptions.use_residue_templates``
     - ``True``
     - Use canonical orders for recognized polymer residues
   * - ``BondOrderOptions.total_charge``
     - ``None``
     - Optional total formal charge; requires one bonded fragment
   * - ``BondOrderOptions.limits``
     - ``SearchLimits()``
     - Per-fragment work limits
   * - ``SearchLimits.max_branches``
     - 5,000,000
     - Maximum visited search nodes per bonded fragment
   * - ``HydrogenOptions.zero_fill_dynamics``
     - ``True``
     - Give new hydrogens zero velocities/forces when dynamics are present
   * - ``PrepareOptions.add_hydrogens``
     - ``None``
     - Set to HydrogenOptions to append explicit hydrogens

``PrepareOptions().bond_orders.input_orders`` defaults to
``InputOrders.PreserveKnown``, as in Rust. This differs from the standalone
``BondOrderOptions()`` default. When supplying a custom BondOrderOptions object,
choose its ``input_orders`` explicitly if preservation is required. Both policies
kekulize aromatic input and keep the resulting aromatic-system bond orders.

Connectivity settings apply only when the topology has no bonds. To replace
existing connectivity, use ``System.perceive_connectivity`` first. Preparation
returns ``None``. Use the separate bond-order perception step when warning
strings must be collected. Preparation is a sequence of changes: if a later
step fails, earlier changes such as connectivity may remain. The system keeps
its topology and state objects after an error.

Pass full options with ``sys.prepare_for_ff(options=options)``. Existing
positional and keyword calls remain valid, such as
``sys.prepare_for_ff(True, True, 0)``. Do not combine full options with
non-default ``infer_hydrogens``, ``add_hydrogens``, or ``total_charge`` keywords;
the binding raises ``ValueError`` rather than silently selecting one set.

Complete example
----------------

.. code-block:: python

   from pathlib import Path
   import pymolar as mol

   # XYZ positions use Angstrom on disk. This creates a carbon-only methane input.
   Path("methane.xyz").write_text("1\nmethane heavy atom\nC 0 0 0\n")
   sys = mol.System("methane.xyz")
   options = mol.PrepareOptions(
       connectivity=mol.ConnectivityOptions(
           tolerance=0.045, minimum_distance=0.04, pbc=mol.PBC_NONE,
           cleanup_overcoordination=True,
       ),
       bond_orders=mol.BondOrderOptions(
           input_orders=mol.InputOrders.PreserveKnown,
           hydrogens=mol.HydrogenPolicy.InferFromGeometry,
           use_functional_groups=True, use_residue_templates=True,
           total_charge=0, limits=mol.SearchLimits(max_branches=100000),
       ),
       add_hydrogens=mol.HydrogenOptions(zero_fill_dynamics=True),
   )
   sys.prepare_for_ff(options=options)
   assert len(sys) == 5
   sys.topology.apply_ff(mol.FFType.Gaff2)
   sys.topology.apply_charges(mol.ChargeModel.Espaloma)
   assert all(atom.type_name for atom in sys.iter_atoms())
   assert abs(sum(atom.charge for atom in sys.iter_atoms())) < 1e-5
   sys().apply_ff("gaff")          # selection entry point also accepts strings
   sys.save("methane.sdf")

Release position views before preparation. Hydrogen addition changes atom
count; recreate selections afterwards to include the new atoms.

Typed errors
------------

``FFError``, ``ChargeError``, and ``BondPerceptionError`` subclass ``ValueError``.
Existing handlers that catch ValueError continue to work. Errors generated
by the library also expose ``kind`` (a Rust variant name) and ``details``
(a dictionary of variant fields). Use these fields instead of parsing the
human-readable message. Invalid argument names and conflicting arguments
remain plain ValueError instances.

.. code-block:: python

   from pathlib import Path
   import pymolar as mol

   Path("pair.xyz").write_text("2\ncarbon pair\nC 0 0 0\nC 1.5 0 0\n")
   sys = mol.System("pair.xyz")
   assert sys.perceive_connectivity() == 1
   try:
       sys.topology.apply_ff(mol.FFType.Gaff)
   except mol.FFError as error:
       assert isinstance(error, ValueError)
       assert error.kind == "MissingBondOrders"
       assert error.details["atoms"] == (0, 1)
       print(error.kind, error.details)

Typing kinds are ``MissingBondOrders``, ``OpenSelection``, ``InvalidAromatic``,
and ``UntypedAtom``. Charge kinds are ``MissingBondOrders``, ``OpenSelection``,
``Kekulize``, ``UnsupportedElement``, and ``Inference``. Common preparation kinds
include ``MissingPeriodicBox``, ``InvalidTolerance``, ``InvalidMinimumDistance``,
``NoValidAssignment``, ``SearchLimitExceeded``,
``TotalChargeWithMultipleComponents``, and ``DynamicsPresent``.
Connectivity, bond-order perception, and hydrogen addition on System also
report the typed BondPerceptionError.

For missing orders, ``details["atoms"]`` is the endpoint pair. For an open
selection, ``details["atom"]`` and ``details["neighbor"]`` are global indices.
An unsupported charge-model element reports ``atomic_number`` and ``model``.
Invalid numeric settings report ``value``. Search and assignment failures
report the relevant ``atom``. Nested engine failures report a ``reason`` string.
``UntypedAtom`` reports ``ff``, ``local`` (index within the typing scope), and
``atomic_number``. Errors without payload use an empty details dictionary.
