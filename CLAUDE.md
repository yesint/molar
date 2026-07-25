# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Code analysis

- Always prefer Rust LSP to grep when exploring the code.
- For “go to definition”, use LSP goToDefinition, not Grep.
- For “find references”, use LSP findReferences, not Grep.
- For type information or docs, use LSP hover.
- Use Grep/Glob only for discovery:
  - finding files
  - searching plain-text patterns
  - locating candidate symbols before LSP
- After identifying the relevant Rust file/symbol, switch to LSP for navigation and understanding.
- Do not use Grep as a substitute for semantic reference search in Rust code when LSP is available.

## Permissions

- Do not ask permissin to run shell commands "with consecutive quote characters at word boundaries or word start"; to run compound commands (piped, chained with `&&`/`;`, or using subshells); any commands involving `git`.

## Commands

```sh
# Build all workspace crates
cargo build

# Build in release mode
cargo build -r

# Run all tests
cargo test

# Run tests for a specific crate
cargo test -p molar
cargo test -p molar_membrane

# Run a single test by name
cargo test -p molar <test_name>

# Run tests with output shown
cargo test -p molar -- --nocapture

# Check compilation without building
cargo check

# Build documentation
cargo doc --open

# Build Python bindings (from molar_python directory)
cd molar_python && maturin build -r && python -m pip install .
```

## Architecture

MolAR is a Cargo workspace (**Rust edition 2024**, MSRV 1.85) with these crates:

| Crate | Purpose |
|---|---|
| `molar` | Core library: SoA atom storage, selections, IO, topology, analysis tasks |
| `molar_gromacs` | Gromacs TPR support via a runtime-`dlopen`ed plugin (built only when Gromacs env vars are set; no compile-time dependency on Gromacs) |
| `molar_ff` | Force-field atom typing (GAFF/GAFF2) and partial charges (espaloma) |
| `molar_membrane` | Lipid membrane analysis (lipid order, curvature, etc.) |
| `molar_bin` | CLI utility (`last`, `rearrange`, `solvate`, `tip3to4` commands) |
| `molar_python` | Python bindings via PyO3/maturin (wheel: `pymolar`) |

All file formats are **pure Rust** (the former `molar_molfile` VMD-plugin crate has been removed).
PowerSASA is an external git dependency, not a workspace crate.

### Core data model (`molar/src/`)

- **`Topology`** (`topology.rs`) — molecules, `atoms: AtomStorage`, and `bonds: BondStorage`; usually read once from file
- **`AtomStorage`** (`atom_storage.rs`) — **Struct-of-Arrays** atom storage: one column per property.
  Ten always-present *core* columns (`name`, `resname`, `resid`, `resindex`, `atomic_number`, `mass`,
  `charge`, `chain`, `bfactor`, `occupancy`) plus four *optional* force-field/chemistry columns
  (`type_name`, `type_id`, `formal_charge`, `flags`) stored as `Option<Vec<T>>` — a `None` column costs
  nothing; when present it is full-length. Atoms are accessed through the borrowed **proxies**
  `AtomRef` / `AtomRefMut` (a two-word `{storage, index}` handle) — there is **no `&Atom`** to borrow.
- **`Atom`** (`atom.rs`) — the owned, densely-packed atom *row*, retained as the detached
  construction/interchange type (builders, IO readers, `From<&AtomLike>`); `AtomStorage::push_row`
  scatters it into the columns. `AtomFlags` holds the ring/aromatic bits (no longer packed into `type_id`).
- **`AtomLike`** (read getters) / **`AtomLikeMut`** (setters) — the atom interface, implemented by
  `Atom`, `AtomRef`, and `AtomRefMut`. Getters for the four *optional* properties return `Option`
  (e.g. `get_type_name() -> Option<&str>`); `charge` is the partial/working charge, `formal_charge`
  is the integer formal charge (kept separate).
- **`BondStorage`** (`bond_storage.rs`) — **Struct-of-Arrays** bond storage, same discipline as
  `AtomStorage`: an always-present pair column (`u32` internally, `usize` at every API boundary —
  caps a system at 4·10⁹ atoms) plus an *optional* `Option<Vec<BondOrder>>` order column, absent for
  connectivity-only sources (PDB CONECT / GRO / TPR) so MD systems allocate nothing for it. Bonds are
  read through the borrowed **`BondRef`** proxy — there is **no `&Bond`** to borrow. The owned
  `Bond` row (`bond.rs`) is the detached construction type; `BondStorage::push` scatters it.
- **`BondAdjacency`** (`bond_storage.rs`) — the per-atom bonded-neighbor index (compressed rows),
  cached inside `BondStorage`. `get_adjacency()` is cheap and parallel-safe; `ensure_adjacency(n_atoms)`
  builds it (a plain `&mut` field, **not** a `OnceCell` — interior mutability would cost `Topology` its
  `Sync`-through-`&` sharing across rayon). Structural change invalidates it, but **`set_order` does
  not** — that asymmetry is the point of splitting the columns. Anything changing the *atom count*
  must call `invalidate_adjacency()`, since `offsets` is sized `n_atoms + 1`.
  Also usable standalone via `BondAdjacency::build(n, pairs)`, which is how `molar_ff` indexes a
  remapped local subgraph. **Neighbor order within an atom's run is a guaranteed invariant**
  (ascending bond index) — the GAFF port indexes neighbors positionally and truncates to the first
  4 or 6. Graph routines (`sssr_rings`, `perception::perceive`, all of `gaff`) take a prebuilt
  adjacency; none builds a throwaway one.
- **`State`** (`state.rs`) — coordinates (`Vec<Pos>`), optional velocities/forces, timestamp, optional `PeriodicBox`
- **`System`** (`selection/system.rs`) — owns `Topology + State`; the primary user-facing container

### Adding a new atom or bond property

**Both storages are designed to be extended** — that is the main reason they are columnar. Adding an
*optional* column costs users who don't set it **nothing** (a `None` column is zero allocation), so
chemistry/force-field properties can be added without taxing MD workloads. Prefer an optional column
unless every atom or bond genuinely always has the property.

Reuse the existing generic helpers; do not hand-roll the materialization logic.

**New optional atom column** `foo: Option<Vec<T>>` (`T: Copy`), mirroring `type_id` /
`formal_charge` / `flags`:

1. `atom.rs` — add `pub foo: Option<T>` to the owned `Atom` row, plus a `with_foo()` builder.
2. `atom_storage.rs` — add the `foo: Option<Vec<T>>` field, then wire it into **every** place the
   other optional columns appear: `push_row` (via `push_opt`), `set_row` (via `set_opt`),
   `remove_by_index` (via `retain_mask`), `reserve`, an `ensure_foo()` materializer, and the
   `invariant_holds()` length check. Missing one of these is how columns silently desynchronize —
   `invariant_holds()` is `debug_assert`ed after every mutation and will catch it in tests.
3. `AtomLike` — `get_foo() -> Option<T>` (optional getters return `Option`); `AtomLikeMut` —
   `set_foo()`. Implement for `Atom`, `AtomRef`, `AtomRefMut`.
4. Extend the `to_atom` / `From<&AtomLike>` round-trip so the property survives a detach/reattach.

An always-present *core* column is the same minus the `Option` plumbing; add a slice accessor
(`foos()`) only if the selection hot path in `ast.rs` actually scans it — that is the sole reason
those accessors exist.

⚠ **Optional-column setters materialize the column and are therefore serial-only.** Materialize
before entering any parallel region (see *Parallel operations* below).

**New optional bond column** `bar: Option<Vec<T>>`, mirroring the order column:

1. `bond.rs` — add the field to the owned `Bond` row.
2. `bond_storage.rs` — add `bar: Option<Vec<T>>`, then: `push` (mirror the `push_order` helper —
   materialize only when the incoming value is *informative*, so a default-valued bond leaves the
   column absent), a `set_bar` (mirror `set_order`), the compact-and-truncate step in
   `remove_by_index`, the manual `Clone` impl, and the `invariant_holds()` length check.
3. `BondRef::bar()` — return the property's default when the column is absent, exactly as
   `order()` yields `Unspecified`. Add `has_bar()` only if a caller must distinguish absent from
   default. Extend `From<&BondRef> for Bond`.

⚠ A new bond column that is **not** connectivity (stereo, wedge direction, rotatable flag, …) must
**not** invalidate the cached `BondAdjacency` — follow `set_order`, not `push`. Only changes to the
pair column or the atom count invalidate it.

Good candidates already identified: bond stereo / E-Z configuration and wedge direction (the SDF
reader currently parses and discards the stereo field), ring-membership and rotatable flags.

### Selection system (`molar/src/selection/`)

The key design: a `Sel` is just a sorted `SVec` of atom indices — it is **detached** from any `System`. To do work, it must be bound:

- `sys.bind(&sel)` → `SelBound<'_>` (read-only, borrows system)
- `sys.bind_mut(&sel)` → `SelBoundMut<'_>` (read-write, mutably borrows system)
- `sys.select_bound("...")` → `SelOwnBound<'_>` (creates and binds in one step)
- `sys.select_bound_mut("...")` → `SelOwnBoundMut<'_>`

Borrow checking is enforced at compile time — you cannot hold a mutable and immutable bound selection simultaneously.

**Empty selections are forbidden** — selection methods return `Err` instead of an empty selection.

Traits that provide behavior live in `selection/traits.rs` and `providers.rs`:
- `AtomProvider` / `AtomMutProvider` (`providers.rs`) — the atom-access layer. The single required
  method is `atom_storage()` / `atom_storage_mut()`; the defaults (`iter_atoms`, `get_atom`,
  `iter_particle`, `par_iter_atoms`, …) hand out `AtomRef` / `AtomRefMut` proxies rather than `&Atom`.
  `PosProvider`/`VelProvider`/… mirror this for coordinates/velocities/forces.
- `Selectable` / `SelectableBound` — creating sub-selections; `split` / `split_par` (blanket-impl
  methods) and `iter_particle*` for particle access.
- Selection keyword evaluation (`ast.rs`) scans **one projected column** at a time (via
  `AtomStorage::names()`/`resids()`/… slice accessors) rather than materializing a proxy per atom.

### IO (`molar/src/io/`)

`FileHandler` dispatches by file extension to format-specific handlers — all **pure Rust**, no
external plugin:
- `.pdb`/`.ent`, `.xyz`, `.gro` → streaming structure/trajectory handlers
- `.dcd`, `.xtc` → binary trajectories (random access / seek supported)
- `.itp` → topology only; `.sdf`/`.mol` → MDL molfile (bonds with order, formal charge)
- `.nc`/`.ncdf` → AMBER NetCDF (optional `netcdf` feature)
- `.tpr` → TPR handler, which `dlopen`s the Gromacs plugin at runtime (see Optional Gromacs linking)

`FileHandler` implements `IntoIterator` yielding `State` for trajectory iteration. IO runs in a background thread with a channel buffer of 10 frames.

### Analysis task framework (`molar/src/analysis_task.rs`)

Implement `AnalysisTask<UserArgs>` trait with three methods:
- `new(context)` — called on first valid frame; create selections here
- `process_frame(context)` — called per frame
- `post_process(context)` — called after all frames

Standard CLI args (`-f files -b begin -e end --log --skip`) are handled automatically by `TrajAnalysisArgs`. Invoke with `TaskType::run()`.

### Parallel operations

Two parallel patterns:
1. **`par_iter_pos()` / `par_iter_atoms()`** — rayon parallel iteration within a single selection
2. **`split_par(closure)`** — produces non-overlapping `ParSplit`; iterate with `sys.iter_par_split_mut(&par)` for parallel processing of distinct fragments (e.g., per-molecule unwrapping)

The mutable parallel path hands out `AtomRefMut` proxies over **disjoint** atom indices. Core-column
setters read the column `Vec` header via a *shared* borrow and write the heap buffer through a raw
element pointer (never forming `&mut [T]`), which is race-free across threads — verified with Miri
under Tree Borrows (`par_atom_column_write_scoped`). **Optional-column setters materialize the column
and are serial-only**: materialize optional columns before entering a parallel region. When adding
`unsafe` to the atom layer, run:
```sh
MIRIFLAGS="-Zmiri-tree-borrows" cargo +nightly miri test -p molar --lib par_atom_column_write_scoped
```

### Optional Gromacs linking

To enable TPR reading, create `.cargo/config.toml` (use `config.toml.template` as a starting point) with:
```toml
[env]
GROMACS_SOURCE_DIR = "<path>"
GROMACS_BUILD_DIR = "<path>"
GROMACS_LIB_DIR = "<path>"
```
Without this, the crate compiles but TPR reading is unavailable.

### Coordinate units

All coordinates and distances use **nanometers** (matching Gromacs convention), not Angstroms.

### Test data

Integration test files live in `molar/tests/`: `protein.pdb`, `protein.xtc`, `membr.gro`, `topol.tpr`, etc.
