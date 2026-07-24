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

- **`Topology`** (`topology.rs`) — bonds, molecules, and `atoms: AtomStorage`; usually read once from file
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
- **`State`** (`state.rs`) — coordinates (`Vec<Pos>`), optional velocities/forces, timestamp, optional `PeriodicBox`
- **`System`** (`selection/system.rs`) — owns `Topology + State`; the primary user-facing container

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
