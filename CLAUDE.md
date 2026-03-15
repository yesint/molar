# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Permissions

- Execute shell commands with consecutive quote characters at word boundaries or word start (e.g. `sed -i 's/'"'"'a/'"'"'b/'`) without asking for confirmation. Never prompt "Command contains consecutive quote characters at word start" — always proceed.
- Execute compound commands (piped, chained with `&&`/`;`, or using subshells) without asking for confirmation.
- Proceed autonomously when fixing compilation errors, resolving merge conflicts, or continuing an in-progress refactoring task — do not stop to ask for permission at each step.
- Fix compilation errors in dependencies (molar_bin, molar_python, molar_membrane) as part of refactoring without asking.

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

MolAR is a Cargo workspace with these crates:

| Crate | Purpose |
|---|---|
| `molar` | Core library: selections, IO, topology, analysis tasks |
| `molar_molfile` | VMD molfile plugin bindings (PDB, DCD, XYZ) |
| `molar_gromacs` | Gromacs bindings for TPR reading (optional, requires Gromacs source) |
| `molar_powersasa` | SASA computation via PowerSASA algorithm |
| `molar_membrane` | Lipid membrane analysis (lipid order, curvature, etc.) |
| `molar_bin` | CLI utility (`last`, `rearrange`, `solvate`, `tip3to4` commands) |
| `molar_python` | Python bindings via PyO3/maturin |

### Core data model (`molar/src/`)

- **`Topology`** (`topology.rs`) — atoms, bonds, molecules; usually read once from file
- **`State`** (`state.rs`) — coordinates (`Vec<Pos>`), timestamp, optional `PeriodicBox`
- **`System`** (`selection/system.rs`) — owns `Topology + State`; the primary user-facing container

### Selection system (`molar/src/selection/`)

The key design: a `Sel` is just a sorted `SVec` of atom indices — it is **detached** from any `System`. To do work, it must be bound:

- `sys.bind(&sel)` → `SelBound<'_>` (read-only, borrows system)
- `sys.bind_mut(&sel)` → `SelBoundMut<'_>` (read-write, mutably borrows system)
- `sys.select_bound("...")` → `SelOwnBound<'_>` (creates and binds in one step)
- `sys.select_bound_mut("...")` → `SelOwnBoundMut<'_>`

Borrow checking is enforced at compile time — you cannot hold a mutable and immutable bound selection simultaneously.

**Empty selections are forbidden** — selection methods return `Err` instead of an empty selection.

Traits that provide behavior live in `selection/traits.rs`:
- `AtomPosAnalysis` — read-only iteration, `split`, `split_par`, particle access
- `Selectable` / `SelectableBound` — creating sub-selections
- Various `*Provider` traits for typed access

### IO (`molar/src/io/`)

`FileHandler` dispatches by file extension to format-specific handlers:
- `.pdb`, `.dcd`, `.xyz` → VMD molfile plugin (C FFI via `molar_molfile`)
- `.xtc` → custom XTC handler (random access supported)
- `.gro` → GRO handler (single frame and multi-frame trajectory)
- `.itp` → ITP handler (topology only)
- `.tpr` → TPR handler (requires Gromacs; optional)

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
