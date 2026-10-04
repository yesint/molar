# pymolar — Python bindings for MolAR

[MolAR](https://github.com/yesint/molar) is a Rust library for molecular
modeling and analysis of MD trajectories. `pymolar` exposes the bindings via
[PyO3](https://pyo3.rs/), with NumPy interop for coordinate access.

## Installation

```sh
pip install pymolar          # default: single precision (float32)
pip install pymolar-f64      # double precision (float64)
```

The two wheels are co-installable in the same virtual environment because
they ship distinct top-level packages and native modules:

```python
import pymolar           # float32 build — `pymolar.molar`
import pymolar_f64       # float64 build — `pymolar_f64.molar`
```

Both expose the same API. The only differences are the dtype of NumPy arrays
returned by the binding (`float32` vs `float64`) and the precision of all
internal computations.

## Building from source

You need a stable Rust toolchain and [maturin](https://www.maturin.rs/).

```sh
pip install maturin
```

The Rust crate (`molar_python`) is shared between the two wheels. Each wheel
has its own project directory containing a `pyproject.toml`. Build by running
`maturin` from the appropriate directory:

```sh
# Single-precision wheel (default)
cd molar_python
maturin build --release

# Double-precision wheel
cd molar_python/pymolar-f64-pkg
maturin build --release
```

The wheels are written to `target/wheels/` at the workspace root.

For local development, replace `build` with `develop` to install the wheel
into the active virtualenv in editable form.

### TPR support (optional)

Reading Gromacs `.tpr` and `.cpt` files requires the runtime plugin
`libmolar_gromacs_plugin.so` to be available. The plugin is built
automatically by `cargo`/`maturin` when the environment variables
`GROMACS_SOURCE_DIR`, `GROMACS_BUILD_DIR` and `GROMACS_LIB_DIR` are set
(see `../config.toml.template`). Without the plugin, `.tpr` and `.cpt`
reading is unavailable. NetCDF also requires its build feature.

### NetCDF support (optional)

AMBER `.nc`/`.ncdf` trajectories are supported when the wheel is built with
the `netcdf` feature:

```sh
maturin build --release --features molar/netcdf            # f32 wheel
cd pymolar-f64-pkg && maturin build --release --features molar/netcdf,f64
```

## Documentation

<https://yesint.github.io/molar/>

### Start here

For agents and new users, read [llms.txt](llms.txt), then use these guides:

- [Agent guide](docs/agent_guide.rst): capability map, units, shared data, and interface limits.
- [Selections](docs/selections.rst): query syntax, index rules, and selection operations.
- [Workflows](docs/workflows.rst): complete examples for trajectories, fitting, contacts, secondary structure, chemistry, IO, and CLI tasks.
- [Python API types](python/pymolar/molar.pyi): names and signatures for static inspection.
- [Coverage audit](docs/coverage-audit.md): gaps found, corrections, and remaining implementation limits.

Coordinates and distances use **nm**; time uses **ps**. Selection coordinate
arrays have shape `(3, n_atoms)` and are copies. Assign a Fortran-contiguous
array with the package's dtype to write coordinates back. Use `rmsd_py()` for
RMSD and `replace_state_deep()` to update existing trajectory selections.
The Python API does not expose every capability of the Rust library.

### Build and check documentation

Install a wheel built from this checkout before generating the reference.
The generator reads docstrings from the installed extension and combines them
with the checked-in guides. Run from the workspace root:

```sh
python -m pip install sphinx
python molar_python/scripts/check_docs.py
python molar_python/scripts/generate_sphinx_docs.py --skip-install --strict
```

The HTML output is `target/pymolar-docs/html/index.html`. The build also writes
`llms.txt` and downloadable type files. Use `--no-build` to generate only the
Sphinx source. To check the double-precision package after installing its wheel:

```sh
python molar_python/scripts/check_docs.py --module pymolar_f64
python molar_python/scripts/generate_sphinx_docs.py --module pymolar_f64 --skip-install --strict \
  --source-dir target/pymolar-f64-docs/sphinx --build-dir target/pymolar-f64-docs/html
```
