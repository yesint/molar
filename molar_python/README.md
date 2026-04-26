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

Reading Gromacs `.tpr` files requires the runtime plugin
`libmolar_gromacs_plugin.so` to be available. The plugin is built
automatically by `cargo`/`maturin` when the environment variables
`GROMACS_SOURCE_DIR`, `GROMACS_BUILD_DIR` and `GROMACS_LIB_DIR` are set
(see `.cargo/config.toml.template`). Without these variables, all formats
except `.tpr` work normally.

### NetCDF support (optional)

AMBER `.nc`/`.ncdf` trajectories are supported when the wheel is built with
the `netcdf` feature:

```sh
maturin build --release --features molar/netcdf            # f32 wheel
cd pymolar-f64-pkg && maturin build --release --features molar/netcdf,f64
```

## Documentation

<https://yesint.github.io/molar/>
