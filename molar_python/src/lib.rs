use std::sync::atomic::{AtomicUsize, Ordering};

use molar::prelude::*;
use numpy::{
    nalgebra::{self},
    PyArrayMethods,
};
use pyo3::{
    exceptions::{PyNotImplementedError, PyTypeError, PyValueError},
    prelude::*,
    IntoPyObjectExt,
};

mod utils;
use utils::*;

mod atom;
use atom::AtomPy;

mod particle;
use particle::ParticlePy;

mod periodic_box;
use periodic_box::PeriodicBoxPy;

mod file_handler;
use file_handler::{FileHandlerPy, FileStatsPy};

mod system;
use system::SystemPy;

mod selection;
use selection::SelPy;

// mod membrane;
// use membrane::*;

mod topology_state;
use topology_state::*;

use crate::{
    selection::{SelAtomIterator, SelPosIterator, TmpSel},
    system::{SysAtomIterator, SysParticleIterator, SysPosIterator},
};
//-------------------------------------------
/// Solvent-accessible surface area and volume measurements for a selection.
///
/// Both a result holder and a persistent calculator: call :meth:`update` on subsequent
/// frames to reuse the power diagram without full reconstruction.
///
/// **Example**
///
/// .. code-block:: python
///
///    sasa = sel.sasa()
///    print(sasa.total_area)    # total SASA in nm²
///    print(sasa.areas[:5])     # per-atom areas

#[pyclass(unsendable, name = "Sasa")]
struct SasaPy(Sasa);

#[pymethods]
impl SasaPy {
    /// Per-atom solvent-accessible areas (nm²).
    #[getter]
    fn areas(&self) -> &[Float] {
        self.0.areas()
    }

    /// Per-atom solvent-excluded volumes (only meaningful when built with ``sasa_vol``).
    #[getter]
    fn volumes(&self) -> &[Float] {
        self.0.volumes()
    }

    /// Total solvent-accessible area (nm²).
    #[getter]
    fn total_area(&self) -> Float {
        self.0.total_area()
    }

    /// Total solvent-excluded volume (only meaningful when built with ``sasa_vol``).
    #[getter]
    fn total_volume(&self) -> Float {
        self.0.total_volume()
    }

    fn __repr__(&self) -> String {
        format!(
            "Sasa(n={}, total_area={:.3}, total_volume={:.3})",
            self.0.areas().len(),
            self.0.total_area(),
            self.0.total_volume()
        )
    }
}

/// Rigid-body isometry (rotation + translation) returned by ``fit_transform``.
///
/// No public constructor; obtained from :func:`pymolar.fit_transform`.
#[pyclass]
struct IsometryTransform(nalgebra::IsometryMatrix3<Float>);

#[pymethods]
impl IsometryTransform {
    fn __repr__(&self) -> String {
        let t = &self.0.translation.vector;
        format!(
            "IsometryTransform(trans=[{:.3}, {:.3}, {:.3}])",
            t[0], t[1], t[2]
        )
    }
}

// Free functions
/// Compute rigid transform that best aligns ``sel1`` onto ``sel2``.
///
/// :param sel1: Mobile selection.
/// :param sel2: Reference selection.
/// :returns: Rigid transform.
/// :rtype: IsometryTransform
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar
///    tr = pymolar.fit_transform(mobile_sel, ref_sel)
///    mobile_sel.apply_transform(tr)

#[pyfunction(name = "fit_transform")]
fn fit_transform_py(sel1: &SelPy, sel2: &SelPy) -> PyResult<IsometryTransform> {
    let tr = molar::prelude::fit_transform(sel1, sel2).map_err(to_py_runtime_err)?;
    Ok(IsometryTransform(tr))
}
/// Align selections by matching atom names before fitting.
///
/// :param sel1: Mobile selection.
/// :param sel2: Reference selection.
/// :returns: Rigid transform computed on matched atoms.
/// :rtype: IsometryTransform
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar
///    tr = pymolar.fit_transform_matching(mobile_sel, ref_sel)
///    mobile_sel.apply_transform(tr)

#[pyfunction(name = "fit_transform_matching")]
fn fit_transform_matching_py(sel1: &SelPy, sel2: &SelPy) -> PyResult<IsometryTransform> {
    let (ind1, ind2) = get_matching_atoms_by_name(sel1, sel2);

    let sub1 = ind1
        .into_sel_index(sel1, Some(&sel1.index()))
        .map_err(to_py_runtime_err)?;
    let sub2 = ind2
        .into_sel_index(sel2, Some(&sel2.index()))
        .map_err(to_py_runtime_err)?;

    let sub_sel1 = TmpSel {
        top: sel1.r_top(),
        st: sel1.r_st(),
        index: &sub1,
    };

    let sub_sel2 = TmpSel {
        top: sel2.r_top(),
        st: sel2.r_st(),
        index: &sub2,
    };

    let tr = fit_transform(&sub_sel1, &sub_sel2).map_err(to_py_runtime_err)?;
    Ok(IsometryTransform(tr))
}
/// Compute RMSD between two selections.
///
/// :param sel1: First selection.
/// :param sel2: Second selection.
/// :returns: RMSD in nm.
/// :rtype: float
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar
///    r = pymolar.rmsd(sel1, sel2)

#[pyfunction]
fn rmsd_py(sel1: &SelPy, sel2: &SelPy) -> PyResult<Float> {
    Ok(rmsd(sel1, sel2).map_err(to_py_runtime_err)?)
}
/// Compute mass-weighted RMSD between two selections.
///
/// :param sel1: First selection.
/// :param sel2: Second selection.
/// :returns: Mass-weighted RMSD in nm.
/// :rtype: float
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar
///    rw = pymolar.rmsd_mw(sel1, sel2)

#[pyfunction(name = "rmsd_mw")]
fn rmsd_mw_py(sel1: &SelPy, sel2: &SelPy) -> PyResult<Float> {
    Ok(rmsd_mw(sel1, sel2).map_err(to_py_runtime_err)?)
}

#[pyclass(frozen)]
struct ParticleIterator {
    sel: Py<SelPy>,
    cur: AtomicUsize,
}

#[pymethods]
impl ParticleIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(slf: &Bound<'_, Self>) -> Option<ParticlePy> {
        let ret = SelPy::__getitem__(slf.get().sel.get(), slf.get().cur.load(Ordering::Relaxed) as isize).ok();
        slf.get().cur.fetch_add(1, Ordering::Relaxed);
        ret
    }
}

#[pyfunction]
#[pyo3(signature = (cutoff,data1,data2=None,dims=None))]
#[pyo3(text_signature = "(cutoff, data1, data2=None, dims=None)")]
/// Find atom pairs within cutoff distance between one or two selections.
///
/// :param cutoff: Distance cutoff in nm, or ``"vdw"`` for van der Waals radii sum.
/// :param data1: First selection.
/// :param data2: Second selection (optional; self-search if omitted).
/// :param dims: Periodic dimensions ``[x, y, z]`` booleans.
/// :returns: Tuple ``(pairs, distances)`` — ``pairs`` is ``[N, 2]`` index array,
///     ``distances`` is length-N float array.
/// :rtype: tuple[numpy.ndarray, numpy.ndarray]
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar
///    pairs, dists = pymolar.distance_search(0.35, sel1, sel2)
///    pairs, dists = pymolar.distance_search(0.35, sel1, dims=[True, True, True])
fn distance_search<'py>(
    py: Python<'py>,
    cutoff: &Bound<'py, PyAny>,
    data1: &Bound<'py, SelPy>,
    data2: Option<&Bound<'py, SelPy>>,
    dims: Option<[bool; 3]>,
) -> PyResult<Bound<'py, PyAny>> {
    let mut res: Vec<(usize, usize, Float)>;
    let dims = dims.unwrap_or([false, false, false]);
    let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
    let sel1 = data1.borrow();

    if let Ok(d) = cutoff.extract::<Float>() {
        // Distance cutoff
        if let Some(d2) = data2 {
            let sel2 = d2.borrow();
            if pbc_dims.any() {
                res = distance_search_double_pbc(
                    d,
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    sel1.iter_index(),
                    sel2.iter_index(),
                    &sel1.require_box().unwrap(),
                    pbc_dims,
                );
            } else {
                res = distance_search_double(
                    d,
                    &sel1 as &SelPy,
                    &sel2 as &SelPy,
                    sel1.iter_index(),
                    sel2.iter_index(),
                );
            }
        } else {
            if pbc_dims.any() {
                res = distance_search_single_pbc(
                    d,
                    sel1.iter_pos(),
                    sel1.iter_index(),
                    &sel1.require_box().unwrap(),
                    pbc_dims,
                );
            } else {
                res = distance_search_single(d, &sel1 as &SelPy, sel1.iter_index());
            }
        }
    } else if let Ok(s) = cutoff.extract::<String>() {
        if s != "vdw" {
            return Err(PyTypeError::new_err(format!("Unknown cutoff type {s}")));
        }

        // VdW cutof
        let vdw1: Vec<Float> = sel1.iter_atoms().map(|a| a.vdw()).collect();

        if sel1.len() != vdw1.len() {
            return Err(PyValueError::new_err(format!(
                "Size mismatch 1: {} {}",
                sel1.len(),
                vdw1.len()
            )));
        }

        if let Some(d2) = data2 {
            let sel2 = d2.borrow();
            let vdw2: Vec<Float> = sel2.iter_atoms().map(|a| a.vdw()).collect();

            if sel2.len() != vdw2.len() {
                return Err(PyValueError::new_err(format!(
                    "Size mismatch 2: {} {}",
                    sel2.len(),
                    vdw2.len()
                )));
            }

            if pbc_dims.any() {
                res = distance_search_double_vdw_pbc(
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    &vdw1,
                    &vdw2,
                    &sel1.require_box().unwrap(),
                    pbc_dims,
                );
            } else {
                res = distance_search_double_vdw(&sel1 as &SelPy, &sel2 as &SelPy, &vdw1, &vdw2);
            }

            // Convert local indices to global
            unsafe {
                for el in &mut res {
                    el.0 = sel1.get_index_unchecked(el.0);
                    el.1 = sel2.get_index_unchecked(el.1);
                }
            }
        } else {
            return Err(PyNotImplementedError::new_err(
                "VdW distance search is not yet supported for single selection",
            ));
        }
    } else {
        return Err(PyTypeError::new_err("cutoff must be a float or 'vdw'"));
    };

    // Subdivide the result into two arrays
    let n = res.len();
    let mut flat_pairs = Vec::with_capacity(n * 2);
    let mut dists = Vec::with_capacity(n);
    for (i, j, d) in res {
        flat_pairs.push(i);
        flat_pairs.push(j);
        dists.push(d);
    }
    let pairs_arr = numpy::PyArray1::from_vec(py, flat_pairs).reshape([n, 2])?;
    let dist_arr = numpy::PyArray1::from_vec(py, dists);
    Ok((pairs_arr, dist_arr).into_bound_py_any(py)?)
}
/// Reader for GROMACS NDX index files.
///
/// **Example**
///
/// .. code-block:: python
///
///    ndx = pymolar.NdxFile("index.ndx")
///    sel = ndx.get_group_as_sel("Protein", system)

#[pyclass(name = "NdxFile")]
struct NdxFilePy(NdxFile);

#[pymethods]
impl NdxFilePy {
    #[new]
    /// Open an NDX file from disk.
    ///
    /// :param fname: Path to ``.ndx`` file.
    /// :returns: Parsed NDX file object.
    /// :rtype: NdxFile
    fn new(fname: &str) -> PyResult<Self> {
        Ok(NdxFilePy(NdxFile::new(fname).map_err(to_py_value_err)?))
    }

    /// Build a selection from named NDX group.
    ///
    /// :param gr_name: Group name in the NDX file.
    /// :param sys: System providing topology/state context.
    /// :returns: Selection with indices from the group.
    /// :rtype: Sel
    fn get_group_as_sel(&self, gr_name: &str, sys: &SystemPy) -> PyResult<SelPy> {
        Python::attach(|py| {
            Ok(SelPy::new(
                sys.py_top().clone_ref(py),
                sys.py_st().clone_ref(py),
                self.0
                    .get_group(gr_name)
                    .map_err(to_py_value_err)?
                    .to_owned(),
            ))
        })
    }
}

//====================================
/// Print library greeting and version banner.
#[pyfunction]
fn greeting() {
    molar::greeting("molar_python");
}

/// A Python module implemented in Rust.
#[pymodule(name = "molar")]
//#[pymodule]
fn molar_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_class::<AtomPy>()?;
    m.add_class::<ParticlePy>()?;
    m.add_class::<TopologyPy>()?;
    m.add_class::<StatePy>()?;
    m.add_class::<PeriodicBoxPy>()?;
    m.add_class::<FileHandlerPy>()?;
    m.add_class::<FileStatsPy>()?;
    m.add_class::<SystemPy>()?;
    m.add_class::<SelPy>()?;
    m.add_class::<SasaPy>()?;
    m.add_class::<NdxFilePy>()?;
    m.add_class::<SysPosIterator>()?;
    m.add_class::<SysAtomIterator>()?;
    m.add_class::<SysParticleIterator>()?;
    m.add_class::<SelPosIterator>()?;
    m.add_class::<SelAtomIterator>()?;
    m.add_function(wrap_pyfunction!(greeting, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_matching_py, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_py, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_mw_py, m)?)?;
    m.add_function(wrap_pyfunction!(distance_search, m)?)?;
    Ok(())
}
