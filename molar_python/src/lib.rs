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
    system::{SysAtomIterator, SysPosIterator},
};
//-------------------------------------------

#[pyclass(unsendable, name = "SasaResults")]
struct SasaResultsPy(SasaResults);

#[pymethods]
impl SasaResultsPy {
    #[getter]
    fn areas(&self) -> &[f32] {
        self.0.areas()
    }

    #[getter]
    fn volumes(&self) -> &[f32] {
        self.0.volumes()
    }

    #[getter]
    fn total_area(&self) -> f32 {
        self.0.total_area()
    }

    #[getter]
    fn total_volume(&self) -> f32 {
        self.0.total_volume()
    }
}

#[pyclass]
struct IsometryTransform(nalgebra::IsometryMatrix3<f32>);

// Free functions

#[pyfunction(name = "fit_transform")]
fn fit_transform_py(sel1: &SelPy, sel2: &SelPy) -> PyResult<IsometryTransform> {
    let tr = molar::prelude::fit_transform(sel1, sel2).map_err(to_py_runtime_err)?;
    Ok(IsometryTransform(tr))
}

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

#[pyfunction]
fn rmsd_py(sel1: &SelPy, sel2: &SelPy) -> PyResult<f32> {
    Ok(rmsd(sel1, sel2).map_err(to_py_runtime_err)?)
}

#[pyfunction(name = "rmsd_mw")]
fn rmsd_mw_py(sel1: &SelPy, sel2: &SelPy) -> PyResult<f32> {
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
#[pyo3(signature = (cutoff,data1,data2=None,dims=[false,false,false]))]
fn distance_search<'py>(
    py: Python<'py>,
    cutoff: &Bound<'py, PyAny>,
    data1: &Bound<'py, SelPy>,
    data2: Option<&Bound<'py, SelPy>>,
    dims: [bool; 3],
) -> PyResult<Bound<'py, PyAny>> {
    let mut res: Vec<(usize, usize, f32)>;
    let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
    let sel1 = data1.borrow();

    if let Ok(d) = cutoff.extract::<f32>() {
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
        let vdw1: Vec<f32> = sel1.iter_atoms().map(|a| a.vdw()).collect();

        if sel1.len() != vdw1.len() {
            return Err(PyValueError::new_err(format!(
                "Size mismatch 1: {} {}",
                sel1.len(),
                vdw1.len()
            )));
        }

        if let Some(d2) = data2 {
            let sel2 = d2.borrow();
            let vdw2: Vec<f32> = sel2.iter_atoms().map(|a| a.vdw()).collect();

            if sel2.len() != vdw2.len() {
                return Err(PyValueError::new_err(format!(
                    "Size mismatch 2: {} {}",
                    sel2.len(),
                    vdw2.len()
                )));
            }

            if pbc_dims.any() {
                res = distance_search_double_vdw(&sel1 as &SelPy, &sel2 as &SelPy, &vdw1, &vdw2);
            } else {
                res = distance_search_double_vdw_pbc(
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    &vdw1,
                    &vdw2,
                    &sel1.require_box().unwrap(),
                    pbc_dims,
                );
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
        unreachable!()
    };

    // Subdivide the result into two arrays
    unsafe {
        // Pairs array
        let pairs_arr = numpy::PyArray2::<usize>::new(py, [res.len(), 2], true);
        for i in 0..res.len() {
            pairs_arr.uget_raw([i, 0]).write(res[i].0);
            pairs_arr.uget_raw([i, 1]).write(res[i].1);
        }

        // Distances array
        let dist_arr = numpy::PyArray1::<f32>::new(py, [res.len()], true);
        for i in 0..res.len() {
            dist_arr.uget_raw(i).write(res[i].2);
        }

        Ok((pairs_arr, dist_arr).into_bound_py_any(py)?)
    }
}

#[pyclass(name = "NdxFile")]
struct NdxFilePy(NdxFile);

#[pymethods]
impl NdxFilePy {
    #[new]
    fn new(fname: &str) -> PyResult<Self> {
        Ok(NdxFilePy(NdxFile::new(fname).map_err(to_py_value_err)?))
    }

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
    m.add_class::<SasaResultsPy>()?;
    m.add_class::<NdxFilePy>()?;
    m.add_class::<SysPosIterator>()?;
    m.add_class::<SysAtomIterator>()?;
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
