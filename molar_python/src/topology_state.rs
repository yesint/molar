use std::ffi::c_void;

use super::periodic_box::PeriodicBoxPy;
use anyhow::{anyhow, bail};
use molar::prelude::*;
use numpy::{
    ndarray::Axis, npyffi, Element, PyArray2, PyArrayDescrMethods, PyArrayMethods,
    PyUntypedArrayMethods, PY_ARRAY_API,
};
use pyo3::{prelude::*, types::PyCapsule};

#[pyclass(name = "State")]
pub(crate) struct StatePy {
    pub(crate) coords: Py<PyArray2<f32>>, // shape (3, n)

    #[pyo3(get, set)]
    pub time: f32,

    pub pbox: Option<PeriodicBox>,
}

#[pymethods]
impl StatePy {
    fn __len__(&self) -> usize {
        Python::attach(|py| self.coords.bind(py).shape()[1])
    }

    #[getter]
    fn get_box(&self) -> anyhow::Result<PeriodicBoxPy> {
        Ok(PeriodicBoxPy(
            self.pbox
                .as_ref()
                .ok_or_else(|| anyhow!("No periodic box"))?
                .clone(),
        ))
    }

    #[setter]
    fn set_box(&mut self, val: Bound<'_, PeriodicBoxPy>) -> anyhow::Result<()> {
        let b = self
            .pbox
            .as_mut()
            .ok_or_else(|| anyhow!("No periodic box"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }
}

impl StatePy {
    pub(crate) fn from_state(st: State) -> anyhow::Result<Self> {
        Python::attach(|py| {
            let n = st.len();
            let data_ptr = st.coords.as_ptr() as *mut f32;
            // Capsule owns the Vec and will drop it when Python frees the capsule
            let capsule =
                PyCapsule::new_with_destructor(py, st.coords, None, |v: Vec<Pos>, _ctx| drop(v))?;
            // dims: (3, n)
            let mut dims: [npyffi::npy_intp; 2] = [3, n as npyffi::npy_intp];

            // strides in BYTES:
            // axis0 step = 1 float, axis1 step = 3 floats  -> xyzxyz... for each column
            let item = std::mem::size_of::<f32>() as npyffi::npy_intp;
            let mut strides: [npyffi::npy_intp; 2] = [item, 3 * item];

            // dtype: PyArray_NewFromDescr steals a reference to descr
            let descr = f32::get_dtype(py).into_dtype_ptr();

            let flags = npyffi::NPY_ARRAY_WRITEABLE | npyffi::NPY_ARRAY_F_CONTIGUOUS;

            let arr_ptr = unsafe {
                PY_ARRAY_API.PyArray_NewFromDescr(
                    py,
                    PY_ARRAY_API.get_type_object(py, npyffi::NpyTypes::PyArray_Type),
                    descr,
                    2,
                    dims.as_mut_ptr(),
                    strides.as_mut_ptr(),
                    data_ptr as *mut c_void,
                    flags,
                    std::ptr::null_mut(),
                )
            };

            // Wrap the returned owned reference
            let arr_any = unsafe { Bound::from_owned_ptr_or_err(py, arr_ptr)? };
            let arr = arr_any.cast_into::<PyArray2<f32>>().unwrap();

            // Link lifetime: ndarray.base = capsule (SetBaseObject steals a reference)
            let cap_ptr = capsule.into_ptr(); // owned ref
            unsafe {
                if PY_ARRAY_API.PyArray_SetBaseObject(py, arr_ptr.cast(), cap_ptr) != 0 {
                    // on failure, we must decref cap_ptr because SetBaseObject did not steal it
                    pyo3::ffi::Py_DECREF(cap_ptr);
                    return Err(PyErr::fetch(py))?;
                }
            }

            Ok(StatePy {
                coords: arr.unbind(),
                time: st.time,
                pbox: st.pbox,
            })
        })
    }

    /// Remove columns (atoms) from `coords` by index and reallocate a fresh NumPy array.
    /// `removed` yields 0-based column indices in `[0, n)`.
    pub fn remove_coords(
        &mut self,
        py: Python<'_>,
        removed: impl Iterator<Item = usize>,
    ) -> anyhow::Result<()> {
        // Read old array
        let ro = self.coords.bind(py).readonly();
        let old = ro.as_array(); // ndarray::ArrayView2<f32>, shape (3, n)

        let n = old.shape()[1];

        // Build a keep mask
        let mut keep = vec![true; n];
        for idx in removed {
            if idx < n {
                keep[idx] = false;
            }
        }

        let kept_indices: Vec<usize> = (0..n).filter(|&i| keep[i]).collect();
        let new_n = kept_indices.len();

        if new_n == 0 {
            bail!("empty system remains after remove");
        }

        // Allocate new array
        let new_py = PyArray2::<f32>::zeros(py, [3, new_n], true);

        // Copy kept columns
        {
            let mut new_rw = new_py.readwrite();
            let mut new = new_rw.as_array_mut(); // ndarray::ArrayViewMut2<f32>

            for (new_j, &old_j) in kept_indices.iter().enumerate() {
                let src = old.column(old_j);
                let mut dst = new.column_mut(new_j);
                dst.assign(&src);
            }
        }

        // Swap in the new allocation
        self.coords = new_py.unbind();

        Ok(())
    }

    pub fn add_coords(&mut self, added: impl Iterator<Item = Pos>) {
        Python::attach(|py| {
            // Collect first so we know how many columns to add
            let added: Vec<Pos> = added.collect();
            if added.is_empty() {
                return;
            }

            // Old array
            let old_b = self.coords.bind(py);
            let old_ro = old_b.readonly();
            let old = old_ro.as_array(); // shape (3, n)

            let n_old = old.shape()[1];
            let n_add = added.len();
            let n_new = n_old + n_add;

            // Preserve memory order if you care
            let fortran = old_b.is_fortran_contiguous();

            // Allocate new array
            let new_b = PyArray2::<f32>::zeros(py, [3, n_new], fortran);

            {
                let mut new_rw = new_b.readwrite();
                let mut new = new_rw.as_array_mut(); // shape (3, n_new)

                // Copy old coords into the left block
                // new
                //     .index_axis_mut(Axis(1), 0) // just to ensure Axis is used; no-op line can be removed
                //     .len(); // no-op line can be removed

                // More efficient block copy:
                new.slice_mut(numpy::ndarray::s![.., 0..n_old]).assign(&old);

                // Fill appended columns
                for (k, p) in added.iter().enumerate() {
                    let j = n_old + k;
                    new[[0, j]] = p.x;
                    new[[1, j]] = p.y;
                    new[[2, j]] = p.z;
                }
            }

            self.coords = new_b.unbind();
        });
    }
}

impl SaveState for StatePy {}

impl LenProvider for StatePy {
    fn len(&self) -> usize {
        self.__len__()
    }
}

impl TimeProvider for StatePy {
    fn get_time(&self) -> f32 {
        self.time
    }
}

impl BoxProvider for StatePy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.pbox.as_ref()
    }
}

impl RandomPosProvider for StatePy {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        Python::attach(|py| {
            let base = self.coords.bind(py).data() as *const Pos;
            unsafe { &*base.add(i) }
        })
    }
}

#[pyclass(unsendable, name = "Topology")]
pub(crate) struct TopologyPy(pub(crate) Topology);
