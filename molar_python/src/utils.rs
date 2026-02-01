use molar::prelude::*;
use numpy::{
    nalgebra::{self, Const},
    npyffi, Element, PyArray1, PyArrayMethods, ToNpyDims, PY_ARRAY_API,
};
use pyo3::prelude::*;
use std::ffi::c_void;

use crate::topology_state::StatePy;

pub(crate) fn to_py_value_err<E: std::fmt::Display>(e: E) -> pyo3::PyErr {
    pyo3::exceptions::PyValueError::new_err(e.to_string())
}

pub(crate) fn to_py_io_err<E: std::fmt::Display>(e: E) -> pyo3::PyErr {
    pyo3::exceptions::PyIOError::new_err(e.to_string())
}

pub(crate) fn to_py_runtime_err<E: std::fmt::Display>(e: E) -> pyo3::PyErr {
    pyo3::exceptions::PyRuntimeError::new_err(e.to_string())
}

// Constructs PyArray backed by existing Pos data.
pub(crate) unsafe fn map_pyarray_to_pos<'py>(
    st: &Bound<'py,StatePy>,
    id: usize,
) -> Bound<'py, PyArray1<f32>> {
    use numpy::Element;
    use numpy::PyArrayDescrMethods;

    let py = st.py();
    let mut dims = numpy::ndarray::Dim(3);

    unsafe {
        let ptr = PY_ARRAY_API.PyArray_NewFromDescr(
            py,
            PY_ARRAY_API.get_type_object(py, npyffi::NpyTypes::PyArray_Type),
            f32::get_dtype(py).into_dtype_ptr(),
            dims.ndim_cint(),
            dims.as_dims_ptr(),
            std::ptr::null_mut(), // no strides
            st.get().inner_mut().coords.as_mut_ptr().add(id) as *mut c_void, // data
            npyffi::NPY_ARRAY_WRITEABLE
                | npyffi::NPY_ARRAY_ALIGNED
                | npyffi::NPY_ARRAY_C_CONTIGUOUS, // flag
            std::ptr::null_mut(), // obj
        );

        // IMPORTANT: SetBaseObject steals a reference to `base`.
        // So we must INCREF first and NOT INCREF after.
        pyo3::ffi::Py_IncRef(st.as_ptr());
        PY_ARRAY_API.PyArray_SetBaseObject(py, ptr.cast(), st.as_ptr());

        // Turn raw pointer into a Bound<PyArray1<f32>>
        let any = Bound::from_owned_ptr(py, ptr.cast::<pyo3::ffi::PyObject>());
        any.cast_into::<PyArray1<f32>>().unwrap()
    }
}

// Constructs read-only PyArray backed by existing Pos data.
pub(crate) fn map_const_pyarray_to_vec3<'py>(
    py: Python<'py>,
    data: &Vector3f,
    parent: &Bound<'py, PyAny>,
) -> *mut npyffi::PyArrayObject {
    use numpy::Element;
    use numpy::PyArrayDescrMethods;

    let mut dims = numpy::ndarray::Dim(3);

    unsafe {
        let ptr = PY_ARRAY_API.PyArray_NewFromDescr(
            py,
            PY_ARRAY_API.get_type_object(py, npyffi::NpyTypes::PyArray_Type),
            f32::get_dtype(py).into_dtype_ptr(),
            dims.ndim_cint(),
            dims.as_dims_ptr(),
            std::ptr::null_mut(),         // no strides
            data.as_ptr() as *mut c_void, // data
            0,                            // no writable flag
            std::ptr::null_mut(),         // obj
        );

        // The following mangling with the ref counting is deduced by
        // tries and errors and seems to work correctly and keeps the parnet object alive
        // until any of the referencing PyArray objects is alive.

        // We set the parent as a base object of the PyArray to link them together.
        PY_ARRAY_API.PyArray_SetBaseObject(py, ptr.cast(), parent.as_ptr());
        // Increase reference count of parent object since
        // our PyArray is now referencing it!
        pyo3::ffi::Py_IncRef(parent.as_ptr());
        ptr.cast()
    }
}

/// numpy built-in [to_pyarray] method always produces 2D arrays,while we need 1D
/// arrays when converting from vectors. This is the fork of original method which
/// does this.
///
/// Note that the NumPy array always has Fortran memory layout
/// matching the [memory layout][memory-layout] used by [`nalgebra`].
///
/// [memory-layout]: https://nalgebra.org/docs/faq/#what-is-the-memory-layout-of-matrices
pub(crate) fn clone_vec_to_pyarray1<'py, N, R, S>(
    v: &nalgebra::Matrix<N, R, Const<1>, S>,
    py: Python<'py>,
) -> Bound<'py, PyArray1<N>>
where
    N: nalgebra::Scalar + Element,
    R: nalgebra::Dim,
    S: nalgebra::Storage<N, R, Const<1>>,
{
    unsafe {
        let array = PyArray1::<N>::new(py, (v.nrows(),), true);
        let mut data_ptr = array.data();
        if v.data.is_contiguous() {
            std::ptr::copy_nonoverlapping(v.data.ptr(), data_ptr, v.len());
        } else {
            for item in v.iter() {
                data_ptr.write(item.clone_ref(py));
                data_ptr = data_ptr.add(1);
            }
        }
        array
    }
}
