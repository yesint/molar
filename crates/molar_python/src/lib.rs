use molar::prelude::*;
use numpy::{
    npyffi::{self},
    PyArray1, PyArrayMethods, ToNpyDims, PY_ARRAY_API,
};
use pyo3::prelude::*;
use std::ffi::c_void;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyclass]
struct Atom(molar::core::Atom);

#[pyclass]
struct Particle {
    //atom:&'static Atom,
    pos: &'static mut molar::core::Pos,

    // Index
    #[pyo3(get)]
    index: usize,
}

#[pymethods]
impl Particle {
    #[getter]
    fn get_pos<'py>(&mut self, py: Python<'py>) -> Bound<'py, PyAny> {
        map_pos_to_pyarray(py, self.pos)
    }

    //pub resname: AsciiString,
    //pub resid: i32, // Could be negative
    //pub resindex: usize,
    // Atom physical properties from topology
    //pub atomic_number: u8,
    //pub mass: f32,
    //pub charge: f32,
    //pub type_name: AsciiString,
    //pub type_id: u32,
    // Specific PDB fields
    //pub chain: AsciiChar,
    //pub bfactor: f32,
    //pub occupancy: f32,
}

#[pyclass]
struct Topology(Option<molar::core::Topology>);

#[pyclass]
struct State(Option<molar::core::State>);

#[pyclass(unsendable)]
struct FileHandler(molar::io::FileHandler);

#[pymethods]
impl FileHandler {
    #[staticmethod]
    fn open(fname: &str) -> PyResult<Self> {
        Ok(FileHandler(
            molar::io::FileHandler::open(fname).map_err(|e| anyhow::anyhow!(e))?,
        ))
    }

    fn read(&mut self) -> PyResult<(Topology, State)> {
        let (top, st) = self.0.read().map_err(|e| anyhow::anyhow!(e))?;
        Ok((Topology(Some(top)), State(Some(st))))
    }
}

#[pyclass(unsendable)]
struct Source(molar::core::Source<molar::core::MutableSerial>);

#[pymethods]
impl Source {
    #[new]
    fn new(topology: &mut Topology, state: &mut State) -> PyResult<Self> {
        Ok(Source(
            molar::core::Source::new_serial(
                topology.0.take().unwrap().into(),
                state.0.take().unwrap().into(),
            )
            .map_err(|e| anyhow::anyhow!(e))?,
        ))
    }

    #[staticmethod]
    fn from_file(fname: &str) -> PyResult<Self> {
        Ok(Source(
            molar::core::Source::serial_from_file(fname).map_err(|e| anyhow::anyhow!(e))?,
        ))
    }

    fn select_all(&mut self) -> PyResult<Sel> {
        Ok(Sel(self.0.select_all().map_err(|e| anyhow::anyhow!(e))?))
    }

    fn select_str(&mut self, sel_str: &str) -> PyResult<Sel> {
        Ok(Sel(self
            .0
            .select_str(sel_str)
            .map_err(|e| anyhow::anyhow!(e))?))
    }
}

//====================================
#[pyclass(unsendable)]
struct Sel(molar::core::Sel<MutableSerial>);

#[pymethods]
impl Sel {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn com<'a>(&self, _py: Python<'a>) -> PyResult<Bound<'a, numpy::PyArray1<f32>>> {
        Ok(copy_pos_to_pyarray(
            _py,
            &self.0.center_of_mass().map_err(|e| anyhow::anyhow!(e))?,
        )) //.into_py(_py))
    }

    fn nth_pos<'a>(&self, _py: Python<'a>, i: usize) -> PyResult<Bound<'a, PyAny>> {
        let pos = self.0.nth_pos_mut(i).ok_or_else(|| anyhow::anyhow!("Out of bounds"))?;
        Ok(map_pos_to_pyarray(_py, pos))
    }

    /*
    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf.0 = 0;
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Py<PyAny>> {
        let py = slf.py();
        let cur = slf.0;
        slf.0 += 1;
        if slf.cur < slf.len() {
            let (i,at,pos) = unsafe{slf.sel.nth_unchecked_mut(cur)};
            //Some((i,map_pos_to_pyarray(py,pos)).into_py(py))
            let p = Particle {
                index: i,
                atom: Py::new(py,Atom(at)).unwrap(),
                pos: map_pos_to_pyarray(py,pos).into_py(py),
            };
            Some(p.into_py(py))
        } else {
            None
        }
    }
    */
}

// Constructs PyArray backed by existing Pos data.
fn map_pos_to_pyarray<'py>(py: Python<'py>, data: &mut molar::core::Pos) -> Bound<'py, PyAny> {
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
            std::ptr::null_mut(),                    // no strides
            data.coords.as_mut_ptr() as *mut c_void, // data
            npyffi::NPY_ARRAY_WRITEABLE | npyffi::NPY_ARRAY_F_CONTIGUOUS, // flag
            std::ptr::null_mut(),                    // obj
        );

        //PyArray1::<f32>::from_borrowed_ptr(py, ptr)
        Bound::from_borrowed_ptr(py, ptr)
    }
}

// Constructs new PyArray that copies Pos
fn copy_pos_to_pyarray<'py>(py: Python<'py>, data: &molar::core::Pos) -> Bound<'py, PyArray1<f32>> {
    unsafe {
        let array = PyArray1::<f32>::new(py, 3, true);
        std::ptr::copy_nonoverlapping(data.coords.as_ptr(), array.data(), 3);
        array
    }
}

//====================================

/// A Python module implemented in Rust.
#[pymodule]
fn molar_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<Particle>()?;
    m.add_class::<Topology>()?;
    m.add_class::<State>()?;
    m.add_class::<FileHandler>()?;
    m.add_class::<Source>()?;
    m.add_class::<Sel>()?;
    Ok(())
}
