use molar::prelude::*;
use numpy::{
    nalgebra, npyffi::{self}, PyArray1, PyArrayMethods, PyReadonlyArray1, ToNpyDims, PY_ARRAY_API
};
use pyo3::{prelude::*, IntoPyObjectExt};
use std::{borrow::Borrow, ffi::c_void};

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyclass]
struct Atom(molar::core::Atom);

#[pyclass(unsendable)]
struct Particle {
    atom: &'static mut molar::core::Atom,
    pos: &'static mut molar::core::Pos,
    // Index is readonly
    #[pyo3(get)]
    id: usize,
}

#[pymethods]
impl Particle {
    #[getter]
    fn get_pos<'py>(&mut self, py: Python<'py>) -> Bound<'py, PyAny> {
        map_pos_to_pyarray(py, self.pos)
    }

    #[setter]
    fn set_pos<'py>(&mut self, value: PyReadonlyArray1<'py,f32>) -> PyResult<()> {
        unsafe {
            std::ptr::copy_nonoverlapping(value.data(), self.pos.coords.as_mut_ptr(), 3);
        }
        Ok(())
    }

    // name
    #[getter(name)]
    fn get_name(&self, _py: Python) -> &str {
        self.atom.name.as_str()
    }

    #[setter(name)]
    fn set_name(&mut self, value: &str) {
        self.atom.name = value.into();
    }

    // resname
    #[getter(resname)]
    fn get_resname(&self, _py: Python) -> &str {
        self.atom.resname.as_str()
    }

    #[setter(resname)]
    fn set_resname(&mut self, value: &str) {
        self.atom.resname = value.into();
    }

    // resid
    #[getter(resname)]
    fn get_resid(&self, _py: Python) -> i32 {
        self.atom.resid
    }

    #[setter(resname)]
    fn set_resid(&mut self, value: i32) {
        self.atom.resid = value;
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self, _py: Python) -> u8 {
        self.atom.atomic_number
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        self.atom.atomic_number = value;
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self, _py: Python) -> f32 {
        self.atom.mass
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        self.atom.mass = value;
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self, _py: Python) -> f32 {
        self.atom.charge
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        self.atom.charge = value;
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self, _py: Python) -> &str {
        self.atom.type_name.as_str()
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, value: &str) {
        self.atom.type_name = value.into();
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self, _py: Python) -> u32 {
        self.atom.type_id
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        self.atom.type_id = value;
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self, _py: Python) -> char {
        self.atom.chain
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        self.atom.chain = value;
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self, _py: Python) -> f32 {
        self.atom.bfactor
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        self.atom.bfactor = value;
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self, _py: Python) -> f32 {
        self.atom.occupancy
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        self.atom.occupancy = value;
    }
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
#[derive(FromPyObject)]
enum SliceOrInt{
    Slice(Py<pyo3::types::PySlice>),
    Int(isize),
}

#[pyclass(unsendable)]
struct Sel(molar::core::Sel<MutableSerial>);

#[pymethods]
impl Sel {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __getitem__(slf: Bound<Self>, ind: SliceOrInt) -> PyResult<Py<PyAny>> {
        let s = slf.borrow();
        match ind {
            SliceOrInt::Int(i) => {
                let ind = if i<0 {
                    s.len() - i.unsigned_abs()
                } else {
                    i as usize
                };
                s.nth_particle(ind).map(|p| p.into_py_any(slf.py()).unwrap())
            },
            SliceOrInt::Slice(s) => {
                todo!("subselect here");
                Ok(s.into_py_any(slf.py()).unwrap())
            }
        }
        
    }

    fn com<'a>(&self, _py: Python<'a>) -> PyResult<Bound<'a, numpy::PyArray1<f32>>> {
        Ok(copy_pos_to_pyarray(
            _py,
            &self.0.center_of_mass().map_err(|e| anyhow::anyhow!(e))?,
        )) //.into_py(_py))
    }

    fn nth_pos<'a>(slf: Bound<Self>, py: Python<'a>, i: usize) -> PyResult<Bound<'a, PyAny>> {
        let s = slf.borrow();
        let pos = s.0.nth_pos_mut(i).ok_or_else(|| anyhow::anyhow!("Out of bounds"))?;
        Ok(map_pos_to_pyarray(py, pos))
    }

    fn nth_particle(&self, i: usize) -> PyResult<Particle> {
        let p = self.0.nth_particle_mut(i).ok_or_else(|| anyhow::anyhow!("Out of bounds"))?;
        let arom_ptr = p.atom as *mut molar::core::Atom;
        let pos_ptr = p.pos as *mut molar::core::Pos;
        Ok(Particle{
            atom: unsafe {&mut *arom_ptr},
            pos: unsafe {&mut *pos_ptr},
            id: p.id,
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, ParticleIterator> {
        Bound::new(slf.py(),ParticleIterator{
            sel: slf.into(),
            cur: 0,
        }).unwrap().borrow()
    }
}

#[pyclass]
struct ParticleIterator {
    sel: Py<Sel>,
    cur: usize,
}

#[pymethods]
impl ParticleIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyObject> {
        let ret = Python::with_gil(|py|{
            slf.sel.borrow(py).nth_particle(slf.cur).ok().map(|p|{
                p.into_py_any(py).unwrap()
            })
        });
        slf.cur += 1;
        ret
    }
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
            std::ptr::null_mut(),                 // no strides
            data.coords.as_mut_ptr() as *mut c_void, // data
            //npyffi::NPY_ARRAY_WRITEABLE | npyffi::NPY_ARRAY_F_CONTIGUOUS, // flag
            npyffi::NPY_ARRAY_F_CONTIGUOUS, // flag
            std::ptr::null_mut(),             // obj
        );
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

// unsafe fn pyarray_from_pos<'py>(
//     py: Python<'py>,
//     data: &mut Pos,
//     container: *mut PyAny,
// ) -> Bound<'py, Self>
// where
//     ID: IntoDimension<Dim = D>,
// {
//     let mut dims = dims.into_dimension();
//     let ptr = PY_ARRAY_API.PyArray_NewFromDescr(
//         py,
//         PY_ARRAY_API.get_type_object(py, npyffi::NpyTypes::PyArray_Type),
//         T::get_dtype(py).into_dtype_ptr(),
//         dims.ndim_cint(),
//         dims.as_dims_ptr(),
//         strides as *mut npy_intp,    // strides
//         data_ptr as *mut c_void,     // data
//         npyffi::NPY_ARRAY_WRITEABLE, // flag
//         ptr::null_mut(),             // obj
//     );

//     PY_ARRAY_API.PyArray_SetBaseObject(
//         py,
//         ptr as *mut npyffi::PyArrayObject,
//         container as *mut ffi::PyObject,
//     );

//     Bound::from_owned_ptr(py, ptr).downcast_into_unchecked()
// }
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
