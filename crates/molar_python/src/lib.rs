use molar::prelude::*;
use numpy::{npyffi, PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods, ToNpyDims, PY_ARRAY_API};
use pyo3::{prelude::*, AsPyPointer, IntoPyObjectExt};
use core::num;
use std::ffi::c_void;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyclass]
struct Atom(molar::core::Atom);

#[pymethods]
impl Atom {
    #[new]
    fn new() -> Self {
        Self {
            0: Default::default(),
        }
    }

    // name
    #[getter(name)]
    fn get_name(&self, _py: Python) -> &str {
        self.0.name.as_str()
    }

    #[setter(name)]
    fn set_name(&mut self, value: &str) {
        self.0.name = value.into();
    }

    // resname
    #[getter(resname)]
    fn get_resname(&self, _py: Python) -> &str {
        self.0.resname.as_str()
    }

    #[setter(resname)]
    fn set_resname(&mut self, value: &str) {
        self.0.resname = value.into();
    }

    // resid
    #[getter(resname)]
    fn get_resid(&self, _py: Python) -> i32 {
        self.0.resid
    }

    #[setter(resname)]
    fn set_resid(&mut self, value: i32) {
        self.0.resid = value;
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self, _py: Python) -> u8 {
        self.0.atomic_number
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        self.0.atomic_number = value;
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self, _py: Python) -> f32 {
        self.0.mass
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        self.0.mass = value;
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self, _py: Python) -> f32 {
        self.0.charge
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        self.0.charge = value;
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self, _py: Python) -> &str {
        self.0.type_name.as_str()
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, value: &str) {
        self.0.type_name = value.into();
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self, _py: Python) -> u32 {
        self.0.type_id
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        self.0.type_id = value;
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self, _py: Python) -> char {
        self.0.chain
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        self.0.chain = value;
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self, _py: Python) -> f32 {
        self.0.bfactor
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        self.0.bfactor = value;
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self, _py: Python) -> f32 {
        self.0.occupancy
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        self.0.occupancy = value;
    }
}

#[pyclass(unsendable)]
struct Particle {
    atom: &'static mut molar::core::Atom,
    // PyArray mapped to pos
    pos: *mut npyffi::PyArrayObject,
    // Index is readonly
    #[pyo3(get)]
    id: usize,
}

#[pymethods]
impl Particle {
    //pos
    #[getter]
    fn get_pos<'py>(&self, py: Python<'py>) -> Bound<'py, PyAny> {
        // Return an owned ptr to avoid incorrect reference count.
        // I still don't quite understand why this is necessary.
        unsafe { Bound::from_owned_ptr(py, self.pos.cast()) }
    }

    #[setter]
    fn set_pos<'py>(&mut self, value: PyReadonlyArray1<'py, f32>) {
        unsafe {
            std::ptr::copy_nonoverlapping(value.data(), (*self.pos).data.cast(), 3);
        }
    }

    #[getter(x)]
    fn get_x(&self) -> f32 {
        unsafe { *(*self.pos).data.cast() }
    }

    #[setter(x)]
    fn set_x(&mut self, value: f32) {
        unsafe { *(*self.pos).data.cast() = value };
    }

    #[getter(y)]
    fn get_y(&self) -> f32 {
        unsafe { *(*self.pos).data.cast::<f32>().offset(1) }
    }

    #[setter(y)]
    fn set_y(&mut self, value: f32) {
        unsafe { *(*self.pos).data.cast::<f32>().offset(1) = value };
    }

    #[getter(z)]
    fn get_z(&self) -> f32 {
        unsafe { *(*self.pos).data.cast::<f32>().offset(2) }
    }

    #[setter(z)]
    fn set_z(&mut self, value: f32) {
        unsafe { *(*self.pos).data.cast::<f32>().offset(2) = value };
    }

    // atom
    #[getter(atom)]
    fn get_atom(&self, _py: Python) -> Atom {
        Atom(self.atom.clone())
    }

    #[setter(atom)]
    fn set_atom(&mut self, value: &Atom) {
        *self.atom = value.0.clone();
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

#[pyclass(sequence, unsendable)]
struct Sel(molar::core::Sel<MutableSerial>);

#[pymethods]
impl Sel {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __call__(&self, sel_str: &str) -> PyResult<Sel> {
        Ok(Sel(self
            .0
            .subsel_str(sel_str)
            .map_err(|e| anyhow::anyhow!(e))?))
    }

    // Indexing
    fn __getitem__(slf: Bound<Self>, i: isize) -> PyResult<Py<PyAny>> {
        let s = slf.borrow();
        let ind = if i < 0 {
            if i.abs() > s.__len__() as isize {
                return Err(anyhow::anyhow!(
                    "Negative index {i} is out of bounds {}:-1",
                    -(s.__len__() as isize)
                )
                .into());
            }
            s.__len__() - i.unsigned_abs()
        } else if i >= s.__len__() as isize {
            return Err(anyhow::anyhow!("Index {} is out of bounds 0:{}", i, s.__len__()).into());
        } else {
            i as usize
        };

        // Call Rust function
        let p = s.0.nth_particle_mut(ind).unwrap();
        Ok(Particle {
            atom: unsafe { &mut *(p.atom as *mut molar::core::Atom) },
            pos: map_pyarray_to_pos(slf.py(), p.pos, slf.clone().into_py_any(slf.py()).unwrap()),
            id: p.id,
        }
        .into_py_any(slf.py())?)
    }

    fn com<'a>(&self, _py: Python<'a>) -> PyResult<Bound<'a, numpy::PyArray1<f32>>> {
        Ok(copy_pos_to_pyarray(
            _py,
            &self.0.center_of_mass().map_err(|e| anyhow::anyhow!(e))?,
        ))
    }

    // Iteration protocol
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, ParticleIterator> {
        Bound::new(
            slf.py(),
            ParticleIterator {
                sel: slf.into(),
                cur: 0,
            },
        )
        .unwrap()
        .borrow()
    }

    fn get_coord<'py>(&self, py: Python<'py>) -> Py<PyAny> {
        use numpy::Element;
        use numpy::PyArrayDescrMethods;

        let mut dims = numpy::ndarray::Dim((3,self.__len__()));

        unsafe {
            let arr = PY_ARRAY_API.PyArray_NewFromDescr(
                py,
                PY_ARRAY_API.get_type_object(py, npyffi::NpyTypes::PyArray_Type),
                f32::get_dtype(py).into_dtype_ptr(),
                dims.ndim_cint(),
                dims.as_dims_ptr(),
                std::ptr::null_mut(), // no strides
                std::ptr::null_mut(), // no data, allocate new buffer
                npyffi::NPY_ARRAY_WRITEABLE | npyffi::NPY_ARRAY_OWNDATA, // flag
                std::ptr::null_mut(), // obj
            ) as *mut npyffi::PyArrayObject;

            let ptr = (*arr).data.cast::<f32>();
            for i in 0..self.__len__() {
                let pos = self.0.nth_pos_unchecked(i);
                std::ptr::copy_nonoverlapping(pos.coords.as_ptr(), ptr.offset(i as isize * 3), 3);
            }

            Py::from_owned_ptr(py, arr.cast())
        }
    }

    fn set_coord(&mut self, arr: PyReadonlyArray2<f32>) {
        // Check if the shape is correct
        if arr.shape() != [3,self.__len__()] {
            panic!("Array shape must be [3, {}], not {:?}",self.__len__(),arr.shape());
        }
        let ptr = arr.data();

        unsafe {
            for i in 0..self.__len__() {
                let pos_ptr = self.0.nth_pos_mut_unchecked(i).coords.as_mut_ptr();
                std::ptr::copy_nonoverlapping(ptr.offset(i as isize * 3), pos_ptr, 3);
            }
        }
    }
}

#[pyclass]
struct ParticleIterator {
    sel: Py<Sel>,
    cur: isize,
}

#[pymethods]
impl ParticleIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyObject> {
        let ret = Python::with_gil(|py| {
            let s = slf.sel.bind(py);
            Sel::__getitem__(s.clone(), slf.cur)
        })
        .ok();
        slf.cur += 1;
        ret
    }
}

// Constructs PyArray backed by existing Pos data.
fn map_pyarray_to_pos<'py>(
    py: Python<'py>,
    data: &mut molar::core::Pos,
    parent: Py<PyAny>,
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
            std::ptr::null_mut(),                    // no strides
            data.coords.as_mut_ptr() as *mut c_void, // data
            npyffi::NPY_ARRAY_WRITEABLE, // flag
            std::ptr::null_mut(), // obj
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
