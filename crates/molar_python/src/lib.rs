use anyhow::bail;
use molar::prelude::*;
use numpy::{
    nalgebra, npyffi, AllowTypeChange, PyArray1, PyArrayLike1, PyArrayLike2, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToNpyDims, PY_ARRAY_API
};
use pyo3::{
    prelude::*,
    types::{PyList, PyString, PyTuple},
    IntoPyObjectExt,
};
use std::ffi::c_void;

pub mod atom;
use atom::Atom;

pub mod particle;
use particle::Particle;

//-------------------------------------------
#[pyclass]
struct PeriodicBox(molar::core::PeriodicBox);

#[pymethods]
impl PeriodicBox {
    #[staticmethod]
    fn from_matrix<'py>(arr: PyArrayLike2<'py, f32, AllowTypeChange>) -> anyhow::Result<Self> {
        if arr.shape() != [3, 3] {
            bail!("Matrix must be 3x3");
        }

        let mut m = Matrix3f::zeros();
        unsafe {
            std::ptr::copy_nonoverlapping(arr.data(), m.as_mut_ptr(), 9);
        }
        Ok(PeriodicBox(molar::core::PeriodicBox::from_matrix(m)?))
    }

    #[staticmethod]
    fn from_vectors_angles<'py>(
        v_arr: PyArrayLike1<'py, f32, AllowTypeChange>,
        a_arr: PyArrayLike1<'py, f32, AllowTypeChange>,
    ) -> anyhow::Result<Self> {
        if v_arr.len() != 3 {
            bail!("Vectors array must be of size 3");
        }
        if a_arr.len() != 3 {
            bail!("Angles array must be of size 3");
        }
        let mut v = Vector3f::zeros();
        let mut a = Vector3f::zeros();
        unsafe {
            std::ptr::copy_nonoverlapping(v_arr.data(), v.as_mut_ptr(), 3);
            std::ptr::copy_nonoverlapping(a_arr.data(), a.as_mut_ptr(), 3);
        }
        Ok(PeriodicBox(molar::core::PeriodicBox::from_vectors_angles(
            v[0], v[1], v[2], a[0], a[1], a[2],
        )?))
    }

    fn to_vectors_angles(slf: Bound<Self>) -> (Bound<PyArray1<f32>>, Bound<PyArray1<f32>>) {
        let (v,a) = slf.borrow().0.to_vectors_angles();
        let v_arr = copy_vec3_to_pyarray(slf.py(), &v);
        let a_arr = copy_vec3_to_pyarray(slf.py(), &a);
        (v_arr, a_arr)
    }

    fn shortest_vector<'py>(&self, py: Python<'py>, arr: PyArrayLike1<'py, f32, AllowTypeChange>) -> Bound<'py,PyArray1<f32>> {
        if arr.len() != 3 {
            panic!("Vector array must be of size 3");
        }
        if arr.strides()[0] != std::mem::size_of::<f32>() as isize {
            panic!("Vector must be contigous");
        };
        let mut inp_v = Vector3f::zeros();
        unsafe {
            std::ptr::copy_nonoverlapping(arr.data(), inp_v.as_mut_ptr(), 3);
        }
        let out_v = self.0.shortest_vector(&inp_v);
        copy_vec3_to_pyarray(py, &out_v)
    }
}

#[pyclass]
struct Topology(Option<molar::core::Topology>);

#[pyclass]
struct State(Option<molar::core::State>);

#[pymethods]
impl State {
    fn __len__(&self) -> anyhow::Result<usize> {
        self.0
            .as_ref()
            .map(|s| s.len())
            .ok_or_else(|| anyhow::anyhow!("State is moved"))
    }

    #[getter]
    fn time(&self) -> anyhow::Result<f32> {
        self.0
            .as_ref()
            .map(|s| s.get_time())
            .ok_or_else(|| anyhow::anyhow!("State is moved"))
    }

    #[getter]
    fn get_box(&self) -> anyhow::Result<PeriodicBox> {
        Ok(PeriodicBox(
            self.0
                .as_ref()
                .map(|s| s.get_box())
                .ok_or_else(|| anyhow::anyhow!("State is moved"))?
                .ok_or_else(|| anyhow::anyhow!("No periodic box"))?
                .clone(),
        ))
    }

    #[setter]
    fn set_box(&mut self, val: Bound<'_, PeriodicBox>) -> anyhow::Result<()> {
        let b = self
            .0
            .as_mut()
            .map(|s| s.get_box_mut())
            .ok_or_else(|| anyhow::anyhow!("State is moved"))?
            .ok_or_else(|| anyhow::anyhow!("No periodic box"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }
}

#[pyclass(unsendable)]
struct FileHandler(molar::io::FileHandler);

#[pymethods]
impl FileHandler {
    #[new]
    fn new(fname: &str) -> anyhow::Result<Self> {
        Ok(FileHandler(molar::io::FileHandler::open(fname)?))
    }

    fn read(&mut self) -> anyhow::Result<(Topology, State)> {
        let (top, st) = self.0.read()?;
        Ok((Topology(Some(top)), State(Some(st))))
    }

    fn read_topology(&mut self) -> anyhow::Result<Topology> {
        let top = self.0.read_topology()?;
        Ok(Topology(Some(top)))
    }

    fn read_state(&mut self) -> anyhow::Result<State> {
        let st = self.0.read_state()?;
        Ok(State(st))
    }

    fn write(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        if let Ok(s) = data.extract::<PyRef<'_, Source>>() {
            self.0.write(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            self.0.write(&s.0)?;
        } else if let Ok(s) = data.downcast::<PyTuple>() {
            if s.len() != 2 {
                return Err(anyhow::anyhow!("Tuple must have two elements"));
            }
            let top = s
                .iter()
                .next()
                .unwrap()
                .extract::<PyRefMut<'_, Topology>>()?
                .0
                .take()
                .ok_or_else(|| anyhow::anyhow!("Topology is moved"))?;
            let st = s
                .iter()
                .next()
                .unwrap()
                .extract::<PyRefMut<'_, State>>()?
                .0
                .take()
                .ok_or_else(|| anyhow::anyhow!("State is moved"))?;
            self.0.write(&(top, st))?;
        } else {
            return Err(anyhow::anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn write_topology(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        if let Ok(s) = data.extract::<PyRef<'_, Source>>() {
            self.0.write_topology(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            self.0.write_topology(&s.0)?;
        } else if let Ok(mut s) = data.extract::<PyRefMut<'_, Topology>>() {
            self.0.write_topology(
                &s.0.take()
                    .ok_or_else(|| anyhow::anyhow!("Topology is moved"))?,
            )?;
        } else {
            return Err(anyhow::anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn write_state(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        if let Ok(s) = data.extract::<PyRef<'_, Source>>() {
            self.0.write_state(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            self.0.write_state(&s.0)?;
        } else if let Ok(mut s) = data.extract::<PyRefMut<'_, State>>() {
            self.0.write_state(
                &s.0.take()
                    .ok_or_else(|| anyhow::anyhow!("State is moved"))?,
            )?;
        } else {
            return Err(anyhow::anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyObject> {
        let st = slf.read_state().expect("expected state");
        if st.0.is_some() {
            Python::with_gil(|py| Some(st.into_py_any(py)))
                .unwrap()
                .ok()
        } else {
            None
        }
    }
}

#[pyclass(unsendable, sequence)]
struct Source(molar::core::Source<molar::core::BuilderSerial>);

#[pymethods]
impl Source {
    #[new]
    fn new(topology: &mut Topology, state: &mut State) -> PyResult<Self> {
        Ok(Source(
            molar::core::Source::new_builder(
                topology.0.take().unwrap().into(),
                state.0.take().unwrap().into(),
            )
            .map_err(|e| anyhow::anyhow!(e))?,
        ))
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    #[staticmethod]
    fn from_file(fname: &str) -> PyResult<Self> {
        Ok(Source(
            molar::core::Source::builder_from_file(fname).map_err(|e| anyhow::anyhow!(e))?,
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

    #[pyo3(signature = (arg=None))]
    fn __call__(&self, arg: Option<&Bound<'_, PyAny>>) -> PyResult<Sel> {
        if let Some(arg) = arg {
            if arg.is_instance_of::<PyString>() {
                Ok(Sel(self
                    .0
                    .select_str(arg.extract::<String>()?)
                    .map_err(|e| anyhow::anyhow!(e))?))
            } else if arg.is_instance_of::<PyTuple>() {
                let (i1, i2) = arg.extract()?;
                Ok(Sel(self
                    .0
                    .select_range(i1..i2)
                    .map_err(|e| anyhow::anyhow!(e))?))
            } else if arg.is_instance_of::<PyList>() {
                let ind = arg.extract()?;
                Ok(Sel(self
                    .0
                    .select_vec(ind)
                    .map_err(|e| anyhow::anyhow!(e))?))
            } else {
                Err(anyhow::anyhow!(
                    "Invalid argument type {} when creating selection",
                    arg.get_type()
                )
                .into())
            }
        } else {
            Ok(Sel(self.0.select_all().map_err(|e| anyhow::anyhow!(e))?))
        }
    }
}

//====================================

#[pyclass(sequence, unsendable)]
struct Sel(molar::core::Sel<BuilderSerial>);

#[pymethods]
impl Sel {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn invert(&mut self) {
        self.0.invert();
    }

    fn exclude(&mut self, other: &Sel) {
        self.0.exclude(&other.0);
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<Sel> {
        if arg.is_instance_of::<PyString>() {
            Ok(Sel(self
                .0
                .subsel_str(arg.extract::<String>()?)
                .map_err(|e| anyhow::anyhow!(e))?))
        } else if arg.is_instance_of::<PyTuple>() {
            let (i1, i2) = arg.extract()?;
            Ok(Sel(self
                .0
                .subsel_local_range(i1..i2)
                .map_err(|e| anyhow::anyhow!(e))?))
        } else if arg.is_instance_of::<PyList>() {
            let ind: Vec<usize> = arg.extract()?;
            Ok(Sel(self
                .0
                .subsel_iter(ind.into_iter())
                .map_err(|e| anyhow::anyhow!(e))?))
        } else {
            Err(anyhow::anyhow!(
                "Invalid argument type {} when creating selection",
                arg.get_type()
            )
            .into())
        }
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
            pos: map_pyarray_to_pos(slf.py(), p.pos, &slf),
            id: p.id,
        }
        .into_py_any(slf.py())?)
    }

    fn com<'a>(&self, _py: Python<'a>) -> PyResult<Bound<'a, numpy::PyArray1<f32>>> {
        Ok(copy_vec3_to_pyarray(
            _py,
            &self.0.center_of_mass().map_err(|e| anyhow::anyhow!(e))?.coords,
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

    fn get_coord<'py>(&self, py: Python<'py>) -> Bound<'py, PyAny> {
        use numpy::Element;
        use numpy::PyArrayDescrMethods;

        let mut dims = numpy::ndarray::Dim((3, self.__len__()));

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

            Bound::from_owned_ptr(py, arr.cast())
        }
    }

    fn set_coord(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        // Check if the shape is correct
        if arr.shape() != [3, self.__len__()] {
            return Err(anyhow::anyhow!(
                "Array shape must be [3, {}], not {:?}",
                self.__len__(),
                arr.shape()
            ))?;
        }
        let ptr = arr.data();

        unsafe {
            for i in 0..self.__len__() {
                let pos_ptr = self.0.nth_pos_mut_unchecked(i).coords.as_mut_ptr();
                std::ptr::copy_nonoverlapping(ptr.offset(i as isize * 3), pos_ptr, 3);
            }
        }

        Ok(())
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
            std::ptr::null_mut(),                    // no strides
            data.coords.as_mut_ptr() as *mut c_void, // data
            npyffi::NPY_ARRAY_WRITEABLE,             // flag
            std::ptr::null_mut(),                    // obj
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

fn copy_vec3_to_pyarray<'py>(py: Python<'py>, data: &Vector3f) -> Bound<'py, PyArray1<f32>> {
    unsafe {
        let array = PyArray1::<f32>::new(py, 3, true);
        std::ptr::copy_nonoverlapping(data.as_ptr(), array.data(), 3);
        array
    }
}

//====================================

/// A Python module implemented in Rust.
#[pymodule]
fn molar_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Atom>()?;
    m.add_class::<Particle>()?;
    m.add_class::<Topology>()?;
    m.add_class::<State>()?;
    m.add_class::<PeriodicBox>()?;
    m.add_class::<FileHandler>()?;
    m.add_class::<Source>()?;
    m.add_class::<Sel>()?;
    Ok(())
}
