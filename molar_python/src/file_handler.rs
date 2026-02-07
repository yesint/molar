use std::path::PathBuf;
use molar::prelude::*;
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::PyValueError;
use pyo3::{exceptions::PyTypeError, prelude::*, types::PyTuple};

use crate::utils::*;
use crate::topology_state::{TopologyPy, StatePy};
use crate::{SystemPy, SelPy};

const ALREADY_TRANDFORMED: &str = "file handler is already transformed to state iterator";

#[pyclass(name = "FileHandler")]
pub struct FileHandlerPy(pub Option<FileHandler>, pub Option<IoStateIterator>);

unsafe impl Send for FileHandlerPy {}
unsafe impl Sync for FileHandlerPy {}

#[pymethods]
impl FileHandlerPy {
    #[new]
    fn new(fname: &str, mode: &str) -> PyResult<Self> {
        match mode {
            "r" => Ok(FileHandlerPy(
                Some(FileHandler::open(fname).map_err(to_py_io_err)?),
                None,
            )),
            "w" => Ok(FileHandlerPy(
                Some(FileHandler::create(fname).map_err(to_py_io_err)?),
                None,
            )),
            _ => Err(PyValueError::new_err("Wrong file open mode")),
        }
    }

    fn read(&mut self) -> PyResult<(TopologyPy, StatePy)> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        let (top, st) = h.read().map_err(to_py_io_err)?;

        Ok((top.into(), st.into()))
    }

    fn read_topology(&mut self) -> PyResult<TopologyPy> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        let top = h.read_topology().map_err(to_py_io_err)?;
        Ok(top.into())
    }

    fn read_state(&mut self) -> PyResult<StatePy> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        let st = h.read_state().map_err(to_py_io_err)?;
        Ok(st.into())
    }

    fn write(&mut self, data: Bound<'_, PyAny>) -> PyResult<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;

        if let Ok(s) = data.extract::<PyRef<'_, SystemPy>>() {
            h.write(&*s).map_err(to_py_io_err)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, SelPy>>() {
            h.write(&*s).map_err(to_py_io_err)?;
        } else if let Ok(s) = data.cast::<PyTuple>() {
            if s.len() != 2 {
                return Err(PyValueError::new_err(format!(
                    "tuple must have two elements, not {}",
                    s.len()
                )));
            }
            let top = s.get_item(0)?.cast::<TopologyPy>()?.as_ptr() as *const TopologyPy;

            let st = s.get_item(1)?.cast::<StatePy>()?.as_ptr() as *const StatePy;
            h.write_topology(unsafe { &*top }).map_err(to_py_io_err)?;
            h.write_state(unsafe { &*st }).map_err(to_py_io_err)?;
        } else {
            return Err(PyTypeError::new_err(format!(
                "Invalid data type {} when writing to file",
                data.get_type().name()?.to_string()
            )));
        }
        Ok(())
    }

    fn write_topology(&mut self, data: Bound<'_, PyAny>) -> PyResult<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        if let Ok(s) = data.cast::<SystemPy>() {
            h.write_topology(s.get().r_top()).map_err(to_py_io_err)?;
        } else if let Ok(s) = data.cast::<SelPy>() {
            h.write_topology(s.get().r_top())
                .map_err(to_py_io_err)?;
        } else if let Ok(s) = data.cast::<TopologyPy>() {
            h.write_topology(&*s.borrow()).map_err(to_py_io_err)?;
        } else {
            return Err(PyTypeError::new_err(format!(
                "Invalid data type {} when writing to file",
                data.get_type().name()?.to_string()
            )));
        }
        Ok(())
    }

    fn write_state(&mut self, data: Bound<'_, PyAny>) -> PyResult<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        if let Ok(s) = data.cast::<SystemPy>() {
            h.write_state(s.get().r_st()).map_err(to_py_io_err)?;
        } else if let Ok(s) = data.cast::<SelPy>() {
            h.write_state(s.get().r_st()).map_err(to_py_io_err)?;
        } else if let Ok(s) = data.cast::<StatePy>() {
            h.write_state(&*s.borrow()).map_err(to_py_io_err)?;
        } else {
            return Err(PyTypeError::new_err(format!(
                "Invalid data type {} when writing to file",
                data.get_type().name()?.to_string()
            )));
        }
        Ok(())
    }

    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        if slf.1.is_none() {
            let h = slf.0.take().unwrap();
            slf.1 = Some(h.into_iter());
        }
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Py<PyAny>> {
        let st = slf
            .1
            .as_mut()
            .expect("Not transformed to state iterator yet")
            .next()
            .map(|st| StatePy::from(st));
        if st.is_some() {
            Some(st.unwrap().into_py_any(slf.py()).unwrap())
        } else {
            None
        }
    }

    fn skip_to_frame(&mut self, fr: usize) -> PyResult<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        h.skip_to_frame(fr).map_err(to_py_io_err)?;
        Ok(())
    }

    fn skip_to_time(&mut self, t: f32) -> PyResult<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        h.skip_to_time(t).map_err(to_py_io_err)?;
        Ok(())
    }

    #[getter]
    fn stats(&self) -> PyResult<FileStatsPy> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        Ok(FileStatsPy(h.stats.clone()))
    }

    #[getter]
    fn file_name(&self) -> PyResult<PathBuf> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| PyTypeError::new_err(ALREADY_TRANDFORMED))?;
        Ok(h.file_path.clone())
    }
}

#[pyclass(name = "FileStats")]
pub struct FileStatsPy(pub FileStats);

#[pymethods]
impl FileStatsPy {
    #[getter]
    fn elapsed_time(&self) -> std::time::Duration {
        self.0.elapsed_time
    }

    #[getter]
    fn frames_processed(&self) -> usize {
        self.0.frames_processed
    }

    #[getter]
    fn cur_t(&self) -> f32 {
        self.0.cur_t
    }

    fn __repr__(&self) -> String {
        format!("{}", self.0)
    }

    fn __str__(&self) -> String {
        format!("{}", self.0)
    }
}
