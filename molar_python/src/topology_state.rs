use pyo3::prelude::*;
use molar::prelude::*;
use super::periodic_box::PeriodicBox;
use anyhow::anyhow;

#[pyclass(unsendable)]
pub(crate) struct State(pub(crate) molar::core::State);

#[pymethods]
impl State {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    #[getter]
    fn get_time(&self) -> f32 {
        self.0.get_time()
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        self.0.set_time(t);
    }

    #[getter]
    fn get_box(&self) -> anyhow::Result<PeriodicBox> {
        Ok(PeriodicBox(
            self.0
                .get_box()
                .ok_or_else(|| anyhow!("No periodic box"))?
                .clone(),
        ))
    }

    #[setter]
    fn set_box(&mut self, val: Bound<'_, PeriodicBox>) -> anyhow::Result<()> {
        let b = self
            .0
            .get_box_mut()
            .ok_or_else(|| anyhow!("No periodic box"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }
}

#[pyclass(unsendable)]
pub(crate) struct Topology(pub(crate) molar::core::Topology);

