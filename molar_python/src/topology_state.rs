use pyo3::prelude::*;
use molar::prelude::*;
use super::periodic_box::PeriodicBoxPy;
use anyhow::anyhow;

#[pyclass(unsendable, name="State")]
pub(crate) struct StatePy(pub(crate) State);

#[pymethods]
impl StatePy {
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
    fn get_box(&self) -> anyhow::Result<PeriodicBoxPy> {
        Ok(PeriodicBoxPy(
            self.0
                .get_box()
                .ok_or_else(|| anyhow!("No periodic box"))?
                .clone(),
        ))
    }

    #[setter]
    fn set_box(&mut self, val: Bound<'_, PeriodicBoxPy>) -> anyhow::Result<()> {
        let b = self
            .0
            .get_box_mut()
            .ok_or_else(|| anyhow!("No periodic box"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }
}

#[pyclass(unsendable, name="Topology")]
pub(crate) struct TopologyPy(pub(crate) Topology);

