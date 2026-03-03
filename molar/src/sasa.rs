use crate::prelude::*;

/// Solvent-Accessible Surface Area calculator backed by the pure-Rust PowerSASA algorithm.
///
/// `Sasa` holds both the power diagram (which is expensive to build) and the computed
/// per-atom areas/volumes. For trajectory analysis, build it once with [`Sasa::new`] or
/// [`Sasa::new_with_volume`] and then call [`Sasa::update`] for subsequent frames to
/// reuse the diagram structure without a full reconstruction.
pub struct Sasa {
    inner: powersasa::PowerSasa<f32>,
    probe_r: f32,
}

impl Sasa {
    /// Standard water probe radius in nm (0.14 nm).
    pub const DEFAULT_PROBE_R: f32 = 0.14;

    /// Build power diagram + compute per-atom SASA (areas only).
    pub fn new<S>(sel: &S) -> Result<Self, MeasureError>
    where
        S: AtomIterProvider + PosIterProvider + LenProvider + ?Sized,
    {
        Self::build(sel, false)
    }

    /// Build power diagram + compute per-atom SASA and volumes.
    pub fn new_with_volume<S>(sel: &S) -> Result<Self, MeasureError>
    where
        S: AtomIterProvider + PosIterProvider + LenProvider + ?Sized,
    {
        Self::build(sel, true)
    }

    fn build<S>(sel: &S, with_vol: bool) -> Result<Self, MeasureError>
    where
        S: AtomIterProvider + PosIterProvider + LenProvider + ?Sized,
    {
        let probe_r = Self::DEFAULT_PROBE_R;
        let mut inner = powersasa::PowerSasa::new(
            sel.iter_pos().map(|p| p.coords),
            sel.iter_atoms().map(|a| a.vdw() + probe_r),
            true,  // with_sasa
            false, // with_dsasa
            with_vol,
            false, // with_dvol
        );
        inner.calc_all()?;
        Ok(Self { inner, probe_r })
    }

    /// Update power diagram with new coordinates and recompute SASA.
    ///
    /// Reuses the existing power diagram structure — significantly cheaper than
    /// rebuilding from scratch. Intended for the inner loop of trajectory analysis.
    pub fn update<S>(&mut self, sel: &S) -> Result<(), MeasureError>
    where
        S: AtomIterProvider + PosIterProvider + LenProvider + ?Sized,
    {
        self.inner.update(
            sel.iter_pos().map(|p| p.coords),
            sel.iter_atoms().map(|a| a.vdw() + self.probe_r),
            sel.len(),
        );
        self.inner.calc_all()?;
        Ok(())
    }

    /// Per-atom solvent-accessible surface areas (nm²).
    pub fn areas(&self) -> &[f32] {
        self.inner.per_atom_sasa()
    }

    /// Sum of per-atom areas (nm²).
    pub fn total_area(&self) -> f32 {
        self.inner.per_atom_sasa().iter().sum()
    }

    /// Per-atom volumes. Only meaningful when constructed with [`Sasa::new_with_volume`].
    pub fn volumes(&self) -> &[f32] {
        self.inner.per_atom_vol()
    }

    /// Sum of per-atom volumes. Only meaningful when constructed with [`Sasa::new_with_volume`].
    pub fn total_volume(&self) -> f32 {
        self.inner.per_atom_vol().iter().sum()
    }
}
