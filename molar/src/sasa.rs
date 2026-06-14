use crate::prelude::*;

/// Triangle mesh of the solvent-accessible surface (vertices/normals in nm, `u32`
/// indices, per-vertex source-atom index in selection order). Re-exported from the
/// `powersasa` backend; produced by [`Sasa::surface_mesh`].
pub use powersasa::SurfaceMesh;

/// Solvent-Accessible Surface Area calculator backed by the pure-Rust PowerSASA algorithm.
///
/// `Sasa` holds both the power diagram (which is expensive to build) and the computed
/// per-atom areas/volumes. For trajectory analysis, build it once with [`Sasa::new`] or
/// [`Sasa::new_with_volume`] and then call [`Sasa::update`] for subsequent frames to
/// reuse the diagram structure without a full reconstruction.
pub struct Sasa {
    inner: powersasa::PowerSasa<Float>,
    probe_r: Float,
}

impl Sasa {
    /// Standard water probe radius in nm (0.14 nm).
    pub const DEFAULT_PROBE_R: Float = 0.14;

    /// Build power diagram + compute per-atom SASA (areas only).
    pub fn new<S>(sel: &S) -> Result<Self, MeasureError>
    where
        S: AtomProvider + PosProvider + LenProvider + ?Sized,
    {
        Self::build(sel, false, Self::DEFAULT_PROBE_R)
    }

    /// Build power diagram + compute per-atom SASA and volumes.
    pub fn new_with_volume<S>(sel: &S) -> Result<Self, MeasureError>
    where
        S: AtomProvider + PosProvider + LenProvider + ?Sized,
    {
        Self::build(sel, true, Self::DEFAULT_PROBE_R)
    }

    /// Build with a custom probe radius (nm). `probe = 0` gives the van-der-Waals
    /// surface (union of vdW spheres); the default (0.14 nm) gives the water SAS.
    pub fn new_with_probe<S>(sel: &S, probe: Float) -> Result<Self, MeasureError>
    where
        S: AtomProvider + PosProvider + LenProvider + ?Sized,
    {
        Self::build(sel, false, probe)
    }

    fn build<S>(sel: &S, with_vol: bool, probe_r: Float) -> Result<Self, MeasureError>
    where
        S: AtomProvider + PosProvider + LenProvider + ?Sized,
    {
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
        S: AtomProvider + PosProvider + LenProvider + ?Sized,
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
    pub fn areas(&self) -> &[Float] {
        self.inner.per_atom_sasa()
    }

    /// Sum of per-atom areas (nm²).
    pub fn total_area(&self) -> Float {
        self.inner.per_atom_sasa().iter().sum()
    }

    /// Per-atom volumes. Only meaningful when constructed with [`Sasa::new_with_volume`].
    pub fn volumes(&self) -> &[Float] {
        self.inner.per_atom_vol()
    }

    /// Sum of per-atom volumes. Only meaningful when constructed with [`Sasa::new_with_volume`].
    pub fn total_volume(&self) -> Float {
        self.inner.per_atom_vol().iter().sum()
    }

    /// Extract the solvent-accessible surface as a triangle mesh (vertices/normals
    /// in nm, in the selection's coordinate frame). `subdiv` is the per-atom
    /// icosphere subdivision level (0 → 20 triangles, 1 → 80, 2 → 320, …): higher
    /// is smoother and heavier. Per-vertex `atom_ids` index into the selection in
    /// iteration order (the same order the areas are reported), so callers can
    /// color the surface per atom. The SAS radius is `vdw + probe` (the probe used
    /// at construction). The diagram + per-atom SASA computed at construction make
    /// fully-buried atoms drop out automatically.
    pub fn surface_mesh(&self, subdiv: usize) -> SurfaceMesh<Float> {
        self.inner.surface_mesh(subdiv)
    }

    /// Extract the **solvent-excluded surface** (SES / Connolly / rolling-probe) as
    /// a smooth triangle mesh (vertices/normals in nm). Unlike [`Self::surface_mesh`]
    /// (the creased SAS union of spheres), this is the smooth surface traced by the
    /// probe rolling over the atoms — convex contact patches bridged by toroidal and
    /// concave reentrant patches. `subdiv` controls tessellation density. Per-vertex
    /// `atom_ids` index into the selection in iteration order, for per-atom coloring.
    pub fn ses_mesh(&self, subdiv: usize) -> SurfaceMesh<Float> {
        self.inner.ses_mesh(self.probe_r, subdiv)
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn surface_mesh_is_nonempty_and_consistent() -> anyhow::Result<()> {
        let sys = System::from_file("tests/albumin.pdb")?;
        let sel = sys.select_bound("resindex 1:20")?;
        let sasa = Sasa::new(&sel)?;
        let mesh = sasa.surface_mesh(2);

        assert!(!mesh.vertices.is_empty(), "surface mesh is empty");
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
        assert_eq!(mesh.vertices.len(), mesh.atom_ids.len());
        assert_eq!(mesh.indices.len() % 3, 0);
        let n = sel.len();
        assert!(
            mesh.atom_ids.iter().all(|&a| (a as usize) < n),
            "atom id out of selection range"
        );
        assert!(
            *mesh.indices.iter().max().unwrap() < mesh.vertices.len() as u32,
            "index out of vertex range"
        );
        Ok(())
    }
}
