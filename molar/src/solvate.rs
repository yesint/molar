//! Filling a periodic box with solvent.
use crate::prelude::*;
use std::collections::HashSet;
use thiserror::Error;

/// Error of [`System::solvate`].
#[derive(Error, Debug)]
pub enum SolvateError {
    #[error("can't solvate a system without a periodic box")]
    NoSoluteBox,
    #[error("the solvent lacks a periodic box")]
    NoSolventBox,
    #[error("triclinic solvent boxes are not supported")]
    TriclinicSolvent,
    #[error(transparent)]
    Selection(#[from] SelectionError),
    #[error(transparent)]
    Builder(#[from] BuilderError),
    #[error(transparent)]
    DistanceSearch(#[from] DistanceSearchError),
}

impl System {
    /// Fill this system's periodic box with copies of a periodic (rectangular) `solvent`
    /// box. Solvent residues are kept whole: a residue is added only if all of its atoms
    /// lie inside the box and none overlaps a solute atom (sum of van der Waals radii,
    /// periodic). The solvent is appended after the existing atoms; returns the number of
    /// residues added.
    ///
    /// The box may be triclinic: the solvent is tiled over the box's whole bounding
    /// region, including parts at negative coordinates (e.g. a hexagonal cell with
    /// b = (−a/2, a·√3/2, 0)), and trimmed to the box.
    pub fn solvate(&mut self, solvent: &System) -> Result<usize, SolvateError> {
        let pbox = self.get_box().ok_or(SolvateError::NoSoluteBox)?.clone();
        let solvent_box = solvent.get_box().ok_or(SolvateError::NoSolventBox)?;
        if solvent_box.is_triclinic() {
            return Err(SolvateError::TriclinicSolvent);
        }
        let ext = solvent_box.get_box_extents();

        // Bounding region of the box in lab coordinates.
        let mut lo = Vector3f::repeat(Float::INFINITY);
        let mut hi = Vector3f::repeat(Float::NEG_INFINITY);
        for corner in 0..8 {
            let f = Vector3f::new((corner & 1) as Float, ((corner >> 1) & 1) as Float, ((corner >> 2) & 1) as Float);
            let p = pbox.to_lab_coords(&f);
            lo = lo.inf(&p);
            hi = hi.sup(&p);
        }
        // Whole solvent cells from below `lo` to above `hi`.
        let first = Vector3f::from_fn(|i, _| (lo[i] / ext[i]).floor());
        let mut nbox = [0usize; 3];
        for i in 0..3 {
            nbox[i] = ((hi[i] / ext[i]).ceil() - first[i]).max(1.) as usize;
        }
        let shift = first.component_mul(&ext);

        let mut tiled = System::new(solvent.topology().clone(), solvent.state().clone())?;
        tiled.multiply_periodically(nbox)?;
        for p in tiled.iter_pos_mut() {
            *p += shift;
        }

        // Whole residues inside the box.
        let mut inside = vec![];
        let all = tiled.select_all_bound();
        'residue: for res in all.split_resindex_bound() {
            for p in res.iter_pos() {
                if !pbox.is_inside(p) {
                    continue 'residue;
                }
            }
            inside.extend(res.iter_index());
        }
        if inside.is_empty() {
            return Ok(0);
        }
        let inside = tiled.select_bound(inside)?;

        // Residues with an atom overlapping the solute.
        let mut clashing = HashSet::new();
        if self.len() > 0 {
            let vdw1: Vec<Float> = inside.iter_atoms().map(|a| a.vdw()).collect();
            let vdw2: Vec<Float> = self.iter_atoms().map(|a| a.vdw()).collect();
            let overlap: Vec<usize> =
                distance_search_double_vdw_pbc(inside.iter_pos(), self.iter_pos(), &vdw1, &vdw2, &pbox, PBC_FULL)?;
            for i in overlap {
                clashing.insert(unsafe { inside.get_atom_unchecked(i) }.get_resindex());
            }
        }
        let mut keep = vec![];
        let mut added = 0;
        for res in inside.split_resindex_bound() {
            if !clashing.contains(&res.first_atom().get_resindex()) {
                keep.extend(res.iter_index());
                added += 1;
            }
        }
        if added > 0 {
            let keep = tiled.select_bound(keep)?;
            self.append(&keep)?;
        }
        Ok(added)
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    /// A 1 nm cubic solvent box with 8 one-atom residues ("O", a water oxygen stand-in).
    fn solvent() -> System {
        let mut top = Topology::default();
        let mut st = State::default();
        let mut k = 0;
        for i in 0..2 {
            for j in 0..2 {
                for l in 0..2 {
                    top.atoms.push(&Atom::new().with_name("OW").with_resname("SOL").with_resid(k).with_resindex(k as usize).guess());
                    st.coords.push(Pos::new(0.25 + 0.5 * i as Float, 0.25 + 0.5 * j as Float, 0.25 + 0.5 * l as Float));
                    k += 1;
                }
            }
        }
        st.pbox = Some(PeriodicBox::from_matrix(Matrix3f::from_diagonal(&Vector3f::new(1., 1., 1.))).unwrap());
        System::new(top, st).unwrap()
    }

    fn solute(pbox: PeriodicBox, at: Pos) -> System {
        let mut top = Topology::default();
        let mut st = State::default();
        top.atoms.push(&Atom::new().with_name("C").with_resname("LIG").guess());
        st.coords.push(at);
        st.pbox = Some(pbox);
        System::new(top, st).unwrap()
    }

    #[test]
    fn fills_rectangular_and_hexagonal_boxes() {
        for pbox in [
            PeriodicBox::from_matrix(Matrix3f::from_diagonal(&Vector3f::new(3., 2., 2.))).unwrap(),
            PeriodicBox::from_vectors_angles(3., 3., 2., 90., 90., 120.).unwrap(),
        ] {
            let centre = pbox.to_lab_coords(&Vector3f::new(0.5, 0.5, 0.5));
            let mut sys = solute(pbox.clone(), Pos::from(centre));
            let added = sys.solvate(&solvent()).unwrap();
            let volume = pbox.get_matrix().determinant().abs();
            // 8 per nm³, less the one overlapping the solute and the boundary losses.
            assert!(added as Float > 0.8 * 8. * volume - 1. && added as Float <= 8. * volume, "{added} in {volume} nm³");
            for p in sys.iter_pos().skip(1) {
                assert!(pbox.is_inside(p), "solvent outside the box");
                assert!(pbox.distance(p, &Pos::from(centre), PBC_FULL) > 0.2, "solvent on the solute");
            }
            // The hexagonal cell reaches x < 0: solvent fills it there too.
            if pbox.is_triclinic() {
                assert!(sys.iter_pos().any(|p| p.x < -0.5));
            }
        }
    }

    #[test]
    fn needs_boxes() {
        let mut sys = solute(PeriodicBox::from_matrix(Matrix3f::from_diagonal(&Vector3f::new(2., 2., 2.))).unwrap(), Pos::origin());
        let mut bare = System::new(solvent().topology().clone(), solvent().state().clone()).unwrap();
        bare.set_box_from(&solute(PeriodicBox::from_vectors_angles(1., 1., 1., 90., 90., 120.).unwrap(), Pos::origin()));
        assert!(matches!(sys.solvate(&bare), Err(super::SolvateError::TriclinicSolvent)));
    }
}
