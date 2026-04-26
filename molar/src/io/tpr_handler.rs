use std::{ffi::CString, path::Path, sync::Arc};

use molar_gromacs::{TprAtom, TprBond, TprHandle, TprMolecule, TprPlugin};
use thiserror::Error;

use crate::atom::{AtomStr, ATOM_NAME_EXPECT, ATOM_RESNAME_EXPECT, ATOM_TYPE_NAME_EXPECT};
use crate::prelude::*;

pub struct TprFileHandler {
    plugin:       Arc<TprPlugin>,
    handle:       *mut TprHandle,
    already_read: bool,
}

// TprHandle is heap-allocated C++ data managed through the plugin functions.
unsafe impl Send for TprFileHandler {}

#[derive(Debug, Error)]
pub enum TprHandlerError {
    #[error("Gromacs plugin not found (is MOLAR_GROMACS_PLUGIN set correctly?): {0}")]
    GromacsNotFound(String),

    #[error("failed to open TPR file: {0}")]
    OpenFailed(String),

    #[error("unexpected null character in path")]
    CStringNull(#[from] std::ffi::NulError),

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),
}

impl Drop for TprFileHandler {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { (self.plugin.fns.tpr_close)(self.handle) };
            self.handle = std::ptr::null_mut();
        }
    }
}

impl FileFormatHandler for TprFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let plugin = TprPlugin::get_cached()
            .map_err(|e| TprHandlerError::GromacsNotFound(e.to_string()))?;

        let c_path = CString::new(fname.as_ref().to_str().unwrap())
            .map_err(TprHandlerError::CStringNull)?;
        let handle = unsafe { (plugin.fns.tpr_open)(c_path.as_ptr()) };

        if handle.is_null() {
            let msg = plugin.last_error();
            return Err(TprHandlerError::OpenFailed(msg).into());
        }

        Ok(TprFileHandler { plugin, handle, already_read: false })
    }

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        if self.already_read {
            return Err(FileFormatError::Eof);
        }
        self.already_read = true;

        let fns = &self.plugin.fns;
        let h   = self.handle;

        let natoms = unsafe { (fns.tpr_natoms)(h) };
        let nbonds = unsafe { (fns.tpr_nbonds)(h) };
        let nmols  = unsafe { (fns.tpr_nmolecules)(h) };

        // Allocate caller-owned buffers without initialization; C++ fills them.
        let mut atoms_buf:  Vec<TprAtom>     = Vec::with_capacity(natoms);
        let mut bonds_buf:  Vec<TprBond>     = Vec::with_capacity(nbonds);
        let mut mols_buf:   Vec<TprMolecule> = Vec::with_capacity(nmols);
        let mut coords_buf: Vec<f32>         = Vec::with_capacity(natoms * 3);
        let mut box_buf                      = std::mem::MaybeUninit::<[f32; 9]>::uninit();

        unsafe {
            atoms_buf.set_len(natoms);
            bonds_buf.set_len(nbonds);
            mols_buf.set_len(nmols);
            coords_buf.set_len(natoms * 3);
            (fns.tpr_fill_atoms)(h, atoms_buf.as_mut_ptr());
            (fns.tpr_fill_bonds)(h, bonds_buf.as_mut_ptr());
            (fns.tpr_fill_molecules)(h, mols_buf.as_mut_ptr());
            (fns.tpr_fill_coords)(h, coords_buf.as_mut_ptr());
            (fns.tpr_fill_box)(h, box_buf.as_mut_ptr() as *mut f32);
        }
        let box_buf = unsafe { box_buf.assume_init() };

        //====================
        // Build Topology
        //====================
        let mut top = Topology::default();
        top.atoms.reserve(natoms);

        for a in &atoms_buf {
            let chain = if a.chain == 0 || a.chain == b'\0' { ' ' } else { a.chain as char };
            top.atoms.push(Atom {
                name:          atom_str_from_buf(&a.name,      ATOM_NAME_EXPECT),
                resname:       atom_str_from_buf(&a.resname,   ATOM_RESNAME_EXPECT),
                type_name:     atom_str_from_buf(&a.type_name, ATOM_TYPE_NAME_EXPECT),
                chain,
                resid:         a.resind as i32,
                type_id:       a.type_id,
                atomic_number: a.atomic_number as u8,
                // TPR plugin C ABI exposes f32; cast at boundary.
                charge:        a.charge   as Float,
                mass:          a.mass     as Float,
                occupancy:     a.occupancy as Float,
                bfactor:       a.bfactor  as Float,
                ..Default::default()
            });
        }

        for b in &bonds_buf {
            top.bonds.push([b.atom1 as usize, b.atom2 as usize]);
        }

        for m in &mols_buf {
            top.molecules.push([m.start as usize, m.end as usize]);
        }

        top.assign_resindex();

        //====================
        // Build State
        //====================
        // TPR plugin returns f32 buffers (see wrapper.hpp). Cast at boundary.
        let mut st = State::default();
        st.coords.reserve(natoms);
        for i in 0..natoms {
            st.coords.push(Pos::new(
                coords_buf[i * 3]     as Float,
                coords_buf[i * 3 + 1] as Float,
                coords_buf[i * 3 + 2] as Float,
            ));
        }

        let m = Matrix3f::from_iterator(box_buf.iter().map(|x| *x as Float));
        st.pbox = Some(PeriodicBox::from_matrix(m).map_err(TprHandlerError::Pbc)?);

        Ok((top, st))
    }
}

fn atom_str_from_buf(buf: &[u8; 8], expect: &str) -> AtomStr {
    let len = buf.iter().position(|&b| b == 0).unwrap_or(8);
    let s = std::str::from_utf8(&buf[..len]).expect("Gromacs atom strings are ASCII");
    AtomStr::try_from_str(s).expect(expect)
}

#[cfg(test)]
mod tests {
    use crate::io::TprFileHandler;
    use crate::prelude::*;

    #[test]
    fn test_tpr() {
        let mut h = match TprFileHandler::open("tests/topol.tpr") {
            Ok(h) => h,
            Err(e) => {
                // Check the error source chain for the "plugin not found" case.
                let is_not_found = std::iter::successors(
                    Some(&e as &dyn std::error::Error),
                    |e| e.source(),
                )
                .any(|e| e.to_string().contains("plugin not found"));
                if is_not_found {
                    eprintln!("Skipping test_tpr: Gromacs plugin not available");
                    return;
                }
                panic!("unexpected error: {e}");
            }
        };
        let (top, st) = h.read().unwrap();
        println!("natoms: {:?}", top.len());
        println!("nbonds: {:?}", BondProvider::num_bonds(&top));
        println!("nmolecules: {:?}", MolProvider::num_molecules(&top));
        println!("state sz: {:?}", st.len());
    }
}
