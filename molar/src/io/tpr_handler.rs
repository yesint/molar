#[cfg(gromacs)]
pub use internal_tpr_enabled::*;
#[cfg(not(gromacs))]
pub use internal_tpr_disabled::*;

#[cfg(gromacs)]
mod internal_tpr_enabled {
    use crate::core::*;
    use molar_gromacs::gromacs_bindings::*;
    use nalgebra::Matrix3;
    use std::{
        ffi::{CStr, CString, NulError},
        ptr::null_mut,
        str::Utf8Error,
    };
    use thiserror::Error;

    pub struct TprFileHandler {
        handle: TprHelper,
    }

    // Allow sending handler between threads
    unsafe impl Send for TprFileHandler {}

    #[derive(Debug, Error)]
    pub enum TprHandlerError {
        #[error("unexpected null characted")]
        CStringNull(#[from] NulError),

        #[error("invalid utf8")]
        CStringUtf8(#[from] Utf8Error),

        #[error("invalid periodic box")]
        Pbc(#[from] PeriodicBoxError),

        #[error("can't read gmx topology")]
        GetTop,
    }

    impl TprFileHandler {
        fn new(fname: &str) -> Result<Self, TprHandlerError> {
            let f_name = CString::new(fname.to_owned())?;
            Ok(TprFileHandler {
                handle: unsafe { TprHelper::new(f_name.as_ptr()) },
            })
        }

        pub fn open(fname: &str) -> Result<Self, TprHandlerError> {
            TprFileHandler::new(fname)
        }

        #[allow(non_snake_case)]
        pub fn read(&mut self) -> Result<(Topology, State), TprHandlerError> {
            //================
            // Read top
            //================
            let gmx_top = unsafe {
                self.handle
                    .get_top()
                    .as_ref()
                    .ok_or_else(|| TprHandlerError::GetTop)?
            };
            let natoms = gmx_top.atoms.nr as usize;
            let nres = gmx_top.atoms.nres as usize;
            println!(">>> {} {}", natoms, nres);

            let mut top = TopologyStorage::default();
            top.atoms.reserve(natoms);

            let gmx_atoms = c_array_to_slice(gmx_top.atoms.atom, natoms);
            let gmx_atomnames = c_array_to_slice(gmx_top.atoms.atomname, natoms);
            let gmx_resinfo = c_array_to_slice(gmx_top.atoms.resinfo, nres);
            let gmx_pdbinfo = if gmx_top.atoms.pdbinfo != null_mut() {
                Some(c_array_to_slice(gmx_top.atoms.pdbinfo, natoms))
            } else {
                None
            };
            let gmx_atomtypes = c_array_to_slice(gmx_top.atoms.atomtype, natoms);

            unsafe {
                for i in 0..natoms {
                    let name = c_ptr_to_str(*gmx_atomnames[i])?;
                    let resi = gmx_atoms[i].resind as usize;
                    let resname = c_ptr_to_str(*gmx_resinfo[resi].name)?;
                    let mut chain = gmx_resinfo[resi].chainid as u8 as char;
                    if chain == '\0' {
                        chain = ' ';
                    }
                    let type_name = c_ptr_to_str(*gmx_atomtypes[i])?;

                    let new_atom = Atom {
                        name,
                        resid: gmx_atoms[i].resind,
                        resname,
                        chain,
                        charge: gmx_atoms[i].q,
                        mass: gmx_atoms[i].m,
                        atomic_number: gmx_atoms[i].atomnumber as u8,
                        type_id: gmx_atoms[i].type_ as u32,
                        type_name,
                        occupancy: gmx_pdbinfo.map(|inf| inf[i].occup).unwrap_or(0.0),
                        bfactor: gmx_pdbinfo.map(|inf| inf[i].bfac).unwrap_or(0.0),
                        ..Default::default()
                    };

                    top.atoms.push(new_atom);
                } //for atoms

                // Parsing idef for bonds
                let functypes = c_array_to_slice(gmx_top.idef.functype, gmx_top.idef.ntypes);
                // Iterate over non-empty interaction lists
                for il in gmx_top.idef.il.iter().filter(|el| el.nr > 0) {
                    let iatoms = c_array_to_slice(il.iatoms, il.nr);
                    // We can check the first type only since we only
                    // need to evaluate the interaction type and not
                    // the concrete parameters for involved atoms
                    match functypes[iatoms[0] as usize] as u32 {
                        F_BONDS | F_G96BONDS | F_HARMONIC | F_FENEBONDS | F_CUBICBONDS
                        | F_CONSTR | F_CONSTRNC => {
                            for el in iatoms.chunks_exact(3) {
                                // el[0] is type and not needed
                                top.bonds.push([el[1] as usize, el[2] as usize]);
                            }
                        }
                        F_SETTLE => {
                            for el in iatoms.chunks_exact(4) {
                                // el[0] is type and not needed
                                // Each settle is 2 bonds
                                top.bonds.push([el[1] as usize, el[2] as usize]);
                                top.bonds.push([el[1] as usize, el[3] as usize]);
                            }
                        }
                        _ => (),
                    } //match
                } //for

                // Read molecules
                let mol_index = c_array_to_slice(gmx_top.mols.index, gmx_top.mols.nr);
                for m in mol_index.chunks_exact(2) {
                    top.molecules.push([m[0] as usize, m[1] as usize - 1]);
                }
            } //unsafe

            // Assign resindexes
            let top: Topology = top.into();
            top.assign_resindex();

            //================
            // Now read state
            //================
            let mut st = StateStorage::default();
            // Gromacs stores coordinates in TPR in internal non-standard vectors
            // So we will need to copy them atom by atom
            st.coords.resize(natoms, Default::default());
            unsafe {
                for i in 0..natoms {
                    // We are passing coords of point
                    st.coords[i]
                        .coords
                        .copy_from_slice(c_array_to_slice(self.handle.get_atom_xyz(i), 3usize));
                }
            }

            // Box is stored as column-major matrix
            let sl = unsafe { std::slice::from_raw_parts(self.handle.get_box(), 9) };
            let m = Matrix3::from_column_slice(sl);
            st.pbox = Some(PeriodicBox::from_matrix(m).map_err(|e| TprHandlerError::Pbc(e))?);

            Ok((top, st.into()))
        }
    }

    impl Drop for TprFileHandler {
        fn drop(&mut self) {
            // Call destructor of C++ helper class
            unsafe { self.handle.destruct() };
        }
    }

    unsafe fn c_ptr_to_str(ptr: *const i8) -> Result<String, TprHandlerError> {
        Ok(CStr::from_ptr(ptr)
            .to_str()
            .map_err(|e| TprHandlerError::CStringUtf8(e))?
            .to_owned())
    }

    fn c_array_to_slice<'a, T, I: TryInto<usize>>(ptr: *mut T, n: I) -> &'a [T] {
        match n.try_into() {
            Ok(sz) => unsafe { std::slice::from_raw_parts(ptr, sz) },
            _ => panic!("Array size is not convertible to usize"),
        }
    }

    #[cfg(test)]
    mod tests {
        use crate::io::TprFileHandler;
        use crate::prelude::*;
        #[test]
        fn test_tpr() {
            let mut h = TprFileHandler::open("tests/topol.tpr").unwrap();
            let (top, st) = h.read().unwrap();
            println!("natoms: {:?}", top.num_atoms());
            println!("nbonds: {:?}", top.num_bonds());
            println!("nmolecules: {:?}", top.num_molecules());
            println!("state sz: {:?}", st.num_coords());
        }
    }
}

#[cfg(not(gromacs))]
mod internal_tpr_disabled {
    use thiserror::Error;
    use crate::core::{State, Topology};

    pub struct TprFileHandler {}

    impl TprFileHandler {
        pub fn open(_fname: &str) -> Result<Self, TprHandlerError> {
            Err(TprHandlerError::GromacsDisabled)
        }

        pub fn read(&mut self) -> Result<(Topology, State), TprHandlerError> {
            Err(TprHandlerError::GromacsDisabled)
        }   
    }

    #[derive(Debug, Error)]
    pub enum TprHandlerError {
        #[error("gromacs support disabled")]
        GromacsDisabled,
    }     
}
