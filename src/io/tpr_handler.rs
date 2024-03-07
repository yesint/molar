use crate::core::*;
use ascii::{AsciiChar, AsciiString};

use anyhow::Result;
use nalgebra::Matrix3;

use molar_gromacs::gromacs_bindings::*;
use std::{
    ffi::{CStr, CString},
    ptr::null_mut,
};

use super::io_splitter::ReadTopAndState;

pub struct TprFileHandler {
    handle: TprHelper,
}

impl TprFileHandler {
    fn new(fname: &str) -> Result<Self> {
        let f_name = CString::new(fname.to_owned())?;
        Ok(TprFileHandler {
            handle: unsafe { TprHelper::new(f_name.as_ptr()) },
        })
    }

    pub fn open(fname: &str) -> Result<Self> {
        TprFileHandler::new(fname)
    }
}

impl ReadTopAndState for TprFileHandler {
    fn read_top_and_state(&mut self) -> Result<(Topology, State)> {
        //================
        // Read top
        //================
        let gmx_top = unsafe { self.handle.get_top().as_ref().unwrap() };
        let natoms = gmx_top.atoms.nr as usize;
        let nres = gmx_top.atoms.nres as usize;

        let mut top = TopologyStorage::default();
        top.atoms.reserve(natoms);

        let gmx_atoms = c_array_to_slice(gmx_top.atoms.atom, natoms);
        let gmx_atomnames = c_array_to_slice(gmx_top.atoms.atomname, natoms);
        let gmx_resinfo = c_array_to_slice(gmx_top.atoms.resinfo, nres);
        let gmx_pdbinfo: Option<&[t_pdbinfo]> = if gmx_top.atoms.pdbinfo != null_mut() {
            Some(c_array_to_slice(gmx_top.atoms.pdbinfo, natoms))
        } else {
            None
        };
        let gmx_atomtypes = c_array_to_slice(gmx_top.atoms.atomtype, natoms);

        unsafe {
            for i in 0..natoms {
                let name = c_ptr_to_ascii_str(*gmx_atomnames[i]);
                let resi = gmx_atoms[i].resind as usize;
                let resname = c_ptr_to_ascii_str(*gmx_resinfo[resi].name);
                let mut chain = AsciiChar::from_ascii(gmx_resinfo[resi].chainid as u8)
                    .unwrap_or(AsciiChar::Space);
                if chain == AsciiChar::Null {
                    chain = AsciiChar::Space;
                }
                let type_name = c_ptr_to_ascii_str(*gmx_atomtypes[i]);

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
                    occupancy: match gmx_pdbinfo {
                        Some(pdbinfo) => pdbinfo[i].occup,
                        None => 0.0,
                    },
                    bfactor: match gmx_pdbinfo {
                        Some(pdbinfo) => pdbinfo[i].bfac,
                        None => 0.0,
                    },
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
                    F_BONDS | F_G96BONDS | F_HARMONIC | F_FENEBONDS | F_CUBICBONDS | F_CONSTR
                    | F_CONSTRNC => {
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
        let mut top: Topology = top.into();
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
                // We are passinh coords of point
                st.coords[i]
                    .coords
                    .copy_from_slice(c_array_to_slice(self.handle.get_atom_xyz(i), 3usize));
            }
        }

        // Box is stored as column-major matrix
        let sl = unsafe { std::slice::from_raw_parts(self.handle.get_box(), 9) };
        let m = Matrix3::from_column_slice(sl);
        st.pbox = Some(PeriodicBox::from_matrix(m)?);

        Ok((top, st.into()))
    }
}

impl Drop for TprFileHandler {
    fn drop(&mut self) {
        // Call destructor of C++ helper class
        unsafe { self.handle.destruct() };
    }
}

unsafe fn c_ptr_to_ascii_str(ptr: *const i8) -> AsciiString {
    let cstr = CStr::from_ptr(ptr).to_bytes();
    AsciiString::from_ascii_unchecked(cstr)
}

fn c_array_to_slice<'a, T, I: TryInto<usize>>(ptr: *mut T, n: I) -> &'a [T] {
    match n.try_into() {
        Ok(sz) => unsafe { std::slice::from_raw_parts(ptr, sz) },
        _ => panic!("Array size is not convertible to usize"),
    }
}

#[cfg(test)]
mod tests {
    use crate::io::{io_splitter::ReadTopAndState, StateProvider, TopologyProvider, TprFileHandler};
    #[test]
    fn test_tpr() {
        let mut h = TprFileHandler::open("tests/topol.tpr").unwrap();
        let (top, st) = h.read_top_and_state().unwrap();
        println!("natoms: {:?}", top.num_atoms());
        //println!("nbonds: {:?}",top.bonds.len());
        //println!("molecules: {:?}",top.molecules.len());
        println!("state sz: {:?}", st.num_coords());
    }
}
