use super::{IoFileOpener, IoStructure, IoState};
use crate::core::*;
use ascii::{AsciiString,AsciiChar};

use crate::io::get_ext;

use anyhow::{anyhow, bail, Result};
use nalgebra::Matrix3;
use std::{ffi::{c_void, CStr, CString}, ptr::null_mut};
use molar_gromacs::gromacs_bindings::*;
use std::collections::HashMap;


pub struct TprFileHandler {
    handle: TprHelper,
}

impl TprFileHandler {
    fn new(fname: &str) -> Result<Self> {
        let f_name = CString::new(fname.clone())?;
        Ok(TprFileHandler { handle: unsafe{ TprHelper::new(f_name.as_ptr()) } })
    }
}

impl Drop for TprFileHandler {
    fn drop(&mut self) {
        // Call destructor of C++ helper class
        unsafe{self.handle.destruct()};
    }
}

impl IoFileOpener for TprFileHandler {
    fn new_reader(fname: &str) -> Result<Self> {
        TprFileHandler::new(fname)        
    }

    fn new_writer(fname: &str) -> Result<Self> {
        bail!("TPR files are not writable!")
    }
}

unsafe fn c_ptr_to_ascii_str(ptr: *const i8) -> AsciiString {
    let cstr = CStr::from_ptr(ptr).to_bytes();
    AsciiString::from_ascii_unchecked(cstr)    
}

fn c_array_to_slice<'a,T,I: TryInto<usize>>(ptr: *mut T, n: I) -> &'a[T] {
    match n.try_into() {
        Ok(sz) => unsafe{ std::slice::from_raw_parts(ptr,sz)  },
        _ => panic!("Array size is not convertible to usize")
    }
}


impl IoStructure for TprFileHandler {
    fn read_structure(&mut self) -> Result<Structure> {
        let top = unsafe{ self.handle.get_top().as_ref().unwrap() };
        let natoms = top.atoms.nr as usize;
        let nres = top.atoms.nres as usize;

        let mut structure: Structure = Default::default();
        structure.atoms.reserve(natoms);

        let gmx_atoms = c_array_to_slice(top.atoms.atom, natoms);
        let gmx_atomnames = c_array_to_slice(top.atoms.atomname, natoms);
        let gmx_resinfo = c_array_to_slice(top.atoms.resinfo, nres);
        let gmx_pdbinfo: Option<&[t_pdbinfo]> =
            if top.atoms.pdbinfo != null_mut() {
                Some(c_array_to_slice(top.atoms.pdbinfo, natoms))
            } else {
                None
            };
        let gmx_atomtypes = c_array_to_slice(top.atoms.atomtype, natoms);
        
        unsafe{
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
                
                let new_atom = Atom{
                    name,
                    resid: gmx_atoms[i].resind,
                    resname,
                    chain,
                    charge: gmx_atoms[i].q,
                    mass: gmx_atoms[i].m,
                    atomic_number: gmx_atoms[i].atomnumber as usize,
                    type_id: gmx_atoms[i].type_ as u32,
                    type_name,
                    occupancy: match gmx_pdbinfo {
                        Some(pdbinfo) => pdbinfo[i].occup,
                        None => 0.0
                    },
                    bfactor: match gmx_pdbinfo {
                        Some(pdbinfo) => pdbinfo[i].bfac,
                        None => 0.0
                    },
                    ..Default::default()
                };
                
                structure.atoms.push(new_atom);
            } //for atoms

            // Parsing idef for bonds
            let functypes = c_array_to_slice(top.idef.functype,top.idef.ntypes);
            // Iterate over non-empty interaction lists 
            for il in top.idef.il.iter().filter(|el| el.nr>0) {
                let iatoms = c_array_to_slice(il.iatoms,il.nr);
                // We can check the first type only since we only
                // need to evaluate the interaction type and not
                // the concrete parameters for involved atoms
                match functypes[iatoms[0] as usize] as u32 {
                    F_BONDS | F_G96BONDS | F_HARMONIC 
                    | F_FENEBONDS | F_CUBICBONDS
                    | F_CONSTR |  F_CONSTRNC => {
                        for el in iatoms.chunks_exact(3) {
                            // el[0] is type and not needed
                            structure.bonds.push([el[1] as usize, el[2] as usize]);
                        }        
                    },
                    F_SETTLE => {
                        for el in iatoms.chunks_exact(4){
                            // el[0] is type and not needed
                            // Each settle is 2 bonds
                            structure.bonds.push([el[1] as usize, el[2] as usize]);
                            structure.bonds.push([el[1] as usize, el[3] as usize]);
                        }
                    }
                    _ => (),
                } //match
                
            } //for

            // Read molecules
            let mol_index = c_array_to_slice(top.mols.index, top.mols.nr);
            for m in mol_index.chunks_exact(2) {
                structure.molecules.push([m[0] as usize, m[1] as usize -1]);
            }


        } //unsafe

        Ok(structure)
    }

    fn write_structure(&mut self, data: &Structure) -> Result<()> {
        bail!("TPR files are not writable!")
    }
}

#[test]
fn test_tpr() {
    let mut h = TprFileHandler::new_reader("tests/topol.tpr").unwrap();
    let structure = h.read_structure().unwrap();
    println!("natoms: {:?}",structure.atoms.len());
    println!("nbonds: {:?}",structure.bonds.len());
    println!("molecules: {:?}",structure.molecules.len());
}