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

fn c_array_to_slice<'a,T>(ptr: *mut T, n: usize) -> &'a[T] {
    unsafe{ std::slice::from_raw_parts(ptr, n) }
}

#[derive(Clone)]
enum IdefType {
    Bond,
    Angle,
    Dihedral,
    None,
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
            Some(unsafe{std::slice::from_raw_parts(top.atoms.pdbinfo, natoms)})
        } else {
            None
        };
        let gmx_atomtypes = unsafe{ std::slice::from_raw_parts(top.atoms.atomtype, natoms) };
        
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

            // Parsing idef
            let ntypes = top.idef.ntypes as usize;
            let mut idef_types = vec![IdefType::None; ntypes];

            for (i,t) in c_array_to_slice(top.idef.functype, ntypes).iter().enumerate() {
                match *t as u32 {
                    F_BONDS | F_G96BONDS | F_HARMONIC | F_FENEBONDS | F_CUBICBONDS
                        | F_CONSTR |  F_CONSTRNC => idef_types[i]=IdefType::Bond,
                    F_ANGLES | F_G96ANGLES | F_UREY_BRADLEY 
                        | F_LINEAR_ANGLES | F_RESTRANGLES => idef_types[i]=IdefType::Angle,
                    F_PDIHS | F_ANG => idef_types[i]=IdefType::Lj14,
                    _ => idef_types[i]=IdefType::None,
                }
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
    let structure = h.read_structure();
    println!("{:?}",structure);
}