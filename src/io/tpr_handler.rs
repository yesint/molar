use super::{IoFileOpener, IoStructure, IoState};
use crate::core::*;
use ascii::AsciiString;

use crate::io::get_ext;

use anyhow::{anyhow, bail, Result};
use nalgebra::Matrix3;
use std::{ffi::{c_void, CStr, CString}, ptr::null_mut};
use molar_gromacs::gromacs_bindings::{t_topology, TprHelper};

pub struct TprFileHandler {
    handle: TprHelper,
}

impl TprFileHandler {
    fn new(fname: &str) -> Result<Self> {
        let f_name = CString::new(fname.clone())?;
        Ok(TprFileHandler { handle: unsafe{ TprHelper::new(f_name.as_ptr()) } })
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

impl IoStructure for TprFileHandler {
    fn read_structure(&mut self) -> Result<Structure> {
        let top = unsafe{ self.handle.get_top().as_ref().unwrap() };
        let natoms = top.atoms.nr as usize;
        let atoms = unsafe{ std::slice::from_raw_parts(top.atoms.atom, natoms) };
        let atomnames = unsafe{ std::slice::from_raw_parts(top.atoms.atomname, natoms) };
        // pointer to c-strings

        for ptr in atomnames.iter() {   
            // Constructs a c-string from *mut i8 
            let cs = unsafe{ CStr::from_ptr(**ptr) };
            // Shoud print a string, but gives rubbish instead
            println!("{:?}",cs);
        }
        Ok(Default::default())
    }

    fn write_structure(&mut self, data: &Structure) -> Result<()> {
        bail!("TPR files are not writable!")
    }
}

#[test]
fn test_tpr() {
    let mut h = TprFileHandler::new_reader("tests/topol.tpr").unwrap();
    let structure = h.read_structure();
}