use super::{IoFileOpener, IoStructure, IoState};
use crate::core::*;
use ascii::AsciiString;

use crate::io::get_ext;

use anyhow::{anyhow, bail, Result};
use nalgebra::Matrix3;
use std::ffi::{c_void, CStr, CString};
use molar_gromacs::gromacs_bindings::{read_tpr,t_topology, t_inputrec};


pub struct TprFileHandler {        
    top: t_topology,
    state: t_state,
    ir: t_inputrec,
    natoms: usize,
}

impl TprFileHandler {
    fn new(fname: &str) -> Result<Self> {
        let f_name = CString::new(fname.clone())?;
        unsafe{ read_tpr(f_name.as_ptr(), ) };
        let natoms = top.atoms.nr as usize;
        Ok(TprFileHandler { top, natoms })
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
        
    }

    fn write_structure(&mut self, data: &Structure) -> Result<()> {
        bail!("TPR files are not writable!")
    }
}