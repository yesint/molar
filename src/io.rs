use crate::core::{Structure,State};
use anyhow::{Result,bail};
use std::path::Path;

#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
mod molfile_bindings;

#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
mod xdrfile_bindings;

pub mod vmd_molfile_handler;
pub use vmd_molfile_handler::VmdMolFileHandler;

pub mod xtc_handler;
pub use xtc_handler::XtcFileHandler;


// Traits for file handler
pub trait FileHandler {
    fn new_reader(fname: &str) -> Result<Self> where Self: Sized;
    fn new_writer(fname: &str) -> Result<Self> where Self: Sized;
}

// Traits for different file types
pub trait IoStructure {
    fn read_structure(&mut self) -> Result<Structure>;
    fn write_structure(&mut self, data: &Structure) -> Result<()>;
}

pub trait IoTraj {
    fn read_next_state(&mut self) -> Result<Option<State>>;
    fn write_next_state(&mut self, data: &State) -> Result<()>;
}

pub trait IoSingleFrame {
    fn read_state(&mut self) -> Result<State>;
    fn write_state(&mut self, data: &State) -> Result<()>;
}

fn get_ext(fname: &str) -> &str {
    // Get extention
    Path::new(fname)
    .extension()
    .expect("File with extension expected!")
    .to_str()
    .unwrap()
}

pub fn get_reader(fname: &str) -> Result<Box<dyn FileHandler>> {   
    let ext = get_ext(fname);
    match ext {
        "pdb"|"xyz"|"dcd" => Ok(Box::new(VmdMolFileHandler::new_reader(fname)?)),
        "xtc" => Ok(Box::new(XtcFileHandler::new_reader(fname)?)),
        _ => bail!("Unrecognized extention {ext}!"),
    }    
}

pub fn get_writer(fname: &str) -> Result<Box<dyn FileHandler>> {   
    let ext = get_ext(fname);
    match ext {
        "pdb"|"xyz"|"dcd" => Ok(Box::new(VmdMolFileHandler::new_writer(fname)?)),
        "xtc" => Ok(Box::new(XtcFileHandler::new_writer(fname)?)),
        _ => bail!("Unrecognized extention {ext}!"),
    }    
}
