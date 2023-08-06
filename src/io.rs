use crate::core::{structure::Structure,state::State};
use anyhow::Result;

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


pub enum OpenMode {
    Read,
    Write,
    None,
}

// Traits for different file types
pub trait MolfileStructure {
    fn read_structure(&mut self) -> Result<Structure>;
    fn write_structure(&mut self, data: &Structure) -> Result<()>;
}

pub trait MolfileMultiFrame {
    fn read_next_state(&mut self) -> Result<Option<State>>;
    fn write_next_state(&mut self, data: &State) -> Result<()>;
}

pub trait MolfileSingleFrame: MolfileMultiFrame {
    fn read_state(&mut self) -> Result<State>;
    fn write_state(&mut self, data: &State) -> Result<()>;
}
