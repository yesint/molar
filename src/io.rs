use crate::core::{Structure,State};
use anyhow::{Result,bail,anyhow};
use std::{path::Path, any};

mod vmd_molfile_handler;
pub use vmd_molfile_handler::VmdMolFileHandler;

mod xtc_handler;
pub use xtc_handler::XtcFileHandler;

mod tpr_handler;
pub use tpr_handler::TprFileHandler;

struct FileContent {
    atoms: bool,
    bonds: bool,
    coords: bool,
    pbox: bool,
    
}
// Traits for file handler
pub trait IoFileOpener {
    fn new_reader(fname: &str) -> Result<Self> where Self: Sized;
    fn new_writer(fname: &str) -> Result<Self> where Self: Sized;
}

// Traits for different file types
pub trait IoStructure {
    fn read_structure(&mut self) -> Result<Structure>;
    fn write_structure(&mut self, data: &Structure) -> Result<()>;
}

pub trait IoState {
    fn read_next_state(&mut self) -> Result<Option<State>>;
    fn write_next_state(&mut self, data: &State) -> Result<()>;
}

pub enum FileHandler<'a> {
    Pdb(VmdMolFileHandler<'a>),
    Dcd(VmdMolFileHandler<'a>),
    Xyz(VmdMolFileHandler<'a>),
    Xtc(XtcFileHandler),
}

pub fn get_ext(fname: &str) -> Result<&str> {
    // Get extention
    Ok(
        Path::new(fname)
        .extension().ok_or(anyhow!("File with extension expected, given {fname}"))?
        .to_str().ok_or(anyhow!("Failed getting file extension from {fname}"))?
    )
}


impl<'a> IoFileOpener for FileHandler<'a> {
    fn new_reader(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::new_reader(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::new_reader(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::new_reader(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::new_reader(fname)?)),
            _ => bail!("Unrecognized extension {ext}"),
        }
    }

    fn new_writer(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::new_writer(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::new_writer(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::new_writer(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::new_writer(fname)?)),
            _ => bail!("Unrecognized extension {ext}"),
        }
    }
}

impl<'a> IoStructure for FileHandler<'a> {
    fn read_structure(&mut self) -> Result<Structure> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) => h.read_structure(),
            _ => bail!("Unable to read structure"),
        }
    }

    fn write_structure(&mut self,data: &Structure) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) => h.write_structure(data),
            _ => bail!("Unable to write structure"),
        }
    }
}

impl<'a> IoState for FileHandler<'a> {
    fn read_next_state(&mut self) -> Result<Option<State>> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) | 
            Self::Dcd(ref mut h) => h.read_next_state(),
            Self::Xtc(ref mut h) => h.read_next_state(),
        }
    }

    fn write_next_state(&mut self,data: &State) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) | 
            Self::Dcd(ref mut h) => h.write_next_state(data),
            Self::Xtc(ref mut h) => h.write_next_state(data),
        }
    }
}

