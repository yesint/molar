use crate::core::{Structure,State};
use anyhow::{Result,bail,anyhow};
use std::path::Path;

mod vmd_molfile_handler;
pub use vmd_molfile_handler::VmdMolFileHandler;

mod xtc_handler;
pub use xtc_handler::XtcFileHandler;

mod tpr_handler;
pub use tpr_handler::TprFileHandler;

// Traits for file handler
pub trait IoReader {
    fn new_reader(fname: &str) -> Result<Self> where Self: Sized;
}

pub trait IoWriter {
    fn new_writer(fname: &str) -> Result<Self> where Self: Sized;
}

// Traits for different file types
pub trait IoStructureReader {
    fn read_structure(&mut self) -> Result<Structure>;
}

pub trait IoStructureWriter {
    fn write_structure(&mut self, data: &Structure) -> Result<()>;
}


pub trait IoStateReader {
    fn read_next_state(&mut self) -> Result<Option<State>>;
}

impl<'a> Iterator for FileHandler<'a> {
    type Item = State;
    fn next(&mut self) -> Option<Self::Item> {
        self.read_next_state().expect("Error reading state")
    }
}

pub trait IoStateWriter {
    fn write_next_state(&mut self, data: &State) -> Result<()>;
}

pub enum FileHandler<'a> {
    Pdb(VmdMolFileHandler<'a>),
    Dcd(VmdMolFileHandler<'a>),
    Xyz(VmdMolFileHandler<'a>),
    Xtc(XtcFileHandler),
    Tpr(TprFileHandler),
}

pub fn get_ext(fname: &str) -> Result<&str> {
    // Get extention
    Ok(
        Path::new(fname)
        .extension().ok_or(anyhow!("File with extension expected, given {fname}"))?
        .to_str().ok_or(anyhow!("Failed getting file extension from {fname}"))?
    )
}


impl<'a> IoReader for FileHandler<'a> {
    fn new_reader(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::new_reader(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::new_reader(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::new_reader(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::new_reader(fname)?)),
            "tpr" => Ok(Self::Tpr(TprFileHandler::new_reader(fname)?)),
            _ => bail!("Unrecognized extension {ext}"),
        }
    }
}

impl<'a> IoWriter for FileHandler<'a> {
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

impl<'a> IoStructureReader for FileHandler<'a> {
    fn read_structure(&mut self) -> Result<Structure> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) => h.read_structure(),
            _ => bail!("Unable to read structure"),
        }
    }
}

impl<'a> IoStructureWriter for FileHandler<'a> {
    fn write_structure(&mut self,data: &Structure) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) => h.write_structure(data),
            _ => bail!("Unable to write structure"),
        }
    }
}

impl<'a> IoStateReader for FileHandler<'a> {
    fn read_next_state(&mut self) -> Result<Option<State>> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) | 
            Self::Dcd(ref mut h) => h.read_next_state(),
            Self::Xtc(ref mut h) => h.read_next_state(),
            Self::Tpr(ref mut h) => h.read_next_state(),
        }
    }
}

impl<'a> IoStateWriter for FileHandler<'a> {
    fn write_next_state(&mut self,data: &State) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) |
            Self::Xyz(ref mut h) | 
            Self::Dcd(ref mut h) => h.write_next_state(data),
            Self::Xtc(ref mut h) => h.write_next_state(data),
            _ => bail!("Unable to write state"),
        }
    }
}

#[test]
fn test_read() {
    use super::io::*;

    let mut h = FileHandler::new_reader("tests/topol.tpr").unwrap();
    let st = h.read_structure().unwrap();
    //let fr = h.read_next_state();
    for fr in h {
        println!("{:?}",fr);
    }
}