use crate::core::{IndexIterator, State, Topology};
use anyhow::{anyhow, bail, Result};
use std::path::Path;

mod vmd_molfile_handler;
mod xtc_handler;

#[cfg(feature = "gromacs")]
mod tpr_handler;

// Reexports
#[cfg(feature = "gromacs")]
pub use tpr_handler::TprFileHandler;
pub use vmd_molfile_handler::VmdMolFileHandler;
pub use xtc_handler::XtcFileHandler;

//===============================
// Traits for file opening
//===============================
pub trait IoReader {
    fn new_reader(fname: &str) -> Result<Self>
    where
        Self: Sized;
}

pub trait IoWriter {
    fn new_writer(fname: &str) -> Result<Self>
    where
        Self: Sized;
}

//===============================
// Traits for Structure IO
//===============================
pub trait IoTopologyReader: IoReader {
    fn read_topology(&mut self) -> Result<Topology>;
}

pub trait IoTopologyWriter: IoWriter {
    fn write_structure_subset(
        &mut self,
        data: &Topology,
        subset_indexes: impl IndexIterator,
    ) -> Result<()>;

    // Default implementation with all indexes
    fn write_structure(&mut self, data: &Topology) -> Result<()> {
        self.write_structure_subset(data, 0..data.atoms.len())
    }
}

//===============================
// Traits for State IO
//===============================
pub trait IoStateReader: IoReader {
    fn read_next_state(&mut self) -> Result<Option<State>>;

    fn into_states_iter(self) -> IoStateIterator<Self>
    where
        Self: Sized,
    {
        IoStateIterator { reader: self }
    }
}

pub trait IoStateWriter: IoWriter {
    fn write_next_state_subset(
        &mut self,
        data: &State,
        subset_indexes: impl IndexIterator,
    ) -> Result<()>;

    // Default implementation with all indexes
    fn write_next_state(&mut self, data: &State) -> Result<()> {
        self.write_next_state_subset(data, 0..data.coords.len())
    }
}

//==================================================================
// Iterator over the frames for any type implementing IoStateReader
//==================================================================
pub struct IoStateIterator<T>
where
    T: IoStateReader,
{
    reader: T,
}

impl<T> Iterator for IoStateIterator<T>
where
    T: IoStateReader,
{
    type Item = State;
    fn next(&mut self) -> Option<Self::Item> {
        self.reader.read_next_state().expect("Error reading state")
    }
}

//================================
// General type for file handlers
//================================

pub enum FileHandler<'a> {
    Pdb(VmdMolFileHandler<'a>),
    Dcd(VmdMolFileHandler<'a>),
    Xyz(VmdMolFileHandler<'a>),
    Xtc(XtcFileHandler),
    #[cfg(feature = "gromacs")]
    Tpr(TprFileHandler),
}

pub fn get_ext(fname: &str) -> Result<&str> {
    // Get extention
    Ok(Path::new(fname)
        .extension()
        .ok_or(anyhow!("File with extension expected, given {fname}"))?
        .to_str()
        .ok_or(anyhow!("Failed getting file extension from {fname}"))?)
}

pub fn get_ext2(fname: &str) -> Result<&str> {
    // Get extention
    Ok(Path::new(fname)
        .extension()
        .ok_or(anyhow!("File with extension expected, given {fname}"))?
        .to_str()
        .ok_or(anyhow!("Failed getting file extension from {fname}"))?)
}

impl<'a> IoReader for FileHandler<'a> {
    fn new_reader(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::new_reader(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::new_reader(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::new_reader(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::new_reader(fname)?)),
            #[cfg(feature = "gromacs")]
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

impl<'a> IoTopologyReader for FileHandler<'a> {
    fn read_topology(&mut self) -> Result<Topology> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => h.read_topology(),
            #[cfg(feature = "gromacs")]
            Self::Tpr(ref mut h) => h.read_topology(),
            _ => bail!("Unable to read structure"),
        }
    }
}

impl<'a> IoTopologyWriter for FileHandler<'a> {
    fn write_structure_subset(
        &mut self,
        data: &Topology,
        subset_indexes: impl IndexIterator,
    ) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => {
                h.write_structure_subset(data, subset_indexes)
            }
            _ => bail!("Unable to write structure"),
        }
    }
}

impl<'a> IoStateReader for FileHandler<'a> {
    fn read_next_state(&mut self) -> Result<Option<State>> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.read_next_state()
            }
            Self::Xtc(ref mut h) => h.read_next_state(),
            #[cfg(feature = "gromacs")]
            Self::Tpr(ref mut h) => h.read_next_state(),
        }
    }
}

impl<'a> IoStateWriter for FileHandler<'a> {
    fn write_next_state_subset(
        &mut self,
        data: &State,
        subset_indexes: impl IndexIterator,
    ) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.write_next_state_subset(data, subset_indexes)
            }
            Self::Xtc(ref mut h) => h.write_next_state_subset(data, subset_indexes),
            _ => bail!("Unable to write state"),
        }
    }
}

#[test]
fn test_read() {
    use super::io::*;

    let mut r = FileHandler::new_reader("tests/topol.tpr").unwrap();
    let mut w = FileHandler::new_writer(concat!(env!("OUT_DIR"), "/1.pdb")).unwrap();

    let st = r.read_topology().unwrap();
    println!("{:?}", st.atoms);

    for fr in r.into_states_iter() {
        //println!("{:?}",fr);
        w.write_structure(&st).unwrap();
        w.write_next_state(&fr).unwrap();
        //w.write_structure(&st).unwrap();
        //w.write_next_state_subset(&fr,0..10).unwrap();
    }
}
