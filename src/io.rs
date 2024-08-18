use crate::prelude::*;
use gro_handler::GroHandlerError;
use triomphe::{Arc, UniqueArc};
use std::path::Path;
use thiserror::Error;
#[cfg(feature = "gromacs")]
use tpr_handler::TprHandlerError;
use vmd_molfile_handler::VmdHandlerError;
use xtc_handler::XtcHandlerError;

mod gro_handler;
#[cfg(feature = "gromacs")]
mod tpr_handler;
mod vmd_molfile_handler;
mod xtc_handler;

// Reexports
pub use gro_handler::GroFileHandler;
#[cfg(feature = "gromacs")]
pub use tpr_handler::TprFileHandler;
pub use vmd_molfile_handler::{VmdMolFileHandler, VmdMolFileType};
pub use xtc_handler::XtcFileHandler;

#[derive(Error, Debug)]
pub enum FileIoError {
    #[error("file {0}")]
    FileHandler(String, #[source] FileHandlerError),

    #[error("file {0} has no extension")]
    NoExtension(String),

    #[error("file {0} is not recognized")]
    NotRecognized(String),

    #[error("file {0} has no states to read")]
    NoStates(String),

    #[error("not a read once format")]
    NotReadOnceFormat,

    #[error("not a write once format")]
    NotWriteOnceFormat,

    #[error("not a topology reading format")]
    NotTopologyReadFormat,

    #[error("not a topology writing format")]
    NotTopologyWriteFormat,

    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizes),

    #[error("not a trajectory reading format")]
    NotTrajectoryReadFormat,

    #[error("not a trajectory write format")]
    NotTrajectoryWriteFormat,

    #[error("not a random access format")]
    NotRandomAccessFormat,

    #[error("Enable gromacs feature in Cargo.toml to get tpr support")]
    TprDisabled,
}

#[derive(Error, Debug)]
pub enum FileHandlerError {
    #[error(transparent)]
    Vmd(#[from] VmdHandlerError),

    #[error(transparent)]
    Gro(#[from] GroHandlerError),

    #[error(transparent)]
    Xtc(#[from] XtcHandlerError),

    #[cfg(feature = "gromacs")]
    #[error(transparent)]
    Tpr(#[from] TprHandlerError),
}

//-------------------------------------------------------
// Machinery for providing a filename context to errors
trait AddFileContext<T> {
    fn with_file<'a>(self, func: impl FnOnce() -> &'a str) -> std::result::Result<T, FileIoError>;
}

impl<T, E> AddFileContext<T> for std::result::Result<T, E>
where
    FileHandlerError: From<E>,
{
    fn with_file<'a>(self, func: impl FnOnce() -> &'a str) -> std::result::Result<T, FileIoError> {
        self.map_err(|e| FileIoError::FileHandler(func().into(), e.into()))
    }
}
//-------------------------------------------------------

//===============================
// Traits for file opening
//===============================

// There are the following types of data file handlers:
// (1)  All-in-once files (GRO, TPR)
//      Topology+State is read and written at once
// (2)  State only multiple times (XTC, TRR)
//      These are classical trajectories
// (3)  Topology once + multiple states (PDB, TNG)

//===============================
// Traits for writing
//===============================
pub trait IndexProvider {
    fn iter_index(&self) -> impl Iterator<Item = usize>;
}

pub trait TopologyProvider: AtomsProvider {
    fn num_atoms(&self) -> usize;
}

pub trait StateProvider: PosProvider + BoxProvider {
    fn get_time(&self) -> f32;
    fn num_coords(&self) -> usize;
}

pub trait WritableToFile: TopologyProvider + StateProvider
where
    Self: Sized,
{
    fn save(&self, fname: &str) -> Result<(), FileIoError> {
        let mut h = FileHandler::create(fname)?;
        h.write(self)
    }
}

//=======================================================================
// Iterator over the frames for any type implementing IoTrajectoryReader
//=======================================================================
pub struct IoStateIterator {
    reader: FileHandler,
}

impl Iterator for IoStateIterator {
    type Item = UniqueArc<State>;
    fn next(&mut self) -> Option<Self::Item> {
        self.reader.read_state()
        .map_err(|e| panic!("error reading next state: {}",e)).unwrap()
    }
}

//================================
// General type for file handlers
//================================

pub enum FileHandler {
    Pdb(VmdMolFileHandler),
    Dcd(VmdMolFileHandler),
    Xyz(VmdMolFileHandler),
    Xtc(XtcFileHandler),
    #[cfg(feature = "gromacs")]
    Tpr(TprFileHandler),
    Gro(GroFileHandler),
}

pub fn get_ext(fname: &str) -> Result<&str, FileIoError> {
    // Get extention
    Ok(Path::new(fname)
        .extension()
        .ok_or_else(|| FileIoError::NoExtension(fname.into()))?
        .to_str()
        .unwrap())
}

impl FileHandler {
    pub fn open(fname: &str) -> Result<Self, FileIoError> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(
                VmdMolFileHandler::open(fname, VmdMolFileType::Pdb).with_file(|| fname)?,
            )),
            "dcd" => Ok(Self::Dcd(
                VmdMolFileHandler::open(fname, VmdMolFileType::Dcd).with_file(|| fname)?,
            )),
            "xyz" => Ok(Self::Xyz(
                VmdMolFileHandler::open(fname, VmdMolFileType::Xyz).with_file(|| fname)?,
            )),

            "xtc" => Ok(Self::Xtc(XtcFileHandler::open(fname).with_file(|| fname)?)),

            "gro" => Ok(Self::Gro(GroFileHandler::open(fname).with_file(|| fname)?)),

            #[cfg(feature = "gromacs")]
            "tpr" => Ok(Self::Tpr(TprFileHandler::open(fname).with_file(|| fname)?)),

            #[cfg(not(feature = "gromacs"))]
            "tpr" => Err(FileIoError::TprDisabled),
            _ => Err(FileIoError::NotRecognized(fname.into())),
        }
    }

    pub fn create(fname: &str) -> Result<Self, FileIoError> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(
                VmdMolFileHandler::create(fname, VmdMolFileType::Pdb).with_file(|| fname)?,
            )),
            "dcd" => Ok(Self::Dcd(
                VmdMolFileHandler::create(fname, VmdMolFileType::Dcd).with_file(|| fname)?,
            )),
            "xyz" => Ok(Self::Xyz(
                VmdMolFileHandler::create(fname, VmdMolFileType::Xyz).with_file(|| fname)?,
            )),
            "xtc" => Ok(Self::Xtc(
                XtcFileHandler::create(fname).with_file(|| fname)?,
            )),
            "gro" => Ok(Self::Gro(
                GroFileHandler::create(fname).with_file(|| fname)?,
            )),
            _ => Err(FileIoError::NotRecognized(fname.into())),
        }
    }

    pub fn read(&mut self) -> Result<(UniqueArc<Topology>, UniqueArc<State>), FileIoError> {
        let (top, st) = match self {
            #[cfg(feature = "gromacs")]
            Self::Tpr(ref mut h) => h.read().with_file(|| h.get_file_name())?,
            Self::Gro(ref mut h) => h.read().with_file(|| h.get_file_name())?,
            Self::Pdb(ref mut h) => {
                let top = h.read_topology().with_file(|| h.get_file_name())?;
                let st = h
                    .read_state()
                    .with_file(|| h.get_file_name())?
                    .ok_or_else(|| FileIoError::NoStates(h.get_file_name().into()))?;
                (top, st)
            }
            _ => return Err(FileIoError::NotReadOnceFormat),
        };
        check_topology_state_sizes(&top, &st)?;
        Ok((top, st))
    }

    pub fn read_shareable(&mut self) -> Result<(Arc<Topology>, Arc<State>), FileIoError> {
        let (top,st) = self.read()?;
        Ok((top.shareable(), st.shareable()))
    }

    pub fn write<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyProvider + StateProvider,
    {
        match self {
            Self::Gro(ref mut h) => Ok(h.write(data).with_file(|| h.get_file_name())?),
            Self::Pdb(ref mut h) => {
                h.write_topology(data).with_file(|| h.get_file_name())?;
                h.write_state(data).with_file(|| h.get_file_name())?;
                Ok(())
            }
            _ => Err(FileIoError::NotWriteOnceFormat),
        }
    }

    pub fn read_topology(&mut self) -> Result<UniqueArc<Topology>, FileIoError> {
        let top = match self {
              Self::Pdb(ref mut h) 
            | Self::Xyz(ref mut h) => {
                h.read_topology().with_file(|| h.get_file_name())?
            },
            _ => return Err(FileIoError::NotTopologyReadFormat),
        };
        Ok(top)
    }

    pub fn read_topology_shareable(&mut self) -> Result<Arc<Topology>, FileIoError> {
        Ok(self.read_topology()?.shareable())
    }

    pub fn write_topology<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyProvider,
    {
        match self {
              Self::Pdb(ref mut h) 
            | Self::Xyz(ref mut h) => {
                h.write_topology(data).with_file(|| h.get_file_name())
            },
            _ => return Err(FileIoError::NotTopologyWriteFormat),
        }
    }

    pub fn read_state(&mut self) -> Result<Option<UniqueArc<State>>, FileIoError> {
        let st = match self {
              Self::Pdb(ref mut h) 
            | Self::Xyz(ref mut h) 
            | Self::Dcd(ref mut h) => {
                h.read_state().with_file(|| h.get_file_name())?
            },
            Self::Xtc(ref mut h) => h.read_state().with_file(|| h.get_file_name())?,
            _ => return Err(FileIoError::NotTrajectoryReadFormat),
        };
        Ok(st)
    }

    pub fn read_state_shareable(&mut self) -> Result<Option<Arc<State>>, FileIoError> {
        Ok(self.read_state()?.map(|v| v.shareable()))
    }

    pub fn write_state<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: StateProvider,
    {
        match self {
              Self::Pdb(ref mut h) 
            | Self::Xyz(ref mut h) 
            | Self::Dcd(ref mut h) => {
                h.write_state(data).with_file(|| h.get_file_name())
            }
            Self::Xtc(ref mut h) => h.write_state(data).with_file(|| h.get_file_name()),
            _ => return Err(FileIoError::NotTrajectoryWriteFormat),
        }
    }

    pub fn seek_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        match self {
            Self::Xtc(ref mut h) => h.seek_frame(fr).with_file(|| h.get_file_name()),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn seek_time(&mut self, t: f32) -> Result<(), FileIoError> {
        match self {
            Self::Xtc(ref mut h) => h.seek_time(t).with_file(|| h.get_file_name()),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn tell_first(&self) -> Result<(usize, f32), FileIoError> {
        match self {
            Self::Xtc(ref h) => h.tell_first().with_file(|| h.get_file_name()),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn tell_current(&self) -> Result<(usize, f32), FileIoError> {
        match self {
            Self::Xtc(ref h) => h.tell_current().with_file(|| h.get_file_name()),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn tell_last(&self) -> Result<(usize, f32), FileIoError> {
        match self {
            Self::Xtc(ref h) => h.tell_last().with_file(|| h.get_file_name()),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }
}

impl IntoIterator for FileHandler {
    type Item = UniqueArc<State>;
    type IntoIter = IoStateIterator;

    fn into_iter(self) -> Self::IntoIter {
        IoStateIterator { reader: self }
    }
}

//----------------------------------------
// Implementation of IO traits for tuples

impl TopologyProvider for (UniqueArc<Topology>,UniqueArc<State>) {
    fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }    
}

impl AtomsProvider for (UniqueArc<Topology>,UniqueArc<State>) {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.0.iter_atoms()
    }
}

impl StateProvider for (UniqueArc<Topology>,UniqueArc<State>) {
    fn num_coords(&self) -> usize {
        self.1.num_coords()
    }

    fn get_time(&self) -> f32 {
        self.1.get_time()
    }
}

impl BoxProvider for (UniqueArc<Topology>,UniqueArc<State>) {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.1.get_box()
    }
}

impl PosProvider for (UniqueArc<Topology>,UniqueArc<State>) {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.1.iter_pos()
    }
}

//----------------------------------------

#[cfg(test)]
mod tests {
    use super::FileHandler;
    use crate::prelude::*;
    use anyhow::Result;
    use triomphe::UniqueArc;

    #[test]
    fn test_read() -> Result<()> {
        let r = FileHandler::open("tests/protein.xtc")?;
        let mut w = FileHandler::create(concat!(env!("OUT_DIR"), "/1.xtc"))?;

        for fr in r {
            w.write_state(&fr)?;
        }

        Ok(())
    }

    #[test]
    fn test_traj() -> Result<()> {
        let mut r = FileHandler::open("tests/protein.xtc")?;
        let (max_fr, max_t) = r.tell_last()?;
        println!("max: {max_fr}:{max_t}");

        let (cur_fr, cur_t) = r.tell_current()?;
        println!("cur: {cur_fr}:{cur_t}");

        r.seek_frame(2000)?;
        let (cur_fr, cur_t) = r.tell_current()?;
        println!("cur after seek to fr 2000: {cur_fr}:{cur_t}");
        Ok(())
    }

    #[test]
    fn test_pdb() -> Result<()> {
        let mut r = FileHandler::open("tests/protein.pdb")?;
        let top1 = r.read_topology()?;
        let st1 = r.read_state()?.unwrap();
        let st2 = UniqueArc::new(st1.clone());
        println!("#1: {}", top1.num_atoms());

        let b = Source::new(top1, st2)?;
        let sel = b.select_all()?;
        sel.rotate(&Vector3f::x_axis(), 45.0_f32.to_radians());

        let outname = concat!(env!("OUT_DIR"), "/2.pdb");
        println!("{outname}");
        Ok(())
    }
}
