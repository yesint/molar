use crate::prelude::*;
use gro_handler::GroHandlerError;
use std::path::Path;
use thiserror::Error;
#[cfg(feature = "gromacs")]
use tpr_handler::TprHandlerError;
use triomphe::{Arc, UniqueArc};
use vmd_molfile_handler::VmdHandlerError;
use xtc_handler::XtcHandlerError;
use log::debug;

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
        self.reader
            .read_state()
            .map_err(|e| panic!("error reading next state: {}", e))
            .unwrap()
    }
}

//================================
// General type for file handlers
//================================

pub struct FileHandler {
    file_name: String,
    format_handler: FileFormat,
    stats: FileStats,
    options: FileOptions,
}

enum FileFormat {
    Pdb(VmdMolFileHandler),
    Dcd(VmdMolFileHandler),
    Xyz(VmdMolFileHandler),
    Xtc(XtcFileHandler),
    #[cfg(feature = "gromacs")]
    Tpr(TprFileHandler),
    Gro(GroFileHandler),
}

#[derive(Default)]
struct FileOptions {
    with_title: bool,
}

#[derive(Default)]
struct FileStats {
    total_time: std::time::Duration,
    num_frames: usize,
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
    fn new(handler: FileFormat, file_name: String) -> Self {
        Self {
            format_handler: handler,
            options: Default::default(),
            file_name,
            stats: Default::default(),
        }
    }

    pub fn open(fname: &str) -> Result<Self, FileIoError> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::new(
                FileFormat::Pdb(
                    VmdMolFileHandler::open(fname, VmdMolFileType::Pdb).with_file(|| fname)?,
                ),
                fname.into(),
            )),
            "dcd" => Ok(Self::new(
                FileFormat::Dcd(
                    VmdMolFileHandler::open(fname, VmdMolFileType::Dcd).with_file(|| fname)?,
                ),
                fname.into(),
            )),
            "xyz" => Ok(Self::new(
                FileFormat::Xyz(
                    VmdMolFileHandler::open(fname, VmdMolFileType::Xyz).with_file(|| fname)?,
                ),
                fname.into(),
            )),

            "xtc" => Ok(Self::new(
                FileFormat::Xtc(XtcFileHandler::open(fname).with_file(|| fname)?),
                fname.into(),
            )),

            "gro" => Ok(Self::new(
                FileFormat::Gro(GroFileHandler::open(fname).with_file(|| fname)?),
                fname.into(),
            )),

            #[cfg(feature = "gromacs")]
            "tpr" => Ok(Self::new(
                FileFormat::Tpr(TprFileHandler::open(fname).with_file(|| fname)?),
                fname.into(),
            )),

            #[cfg(not(feature = "gromacs"))]
            "tpr" => Err(FileIoError::TprDisabled),
            _ => Err(FileIoError::NotRecognized(fname.into())),
        }
    }

    pub fn create(fname: &str) -> Result<Self, FileIoError> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::new(
                FileFormat::Pdb(
                    VmdMolFileHandler::create(fname, VmdMolFileType::Pdb).with_file(|| fname)?,
                ),
                fname.into(),
            )),
            "dcd" => Ok(Self::new(
                FileFormat::Dcd(
                    VmdMolFileHandler::create(fname, VmdMolFileType::Dcd).with_file(|| fname)?,
                ),
                fname.into(),
            )),
            "xyz" => Ok(Self::new(
                FileFormat::Xyz(
                    VmdMolFileHandler::create(fname, VmdMolFileType::Xyz).with_file(|| fname)?,
                ),
                fname.into(),
            )),
            "xtc" => Ok(Self::new(
                FileFormat::Xtc(XtcFileHandler::create(fname).with_file(|| fname)?),
                fname.into(),
            )),
            "gro" => Ok(Self::new(
                FileFormat::Gro(GroFileHandler::create(fname).with_file(|| fname)?),
                fname.into(),
            )),
            _ => Err(FileIoError::NotRecognized(fname.into())),
        }
    }

    pub fn read(&mut self) -> Result<(UniqueArc<Topology>, UniqueArc<State>), FileIoError> {
        let t = std::time::Instant::now();

        let (top, st) = match self.format_handler {
            #[cfg(feature = "gromacs")]
            FileFormat::Tpr(ref mut h) => h.read().with_file(|| &self.file_name)?,
            FileFormat::Gro(ref mut h) => h.read().with_file(|| &self.file_name)?,
            FileFormat::Pdb(ref mut h) => {
                let top = h.read_topology().with_file(|| &self.file_name)?;
                let st = h
                    .read_state()
                    .with_file(|| &self.file_name)?
                    .ok_or_else(|| FileIoError::NoStates(self.file_name.clone()))?;
                (top, st)
            }
            _ => return Err(FileIoError::NotReadOnceFormat),
        };
        check_topology_state_sizes(&top, &st)?;

        self.stats.total_time += t.elapsed();
        self.stats.num_frames += 1;

        Ok((top, st))
    }

    pub fn read_shareable(&mut self) -> Result<(Arc<Topology>, Arc<State>), FileIoError> {
        let (top, st) = self.read()?;
        Ok((top.shareable(), st.shareable()))
    }

    pub fn write<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyProvider + StateProvider,
    {
        let t = std::time::Instant::now();

        match self.format_handler {
            FileFormat::Gro(ref mut h) => h.write(data).with_file(|| &self.file_name)?,
            FileFormat::Pdb(ref mut h) => {
                h.write_topology(data).with_file(|| &self.file_name)?;
                h.write_state(data).with_file(|| &self.file_name)?;
            }
            _ => return Err(FileIoError::NotWriteOnceFormat),
        }

        self.stats.total_time += t.elapsed();
        self.stats.num_frames += 1;

        Ok(())
    }

    pub fn read_topology(&mut self) -> Result<UniqueArc<Topology>, FileIoError> {
        let t = std::time::Instant::now();

        let top = match self.format_handler {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) => {
                h.read_topology().with_file(|| &self.file_name)?
            }
            _ => return Err(FileIoError::NotTopologyReadFormat),
        };

        self.stats.total_time += t.elapsed();
        self.stats.num_frames += 1;

        Ok(top)
    }

    pub fn read_topology_shareable(&mut self) -> Result<Arc<Topology>, FileIoError> {
        Ok(self.read_topology()?.shareable())
    }

    pub fn write_topology<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyProvider,
    {
        let t = std::time::Instant::now();

        match self.format_handler {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) => {
                h.write_topology(data).with_file(|| &self.file_name)?
            }
            _ => return Err(FileIoError::NotTopologyWriteFormat),
        }

        self.stats.total_time += t.elapsed();
        self.stats.num_frames += 1;

        Ok(())
    }

    pub fn read_state(&mut self) -> Result<Option<UniqueArc<State>>, FileIoError> {
        let t = std::time::Instant::now();

        let st = match self.format_handler {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) | FileFormat::Dcd(ref mut h) => {
                h.read_state().with_file(|| &self.file_name)?
            }
            FileFormat::Xtc(ref mut h) => h.read_state().with_file(|| &self.file_name)?,
            _ => return Err(FileIoError::NotTrajectoryReadFormat),
        };

        self.stats.total_time += t.elapsed();
        self.stats.num_frames += 1;

        Ok(st)
    }

    pub fn read_state_shareable(&mut self) -> Result<Option<Arc<State>>, FileIoError> {
        Ok(self.read_state()?.map(|v| v.shareable()))
    }

    pub fn write_state<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: StateProvider,
    {
        let t = std::time::Instant::now();

        match self.format_handler {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) | FileFormat::Dcd(ref mut h) => {
                h.write_state(data).with_file(|| &self.file_name)?
            }
            FileFormat::Xtc(ref mut h) => h.write_state(data).with_file(|| &self.file_name)?,
            _ => return Err(FileIoError::NotTrajectoryWriteFormat),
        }

        self.stats.total_time += t.elapsed();
        self.stats.num_frames += 1;

        Ok(())
    }

    pub fn seek_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref mut h) => h.seek_frame(fr).with_file(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn seek_time(&mut self, t: f32) -> Result<(), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref mut h) => h.seek_time(t).with_file(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn tell_first(&self) -> Result<(usize, f32), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref h) => h.tell_first().with_file(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn tell_current(&self) -> Result<(usize, f32), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref h) => h.tell_current().with_file(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }

    pub fn tell_last(&self) -> Result<(usize, f32), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref h) => h.tell_last().with_file(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat),
        }
    }
}

impl Drop for FileHandler {
    fn drop(&mut self) {
        debug!(//target: "FileHandler",
            "Done with file '{}': total IO time {:.4}s, {} frames, {:.4}s per frame",
            self.file_name,
            self.stats.total_time.as_secs_f32(),
            self.stats.num_frames,
            self.stats.total_time.as_secs_f32() / self.stats.num_frames as f32
        );
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

impl TopologyProvider for (UniqueArc<Topology>, UniqueArc<State>) {
    fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }
}

impl AtomsProvider for (UniqueArc<Topology>, UniqueArc<State>) {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.0.iter_atoms()
    }
}

impl StateProvider for (UniqueArc<Topology>, UniqueArc<State>) {
    fn num_coords(&self) -> usize {
        self.1.num_coords()
    }

    fn get_time(&self) -> f32 {
        self.1.get_time()
    }
}

impl BoxProvider for (UniqueArc<Topology>, UniqueArc<State>) {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.1.get_box()
    }
}

impl PosProvider for (UniqueArc<Topology>, UniqueArc<State>) {
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
