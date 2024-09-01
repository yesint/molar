use crate::prelude::*;
use gro_handler::GroHandlerError;
use log::debug;
use std::{fmt::Display, path::Path};
use thiserror::Error;
#[cfg(feature = "gromacs")]
use tpr_handler::TprHandlerError;
use triomphe::{Arc, UniqueArc};
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

use thiserror_string_context::*;

//--------------------------------------------------------
#[string_context("in file '{0}'")]
#[derive(Error, Debug)]
pub enum FileIoError {
    #[error("in vmd format handler")]
    Vmd(#[from] VmdHandlerError),

    #[error("in gro format handler")]
    Gro(#[from] GroHandlerError),

    #[error("in xtc format handler")]
    Xtc(#[from] XtcHandlerError),

    #[cfg(feature = "gromacs")]
    #[error("in tpr format handler")]
    Tpr(#[from] TprHandlerError),

    #[error("file has no extension")]
    NoExtension,

    #[error("file is not recognized")]
    NotRecognized,

    #[error("file has no states to read")]
    NoStates,

    #[error("not a format containing both structure and coordinates")]
    NotTopologyAndStateFormat,

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

    #[error("enable gromacs feature in Cargo.toml to get tpr support")]
    TprDisabled,
    
    #[error("unexpected end of file")]
    UnexpectedEof,

    #[error("can't seek to frame {0}")]
    SeekFrameError(usize),

    #[error("can't seek to time {0}")]
    SeekTimeError(f32),
}

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
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize>;
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

#[derive(Default, Debug)]
struct FileStats {
    elapsed_time: std::time::Duration,
    frames_processed: usize,
    cur_t: f32,
}

impl Display for FileStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "total IO time {:.4}s, {} frames, {:.4}s per frame",
            self.elapsed_time.as_secs_f32(),
            self.frames_processed,
            self.elapsed_time.as_secs_f32() / self.frames_processed as f32
        )
    }
}

pub fn get_ext(fname: &str) -> Result<&str, FileIoError> {
    // Get extention
    Ok(Path::new(fname)
        .extension()
        .ok_or_else(|| FileIoError::NoExtension)?
        .to_str()
        .unwrap())
}

impl FileHandler {
    fn new(handler: FileFormat, file_name: String) -> Self {
        Self {
            file_name,
            format_handler: handler,
            stats: Default::default(),
        }
    }

    pub fn open(fname: &str) -> Result<Self, FileIoError> {
        let ext = get_ext(fname).with_context(|| fname)?;
        match ext {
            "pdb" => Ok(Self::new(
                FileFormat::Pdb(
                    VmdMolFileHandler::open(fname, VmdMolFileType::Pdb).with_context(|| fname)?,
                ),
                fname.into(),
            )),
            "dcd" => Ok(Self::new(
                FileFormat::Dcd(
                    VmdMolFileHandler::open(fname, VmdMolFileType::Dcd).with_context(|| fname)?,
                ),
                fname.into(),
            )),
            "xyz" => Ok(Self::new(
                FileFormat::Xyz(
                    VmdMolFileHandler::open(fname, VmdMolFileType::Xyz).with_context(|| fname)?,
                ),
                fname.into(),
            )),

            "xtc" => Ok(Self::new(
                FileFormat::Xtc(XtcFileHandler::open(fname).with_context(|| fname)?),
                fname.into(),
            )),

            "gro" => Ok(Self::new(
                FileFormat::Gro(GroFileHandler::open(fname).with_context(|| fname)?),
                fname.into(),
            )),

            #[cfg(feature = "gromacs")]
            "tpr" => Ok(Self::new(
                FileFormat::Tpr(TprFileHandler::open(fname).with_context(|| fname)?),
                fname.into(),
            )),

            #[cfg(not(feature = "gromacs"))]
            "tpr" => Err(FileIoError::TprDisabled),
            _ => Err(FileIoError::NotRecognized).with_context(|| fname),
        }
    }

    pub fn create(fname: &str) -> Result<Self, FileIoError> {
        let ext = get_ext(fname).with_context(|| fname)?;
        match ext {
            "pdb" => Ok(Self::new(
                FileFormat::Pdb(
                    VmdMolFileHandler::create(fname, VmdMolFileType::Pdb).with_context(|| fname)?,
                ),
                fname.into(),
            )),
            "dcd" => Ok(Self::new(
                FileFormat::Dcd(
                    VmdMolFileHandler::create(fname, VmdMolFileType::Dcd).with_context(|| fname)?,
                ),
                fname.into(),
            )),
            "xyz" => Ok(Self::new(
                FileFormat::Xyz(
                    VmdMolFileHandler::create(fname, VmdMolFileType::Xyz).with_context(|| fname)?,
                ),
                fname.into(),
            )),
            "xtc" => Ok(Self::new(
                FileFormat::Xtc(XtcFileHandler::create(fname).with_context(|| fname)?),
                fname.into(),
            )),
            "gro" => Ok(Self::new(
                FileFormat::Gro(GroFileHandler::create(fname).with_context(|| fname)?),
                fname.into(),
            )),
            _ => Err(FileIoError::NotRecognized).with_context(|| fname),
        }
    }

    pub fn read(&mut self) -> Result<(UniqueArc<Topology>, UniqueArc<State>), FileIoError> {
        let t = std::time::Instant::now();

        let (top, st) = match self.format_handler {
            #[cfg(feature = "gromacs")]
            FileFormat::Tpr(ref mut h) => h.read().with_context(|| &self.file_name)?,
            FileFormat::Gro(ref mut h) => h.read().with_context(|| &self.file_name)?,
            FileFormat::Pdb(ref mut h) => {
                let top = h.read_topology().with_context(|| &self.file_name)?;
                let st = h
                    .read_state()
                    .with_context(|| &self.file_name)?
                    .ok_or_else(|| FileIoError::NoStates).with_context(|| &self.file_name)?;
                (top, st)
            }
            _ => return Err(FileIoError::NotTopologyAndStateFormat).with_context(|| &self.file_name),
        };
        check_topology_state_sizes(&top, &st).with_context(|| &self.file_name)?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

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
            FileFormat::Gro(ref mut h) => h.write(data).with_context(|| &self.file_name)?,
            FileFormat::Pdb(ref mut h) => {
                h.write_topology(data).with_context(|| &self.file_name)?;
                h.write_state(data).with_context(|| &self.file_name)?;
            }
            _ => return Err(FileIoError::NotWriteOnceFormat).with_context(|| &self.file_name),
        }

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    pub fn read_topology(&mut self) -> Result<UniqueArc<Topology>, FileIoError> {
        let t = std::time::Instant::now();

        let top = match self.format_handler {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) => {
                h.read_topology().with_context(|| &self.file_name)?
            }
            _ => return Err(FileIoError::NotTopologyReadFormat).with_context(|| &self.file_name),
        };

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

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
                h.write_topology(data).with_context(|| &self.file_name)?
            }
            _ => return Err(FileIoError::NotTopologyWriteFormat).with_context(|| &self.file_name)
        }

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    pub fn read_state(&mut self) -> Result<Option<UniqueArc<State>>, FileIoError> {
        let t = std::time::Instant::now();

        let st = match self.format_handler {
            FileFormat::Pdb(ref mut h)
            | FileFormat::Xyz(ref mut h)
            | FileFormat::Dcd(ref mut h) => h.read_state().with_context(|| &self.file_name)?,
            FileFormat::Xtc(ref mut h) => h.read_state().with_context(|| &self.file_name)?,
            _ => return Err(FileIoError::NotTrajectoryReadFormat).with_context(|| &self.file_name),
        };

        self.stats.elapsed_time += t.elapsed();        
        // Update stats if not None
        if st.is_some() {
            self.stats.frames_processed += 1;
            self.stats.cur_t = st.as_ref().unwrap().get_time();
        }

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
            FileFormat::Pdb(ref mut h)
            | FileFormat::Xyz(ref mut h)
            | FileFormat::Dcd(ref mut h) => h.write_state(data).with_context(|| &self.file_name)?,
            FileFormat::Xtc(ref mut h) => h.write_state(data).with_context(|| &self.file_name)?,
            _ => return Err(FileIoError::NotTrajectoryWriteFormat).with_context(|| &self.file_name),
        }

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    pub fn seek_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref mut h) => h.seek_frame(fr).with_context(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat).with_context(|| &self.file_name),
        }
    }

    pub fn seek_time(&mut self, t: f32) -> Result<(), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref mut h) => h.seek_time(t).with_context(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat).with_context(|| &self.file_name),
        }
    }

    pub fn tell_first(&self) -> Result<(usize, f32), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref h) => h.tell_first().with_context(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat).with_context(|| &self.file_name),
        }
    }

    pub fn tell_current(&self) -> Result<(usize, f32), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref h) => h.tell_current().with_context(|| &self.file_name),
            // For non-random-access trajectories report FileHandler stats
            _ => Ok((self.stats.frames_processed,self.stats.cur_t)),
        }
    }

    pub fn tell_last(&self) -> Result<(usize, f32), FileIoError> {
        match self.format_handler {
            FileFormat::Xtc(ref h) => h.tell_last().with_context(|| &self.file_name),
            _ => return Err(FileIoError::NotRandomAccessFormat).with_context(|| &self.file_name),
        }
    }

    /// Consumes frames until reaching serial frame number `fr` (which is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.    
    pub fn skip_to_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        // Try random-access first
        match self.seek_frame(fr) {
            Ok(_) => { return Ok(()) },
            Err(_) => {                
                // Not a random access trajectory
                if self.stats.frames_processed == fr {
                    // We are at needed frame, nothing to do
                    return Ok(())
                } else if self.stats.frames_processed > fr {
                    // We are past the needed frame, return error
                    return Err(FileIoError::SeekFrameError(fr)).with_context(|| &self.file_name)
                } else {
                    // Do serial read until reaching needed frame or EOF
                    while self.stats.frames_processed < fr {
                        self.read_state()?.ok_or_else(|| FileIoError::SeekFrameError(fr)).with_context(|| &self.file_name)?;
                    }
                }
            },
        }
        Ok(())
    }

    /// Consumes frames until reaching beyond time `t` (frame with time exactly equal to `t` is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.    
    pub fn skip_to_time(&mut self, t: f32) -> Result<(), FileIoError> {
        // Try random-access first
        match self.seek_time(t) {
            Ok(_) => { return Ok(()) },
            Err(_) => {                
                // Not a random access trajectory
                if self.stats.cur_t > t {                    
                    // We are past the needed time, return error
                    return Err(FileIoError::SeekTimeError(t)).with_context(|| &self.file_name)
                } else {
                    // Do serial read until reaching needed time or EOF
                    while self.stats.cur_t < t {
                        self.read_state()?.ok_or_else(|| FileIoError::SeekTimeError(t)).with_context(|| &self.file_name)?;
                    }
                }
            },
        }
        Ok(())
    }
}

impl Drop for FileHandler {
    fn drop(&mut self) {
        debug!("Done with file '{}': {}", self.file_name, self.stats);
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

        let mut b = Source::new_serial(top1, st2)?;
        let sel = b.select_all()?;
        sel.rotate(&Vector3f::x_axis(), 45.0_f32.to_radians());

        let outname = concat!(env!("OUT_DIR"), "/2.pdb");
        println!("{outname}");
        Ok(())
    }
}
