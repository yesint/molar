use crate::prelude::*;
use gro_handler::GroHandlerError;
use itp_handler::ItpHandlerError;
use log::debug;
use sorted_vec::SortedSet;
use std::{fmt::Display, path::Path};
use thiserror::Error;
#[cfg(feature = "gromacs")]
use tpr_handler::TprHandlerError;
use vmd_molfile_handler::VmdHandlerError;
use xtc_handler::XtcHandlerError;

mod gro_handler;
mod itp_handler;
#[cfg(feature = "gromacs")]
mod tpr_handler;
mod vmd_molfile_handler;
mod xtc_handler;

// Reexports
pub use gro_handler::GroFileHandler;
pub use itp_handler::ItpFileHandler;
#[cfg(feature = "gromacs")]
pub use tpr_handler::TprFileHandler;
pub use vmd_molfile_handler::{VmdMolFileHandler, VmdMolFileType};
pub use xtc_handler::XtcFileHandler;

//--------------------------------------------------------
#[derive(Error, Debug)]
#[error("in file {0}:")]
pub struct FileIoError(String, #[source] FileFormatError);

#[derive(Error, Debug)]
pub enum FileFormatError {
    #[error("in vmd format handler")]
    Vmd(#[from] VmdHandlerError),

    #[error("in gro format handler")]
    Gro(#[from] GroHandlerError),

    #[error("in xtc format handler")]
    Xtc(#[from] XtcHandlerError),

    #[cfg(feature = "gromacs")]
    #[error("in tpr format handler")]
    Tpr(#[from] TprHandlerError),

    #[error("in itp format handler")]
    Itp(#[from] ItpHandlerError),

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

impl IndexProvider for SortedSet<usize> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.iter().cloned()
    }
}

impl IndexProvider for Vec<usize> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.iter().cloned()
    }
}

//=======================================================================
// Iterator over the frames for any type implementing IoTrajectoryReader
//=======================================================================
pub struct IoStateIterator {
    //reader: FileHandler,
    receiver: std::sync::mpsc::Receiver<Option<State>>,
}

impl IoStateIterator {
    fn new(mut fh: FileHandler) -> Self {
        use std::sync::mpsc::sync_channel;
        use std::thread;

        let (sender, receiver) = sync_channel(10);

        // Spawn reading thread
        thread::spawn(move || {
            let mut ok = true;
            while ok {
                let res = fh
                    .read_state()
                    .map_err(|e| panic!("error reading next state: {}", e))
                    .unwrap();
                // If None returned exit reading loop
                if res.is_none() {
                    ok = false;
                }
                // Send state to channel
                sender.send(res).unwrap();
            }
        });

        IoStateIterator { receiver }
    }
}

impl Iterator for IoStateIterator {
    type Item = State;
    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.recv().unwrap()
    }
}

//================================
// General type for file handlers
//================================

pub struct FileHandler {
    pub file_name: String,
    format_handler: FileFormat,
    pub stats: FileStats,
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

impl FileFormat {
    pub fn open(fname: &str) -> Result<Self, FileFormatError> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(FileFormat::Pdb(VmdMolFileHandler::open(
                fname,
                VmdMolFileType::Pdb,
            )?)),
            "dcd" => Ok(FileFormat::Dcd(VmdMolFileHandler::open(
                fname,
                VmdMolFileType::Dcd,
            )?)),
            "xyz" => Ok(FileFormat::Xyz(VmdMolFileHandler::open(
                fname,
                VmdMolFileType::Xyz,
            )?)),

            "xtc" => Ok(FileFormat::Xtc(XtcFileHandler::open(fname)?)),

            "gro" => Ok(FileFormat::Gro(
                GroFileHandler::open(fname)?,
            )),

            #[cfg(feature = "gromacs")]
            "tpr" => Ok(FileFormat::Tpr(TprFileHandler::open(fname)?)),

            #[cfg(not(feature = "gromacs"))]
            "tpr" => Err(FileFormatError::TprDisabled),
            _ => Err(FileFormatError::NotRecognized),
        }
    }

    pub fn create(fname: &str) -> Result<Self, FileFormatError> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(FileFormat::Pdb(VmdMolFileHandler::create(
                fname,
                VmdMolFileType::Pdb,
            )?)),
            "dcd" => Ok(FileFormat::Dcd(VmdMolFileHandler::create(
                fname,
                VmdMolFileType::Dcd,
            )?)),
            "xyz" => Ok(FileFormat::Xyz(VmdMolFileHandler::create(
                fname,
                VmdMolFileType::Xyz,
            )?)),
            "xtc" => Ok(FileFormat::Xtc(XtcFileHandler::create(fname)?)),
            "gro" => Ok(FileFormat::Gro(
                GroFileHandler::create(fname)?,
            )),
            _ => Err(FileFormatError::NotRecognized),
        }
    }

    pub fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        let (top, st) = match self {
            #[cfg(feature = "gromacs")]
            FileFormat::Tpr(ref mut h) => h.read()?,
            FileFormat::Gro(ref mut h) => h.read()?,
            FileFormat::Pdb(ref mut h) => {
                let top = h.read_topology()?;
                let st = h.read_state()?.ok_or_else(|| FileFormatError::NoStates)?;
                (top, st)
            }
            _ => return Err(FileFormatError::NotTopologyAndStateFormat),
        };
        check_topology_state_sizes(&top, &st)?;
        Ok((top, st))
    }

    pub fn write<T>(&mut self, data: &T) -> Result<(), FileFormatError>
    where
        T: TopologyProvider + StateProvider,
    {
        match self {
            FileFormat::Gro(ref mut h) => {
                h.write(data)?;
            }
            FileFormat::Pdb(ref mut h) => {
                h.write_topology(data)?;
                h.write_state(data)?;
            }
            _ => return Err(FileFormatError::NotWriteOnceFormat),
        }
        Ok(())
    }

    pub fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        let top = match self {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) => h.read_topology()?,
            FileFormat::Gro(ref mut h) => {
                let (top, _) = h.read()?;
                top
            }
            _ => return Err(FileFormatError::NotTopologyReadFormat),
        };
        Ok(top)
    }

    pub fn write_topology<T>(&mut self, data: &T) -> Result<(), FileFormatError>
    where
        T: TopologyProvider,
    {
        match self {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) => h.write_topology(data)?,
            _ => return Err(FileFormatError::NotTopologyWriteFormat),
        }
        Ok(())
    }

    pub fn read_state(&mut self) -> Result<Option<State>, FileFormatError> {
        let st = match self {
            FileFormat::Pdb(ref mut h)
            | FileFormat::Xyz(ref mut h)
            | FileFormat::Dcd(ref mut h) => h.read_state()?,
            FileFormat::Xtc(ref mut h) => h.read_state()?,
            FileFormat::Gro(ref mut h) => {
                let (_, st) = h.read()?;
                Some(st)
            }
        };
        Ok(st)
    }

    pub fn write_state<T>(&mut self, data: &T) -> Result<(), FileFormatError>
    where
        T: StateProvider,
    {
        match self {
            FileFormat::Pdb(ref mut h)
            | FileFormat::Xyz(ref mut h)
            | FileFormat::Dcd(ref mut h) => h.write_state(data)?,
            FileFormat::Xtc(ref mut h) => h.write_state(data)?,
            _ => return Err(FileFormatError::NotTrajectoryWriteFormat),
        }
        Ok(())
    }

    pub fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
        match self {
            FileFormat::Xtc(ref mut h) => Ok(h.seek_frame(fr)?),
            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    pub fn seek_time(&mut self, t: f32) -> Result<(), FileFormatError> {
        match self {
            FileFormat::Xtc(ref mut h) => Ok(h.seek_time(t)?),
            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    pub fn tell_first(&self) -> Result<(usize, f32), FileFormatError> {
        match self {
            FileFormat::Xtc(ref h) => Ok(h.tell_first()?),
            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    pub fn tell_current(&self, stats: &FileStats) -> Result<(usize, f32), FileFormatError> {
        match self {
            FileFormat::Xtc(ref h) => Ok(h.tell_current()?),
            // For non-random-access trajectories report FileHandler stats
            _ => Ok((stats.frames_processed, stats.cur_t)),
        }
    }

    pub fn tell_last(&self) -> Result<(usize, f32), FileFormatError> {
        match self {
            FileFormat::Xtc(ref h) => Ok(h.tell_last()?),
            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    /// Consumes frames until reaching serial frame number `fr` (which is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.    
    pub fn skip_to_frame(&mut self, fr: usize, stats: &FileStats) -> Result<(), FileFormatError> {
        // Try random-access first
        match self.seek_frame(fr) {
            Ok(_) => return Ok(()),
            Err(_) => {
                // Not a random access trajectory
                if stats.frames_processed == fr {
                    // We are at needed frame, nothing to do
                    return Ok(());
                } else if stats.frames_processed > fr {
                    // We are past the needed frame, return error
                    return Err(FileFormatError::SeekFrameError(fr));
                } else {
                    // Do serial read until reaching needed frame or EOF
                    while stats.frames_processed < fr {
                        self.read_state()?
                            .ok_or_else(|| FileFormatError::SeekFrameError(fr))?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Consumes frames until reaching beyond time `t` (frame with time exactly equal to `t` is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.    
    pub fn skip_to_time(&mut self, t: f32, stats: &FileStats) -> Result<(), FileFormatError> {
        // Try random-access first
        match self.seek_time(t) {
            Ok(_) => return Ok(()),
            Err(_) => {
                // Not a random access trajectory
                if stats.cur_t > t {
                    // We are past the needed time, return error
                    return Err(FileFormatError::SeekTimeError(t));
                } else {
                    // Do serial read until reaching needed time or EOF
                    while stats.cur_t < t {
                        self.read_state()?
                            .ok_or_else(|| FileFormatError::SeekTimeError(t))?;
                    }
                }
            }
        }
        Ok(())
    }
}

#[derive(Default, Debug, Clone)]
pub struct FileStats {
    pub elapsed_time: std::time::Duration,
    pub frames_processed: usize,
    pub cur_t: f32,
}

impl Display for FileStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "IO time {:.4}s, {} frames, {:.4}s per frame",
            self.elapsed_time.as_secs_f32(),
            self.frames_processed,
            self.elapsed_time.as_secs_f32() / self.frames_processed as f32
        )
    }
}

pub fn get_ext(fname: &str) -> Result<&str, FileFormatError> {
    // Get extention
    Ok(Path::new(fname)
        .extension()
        .ok_or_else(|| FileFormatError::NoExtension)?
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
        Ok(Self::new(
            FileFormat::open(fname).map_err(|e| FileIoError(fname.to_owned(), e))?,
            fname.into(),
        ))
    }

    pub fn create(fname: &str) -> Result<Self, FileIoError> {
        Ok(Self::new(
            FileFormat::create(fname).map_err(|e| FileIoError(fname.to_owned(), e))?,
            fname.into(),
        ))
    }

    pub fn read(&mut self) -> Result<(Topology, State), FileIoError> {
        let t = std::time::Instant::now();

        let (top, st) = self
            .format_handler
            .read()
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok((top, st))
    }

    pub fn write<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyProvider + StateProvider,
    {
        let t = std::time::Instant::now();

        self.format_handler
            .write(data)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    pub fn read_topology(&mut self) -> Result<Topology, FileIoError> {
        let t = std::time::Instant::now();

        let top = self
            .format_handler
            .read_topology()
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(top)
    }

    pub fn write_topology<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyProvider,
    {
        let t = std::time::Instant::now();

        self.format_handler
            .write_topology(data)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    pub fn read_state(&mut self) -> Result<Option<State>, FileIoError> {
        let t = std::time::Instant::now();

        let st = self
            .format_handler
            .read_state()
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        // Update stats if not None
        if st.is_some() {
            self.stats.frames_processed += 1;
            self.stats.cur_t = st.as_ref().unwrap().get_time();
        }

        Ok(st)
    }

    pub fn write_state<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: StateProvider,
    {
        let t = std::time::Instant::now();

        self.format_handler
            .write_state(data)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    pub fn seek_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        Ok(self
            .format_handler
            .seek_frame(fr)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?)
    }

    pub fn seek_time(&mut self, t: f32) -> Result<(), FileIoError> {
        Ok(self
            .format_handler
            .seek_time(t)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?)
    }

    pub fn tell_first(&self) -> Result<(usize, f32), FileIoError> {
        Ok(self
            .format_handler
            .tell_first()
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?)
    }

    pub fn tell_current(&self) -> Result<(usize, f32), FileIoError> {
        Ok(self
            .format_handler
            .tell_current(&self.stats)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?)
    }

    pub fn tell_last(&self) -> Result<(usize, f32), FileIoError> {
        Ok(self
            .format_handler
            .tell_last()
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?)
    }

    /// Consumes frames until reaching serial frame number `fr` (which is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.    
    pub fn skip_to_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        // Try random-access first
        Ok(self
            .format_handler
            .skip_to_frame(fr, &self.stats)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?)
    }

    /// Consumes frames until reaching beyond time `t` (frame with time exactly equal to `t` is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.    
    pub fn skip_to_time(&mut self, t: f32) -> Result<(), FileIoError> {
        // Try random-access first
        Ok(self
            .format_handler
            .skip_to_time(t, &self.stats)
            .map_err(|e| FileIoError(self.file_name.to_owned(), e))?)
    }
}

impl Drop for FileHandler {
    fn drop(&mut self) {
        debug!("Done with file '{}': {}", self.file_name, self.stats);
    }
}

impl IntoIterator for FileHandler {
    type Item = State;
    type IntoIter = IoStateIterator;

    fn into_iter(self) -> Self::IntoIter {
        IoStateIterator::new(self)
    }
}

//----------------------------------------
// Implementation of IO traits for tuples
macro_rules! impl_io_traits_for_tuples {
    ( $t:ty, $s:ty ) => {
        impl TopologyProvider for ($t, $s) {
            fn num_atoms(&self) -> usize {
                self.0.num_atoms()
            }
        }

        impl AtomsProvider for ($t, $s) {
            fn iter_atoms(&self) -> impl AtomIterator<'_> {
                self.0.iter_atoms()
            }
        }

        impl StateProvider for ($t, $s) {
            fn num_coords(&self) -> usize {
                self.1.num_coords()
            }

            fn get_time(&self) -> f32 {
                self.1.get_time()
            }
        }

        impl BoxProvider for ($t, $s) {
            fn get_box(&self) -> Option<&PeriodicBox> {
                self.1.get_box()
            }
        }

        impl PosProvider for ($t, $s) {
            fn iter_pos(&self) -> impl PosIterator<'_> {
                self.1.iter_pos()
            }
        }
    };
}

impl_io_traits_for_tuples!(Topology, State);
impl_io_traits_for_tuples!(Holder<Topology,BuilderSerial>,Holder<State,BuilderSerial>);

//----------------------------------------

#[cfg(test)]
mod tests {
    use super::FileHandler;
    use crate::prelude::*;
    use anyhow::Result;

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
    fn test_into_iter() -> Result<()> {
        let it = FileHandler::open("tests/protein.xtc")?.into_iter();
        for st in it {
            println!("{}", st.get_time());
        }
        Ok(())
    }

    #[test]
    fn test_pdb() -> Result<()> {
        let mut r = FileHandler::open("tests/protein.pdb")?;
        let top1 = r.read_topology()?;
        let st1 = r.read_state()?.unwrap();
        let st2 = st1.clone();
        println!("#1: {}", top1.num_atoms());

        let b = Source::new_serial(top1.into(), st2.into())?;
        let sel = b.select_all()?;
        sel.rotate(&Vector3f::x_axis(), 45.0_f32.to_radians());

        let outname = concat!(env!("OUT_DIR"), "/2.pdb");
        println!("{outname}");
        Ok(())
    }
}
