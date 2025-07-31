//! IO handlers for different file formats and associated traits

use crate::prelude::*;
use gro_handler::GroHandlerError;
use itp_handler::ItpHandlerError;
use log::{debug, error, warn};
use triomphe::Arc;
use std::{
    fmt::Display,
    path::{Path, PathBuf},
};
use thiserror::Error;

mod gro_handler;
mod itp_handler;
mod tpr_handler;
mod vmd_molfile_handler;
mod xtc_handler;

use tpr_handler::TprHandlerError;
use vmd_molfile_handler::VmdHandlerError;
use xtc_handler::XtcHandlerError;

use gro_handler::GroFileHandler;
use itp_handler::ItpFileHandler;
use tpr_handler::TprFileHandler;
use vmd_molfile_handler::{VmdMolFileHandler, VmdMolFileType};
use xtc_handler::XtcFileHandler;

//===============================
// Traits for file opening
//===============================

// There are the following types of data file handlers:
// (1)  All-in-once files (GRO, TPR)
//      Topology+State is read and written at once
// (2)  State only multiple times (XTC, TRR)
//      These are classical trajectories
// (3)  Topology once + multiple states (PDB, TNG)

//=======================================================================
// Iterator over the frames for any type implementing IoTrajectoryReader
//=======================================================================
/// Iterator over states in a trajectory file
pub struct IoStateIterator {
    //reader: FileHandler,
    receiver: std::sync::mpsc::Receiver<Result<Option<State>, FileIoError>>,
}

impl IoStateIterator {
    fn new(mut fh: FileHandler) -> Self {
        use std::sync::mpsc::sync_channel;
        let (sender, receiver) = sync_channel(10);

        // Spawn reading thread
        std::thread::spawn(move || {
            let mut terminate = false;
            while !terminate {
                let res = fh.read_state();
                terminate = match res.as_ref() {
                    Err(_) => true,   // terminate if reading failed
                    Ok(None) => true, // terminate if none is returned (trajectory done)
                    _ => false,       // otherwise continut reading
                };

                // Send to channel.
                // An error means that the reciever has closed the channel already.
                // This is fine, just exit in this case.
                if sender.send(res).is_err() {
                    break;
                }
            }
        });

        IoStateIterator { receiver }
    }
}

impl Iterator for IoStateIterator {
    type Item = State;
    fn next(&mut self) -> Option<Self::Item> {
        // Reader thread should never crash since it catches errors and ends on them.
        // If it does then something is horrible anyway, so unwrap is fine here.
        match self.receiver.recv().expect("reader thread shouldn't crash") {
            Ok(opt_st) => opt_st,
            Err(e) => {
                warn!("reader thread can't read state from '{}'. File is likely corrupted, reading stopped.",e.0.display());
                None
            }
        }
    }
}

//================================
// General type for file handlers
//================================

/// A file handler for reading molecular structure and trajectory files.
/// Supports various formats including PDB, DCD, XYZ, XTC, GRO, ITP and TPR.
pub struct FileHandler {
    pub file_path: PathBuf,
    format_handler: FileFormat,
    pub stats: FileStats,
}

enum FileFormat {
    Pdb(VmdMolFileHandler),
    Dcd(VmdMolFileHandler),
    Xyz(VmdMolFileHandler),
    Xtc(XtcFileHandler),
    Tpr(TprFileHandler),
    Gro(GroFileHandler),
    Itp(ItpFileHandler),
}

impl FileFormat {
    pub(crate) fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError> {
        let ext = get_ext(fname.as_ref())?;
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

            "gro" => Ok(FileFormat::Gro(GroFileHandler::open(fname)?)),

            "itp" => Ok(FileFormat::Itp(ItpFileHandler::open(fname)?)),

            "tpr" => Ok(FileFormat::Tpr(TprFileHandler::open(fname)?)),

            _ => Err(FileFormatError::NotReadable),
        }
    }

    pub(crate) fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError> {
        let ext = get_ext(fname.as_ref())?;
        match ext {
            "pdb" | "ent" => Ok(FileFormat::Pdb(VmdMolFileHandler::create(
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

            "gro" => Ok(FileFormat::Gro(GroFileHandler::create(fname)?)),

            _ => Err(FileFormatError::NotWritable),
        }
    }

    pub(crate) fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        let (top, st) = match self {
            FileFormat::Tpr(ref mut h) => h.read()?,

            FileFormat::Gro(ref mut h) => h.read()?.ok_or(FileFormatError::NoStates)?,

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

    pub(crate) fn write<T>(&mut self, data: &T) -> Result<(), FileFormatError>
    where
        T: TopologyIoProvider + StateIoProvider,
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

    pub(crate) fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        let top = match self {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) => h.read_topology()?,

            FileFormat::Gro(ref mut h) => h.read_topology()?,

            FileFormat::Tpr(ref mut h) => {
                let (top, _) = h.read()?;
                top
            }

            FileFormat::Itp(ref mut h) => h.read_topology()?,

            _ => return Err(FileFormatError::NotTopologyReadFormat),
        };
        Ok(top)
    }

    pub(crate) fn write_topology<T>(&mut self, data: &T) -> Result<(), FileFormatError>
    where
        T: TopologyIoProvider,
    {
        match self {
            FileFormat::Pdb(ref mut h) | FileFormat::Xyz(ref mut h) => h.write_topology(data)?,

            _ => return Err(FileFormatError::NotTopologyWriteFormat),
        }
        Ok(())
    }

    pub(crate) fn read_state(&mut self) -> Result<Option<State>, FileFormatError> {
        let st = match self {
            FileFormat::Pdb(ref mut h)
            | FileFormat::Xyz(ref mut h)
            | FileFormat::Dcd(ref mut h) => h.read_state()?,

            FileFormat::Xtc(ref mut h) => h.read_state()?,

            FileFormat::Gro(ref mut h) => h.read_state()?,

            FileFormat::Tpr(ref mut h) => {
                let (_, st) = h.read()?;
                Some(st)
            }
            
            _ => return Err(FileFormatError::NotStateReadFormat),
        };
        Ok(st)
    }

    pub(crate) fn write_state<T>(&mut self, data: &T) -> Result<(), FileFormatError>
    where
        T: StateIoProvider,
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

    pub(crate) fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
        match self {
            FileFormat::Xtc(ref mut h) => Ok(h.seek_frame(fr)?),

            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    pub(crate) fn seek_time(&mut self, t: f32) -> Result<(), FileFormatError> {
        match self {
            FileFormat::Xtc(ref mut h) => Ok(h.seek_time(t)?),

            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    pub(crate) fn tell_first(&self) -> Result<(usize, f32), FileFormatError> {
        match self {
            FileFormat::Xtc(ref h) => Ok(h.tell_first()?),

            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    pub(crate) fn tell_current(&self, stats: &FileStats) -> Result<(usize, f32), FileFormatError> {
        match self {
            FileFormat::Xtc(ref h) => Ok(h.tell_current()?),

            // For non-random-access trajectories report FileHandler stats
            _ => Ok((stats.frames_processed, stats.cur_t)),
        }
    }

    pub(crate) fn tell_last(&self) -> Result<(usize, f32), FileFormatError> {
        match self {
            FileFormat::Xtc(ref h) => Ok(h.tell_last()?),

            _ => return Err(FileFormatError::NotRandomAccessFormat),
        }
    }

    /// Consumes frames until reaching serial frame number `fr` (which is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.
    pub(crate) fn skip_to_frame(
        &mut self,
        fr: usize,
        stats: &FileStats,
    ) -> Result<(), FileFormatError> {
        // Try random-access first
        self.seek_frame(fr).or_else(|_| {
            // Not a random access trajectory
            if stats.frames_processed == fr {
                // We are at needed frame, nothing to do
                Ok(())
            } else if stats.frames_processed > fr {
                // We are past the needed frame, return error
                Err(FileFormatError::SeekFrameError(fr))
            } else {
                // Do serial read until reaching needed frame or EOF
                while stats.frames_processed < fr {
                    self.read_state()?
                        .ok_or_else(|| FileFormatError::SeekFrameError(fr))?;
                }
                Ok(())
            }
        })
    }

    /// Consumes frames until reaching beyond time `t` (frame with time exactly equal to `t` is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.
    pub(crate) fn skip_to_time(
        &mut self,
        t: f32,
        stats: &FileStats,
    ) -> Result<(), FileFormatError> {
        // Try random-access first
        self.seek_time(t).or_else(|_| {
            // Not a random access trajectory
            if stats.cur_t > t {
                // We are past the needed time, return error
                Err(FileFormatError::SeekTimeError(t))
            } else {
                // Do serial read until reaching needed time or EOF
                while stats.cur_t < t {
                    self.read_state()?
                        .ok_or_else(|| FileFormatError::SeekTimeError(t))?;
                }
                Ok(())
            }
        })
    }
}

/// Statistics about file operations including elapsed time and number of processed frames
#[derive(Default, Debug, Clone)]
pub struct FileStats {
    /// Total time spent on IO operations
    pub elapsed_time: std::time::Duration,
    /// Number of frames processed
    pub frames_processed: usize,
    /// Current time in the trajectory
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

fn get_ext(fname: &Path) -> Result<&str, FileFormatError> {
    // Get extention
    Ok(fname
        .extension()
        .ok_or_else(|| FileFormatError::NoExtension)?
        .to_str()
        .unwrap())
}

//------------------------------------------------------------------

impl FileHandler {
    /// Opens a file for reading. Format is determined by extension.
    ///
    /// Supported formats:
    /// - pdb,ent: Protein Data Bank structure format
    /// - dcd: DCD trajectory format
    /// - xyz: XYZ format
    /// - xtc: GROMACS compressed trajectory format
    /// - gro: GROMACS structure format
    /// - itp: GROMACS topology format (read-only)
    /// - tpr: GROMACS run input format (read-only)
    ///
    /// # Errors
    /// Returns [FileIoError] if:
    /// - File extension is not recognized
    /// - File cannot be opened
    pub fn open(fname: impl AsRef<Path>) -> Result<Self, FileIoError> {
        let fname = fname.as_ref();
        Ok(Self {
            file_path: fname.to_path_buf(),
            format_handler: FileFormat::open(fname)
                .map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            stats: Default::default(),
        })
    }

    /// Creates a new file for writing. Format is determined by extension.
    ///
    /// Supported formats:
    /// - pdb: Protein Data Bank structure format
    /// - dcd: DCD trajectory format
    /// - xyz: XYZ format
    /// - xtc: GROMACS compressed trajectory format
    /// - gro: GROMACS structure format
    ///
    /// # Errors
    /// Returns [FileIoError] if:
    /// - File extension is not recognized as writable format
    /// - File cannot be created
    pub fn create(fname: impl AsRef<Path>) -> Result<Self, FileIoError> {
        Ok(Self {
            file_path: fname.as_ref().to_owned(),
            format_handler: FileFormat::create(fname.as_ref())
                .map_err(|e| FileIoError(fname.as_ref().to_owned(), e))?,
            stats: Default::default(),
        })
    }

    /// Reads both topology and state from formats that contain both.
    ///
    /// Only works for formats that contain both topology and coordinates in a single file:
    /// - pdb
    /// - gro
    /// - tpr
    ///
    /// # Errors
    /// Returns [FileIoError] if:
    /// - Format doesn't support combined topology+state reading
    /// - File format is invalid
    /// - State doesn't match topology size
    pub fn read(&mut self) -> Result<(Topology, State), FileIoError> {
        let t = std::time::Instant::now();

        let (top, st) = self
            .format_handler
            .read()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok((top, st))
    }

    /// Writes both topology and state to formats that support it.
    ///
    /// Only works for formats that can write both topology and coordinates:
    /// - pdb
    /// - gro
    ///
    /// # Errors
    /// Returns [FileIoError] if format doesn't support writing both topology and state
    pub fn write<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyIoProvider + StateIoProvider,
    {
        let t = std::time::Instant::now();

        self.format_handler
            .write(data)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    /// Reads topology from supported formats.
    ///
    /// Works for formats containing topology information:
    /// - pdb
    /// - gro
    /// - tpr
    /// - itp
    /// - xyz
    ///
    /// # Errors
    /// Returns [FileIoError] if:
    /// - Format doesn't support topology reading
    /// - File format is invalid
    pub fn read_topology(&mut self) -> Result<Topology, FileIoError> {
        let t = std::time::Instant::now();

        let top = self
            .format_handler
            .read_topology()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(top)
    }

    /// Writes topology to supported formats.
    ///
    /// Currently supported formats:
    /// - pdb
    /// - xyz
    ///
    /// # Errors
    /// Returns [FileIoError] if format doesn't support topology writing
    pub fn write_topology<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: TopologyIoProvider,
    {
        let t = std::time::Instant::now();

        self.format_handler
            .write_topology(data)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    /// Reads next state from the file.
    ///
    /// For trajectory formats (XTC, DCD) returns sequential frames.
    /// For single-frame formats (TPR) returns that frame once
    /// then None afterwards.
    ///
    /// # Errors
    /// Returns [FileIoError] if:
    /// - Format doesn't support state reading
    /// - File format is invalid
    pub fn read_state(&mut self) -> Result<Option<State>, FileIoError> {
        let t = std::time::Instant::now();

        let st = self
            .format_handler
            .read_state()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        // Update stats if not None
        if st.is_some() {
            self.stats.frames_processed += 1;
            self.stats.cur_t = st.as_ref().unwrap().get_time();
        }

        Ok(st)
    }

    /// Writes state to trajectory or structure formats.
    ///
    /// Supported formats:
    /// - pdb
    /// - xyz
    /// - dcd
    /// - xtc
    ///
    /// # Errors
    /// Returns [FileIoError] if format doesn't support state writing
    pub fn write_state<T>(&mut self, data: &T) -> Result<(), FileIoError>
    where
        T: StateIoProvider,
    {
        let t = std::time::Instant::now();

        self.format_handler
            .write_state(data)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(())
    }

    /// Seeks to specified frame number in random-access trajectories.
    ///
    /// Currently only works for XTC format.
    ///
    /// # Errors
    /// Returns [FileIoError] if:
    /// - Format doesn't support random access
    /// - Frame number is invalid
    pub fn seek_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        Ok(self
            .format_handler
            .seek_frame(fr)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }

    /// Seeks to specified time in random-access trajectories.
    ///
    /// Currently only works for XTC format.
    ///
    /// # Errors
    /// Returns [FileIoError] if:
    /// - Format doesn't support random access  
    /// - Time value is invalid
    pub fn seek_time(&mut self, t: f32) -> Result<(), FileIoError> {
        Ok(self
            .format_handler
            .seek_time(t)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }

    /// Returns frame number and time of first frame in trajectory.
    ///
    /// Currently only works for XTC format.
    ///
    /// # Errors
    /// Returns [FileIoError] if format doesn't support random access
    pub fn tell_first(&self) -> Result<(usize, f32), FileIoError> {
        Ok(self
            .format_handler
            .tell_first()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }

    /// Returns current frame number and time in trajectory.
    ///
    /// For random access formats (XTC) returns actual frame info.
    /// For other formats returns number of processed frames and current time.
    ///
    /// # Errors
    /// Returns [FileIoError] if format info cannot be obtained
    pub fn tell_current(&self) -> Result<(usize, f32), FileIoError> {
        Ok(self
            .format_handler
            .tell_current(&self.stats)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }

    /// Returns last frame number and time in trajectory.
    ///
    /// Currently only works for XTC format.
    ///
    /// # Errors
    /// Returns [FileIoError] if format doesn't support random access
    pub fn tell_last(&self) -> Result<(usize, f32), FileIoError> {
        Ok(self
            .format_handler
            .tell_last()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }

    /// Consumes frames until reaching serial frame number `fr` (which is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.
    pub fn skip_to_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        // Try random-access first
        Ok(self
            .format_handler
            .skip_to_frame(fr, &self.stats)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }

    /// Consumes frames until reaching beyond time `t` (frame with time exactly equal to `t` is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.
    pub fn skip_to_time(&mut self, t: f32) -> Result<(), FileIoError> {
        // Try random-access first
        Ok(self
            .format_handler
            .skip_to_time(t, &self.stats)
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }
}

impl Drop for FileHandler {
    fn drop(&mut self) {
        debug!(
            "Done with file '{}': {}",
            self.file_path.display(),
            self.stats
        );
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
        impl LenProvider for ($t,$s) {
            fn len(&self) -> usize {
                self.0.len()
            }
        }

        impl TopologyIoProvider for ($t, $s) {}

        impl RandomAtomProvider for ($t, $s) {
            unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
                self.0.get_atom_unchecked(i)
            }
        }

        impl AtomIterProvider for ($t, $s) {
            fn iter_atoms(&self) -> impl AtomIterator<'_> {
                self.0.iter_atoms()
            }
        }

        impl MoleculesProvider for ($t, $s) {
            fn num_molecules(&self) -> usize {
                self.0.num_molecules()
            }

            fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.0.iter_molecules()
            }

            unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
                self.0.get_molecule_unchecked(i)
            }
        }

        impl BondsProvider for ($t, $s) {
            fn num_bonds(&self) -> usize {
                self.0.num_bonds()
            }

            fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.0.iter_bonds()
            }

            unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
                self.0.get_bond_unchecked(i)
            }
        }

        impl StateIoProvider for ($t, $s) {}

        impl TimeProvider for ($t, $s) {
            fn get_time(&self) -> f32 {
                self.1.get_time()
            }
        }

        impl BoxProvider for ($t, $s) {
            fn get_box(&self) -> Option<&PeriodicBox> {
                self.1.get_box()
            }
        }

        impl PosIterProvider for ($t, $s) {
            fn iter_pos(&self) -> impl PosIterator {
                self.1.iter_pos()
            }
        }

        impl RandomPosProvider for ($t, $s) {
            unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
                self.1.get_pos_unchecked(i)
            }
        }
    };
}

impl_io_traits_for_tuples!(Topology, State);
impl_io_traits_for_tuples!(Arc<Topology>,Arc<State>);

//--------------------------------------------------------
/// An error that occurred during file IO operations
#[derive(Error, Debug)]
#[error("in file {0}:")]
pub struct FileIoError(PathBuf, #[source] FileFormatError);

/// Errors that can occur when handling different file formats
#[derive(Error, Debug)]
enum FileFormatError {
    /// VMD molfile format handler error
    #[error("in vmd format handler")]
    Vmd(#[from] VmdHandlerError),

    /// GRO format handler error
    #[error("in gro format handler")]
    Gro(#[from] GroHandlerError),

    /// XTC format handler error
    #[error("in xtc format handler")]
    Xtc(#[from] XtcHandlerError),

    /// TPR format handler error
    #[error("in tpr format handler")]
    Tpr(#[from] TprHandlerError),

    /// ITP format handler error
    #[error("in itp format handler")]
    Itp(#[from] ItpHandlerError),

    #[error("file has no extension")]
    NoExtension,

    #[error("format is not readable")]
    NotReadable,

    #[error("format is not writable")]
    NotWritable,

    #[error("file has no more states to read")]
    NoStates,

    #[error("not a format containing both topology and state")]
    NotTopologyAndStateFormat,

    #[error("not a write once format")]
    NotWriteOnceFormat,

    #[error("not a topology reading format")]
    NotTopologyReadFormat,

    #[error("not a topology writing format")]
    NotTopologyWriteFormat,

    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizes),

    #[error("not a trajectory write format")]
    NotTrajectoryWriteFormat,

    #[error("not a random access format")]
    NotRandomAccessFormat,

    // #[error("unexpected end of file")]
    // UnexpectedEof,
    #[error("can't seek to frame {0}")]
    SeekFrameError(usize),

    #[error("can't seek to time {0}")]
    SeekTimeError(f32),

    #[error("not a state read format")]
    NotStateReadFormat,
}
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
        println!("#1: {}", top1.len());

        let b = System::new(top1, st2)?;
        let sel = b.select_all()?;
        sel.rotate(&Vector3f::x_axis(), 45.0_f32.to_radians());

        let outname = concat!(env!("OUT_DIR"), "/2.pdb");
        println!("{outname}");
        Ok(())
    }

    #[test]
    fn test_itp() -> Result<()> {
        let mut h = FileHandler::open("../molar_membrane/tests/POPE.itp")?;
        let top = h.read_topology()?;
        for a in top.iter_atoms() {
            println!("{:?}", a);
        }
        Ok(())
    }

    #[test]
    fn test_io_gro_traj() -> Result<()> {
        let src = System::from_file("tests/protein.pdb")?;
        println!(env!("OUT_DIR"));
        let groname = concat!(env!("OUT_DIR"), "/multi.gro");
        let mut w = FileHandler::create(groname)?;
        src.set_time(1.0);
        w.write(&src)?;
        src.set_time(2.0);
        w.write(&src)?;
        drop(w); // To flush buffer

        // Read it back
        let mut h = FileHandler::open(groname)?;
        let st = h.read_state()?.unwrap();
        println!("t1={}", st.get_time());
        let st = h.read_state()?.unwrap();
        println!("t2={}", st.get_time());
        Ok(())
    }

    #[test]
    fn test_iter_gro_traj() -> Result<()> {
        let mut h = FileHandler::open("tests/multi.gro")?;
        let _ = h.read_topology()?;
        for st in h.into_iter() {
            println!("{}", st.get_time());
        }
        Ok(())
    }
}
