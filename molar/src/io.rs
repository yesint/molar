//! IO handlers for different file formats and associated traits

use crate::prelude::*;
use log::{debug, warn};
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

use gro_handler::{GroFileHandler, GroHandlerError};
use itp_handler::{ItpFileHandler, ItpHandlerError};
use tpr_handler::{TprFileHandler, TprHandlerError};
use vmd_molfile_handler::{VmdHandlerError, VmdMolFileHandler};
use xtc_handler::{XtcFileHandler, XtcHandlerError};

/// Trait for saving [Topology] to file
pub trait SaveTopology: RandomBondProvider + LenProvider {
    fn iter_atoms_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Atom> + 'a>;
}

/// Trait for saving [State] to file
pub trait SaveState: BoxProvider + TimeProvider + LenProvider {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Pos> + 'a>;
}

/// Trait for saving both [Topology] and [State] to file
pub trait SaveTopologyState: SaveTopology + SaveState {
    fn save(&self, fname: &str) -> Result<(), FileIoError>
    where
        Self: Sized,
    {
        let mut h = FileHandler::create(fname)?;
        h.write(self)?;
        Ok(())
    }
}

/// Trait for file format handlers.
/// Concrete handlers implement only methods which are
/// relevant for a particular format.
#[allow(unused_variables)]
pub(crate) trait FileFormatHandler: Send {
    /// Open file for reading
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Err(FileFormatError::NotReadable)
    }

    /// Open file for writing (overwrites if it exists)
    fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Err(FileFormatError::NotWritable)
    }

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        Err(FileFormatError::NotTopologyStateReadFormat)
    }

    fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        Err(FileFormatError::NotTopologyReadFormat)
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        Err(FileFormatError::NotStateReadFormat)
    }

    fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileFormatError> {
        Err(FileFormatError::NotTopologyStateWriteFormat)
    }

    fn write_topology(&mut self, data: &dyn SaveTopology) -> Result<(), FileFormatError> {
        Err(FileFormatError::NotTopologyWriteFormat)
    }

    fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileFormatError> {
        Err(FileFormatError::NotStateWriteFormat)
    }

    fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
        Err(FileFormatError::NotRandomAccessFormat)
    }

    fn seek_time(&mut self, t: f32) -> Result<(), FileFormatError> {
        Err(FileFormatError::NotRandomAccessFormat)
    }

    fn seek_last(&mut self) -> Result<(), FileFormatError> {
        Err(FileFormatError::NotRandomAccessFormat)
    }
}

/// Iterator over states in a trajectory file
pub struct IoStateIterator {
    receiver: std::sync::mpsc::Receiver<Result<State, FileIoError>>,
}

impl IoStateIterator {
    /// Create new state iterator by consuming the file handler
    fn new(mut fh: FileHandler) -> Self {
        use std::sync::mpsc::sync_channel;
        let (sender, receiver) = sync_channel(10);

        //Spawn reading thread
        std::thread::spawn(move || {
            let mut terminate = false;
            while !terminate {
                let res = fh.read_state();
                terminate = match res.as_ref() {
                    Err(_) => true, // terminate if reading failed or Eof returned
                    _ => false,     // otherwise continut reading
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
        // If it does then something is horrible wrong anyway, so panicing is fine here.
        match self.receiver.recv().expect("reader thread shouldn't crash") {
            Ok(opt_st) => Some(opt_st),
            Err(FileIoError(f, e)) => {
                match e {
                    // Do nothing on EOF, just finish normally
                    FileFormatError::Eof => {}
                    // On any other error complain
                    _ => warn!(
                        "file '{}' is likely corrupted, reading stopped: {}",
                        f.display(),
                        e
                    ),
                }
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
    format_handler: Box<dyn FileFormatHandler>,
    pub stats: FileStats,
}

/// Statistics of file operations including elapsed time and number of processed frames
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
        let ext = get_ext(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?;
        let format_handler: Box<dyn FileFormatHandler> = match ext {
            "pdb" | "dcd" | "xyz" => Box::new(
                VmdMolFileHandler::open(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            "xtc" => Box::new(
                XtcFileHandler::open(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            "gro" => Box::new(
                GroFileHandler::open(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            "itp" => Box::new(
                ItpFileHandler::open(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            "tpr" => Box::new(
                TprFileHandler::open(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            _ => Err(FileFormatError::NotRecognized)
                .map_err(|e| FileIoError(fname.to_path_buf(), e))?,
        };
        Ok(Self {
            file_path: fname.to_path_buf(),
            format_handler,
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
        let fname = fname.as_ref();
        let ext = get_ext(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?;
        let format_handler: Box<dyn FileFormatHandler> = match ext {
            "pdb" | "dcd" | "xyz" => Box::new(
                VmdMolFileHandler::create(fname)
                    .map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            "xtc" => Box::new(
                XtcFileHandler::create(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            "gro" => Box::new(
                GroFileHandler::create(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            "itp" => Box::new(
                ItpFileHandler::create(fname).map_err(|e| FileIoError(fname.to_path_buf(), e))?,
            ),
            _ => Err(FileFormatError::NotRecognized)
                .map_err(|e| FileIoError(fname.to_path_buf(), e))?,
        };
        Ok(Self {
            file_path: fname.to_path_buf(),
            format_handler,
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

        let ret = self
            .format_handler
            .read()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?;

        self.stats.elapsed_time += t.elapsed();
        self.stats.frames_processed += 1;

        Ok(ret)
    }

    /// Writes both topology and state to formats that support it.
    ///
    /// Only works for formats that can write both topology and coordinates:
    /// - pdb
    /// - gro
    ///
    /// # Errors
    /// Returns [FileIoError] if format doesn't support writing both topology and state
    pub fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileIoError> {
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
    pub fn write_topology(&mut self, data: &dyn SaveTopology) -> Result<(), FileIoError> {
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
    pub fn read_state(&mut self) -> Result<State, FileIoError> {
        let t = std::time::Instant::now();

        let res = self
            .format_handler
            .read_state()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e));

        // Update stats if not error
        if res.is_ok() {
            self.stats.frames_processed += 1;
            self.stats.cur_t = res.as_ref().unwrap().get_time();
        }

        self.stats.elapsed_time += t.elapsed();
        Ok(res?)
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
    pub fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileIoError> {
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

    /// Jumps to the last frame in trajectory.
    ///
    /// # Errors
    /// Returns [FileIoError] if format doesn't support random access
    pub fn seek_last(&mut self) -> Result<(), FileIoError> {
        Ok(self
            .format_handler
            .seek_last()
            .map_err(|e| FileIoError(self.file_path.to_owned(), e))?)
    }

    /// Skips frames until reaching serial frame number `fr` (which is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.
    pub fn skip_to_frame(&mut self, fr: usize) -> Result<(), FileIoError> {
        // Try random-access first
        self.seek_frame(fr).or_else(|_| {
            // Not a random access trajectory
            if self.stats.frames_processed == fr {
                // We are at needed frame, nothing to do
                Ok(())
            } else if self.stats.frames_processed > fr {
                // We are past the needed frame, return error
                Err(FileIoError(
                    self.file_path.to_path_buf(),
                    FileFormatError::SeekFrame(fr),
                ))
            } else {
                // Do serial read until reaching needed frame or EOF
                while self.stats.frames_processed < fr {
                    self.read_state()?;
                }
                Ok(())
            }
        })
    }

    /// Skips frames until reaching beyond time `t` (frame with time exactly equal to `t` is not consumed)
    /// This uses random-access if available and falls back to serial reading if it is not.
    pub fn skip_to_time(&mut self, t: f32) -> Result<(), FileIoError> {
        // Try random-access first
        self.seek_time(t).or_else(|_| {
            // Not a random access trajectory
            if self.stats.cur_t > t {
                // We are past the needed time, return error
                Err(FileIoError(
                    self.file_path.to_path_buf(),
                    FileFormatError::SeekTime(t),
                ))
            } else {
                // Do serial read until reaching needed time or EOF
                while self.stats.cur_t < t {
                    self.read_state()?;
                }
                Ok(())
            }
        })
    }

    // pub fn skip_to_last(&mut self) -> Result<(), FileIoError> {
    //     // Try random-access first
    //     self.seek_last().or_else(|_| {
    //         // Not a random access trajectory
    //         // Do serial read until reaching EOF
    //         while let Ok(_) = self.read_state() {}
            
    //         Ok(())
            
    //     })
    // }
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

//--------------------------------------------------------
/// An error that occurred during file IO operations
#[derive(Error, Debug)]
#[error("in file {0}:")]
pub struct FileIoError(pub(crate) PathBuf, #[source] pub(crate) FileFormatError);

/// Errors that can occur when handling different file formats
#[derive(Error, Debug)]
pub(crate) enum FileFormatError {
    #[error("end of file reached")]
    Eof,

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

    #[error("not a format able to read topology and state at once")]
    NotTopologyStateReadFormat,

    #[error("not a format able to write topology and state at once")]
    NotTopologyStateWriteFormat,

    #[error("not a topology reading format")]
    NotTopologyReadFormat,

    #[error("not a topology writing format")]
    NotTopologyWriteFormat,

    #[error("not a state reading format")]
    NotStateReadFormat,

    #[error("not a state writing format")]
    NotStateWriteFormat,

    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizesError),

    #[error("not a random access format")]
    NotRandomAccessFormat,

    #[error("can't seek to frame {0}")]
    SeekFrame(usize),

    #[error("can't seek to time {0}")]
    SeekTime(f32),

    #[error("unexpected io error")]
    Io(#[from] std::io::Error),

    #[error("file extension is not recognized")]
    NotRecognized,
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
    fn test_seek_frame() -> Result<()> {
        let mut r = FileHandler::open("tests/protein.xtc")?;
        r.seek_frame(2)?;
        let st = r.read_state()?;
        println!("{} {}",r.stats.frames_processed, st.time);
        r.seek_frame(1)?;
        let st = r.read_state()?;
        println!("{} {}",r.stats.frames_processed, st.time);
        r.seek_frame(3)?;
        let st = r.read_state()?;
        println!("{} {}",r.stats.frames_processed, st.time);
        Ok(())
    }

    #[test]
    fn test_seek_time() -> Result<()> {
        let mut r = FileHandler::open("tests/protein.xtc")?;
        r.seek_time(200_001.0)?;
        let st = r.read_state()?;
        println!("{} {}",r.stats.frames_processed, st.time);
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
        let st1 = r.read_state()?;
        let st2 = st1.clone();
        println!("#1: {}", top1.len());

        let mut sys = System::new(top1, st2)?;
        let mut sel = sys.select_all_bound_mut();
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
        let mut sys = System::from_file("tests/protein.pdb")?;
        println!(env!("OUT_DIR"));
        let groname = concat!(env!("OUT_DIR"), "/multi.gro");
        let mut w = FileHandler::create(groname)?;
        sys.set_time(1.0);
        w.write(&sys)?;
        sys.set_time(2.0);
        w.write(&sys)?;
        drop(w); // To flush buffer

        // Read it back
        let mut h = FileHandler::open(groname)?;
        let st = h.read_state()?;
        println!("t1={}", st.get_time());
        let st = h.read_state()?;
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

    #[test]
    fn xyz_test() -> anyhow::Result<()> {
        let sys = System::from_file("tests/test.xyz")?;
        for atom in sys.iter_atoms() {
            println!("mass = {}", atom.mass);
        }
        Ok(())
    }
}
