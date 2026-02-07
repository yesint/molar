use super::{PeriodicBoxError, State};
use nalgebra::Matrix3;
use thiserror::Error;

use crate::io::FileFormatHandler;
use crate::periodic_box::PeriodicBox;
use crate::{FileFormatError, Matrix3f, Pos};

use std::io::{self, Seek, SeekFrom};
use std::mem::ManuallyDrop;
use std::path::Path;

pub(crate) struct XtcReader {
    reader: molly::XTCReader<std::fs::File>,
    cur_fr: usize,
}

pub(crate) struct XtcWriter {
    writer: molly::XTCWriter<std::fs::File>,
    cur_fr: usize,
}

pub enum XtcFileHandler {
    Reader(XtcReader),
    Writer(XtcWriter),
}

impl XtcReader {
    // Seek the start of the previous frame and read its header.
    // The file position remains at the beginning of the frame.
    fn goto_prev_and_inspect_header(&mut self) -> io::Result<molly::Header> {
        // Save where we are initially
        let initial_pos = self.reader.file.stream_position()?;
        // Logic is taken from xdrfile extention for mdtraj, no idea how it works :)
        const XDR_INT_SIZE: i64 = 4;
        self.reader
            .file
            .seek(SeekFrom::Current(-3 * XDR_INT_SIZE))?;
        // In a loop try reading the header until we succeed
        // or beginning of the file reached
        loop {
            // Remember current offset to return to it if header will be read successfully
            let pos = self.reader.file.stream_position()?;
            let h_res = self.reader.read_header();
            
            if let Ok(h) = h_res {
                // Success! Update frame count and return header
                self.cur_fr -= 1;
                self.reader.file.seek(SeekFrom::Start(pos))?;
                return Ok(h);
            } else {
                // Header not found, seek backwards
                self.reader
                    .file
                    .seek(SeekFrom::Current(-2 * XDR_INT_SIZE))
                    .or_else(|e| {
                        // This will fail if passing beginning of the file
                        // Go back to initial file position
                        self.reader.file.seek(SeekFrom::Start(initial_pos))?;
                        return Err(e);
                    })?;
            }
        }
    }

    /// Skip positions of the current frame and jump to the beginning of the next one
    /// Assumes that we have just read the header of the current frame and it is passed
    /// as an argument.
    /// Returns the file offset of the beginning of new frame.
    fn skip_positions(&mut self, header: &molly::Header) -> io::Result<u64> {
        let skip = if header.natoms <= 9 {
            // Know how many bytes are in this frame until the next header since the positions
            // are uncompressed.
            header.natoms as u64 * 3 * 4
        } else {
            // We need to read the nbytes value to get the offset until the next header.
            self.reader.file.seek(SeekFrom::Current(32))?;
            // The size of the buffer is stored either as a 64 or 32-bit integer, depending on
            // the magic number in the header.
            let nbytes = molly::reader::read_nbytes(&mut self.reader.file, header.magic)? as u64;
            nbytes + molly::padding(nbytes as usize) as u64
        };
        let offset = self.reader.file.seek(SeekFrom::Current(skip as i64))?;
        Ok(offset)
    }

    // Seek the start of the next frame and read its header.
    // The file position remains at the beginning of the frame.
    fn goto_next_and_inspect_header(&mut self) -> io::Result<molly::Header> {
        // Remember where we start so we can return to it later.
        let initial_pos = self.reader.file.stream_position()?;
        // Read current header
        let mut header = self.reader.read_header()?;
        let frame_offset = self.skip_positions(&header)?;
        // Read next header. This may fail and in this case we should
        // restore initial file position
        header = self.reader.read_header().or_else(|e| {
            self.reader.file.seek(SeekFrom::Start(initial_pos))?;
            Err(e)
        })?;
        // Update frame counter
        self.cur_fr += 1;
        // Seek back to frame start after reading header
        self.reader.file.seek(SeekFrom::Start(frame_offset))?;
        Ok(header)
    }
}

impl FileFormatHandler for XtcFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, super::FileFormatError>
    where
        Self: Sized,
    {
        Ok(XtcFileHandler::Reader(XtcReader {
            reader: molly::XTCReader::open(fname)?,
            cur_fr: 0,
        }))
    }

    fn create(fname: impl AsRef<Path>) -> Result<Self, super::FileFormatError>
    where
        Self: Sized,
    {
        Ok(XtcFileHandler::Writer(XtcWriter {
            writer: molly::XTCWriter::create(fname)?,
            cur_fr: 0,
        }))
    }

    fn read_state(&mut self) -> Result<State, super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                let mut molly_fr = molly::Frame::default();
                r.reader
                    .read_frame(&mut molly_fr)
                    .map_err(|e| XtcHandlerError::ReadFrame(e))?;
                // Destructure the frame
                let molly::Frame {
                    positions,
                    time,
                    step: _,
                    boxvec,
                    precision: _,
                } = molly_fr;

                let mut positions = ManuallyDrop::new(positions);
                // Destructure it into parts
                let p = positions.as_mut_ptr() as *mut Pos;
                let len = positions.len() / 3;
                let cap = positions.capacity() / 3;
                // Assamble it back to our coord vector
                let coords = unsafe { Vec::from_raw_parts(p, len, cap) };
                // Create a periodic box matrix
                let m = Matrix3f::from_iterator(boxvec.to_cols_array().into_iter());
                // Update frame counter
                r.cur_fr += 1;
                // Create and return our State
                Ok(State {
                    coords,
                    time,
                    pbox: Some(PeriodicBox::from_matrix(m).map_err(|e| XtcHandlerError::Pbc(e))?),
                })
            }
            Self::Writer(_) => unreachable!(),
        }
    }

    fn write_state(&mut self, data: &dyn super::SaveState) -> Result<(), super::FileFormatError> {
        match self {
            Self::Writer(w) => {
                // Coordinate buffer with ownership to be transfered to molly frame
                let mut buf =
                    ManuallyDrop::new(Vec::<Pos>::from_iter(data.iter_pos_dyn().cloned()));
                // Decompose the buffer
                let p = buf.as_mut_ptr() as *mut f32;
                let len = buf.len() * 3;
                let cap = buf.capacity() * 3;
                // Assamble it back to molly positions
                let positions = unsafe { Vec::from_raw_parts(p, len, cap) };
                // Construct molly box
                let m = match data.get_box() {
                    Some(b) => b.get_matrix(),
                    None => Matrix3::<f32>::zeros(),
                };
                let boxvec = glam::Mat3::from_cols_slice(m.as_slice());
                let fr = molly::Frame {
                    precision: 1000.0,
                    time: data.get_time(),
                    step: w.cur_fr as u32,
                    positions,
                    boxvec,
                };
                // Update frame counter
                w.cur_fr += 1;
                w.writer
                    .write_frame(&fr)
                    .map_err(|e| XtcHandlerError::WriteFrame(e))?
            }

            Self::Reader(_) => unreachable!(),
        }

        Ok(())
    }

    fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
        match self {
            Self::Reader(r) => {
                if fr == r.cur_fr {
                    return Ok(());
                }
                if fr < r.cur_fr {
                    for _ in 0..r.cur_fr - fr {
                        r.goto_prev_and_inspect_header()
                            .map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                    }
                } else {
                    for _ in 0..fr - r.cur_fr {
                        r.goto_next_and_inspect_header()
                            .map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                    }
                }
            }
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }

    fn seek_time(&mut self, t: f32) -> Result<(), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                loop {
                    let pos = r.reader.file.stream_position()?;
                    let header = r.reader.read_header()?;
                    if header.time >= t {
                        r.reader.file.seek(SeekFrom::Start(pos))?;
                        break;
                    }
                    r.skip_positions(&header)
                        .map_err(|e| XtcHandlerError::SeekTime(t, e))?;
                }
            }
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }
}

#[derive(Error, Debug)]
pub enum XtcHandlerError {
    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),

    #[error("failed to read frame")]
    ReadFrame(#[source] std::io::Error),

    #[error("failed to write frame")]
    WriteFrame(#[source] std::io::Error),

    #[error("seek to frame {0} failed")]
    SeekFrame(usize, #[source] std::io::Error),

    #[error("seek to time {0} failed")]
    SeekTime(f32, #[source] std::io::Error),

    #[error("unexpected io error in xtc file")]
    Io(#[from] std::io::Error),
}
