use super::{PeriodicBoxError, State};
use nalgebra::Matrix3;
use thiserror::Error;

use crate::io::FileFormatHandler;
use crate::periodic_box::PeriodicBox;
use crate::{Matrix3f, Pos};

use std::io::{Seek, SeekFrom};
use std::mem::ManuallyDrop;
use std::path::Path;

pub(crate) struct XtcReader {
    reader: molly::XTCReader<std::fs::File>,
}

pub(crate) struct XtcWriter {
    writer: molly::XTCWriter<std::fs::File>,
    cur_fr: u32,
}

pub enum XtcFileHandler {
    Reader(XtcReader),
    Writer(XtcWriter),
}

impl XtcReader {
    // Reads a header and computes an offset to the next frame
    // The file is not advanced
    fn get_cur_header_and_offset(&mut self) -> std::io::Result<(molly::Header, u64)> {
        // Remember where we start so we can return to it later.
        let start_pos = self.reader.file.stream_position()?;
        // Read header
        let header = self.reader.read_header()?;
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
        // Return back to where we started.
        self.reader.file.seek(SeekFrom::Start(start_pos))?;
        Ok((header, offset))
    }

    fn get_last_header_and_offset(&mut self) -> Result<(molly::Header, u64), XtcHandlerError> {
        // We start seeking from the end
        // Save where we are
        let saved_pos = self.reader.file.stream_position()?;
        // Logic is taken from xdrfile extention for mdtraj, no idea how it works :)
        const XDR_INT_SIZE: i64 = 4;
        let mut offset = self.reader.file.seek(SeekFrom::End(-3 * XDR_INT_SIZE))?; // 4 is an XDR_INT_SIZE
        // In a loop try reading the header until we can
        loop {
            if let Ok(h) = self.reader.read_header() {
                // Success!
                // Go back to saved position
                self.reader.file.seek(SeekFrom::Start(saved_pos))?;
                // Return what we found
                return Ok((h, offset));
            } else if offset > 0 {
                // Header not found, seek backwards
                offset = self
                    .reader
                    .file
                    .seek(SeekFrom::Current(-2 * XDR_INT_SIZE))?;
            } else {
                // We reached beginning of the file, this is an error
                // Go back to saved position
                self.reader.file.seek(SeekFrom::Start(saved_pos))?;
                return Err(XtcHandlerError::SeekLastFrame);
            }
        }
    }
}

impl FileFormatHandler for XtcFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, super::FileFormatError>
    where
        Self: Sized,
    {
        Ok(XtcFileHandler::Reader(XtcReader {
            reader: molly::XTCReader::open(fname)?,
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
                r.reader.read_frame(&mut molly_fr).map_err(|e| XtcHandlerError::ReadFrame(e))?;
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
                // Create our State
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
                    step: w.cur_fr,
                    positions,
                    boxvec,
                };
                w.writer.write_frame(&fr).map_err(|e| XtcHandlerError::WriteFrame(e))?
            }

            Self::Reader(_) => unreachable!(),
        }

        Ok(())
    }

    // Skip fr frames starting from the current position
    fn seek_frame(&mut self, fr: usize) -> Result<(), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                for _ in 0..fr {
                    let (_, offset) = r
                        .get_cur_header_and_offset()
                        .map_err(|e| XtcHandlerError::SeekFrame(fr,e))?;
                    r.reader.file.seek(SeekFrom::Start(offset))?;
                }
            }
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }

    fn seek_time(&mut self, t: f32) -> Result<(), super::FileFormatError> {
        match self {
            Self::Reader(r) => loop {
                let (header, offset) = r
                    .get_cur_header_and_offset()
                    .map_err(|e| XtcHandlerError::SeekTime(t,e))?;
                if header.time >= t {
                    break;
                }
                r.reader.file.seek(SeekFrom::Start(offset))?;
            },
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }

    fn tell_last(&mut self) -> Result<(usize, f32), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                let (h, _) = r.get_last_header_and_offset()?;
                Ok((h.step as usize, h.time))
            }
            Self::Writer(_) => unreachable!(),
        }
    }

    fn tell_current(&mut self) -> Result<(usize, f32), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                let saved_pos = r.reader.file.stream_position()?;
                let h = r.reader.read_header()?;
                r.reader.file.seek(SeekFrom::Start(saved_pos))?;
                Ok((h.step as usize, h.time))
            }
            Self::Writer(_) => unreachable!(),
        }
    }

    fn tell_first(&mut self) -> Result<(usize, f32), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                // Save where we are
                let saved_pos = r.reader.file.stream_position()?;
                r.reader.file.seek(SeekFrom::Start(0))?;
                let h = r.reader.read_header()?;
                // Go back to saved position
                r.reader.file.seek(SeekFrom::Start(saved_pos))?;
                Ok((h.step as usize, h.time))
            }
            Self::Writer(_) => unreachable!(),
        }
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

    #[error("failed to find the header of last frame")]
    SeekLastFrame,

    #[error("seek to frame {0} failed")]
    SeekFrame(usize, #[source] std::io::Error),

    #[error("seek to time {0} failed")]
    SeekTime(f32, #[source] std::io::Error),

    #[error("unexpected io error in xtc file")]
    Io(#[from] std::io::Error),
}
