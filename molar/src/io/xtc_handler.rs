use super::{PeriodicBoxError, State};
use nalgebra::Matrix3;
use thiserror::Error;

use crate::io::FileFormatHandler;
use crate::periodic_box::PeriodicBox;
use crate::{FileFormatError, Matrix3f, Pos};

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
                let m = Matrix3f::from_iterator(boxvec.into_iter());
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
                let fr = molly::Frame {
                    precision: 1000.0,
                    time: data.get_time(),
                    step: w.cur_fr as u32,
                    positions,
                    boxvec: m.as_slice().try_into().unwrap(),
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
                        r.reader.seek_prev()
                            .map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                    }
                } else {
                    r.reader.skip_frames((fr - r.cur_fr) as u64)
                        .map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                }
            }
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }

    fn seek_time(&mut self, t: f32) -> Result<(), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                r.reader.skip_to_time(t)
                    .map_err(|e| XtcHandlerError::SeekTime(t, e))?;
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
