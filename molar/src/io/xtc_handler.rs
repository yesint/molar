use super::{PeriodicBoxError, State};
use nalgebra::Matrix3;
use thiserror::Error;

use crate::io::FileFormatHandler;
use crate::periodic_box::PeriodicBox;
use crate::{FileFormatError, Float, Matrix3f, Pos};

#[cfg(not(feature = "f64"))]
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
                    .map_err(|e| xtc_to_ff_err(XtcHandlerError::ReadFrame(e)))?;
                // Destructure the frame
                let molly::Frame {
                    positions,
                    time,
                    step: _,
                    boxvec,
                    precision: _,
                } = molly_fr;

                // XTC is f32-on-disk. In f32 builds we reuse the buffer as
                // Vec<Pos> via raw-parts (zero-copy). In f64 builds we must
                // allocate a fresh Vec<Pos> and upcast triplet-by-triplet.
                #[cfg(not(feature = "f64"))]
                let coords = {
                    const _: () = assert!(
                        std::mem::size_of::<Pos>() == std::mem::size_of::<[f32; 3]>(),
                        "Pos must be layout-compatible with [f32; 3] for the f32 raw-parts fast path",
                    );
                    let mut positions = ManuallyDrop::new(positions);
                    let p = positions.as_mut_ptr() as *mut Pos;
                    let len = positions.len() / 3;
                    let cap = positions.capacity() / 3;
                    unsafe { Vec::from_raw_parts(p, len, cap) }
                };
                #[cfg(feature = "f64")]
                let coords: Vec<Pos> = positions
                    .chunks_exact(3)
                    .map(|c| Pos::new(c[0] as Float, c[1] as Float, c[2] as Float))
                    .collect();

                // Create a periodic box matrix (boxvec is [f32; 9])
                let m = Matrix3f::from_iterator(boxvec.into_iter().map(|x| x as Float));
                // Update frame counter
                r.cur_fr += 1;
                // Create and return our State
                Ok(State {
                    coords,
                    time: time as Float,
                    pbox: Some(PeriodicBox::from_matrix(m).map_err(|e| XtcHandlerError::Pbc(e))?),
                    ..Default::default()
                })
            }
            Self::Writer(_) => unreachable!(),
        }
    }

    fn write_state(&mut self, data: &dyn super::SaveState) -> Result<(), super::FileFormatError> {
        match self {
            Self::Writer(w) => {
                // Construct molly box (always f32 on the wire). In f64 builds the
                // matrix is downcast component-by-component.
                let m: Matrix3<f32> = match data.get_box() {
                    #[cfg(not(feature = "f64"))]
                    Some(b) => b.get_matrix(),
                    #[cfg(feature = "f64")]
                    Some(b) => Matrix3::<f32>::from_iterator(
                        b.get_matrix().iter().map(|x| *x as f32),
                    ),
                    None => Matrix3::<f32>::zeros(),
                };

                // Update frame counter
                w.cur_fr += 1;

                // f32: zero-copy slice view of `Pos` (which is `Point3<f32>`).
                // f64: allocate triplet array per atom — XTC is f32-on-disk.
                #[cfg(not(feature = "f64"))]
                {
                    let it = data.iter_pos_dyn().map(|p| p.coords.as_ref());
                    w.writer.write_frame_parts(
                        w.cur_fr as u32,
                        data.get_time(),
                        m.as_slice().try_into().unwrap(),
                        it,
                        1000.0,
                    ).map_err(|e| XtcHandlerError::WriteFrame(e))?;
                }
                #[cfg(feature = "f64")]
                {
                    let buf: Vec<[f32; 3]> = data.iter_pos_dyn()
                        .map(|p| [p.x as f32, p.y as f32, p.z as f32])
                        .collect();
                    let it = buf.iter().map(|a| a as &[f32; 3]);
                    w.writer.write_frame_parts(
                        w.cur_fr as u32,
                        data.get_time() as f32,
                        m.as_slice().try_into().unwrap(),
                        it,
                        1000.0,
                    ).map_err(|e| XtcHandlerError::WriteFrame(e))?;
                }
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
                        r.reader
                            .seek_prev()
                            .map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                    }
                } else {
                    r.reader
                        .skip_frames((fr - r.cur_fr) as u64)
                        .map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                }
            }
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }

    fn seek_time(&mut self, t: Float) -> Result<(), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                r.reader
                    .skip_to_time(t as f32)
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
    SeekTime(Float, #[source] std::io::Error),

    #[error("unexpected io error in xtc file")]
    Io(#[from] std::io::Error),
}

fn xtc_to_ff_err(e: XtcHandlerError) -> FileFormatError {
    if let XtcHandlerError::ReadFrame(ref io_err) = e {
        if io_err.kind() == std::io::ErrorKind::UnexpectedEof {
            return FileFormatError::Eof;
        }
    }
    FileFormatError::from(e)
}

