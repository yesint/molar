use super::{PeriodicBoxError, State};
use nalgebra::Matrix3;
use thiserror::Error;

use crate::io::{DynSource, FileFormatHandler};
use crate::periodic_box::PeriodicBox;
use crate::{FileFormatError, Float, Matrix3f, Pos};

#[cfg(not(feature = "f64"))]
use std::mem::ManuallyDrop;
use std::io::{Seek, SeekFrom};
use std::path::Path;

// XTC reads from any `Read + Seek` source via molly's generic `XTCReader<R>`
// sequential path; random-access seeking (molly's fast path is `File`-only,
// using its internal `Buffer`) is reimplemented here over molly's *public* API
// (`XTCReader::{file, step, read_header}`, `Header`, `molly::reader::read_nbytes`,
// `molly::padding`). This is a faithful port of molly's File-bound seek methods,
// minus the `Buffer` optimization, so it works for files and browser blobs alike.
pub(crate) struct XtcReader {
    reader: molly::XTCReader<DynSource>,
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

impl XtcFileHandler {
    /// Build a reading handler from an arbitrary byte source.
    pub(crate) fn from_source(src: DynSource) -> Result<Self, super::FileFormatError> {
        Ok(XtcFileHandler::Reader(XtcReader {
            reader: molly::XTCReader::new(src),
            cur_fr: 0,
        }))
    }
}

impl FileFormatHandler for XtcFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, super::FileFormatError>
    where
        Self: Sized,
    {
        let file = std::fs::File::open(fname).map_err(XtcHandlerError::Io)?;
        Self::from_source(DynSource(Box::new(file)))
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
                        seek_prev(&mut r.reader).map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                    }
                } else {
                    skip_frames(&mut r.reader, (fr - r.cur_fr) as u64)
                        .map_err(|e| XtcHandlerError::SeekFrame(fr, e))?;
                }
                r.cur_fr = fr;
            }
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }

    fn seek_time(&mut self, t: Float) -> Result<(), super::FileFormatError> {
        match self {
            Self::Reader(r) => {
                skip_to_time(&mut r.reader, t as f32)
                    .map_err(|e| XtcHandlerError::SeekTime(t, e))?;
            }
            Self::Writer(_) => unreachable!(),
        }
        Ok(())
    }
}

// ── Generic XTC seeking over any `Read + Seek` source ─────────────────────────
// Faithful ports of molly's File-bound seek helpers, written against molly's
// public `XTCReader` fields and `Header`/`read_nbytes`/`padding`. No `Buffer`
// optimization (that one is genuinely `File`-only), so these read headers
// directly; correct for files and browser blobs alike.

/// Skip the positions block of the current frame, leaving the reader at the
/// start of the next frame. Assumes the position is right after `header`.
fn skip_positions(
    r: &mut molly::XTCReader<DynSource>,
    header: &molly::Header,
) -> std::io::Result<u64> {
    let skip = if header.natoms <= 9 {
        // Small frames store uncompressed positions: 3 f32 per atom.
        header.natoms as u64 * 3 * 4
    } else {
        // Skip the compression preamble, read the byte count, then skip the
        // compressed block plus its XDR padding.
        r.file.seek(SeekFrom::Current(32))?;
        let nbytes = molly::reader::read_nbytes(&mut r.file, header.magic)? as u64;
        nbytes + molly::padding(nbytes as usize) as u64
    };
    r.file.seek(SeekFrom::Current(skip as i64))
}

/// Advance one frame (assumes position at a frame start).
fn seek_next(r: &mut molly::XTCReader<DynSource>) -> std::io::Result<()> {
    let header = r.read_header()?;
    skip_positions(r, &header)?;
    Ok(())
}

/// Advance `n` frames.
fn skip_frames(r: &mut molly::XTCReader<DynSource>, n: u64) -> std::io::Result<()> {
    for _ in 0..n {
        seek_next(r)?;
    }
    Ok(())
}

/// Step back one frame: scan backwards for the previous valid header and leave
/// the reader at its start. Mirrors molly's `seek_prev`.
fn seek_prev(r: &mut molly::XTCReader<DynSource>) -> std::io::Result<molly::Header> {
    const XDR_INT_SIZE: i64 = 4;
    r.file
        .seek(SeekFrom::Current(-(molly::Header::SIZE as i64 + XDR_INT_SIZE * 2)))?;
    loop {
        let pos = r.file.stream_position()?;
        if let Ok(header) = r.read_header() {
            r.file.seek(SeekFrom::Start(pos))?;
            return Ok(header);
        }
        r.file.seek(SeekFrom::Current(-2 * XDR_INT_SIZE))?;
    }
}

/// Skip forward until a frame with `time >= t`, leaving the reader at its start.
fn skip_to_time(
    r: &mut molly::XTCReader<DynSource>,
    t: f32,
) -> std::io::Result<molly::Header> {
    loop {
        let pos = r.file.stream_position()?;
        let header = r.read_header()?;
        if header.time >= t {
            r.file.seek(SeekFrom::Start(pos))?;
            return Ok(header);
        }
        skip_positions(r, &header)?;
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

