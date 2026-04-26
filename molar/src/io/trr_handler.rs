use crate::prelude::*;
use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
};
use thiserror::Error;

use super::{FileFormatError, FileFormatHandler, SaveState};

//=============================== Error type ===============================

#[derive(Error, Debug)]
pub enum TrrHandlerError {
    #[error("invalid TRR magic number {0}")]
    MagicNumber(i32),

    #[error("invalid TRR header")]
    Header,

    #[error("TRR seek to frame {0} failed")]
    SeekFrame(usize, #[source] io::Error),

    #[error("TRR seek to time {0} failed")]
    SeekTime(Float, #[source] io::Error),

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),

    #[error("unexpected io error")]
    Io(#[from] io::Error),
}

//=============================== XDR primitives ===========================

fn xdr_read_i32(r: &mut impl Read) -> io::Result<i32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(i32::from_be_bytes(buf))
}

fn xdr_read_f32(r: &mut impl Read) -> io::Result<f32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(f32::from_bits(u32::from_be_bytes(buf)))
}

fn xdr_read_f64(r: &mut impl Read) -> io::Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(f64::from_bits(u64::from_be_bytes(buf)))
}

fn xdr_write_i32(w: &mut impl Write, v: i32) -> io::Result<()> {
    w.write_all(&v.to_be_bytes())
}

fn xdr_write_f32(w: &mut impl Write, v: f32) -> io::Result<()> {
    w.write_all(&v.to_bits().to_be_bytes())
}

fn xdr_skip(r: &mut impl Read, n_bytes: usize) -> io::Result<()> {
    let mut buf = vec![0u8; n_bytes];
    r.read_exact(&mut buf)
}

/// Read `n` XDR-encoded 3D vectors (velocities or forces). Returns `Vec<Vel>`.
/// Since `Vel` and `Force` share the same alias, the same function serves both.
/// TRR is f32 or f64 on disk depending on the writer; we cast to `Float` at the boundary.
fn read_xvf(r: &mut impl Read, n: usize, b_double: bool) -> io::Result<Vec<Vel>> {
    if b_double {
        (0..n)
            .map(|_| {
                let x = xdr_read_f64(r)? as Float;
                let y = xdr_read_f64(r)? as Float;
                let z = xdr_read_f64(r)? as Float;
                Ok(Vel::new(x, y, z))
            })
            .collect()
    } else {
        (0..n)
            .map(|_| {
                let x = xdr_read_f32(r)? as Float;
                let y = xdr_read_f32(r)? as Float;
                let z = xdr_read_f32(r)? as Float;
                Ok(Vel::new(x, y, z))
            })
            .collect()
    }
}

//=============================== Header struct ============================

const GROMACS_MAGIC: i32 = 1993;
const TRR_VERSION: &str = "GMX_trn_file";

struct TrrHeader {
    box_size: i32,
    vir_size: i32,
    pres_size: i32,
    x_size: i32,
    v_size: i32,
    f_size: i32,
    natoms: i32,
    step: i32,
    time: f32,
    b_double: bool,
}

fn read_trr_header(r: &mut impl Read) -> Result<TrrHeader, TrrHandlerError> {
    // magic
    let magic = xdr_read_i32(r)?;
    if magic != GROMACS_MAGIC {
        return Err(TrrHandlerError::MagicNumber(magic));
    }

    // slen (= 13)
    let _slen = xdr_read_i32(r)?;

    // XDR string: i32 length (12) + 12 bytes "GMX_trn_file"
    let str_len = xdr_read_i32(r)? as usize;
    // XDR strings are padded to 4-byte boundary
    let padded = (str_len + 3) & !3;
    let mut str_buf = vec![0u8; padded];
    r.read_exact(&mut str_buf)?;
    let version = std::str::from_utf8(&str_buf[..str_len])
        .map_err(|_| TrrHandlerError::Header)?;
    if version != TRR_VERSION {
        return Err(TrrHandlerError::Header);
    }

    let ir_size = xdr_read_i32(r)?;
    let e_size = xdr_read_i32(r)?;
    let box_size = xdr_read_i32(r)?;
    let vir_size = xdr_read_i32(r)?;
    let pres_size = xdr_read_i32(r)?;
    let top_size = xdr_read_i32(r)?;
    let sym_size = xdr_read_i32(r)?;
    let x_size = xdr_read_i32(r)?;
    let v_size = xdr_read_i32(r)?;
    let f_size = xdr_read_i32(r)?;
    let natoms = xdr_read_i32(r)?;
    let step = xdr_read_i32(r)?;
    let _nre = xdr_read_i32(r)?;

    // Detect precision: if any size == natoms*3*8 or box == 9*8, it's double
    let n3 = natoms as i32 * 3;
    let b_double = (box_size == 72)
        || (x_size == n3 * 8)
        || (v_size != 0 && v_size == n3 * 8)
        || (f_size != 0 && f_size == n3 * 8);

    let elem_size = if b_double { 8 } else { 4 };

    // time and lambda
    let time = if b_double {
        xdr_read_f64(r)? as f32
    } else {
        xdr_read_f32(r)?
    };
    // lambda (skip)
    xdr_skip(r, elem_size)?;

    // Validate sizes (ignore ir_size, e_size, top_size, sym_size — should all be 0)
    let _ = (ir_size, e_size, top_size, sym_size);

    Ok(TrrHeader {
        box_size,
        vir_size,
        pres_size,
        x_size,
        v_size,
        f_size,
        natoms,
        step,
        time,
        b_double,
    })
}

fn write_trr_header(w: &mut impl Write, h: &TrrHeader) -> Result<(), TrrHandlerError> {
    xdr_write_i32(w, GROMACS_MAGIC)?;

    // slen = strlen("GMX_trn_file") + 1 = 13
    xdr_write_i32(w, 13)?;

    // XDR string: length field (12) + 12 bytes (no padding needed since 12 % 4 == 0)
    xdr_write_i32(w, TRR_VERSION.len() as i32)?;
    w.write_all(TRR_VERSION.as_bytes())?;

    xdr_write_i32(w, 0)?; // ir_size
    xdr_write_i32(w, 0)?; // e_size
    xdr_write_i32(w, h.box_size)?;
    xdr_write_i32(w, 0)?; // vir_size
    xdr_write_i32(w, 0)?; // pres_size
    xdr_write_i32(w, 0)?; // top_size
    xdr_write_i32(w, 0)?; // sym_size
    xdr_write_i32(w, h.x_size)?;
    xdr_write_i32(w, h.v_size)?;
    xdr_write_i32(w, h.f_size)?;
    xdr_write_i32(w, h.natoms)?;
    xdr_write_i32(w, h.step)?;
    xdr_write_i32(w, 0)?; // nre

    // time and lambda (always f32 when writing)
    xdr_write_f32(w, h.time)?;
    xdr_write_f32(w, 0.0)?; // lambda

    Ok(())
}

/// Compute data size for a single frame (without reading data).
/// Returns None if we can't determine size (shouldn't happen for well-formed files).
fn frame_data_size(h: &TrrHeader) -> usize {
    let elem = if h.b_double { 8usize } else { 4usize };
    let mut size = 0usize;
    if h.box_size != 0 { size += 9 * elem; }
    if h.vir_size != 0 { size += 9 * elem; }
    if h.pres_size != 0 { size += 9 * elem; }
    let n3 = h.natoms as usize * 3;
    if h.x_size != 0 { size += n3 * elem; }
    if h.v_size != 0 { size += n3 * elem; }
    if h.f_size != 0 { size += n3 * elem; }
    size
}

//=============================== Structs ==================================

pub(crate) struct TrrReader {
    file: BufReader<File>,
    cur_frame: usize,
}

pub(crate) struct TrrWriter {
    file: BufWriter<File>,
    cur_frame: usize,
    natoms: Option<usize>,
}

pub enum TrrFileHandler {
    Reader(TrrReader),
    Writer(TrrWriter),
}

//=============================== Helper to read a frame's data ============

/// Read a frame's data from current position (after header has been read).
/// The bool flags control which fields are populated; unwanted data is skipped at the I/O level.
fn read_frame_data(
    r: &mut impl Read,
    h: &TrrHeader,
    coords: bool,
    velocities: bool,
    forces: bool,
) -> Result<State, TrrHandlerError> {
    let elem = if h.b_double { 8usize } else { 4usize };
    let n = h.natoms as usize;

    // Box (TRR is f32 or f64 on disk; cast to Float at the boundary)
    let pbox = if h.box_size != 0 {
        let vals: Vec<Float> = if h.b_double {
            (0..9).map(|_| xdr_read_f64(r).map(|v| v as Float)).collect::<io::Result<_>>()?
        } else {
            (0..9).map(|_| xdr_read_f32(r).map(|v| v as Float)).collect::<io::Result<_>>()?
        };
        let m = Matrix3f::from_iterator(vals.into_iter());
        Some(PeriodicBox::from_matrix(m)?)
    } else {
        None
    };

    // vir (skip)
    if h.vir_size != 0 {
        xdr_skip(r, 9 * elem)?;
    }

    // pres (skip)
    if h.pres_size != 0 {
        xdr_skip(r, 9 * elem)?;
    }

    // Coordinates (cast to Float at the boundary)
    let coord_data: Vec<Pos> = if h.x_size != 0 {
        if coords {
            if h.b_double {
                (0..n).map(|_| {
                    let x = xdr_read_f64(r)? as Float;
                    let y = xdr_read_f64(r)? as Float;
                    let z = xdr_read_f64(r)? as Float;
                    Ok(Pos::new(x, y, z))
                }).collect::<io::Result<_>>()?
            } else {
                (0..n).map(|_| {
                    let x = xdr_read_f32(r)? as Float;
                    let y = xdr_read_f32(r)? as Float;
                    let z = xdr_read_f32(r)? as Float;
                    Ok(Pos::new(x, y, z))
                }).collect::<io::Result<_>>()?
            }
        } else {
            xdr_skip(r, n * 3 * elem)?;
            Vec::new()
        }
    } else {
        Vec::new()
    };

    // Velocities
    let vel_data: Vec<Vel> = if h.v_size != 0 {
        if velocities {
            read_xvf(r, n, h.b_double)?
        } else {
            xdr_skip(r, n * 3 * elem)?;
            Vec::new()
        }
    } else {
        Vec::new()
    };

    // Forces
    let force_data: Vec<Force> = if h.f_size != 0 {
        if forces {
            read_xvf(r, n, h.b_double)?
        } else {
            xdr_skip(r, n * 3 * elem)?;
            Vec::new()
        }
    } else {
        Vec::new()
    };

    Ok(State {
        coords: coord_data,
        velocities: vel_data,
        forces: force_data,
        pbox,
        time: h.time as Float,
    })
}

/// Skip a frame's data section (already read header).
fn skip_frame_data(r: &mut (impl Read + Seek), h: &TrrHeader) -> io::Result<()> {
    let size = frame_data_size(h) as i64;
    r.seek(SeekFrom::Current(size))?;
    Ok(())
}

//=============================== FileFormatHandler impl ===================

impl FileFormatHandler for TrrFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let fname = fname.as_ref();
        let file = File::open(fname).map_err(TrrHandlerError::Io)?;
        let mut reader = BufReader::new(file);

        // Read first header to get natoms, then rewind
        let h = read_trr_header(&mut reader).map_err(trr_to_ff_err)?;
        let natoms = h.natoms as usize;
        reader.seek(SeekFrom::Start(0)).map_err(TrrHandlerError::Io)?;

        let _ = natoms; // validated but not stored
        Ok(TrrFileHandler::Reader(TrrReader {
            file: reader,
            cur_frame: 0,
        }))
    }

    fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let fname = fname.as_ref();
        let file = File::create(fname).map_err(TrrHandlerError::Io)?;
        Ok(TrrFileHandler::Writer(TrrWriter {
            file: BufWriter::new(file),
            cur_frame: 0,
            natoms: None,
        }))
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        let TrrFileHandler::Reader(ref mut r) = self else {
            return Err(FileFormatError::NotStateReadFormat);
        };

        let h = read_trr_header(&mut r.file).map_err(trr_to_ff_err)?;
        let st = read_frame_data(&mut r.file, &h, true, true, true).map_err(trr_to_ff_err)?;
        r.cur_frame += 1;
        Ok(st)
    }

    fn read_state_pick(&mut self, coords: bool, velocities: bool, forces: bool) -> Result<State, FileFormatError> {
        if !coords {
            return Err(FileFormatError::NoCoords);
        }
        let TrrFileHandler::Reader(ref mut r) = self else {
            return Err(FileFormatError::NotStateReadFormat);
        };

        let h = read_trr_header(&mut r.file).map_err(trr_to_ff_err)?;
        let st = read_frame_data(&mut r.file, &h, coords, velocities, forces).map_err(trr_to_ff_err)?;
        r.cur_frame += 1;
        Ok(st)
    }

    fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileFormatError> {
        self.write_state_pick(data, true, true, true)
    }

    fn write_state_pick(&mut self, data: &dyn SaveState, coords: bool, velocities: bool, forces: bool) -> Result<(), FileFormatError> {
        let TrrFileHandler::Writer(ref mut w) = self else {
            return Err(FileFormatError::NotStateWriteFormat);
        };

        let natoms = data.len();
        if w.natoms.is_none() {
            w.natoms = Some(natoms);
        }

        let has_box = data.get_box().is_some();
        let box_size = if has_box { 9 * 4i32 } else { 0i32 };
        let x_size = if coords     { natoms as i32 * 3 * 4 } else { 0 };

        let vel_it   = if velocities { data.iter_vel_dyn()   } else { Box::new(std::iter::empty()) };
        let force_it = if forces     { data.iter_force_dyn() } else { Box::new(std::iter::empty()) };

        let v_size = if vel_it.len() > 0   { natoms as i32 * 3 * 4 } else { 0 };
        let f_size = if force_it.len() > 0 { natoms as i32 * 3 * 4 } else { 0 };

        let h = TrrHeader {
            box_size,
            vir_size: 0,
            pres_size: 0,
            x_size,
            v_size,
            f_size,
            natoms: natoms as i32,
            step: w.cur_frame as i32,
            time: data.get_time() as f32,
            b_double: false,
        };

        write_trr_header(&mut w.file, &h)?;

        // Write box (row-major, 9 f32). f64 builds downcast at the boundary.
        if let Some(b) = data.get_box() {
            for &v in b.get_matrix().as_slice() {
                xdr_write_f32(&mut w.file, v as f32)?;
            }
        }

        // Write coordinates
        if coords {
            for p in data.iter_pos_dyn() {
                xdr_write_f32(&mut w.file, p.x as f32)?;
                xdr_write_f32(&mut w.file, p.y as f32)?;
                xdr_write_f32(&mut w.file, p.z as f32)?;
            }
        }

        // Write velocities
        for v in vel_it {
            xdr_write_f32(&mut w.file, v.x as f32)?;
            xdr_write_f32(&mut w.file, v.y as f32)?;
            xdr_write_f32(&mut w.file, v.z as f32)?;
        }

        // Write forces
        for f in force_it {
            xdr_write_f32(&mut w.file, f.x as f32)?;
            xdr_write_f32(&mut w.file, f.y as f32)?;
            xdr_write_f32(&mut w.file, f.z as f32)?;
        }

        w.cur_frame += 1;
        Ok(())
    }

    fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
        let TrrFileHandler::Reader(ref mut r) = self else {
            return Err(FileFormatError::NotRandomAccessFormat);
        };

        if fr == r.cur_frame {
            return Ok(());
        }

        // If seeking backwards, rewind to start
        if fr < r.cur_frame {
            r.file
                .seek(SeekFrom::Start(0))
                .map_err(|e| TrrHandlerError::SeekFrame(fr, e))?;
            r.cur_frame = 0;
        }

        // Skip forward fr - cur_frame frames
        let skip = fr - r.cur_frame;
        for _ in 0..skip {
            let h = read_trr_header(&mut r.file)
                .map_err(|e| match e {
                    TrrHandlerError::Io(io_e) => TrrHandlerError::Io(io_e),
                    other => TrrHandlerError::Io(io::Error::new(io::ErrorKind::Other, other.to_string())),
                })?;
            skip_frame_data(&mut r.file, &h)
                .map_err(|e| TrrHandlerError::SeekFrame(fr, e))?;
            r.cur_frame += 1;
        }

        Ok(())
    }

    fn seek_last(&mut self) -> Result<(), FileFormatError> {
        let TrrFileHandler::Reader(ref mut r) = self else {
            return Err(FileFormatError::NotRandomAccessFormat);
        };

        // Rewind and scan all frames, track last valid frame start position
        r.file
            .seek(SeekFrom::Start(0))
            .map_err(TrrHandlerError::Io)?;
        r.cur_frame = 0;

        let mut last_pos: Option<u64> = None;
        let mut last_frame: usize = 0;
        let mut frame_idx: usize = 0;

        loop {
            let pos = r.file.stream_position().map_err(TrrHandlerError::Io)?;
            match read_trr_header(&mut r.file) {
                Ok(h) => {
                    last_pos = Some(pos);
                    last_frame = frame_idx;
                    match skip_frame_data(&mut r.file, &h) {
                        Ok(()) => {}
                        Err(_) => break,
                    }
                    frame_idx += 1;
                }
                Err(TrrHandlerError::Io(e)) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(FileFormatError::from(e)),
            }
        }

        match last_pos {
            Some(pos) => {
                r.file
                    .seek(SeekFrom::Start(pos))
                    .map_err(TrrHandlerError::Io)?;
                r.cur_frame = last_frame;
                Ok(())
            }
            None => Err(FileFormatError::Eof),
        }
    }

    fn seek_time(&mut self, t: Float) -> Result<(), FileFormatError> {
        let TrrFileHandler::Reader(ref mut r) = self else {
            return Err(FileFormatError::NotRandomAccessFormat);
        };

        // Rewind then scan until frame.time >= t
        r.file
            .seek(SeekFrom::Start(0))
            .map_err(TrrHandlerError::Io)?;

        let mut frame_idx: usize = 0;
        // TRR header time is on-disk f32; compare in f32 space to match.
        let t_disk = t as f32;

        loop {
            let pos = r.file.stream_position().map_err(TrrHandlerError::Io)?;
            match read_trr_header(&mut r.file) {
                Ok(h) => {
                    if h.time >= t_disk {
                        // Seek back to this frame's start
                        r.file
                            .seek(SeekFrom::Start(pos))
                            .map_err(|e| TrrHandlerError::SeekTime(t, e))?;
                        r.cur_frame = frame_idx;
                        return Ok(());
                    }
                    skip_frame_data(&mut r.file, &h)
                        .map_err(|e| TrrHandlerError::SeekTime(t, e))?;
                    frame_idx += 1;
                }
                Err(TrrHandlerError::Io(e)) if e.kind() == io::ErrorKind::UnexpectedEof => {
                    return Err(FileFormatError::Eof);
                }
                Err(e) => return Err(FileFormatError::from(e)),
            }
        }
    }
}

/// Map TrrHandlerError to FileFormatError, converting UnexpectedEof to FileFormatError::Eof.
fn trr_to_ff_err(e: TrrHandlerError) -> FileFormatError {
    if let TrrHandlerError::Io(ref io_err) = e {
        if io_err.kind() == io::ErrorKind::UnexpectedEof {
            return FileFormatError::Eof;
        }
    }
    FileFormatError::from(e)
}
