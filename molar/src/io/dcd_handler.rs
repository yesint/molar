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
pub enum DcdHandlerError {
    #[error("io error")]
    Io(#[from] io::Error),

    #[error("not a DCD file (bad record length or magic)")]
    BadMagic,

    #[error("unexpected record length: expected {0}, got {1}")]
    BadRecord(u32, u32),

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),
}

//=============================== Structs ==================================

pub struct DcdReader {
    reader: BufReader<File>,
    n_atoms: usize,
    n_frames: usize,
    cur_frame: usize,
    has_extra_block: bool,
    has_4dims: bool,
    reverse_endian: bool,
    n_fixed: usize,
    free_indices: Vec<usize>,    // 0-indexed free atom indices
    fixed_coords: Vec<[f32; 3]>, // fixed atom coords in nm (set after first frame)
    frame_offset: u64,           // file offset of first frame
    frame_size: u64,             // bytes per frame (0 = can't seek)
    // For physical_time computation
    istart: i32,
    nsavc: i32,
    delta: f32,
}

pub struct DcdWriter {
    writer: BufWriter<File>,
    n_atoms: usize,          // 0 = uninitialized (header not yet written)
    n_frames_written: usize,
}

pub enum DcdFileHandler {
    Reader(DcdReader),
    Writer(DcdWriter),
}

//=============================== Low-level helpers ========================

fn read_marker(r: &mut impl Read, swap: bool) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    let val = u32::from_le_bytes(buf);
    Ok(if swap { val.swap_bytes() } else { val })
}

fn read_byte_record(r: &mut impl Read, swap: bool) -> Result<Vec<u8>, DcdHandlerError> {
    let len = read_marker(r, swap)?;
    let mut bytes = vec![0u8; len as usize];
    r.read_exact(&mut bytes)?;
    let end = read_marker(r, swap)?;
    if end != len {
        return Err(DcdHandlerError::BadRecord(len, end));
    }
    Ok(bytes)
}

fn read_f32_record(r: &mut impl Read, n: usize, swap: bool) -> Result<Vec<f32>, DcdHandlerError> {
    let len = read_marker(r, swap)?;
    let expected = (n * 4) as u32;
    if len != expected {
        return Err(DcdHandlerError::BadRecord(expected, len));
    }
    let mut bytes = vec![0u8; n * 4];
    r.read_exact(&mut bytes)?;
    let data: Vec<f32> = bytes
        .chunks_exact(4)
        .map(|c| {
            let arr: [u8; 4] = c.try_into().unwrap();
            if swap {
                f32::from_bits(u32::from_be_bytes(arr))
            } else {
                f32::from_le_bytes(arr)
            }
        })
        .collect();
    let end = read_marker(r, swap)?;
    if end != len {
        return Err(DcdHandlerError::BadRecord(len, end));
    }
    Ok(data)
}

fn read_f64_record(r: &mut impl Read, n: usize, swap: bool) -> Result<Vec<f64>, DcdHandlerError> {
    let len = read_marker(r, swap)?;
    let expected = (n * 8) as u32;
    if len != expected {
        return Err(DcdHandlerError::BadRecord(expected, len));
    }
    let mut bytes = vec![0u8; n * 8];
    r.read_exact(&mut bytes)?;
    let data: Vec<f64> = bytes
        .chunks_exact(8)
        .map(|c| {
            let arr: [u8; 8] = c.try_into().unwrap();
            if swap {
                f64::from_bits(u64::from_be_bytes(arr))
            } else {
                f64::from_le_bytes(arr)
            }
        })
        .collect();
    let end = read_marker(r, swap)?;
    if end != len {
        return Err(DcdHandlerError::BadRecord(len, end));
    }
    Ok(data)
}

fn write_record(w: &mut impl Write, data: &[u8]) -> io::Result<()> {
    let len = data.len() as u32;
    w.write_all(&len.to_le_bytes())?;
    w.write_all(data)?;
    w.write_all(&len.to_le_bytes())?;
    Ok(())
}

fn i32_from_header(buf: &[u8], offset: usize, swap: bool) -> i32 {
    let arr: [u8; 4] = buf[offset..offset + 4].try_into().unwrap();
    let val = u32::from_le_bytes(arr);
    (if swap { val.swap_bytes() } else { val }) as i32
}

fn f32_from_header(buf: &[u8], offset: usize, swap: bool) -> f32 {
    let arr: [u8; 4] = buf[offset..offset + 4].try_into().unwrap();
    if swap {
        f32::from_bits(u32::from_be_bytes(arr))
    } else {
        f32::from_le_bytes(arr)
    }
}

fn f64_from_header(buf: &[u8], offset: usize, swap: bool) -> f64 {
    let arr: [u8; 8] = buf[offset..offset + 8].try_into().unwrap();
    if swap {
        f64::from_bits(u64::from_be_bytes(arr))
    } else {
        f64::from_le_bytes(arr)
    }
}

fn f32_slice_to_bytes(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|&f| f.to_le_bytes()).collect()
}

fn encode_unit_cell(pbox: Option<&PeriodicBox>) -> [u8; 48] {
    let mut buf = [0u8; 48];
    if let Some(b) = pbox {
        let (lengths, angles) = b.to_vectors_angles();
        // DCD unit cell: [A, cos(γ), B, cos(β), cos(α), C] in Ångström + cosines
        let cell = [
            (lengths[0] * 10.0) as f64,           // A in Å
            (angles[2] as f64).to_radians().cos(), // cos(γ)
            (lengths[1] * 10.0) as f64,            // B in Å
            (angles[1] as f64).to_radians().cos(), // cos(β)
            (angles[0] as f64).to_radians().cos(), // cos(α)
            (lengths[2] * 10.0) as f64,            // C in Å
        ];
        for (i, &v) in cell.iter().enumerate() {
            buf[i * 8..(i + 1) * 8].copy_from_slice(&v.to_le_bytes());
        }
    }
    buf
}

fn parse_unit_cell(cell: &[f64]) -> Result<Option<PeriodicBox>, DcdHandlerError> {
    let a = cell[0];
    let b = cell[2];
    let c = cell[5];
    if a == 0.0 {
        return Ok(None);
    }
    // cell[1]=cos(γ) or γ in degrees, cell[3]=cos(β) or β, cell[4]=cos(α) or α
    // Detect cosine vs degrees: if |value| <= 1.0 it's a cosine
    let (alpha, beta, gamma) = if cell[4].abs() <= 1.0 {
        (
            cell[4].acos().to_degrees(),
            cell[3].acos().to_degrees(),
            cell[1].acos().to_degrees(),
        )
    } else {
        // Legacy NAMD: stored as degrees
        (cell[4], cell[3], cell[1])
    };
    Ok(Some(PeriodicBox::from_vectors_angles(
        (a * 0.1) as f32,
        (b * 0.1) as f32,
        (c * 0.1) as f32,
        alpha as f32,
        beta as f32,
        gamma as f32,
    )?))
}

fn write_dcd_header(w: &mut impl Write, n_atoms: usize) -> io::Result<()> {
    // Block 1: 84-byte CORD header (CHARMM format)
    let mut header = [0u8; 84];
    header[0..4].copy_from_slice(b"CORD");
    // NSET = 0 at [4..8] initially (updated on each frame write)
    // ISTART = 0 at [8..12]
    // NSAVC = 1 at [12..16]
    header[12..16].copy_from_slice(&1i32.to_le_bytes());
    // NAMNF = 0 at [32..36] (no fixed atoms)
    // DELTA = 0.0f32 at [36..40]
    // extra_block = 1 at [40..44] (always write unit cell)
    header[40..44].copy_from_slice(&1i32.to_le_bytes());
    // has_4dims = 0 at [44..48]
    // CHARMM version = 24 at [76..80]
    header[76..80].copy_from_slice(&24i32.to_le_bytes());
    write_record(w, &header)?;

    // Block 2: title block (NTITLE=1 + 80 bytes of empty title)
    let mut title_block = [0u8; 84];
    title_block[0..4].copy_from_slice(&1i32.to_le_bytes()); // NTITLE=1
    write_record(w, &title_block)?;

    // Block 3: atom count
    write_record(w, &(n_atoms as i32).to_le_bytes())?;
    Ok(())
}

/// Map DcdHandlerError to FileFormatError, converting UnexpectedEof to FileFormatError::Eof.
fn dcd_to_ff_err(e: DcdHandlerError) -> FileFormatError {
    if let DcdHandlerError::Io(ref io_err) = e {
        if io_err.kind() == io::ErrorKind::UnexpectedEof {
            return FileFormatError::Eof;
        }
    }
    FileFormatError::from(e)
}

//=============================== FileFormatHandler impl ===================

impl FileFormatHandler for DcdFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let fname = fname.as_ref();
        let file = File::open(fname).map_err(|e| DcdHandlerError::Io(e))?;
        let mut reader = BufReader::new(file);

        // Endianness detection via first 4-byte record marker (must equal 84)
        let mut marker_buf = [0u8; 4];
        reader
            .read_exact(&mut marker_buf)
            .map_err(|e| DcdHandlerError::Io(e))?;
        let le_val = u32::from_le_bytes(marker_buf);
        let be_val = u32::from_be_bytes(marker_buf);
        let reverse_endian = if le_val == 84 {
            false
        } else if be_val == 84 {
            true
        } else {
            return Err(DcdHandlerError::BadMagic.into());
        };

        // Read the 84-byte header payload
        let mut header = [0u8; 84];
        reader
            .read_exact(&mut header)
            .map_err(|e| DcdHandlerError::Io(e))?;
        let end = read_marker(&mut reader, reverse_endian).map_err(|e| DcdHandlerError::Io(e))?;
        if end != 84 {
            return Err(DcdHandlerError::BadRecord(84, end).into());
        }

        if &header[0..4] != b"CORD" {
            return Err(DcdHandlerError::BadMagic.into());
        }

        let n_frames = i32_from_header(&header, 4, reverse_endian) as usize;
        let istart = i32_from_header(&header, 8, reverse_endian);
        let nsavc = i32_from_header(&header, 12, reverse_endian);
        let n_fixed = i32_from_header(&header, 32, reverse_endian).max(0) as usize;
        let is_charmm = i32_from_header(&header, 76, reverse_endian) != 0;
        let has_extra_block = i32_from_header(&header, 40, reverse_endian) != 0;
        let has_4dims = i32_from_header(&header, 44, reverse_endian) != 0;

        let delta: f32 = if is_charmm {
            f32_from_header(&header, 36, reverse_endian)
        } else {
            f64_from_header(&header, 36, reverse_endian) as f32
        };

        // Skip title block
        read_byte_record(&mut reader, reverse_endian)?;

        // Read atom count block
        let atom_bytes = read_byte_record(&mut reader, reverse_endian)?;
        if atom_bytes.len() < 4 {
            return Err(DcdHandlerError::BadRecord(4, atom_bytes.len() as u32).into());
        }
        let n_atoms = i32_from_header(&atom_bytes, 0, reverse_endian).max(0) as usize;

        // Read fixed atom index block if n_fixed > 0
        let mut free_indices: Vec<usize> = Vec::new();
        if n_fixed > 0 {
            let idx_bytes = read_byte_record(&mut reader, reverse_endian)?;
            let n_free = n_atoms - n_fixed;
            free_indices = idx_bytes
                .chunks_exact(4)
                .take(n_free)
                .map(|c| {
                    let arr: [u8; 4] = c.try_into().unwrap();
                    let val = if reverse_endian {
                        u32::from_be_bytes(arr)
                    } else {
                        u32::from_le_bytes(arr)
                    };
                    val.saturating_sub(1) as usize // convert 1-indexed to 0-indexed
                })
                .collect();
        }

        // Record offset of first frame
        let frame_offset = reader
            .stream_position()
            .map_err(|e| DcdHandlerError::Io(e))?;

        // Compute frame_size for seek support (only when n_fixed==0 and !has_4dims)
        let frame_size: u64 = if n_fixed == 0 && !has_4dims {
            let extra = if has_extra_block { 56u64 } else { 0u64 }; // 4+48+4
            extra + 3 * (n_atoms as u64 * 4 + 8)
        } else {
            0 // can't seek
        };

        Ok(DcdFileHandler::Reader(DcdReader {
            reader,
            n_atoms,
            n_frames,
            cur_frame: 0,
            has_extra_block,
            has_4dims,
            reverse_endian,
            n_fixed,
            free_indices,
            fixed_coords: Vec::new(),
            frame_offset,
            frame_size,
            istart,
            nsavc,
            delta,
        }))
    }

    fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let file = File::create(fname.as_ref()).map_err(|e| DcdHandlerError::Io(e))?;
        Ok(DcdFileHandler::Writer(DcdWriter {
            writer: BufWriter::new(file),
            n_atoms: 0,
            n_frames_written: 0,
        }))
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        let DcdFileHandler::Reader(ref mut r) = self else {
            return Err(FileFormatError::NotStateReadFormat);
        };

        let swap = r.reverse_endian;
        let is_first_frame = r.cur_frame == 0;
        let n_read = if r.n_fixed > 0 && !is_first_frame {
            r.n_atoms - r.n_fixed
        } else {
            r.n_atoms
        };

        // Unit cell block — catch EOF at start of frame
        let pbox = if r.has_extra_block {
            let cell = read_f64_record(&mut r.reader, 6, swap).map_err(dcd_to_ff_err)?;
            parse_unit_cell(&cell)?
        } else {
            None
        };

        // X block — catch EOF here when no extra block
        let x = read_f32_record(&mut r.reader, n_read, swap).map_err(|e| {
            if r.has_extra_block {
                FileFormatError::from(e)
            } else {
                dcd_to_ff_err(e)
            }
        })?;
        let y = read_f32_record(&mut r.reader, n_read, swap)?;
        let z = read_f32_record(&mut r.reader, n_read, swap)?;

        // Skip 4D block if present
        if r.has_4dims {
            read_byte_record(&mut r.reader, swap)?;
        }

        // Assemble coordinates in nm (Å × 0.1)
        let coords: Vec<Pos> = if r.n_fixed == 0 || is_first_frame {
            // All atoms present in the coordinate arrays
            (0..r.n_atoms)
                .map(|i| Pos::new(x[i] * 0.1, y[i] * 0.1, z[i] * 0.1))
                .collect()
        } else {
            // Merge free (from file) and fixed (stored) atoms
            let free_xyz: Vec<(f32, f32, f32)> = x
                .iter()
                .zip(y.iter())
                .zip(z.iter())
                .map(|((&xi, &yi), &zi)| (xi, yi, zi))
                .collect();
            let mut free_iter = free_xyz.into_iter();
            let mut fi = 0usize; // index into free_indices
            (0..r.n_atoms)
                .map(|i| {
                    if fi < r.free_indices.len() && r.free_indices[fi] == i {
                        fi += 1;
                        let (xi, yi, zi) = free_iter.next().unwrap();
                        Pos::new(xi * 0.1, yi * 0.1, zi * 0.1)
                    } else {
                        let fc = r.fixed_coords[i];
                        Pos::new(fc[0], fc[1], fc[2]) // already in nm
                    }
                })
                .collect()
        };

        // Save fixed_coords after first frame
        if is_first_frame && r.n_fixed > 0 {
            r.fixed_coords = coords.iter().map(|p| [p.x, p.y, p.z]).collect();
        }

        let time = (r.istart + r.cur_frame as i32 * r.nsavc) as f32 * r.delta;
        r.cur_frame += 1;

        Ok(State { coords, pbox, time })
    }

    fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileFormatError> {
        let DcdFileHandler::Writer(ref mut w) = self else {
            return Err(FileFormatError::NotStateWriteFormat);
        };

        // Write header on first call (lazy init)
        if w.n_atoms == 0 {
            w.n_atoms = data.len();
            write_dcd_header(&mut w.writer, w.n_atoms).map_err(DcdHandlerError::Io)?;
        }

        let n = w.n_atoms;

        // Unit cell block (always written, zeros if no box)
        let cell_bytes = encode_unit_cell(data.get_box());
        write_record(&mut w.writer, &cell_bytes).map_err(DcdHandlerError::Io)?;

        // Collect X, Y, Z arrays in Å (nm × 10)
        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);
        let mut zs = Vec::with_capacity(n);
        for p in data.iter_pos_dyn() {
            xs.push(p.x * 10.0);
            ys.push(p.y * 10.0);
            zs.push(p.z * 10.0);
        }

        write_record(&mut w.writer, &f32_slice_to_bytes(&xs)).map_err(DcdHandlerError::Io)?;
        write_record(&mut w.writer, &f32_slice_to_bytes(&ys)).map_err(DcdHandlerError::Io)?;
        write_record(&mut w.writer, &f32_slice_to_bytes(&zs)).map_err(DcdHandlerError::Io)?;

        // Update NSET counter at file offset 8 (inside the CORD header payload)
        w.n_frames_written += 1;
        w.writer.seek(SeekFrom::Start(8)).map_err(DcdHandlerError::Io)?;
        w.writer
            .write_all(&(w.n_frames_written as i32).to_le_bytes())
            .map_err(DcdHandlerError::Io)?;
        w.writer.seek(SeekFrom::End(0)).map_err(DcdHandlerError::Io)?;

        Ok(())
    }

    fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
        let DcdFileHandler::Reader(ref mut r) = self else {
            return Err(FileFormatError::NotRandomAccessFormat);
        };
        if r.frame_size == 0 {
            // n_fixed > 0 or has_4dims — can't seek
            return Err(FileFormatError::NotRandomAccessFormat);
        }
        if fr >= r.n_frames {
            return Err(FileFormatError::SeekFrame(fr));
        }
        let offset = r.frame_offset + fr as u64 * r.frame_size;
        r.reader
            .seek(SeekFrom::Start(offset))
            .map_err(DcdHandlerError::Io)?;
        r.cur_frame = fr;
        Ok(())
    }

    fn seek_last(&mut self) -> Result<(), FileFormatError> {
        let n_frames = match self {
            DcdFileHandler::Reader(ref r) => r.n_frames,
            _ => return Err(FileFormatError::NotRandomAccessFormat),
        };
        if n_frames == 0 {
            return Err(FileFormatError::Eof);
        }
        self.seek_frame(n_frames - 1)
    }
}
