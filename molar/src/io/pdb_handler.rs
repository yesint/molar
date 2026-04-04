use crate::atom::element_symbol;
use crate::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
    path::Path,
};
use thiserror::Error;

pub struct PdbFileHandler {
    // Reading
    reader: Option<BufReader<File>>,
    bonds: Vec<[usize; 2]>, // 0-indexed, collected from CONECT records
    stored_topology: Option<Topology>,
    stored_state: Option<State>,
    at_least_one_state_read: bool,
    // Writing
    writer: Option<BufWriter<File>>,
}

#[derive(Debug, Error)]
pub enum PdbHandlerError {
    #[error("can't open pdb file for reading")]
    OpenRead(#[source] std::io::Error),

    #[error("can't open pdb file for writing")]
    OpenWrite(#[source] std::io::Error),

    #[error("pdb file is empty or contains no ATOM/HETATM records")]
    Empty,

    #[error("unexpected end of file reached")]
    Eof,

    #[error(transparent)]
    ParseInt(#[from] ParseIntError),

    #[error(transparent)]
    ParseFloat(#[from] ParseFloatError),

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),

    #[error("io error")]
    Io(#[from] std::io::Error),
}

// ── helper parsers ──────────────────────────────────────────────────────────

fn parse_f32(line: &str, start: usize, end: usize) -> Result<f32, PdbHandlerError> {
    line.get(start..end)
        .unwrap_or("")
        .trim()
        .parse::<f32>()
        .map_err(PdbHandlerError::ParseFloat)
}

fn parse_f32_opt(line: &str, start: usize, end: usize) -> Option<f32> {
    line.get(start..end)?.trim().parse::<f32>().ok()
}

fn parse_i32_opt(line: &str, start: usize, end: usize) -> Option<i32> {
    line.get(start..end)?.trim().parse::<i32>().ok()
}

/// Parse a 1-indexed PDB atom serial from fixed-width column, convert to 0-indexed.
fn parse_serial_opt(line: &str, start: usize, end: usize) -> Option<usize> {
    let s = line.get(start..end)?.trim();
    if s.is_empty() {
        return None;
    }
    s.parse::<usize>().ok().map(|n| n.saturating_sub(1))
}

/// Format PDB atom name:
/// - ≤3 chars: prefix with a space → " CA "
/// - 4 chars: use as-is
fn format_atom_name(name: &str) -> String {
    if name.len() >= 4 {
        format!("{:<4}", &name[..4])
    } else {
        format!(" {:<3}", name)
    }
}

// ── FileFormatHandler ────────────────────────────────────────────────────────

impl FileFormatHandler for PdbFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Ok(Self {
            reader: Some(BufReader::new(
                File::open(fname).map_err(PdbHandlerError::OpenRead)?,
            )),
            writer: None,
            bonds: Vec::new(),
            stored_topology: None,
            stored_state: None,
            at_least_one_state_read: false,
        })
    }

    fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Ok(Self {
            reader: None,
            writer: Some(BufWriter::new(
                File::create(fname).map_err(PdbHandlerError::OpenWrite)?,
            )),
            bonds: Vec::new(),
            stored_topology: None,
            stored_state: None,
            at_least_one_state_read: false,
        })
    }

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        let buf = self.reader.as_mut().unwrap();

        // Check if we are at EOF
        if buf.fill_buf()?.is_empty() {
            return if self.at_least_one_state_read {
                Err(PdbHandlerError::Eof)?
            } else {
                Err(PdbHandlerError::Empty)?
            };
        }

        // Fresh bond list for this model
        self.bonds.clear();

        let mut atoms: Vec<Atom> = Vec::new();
        let mut coords: Vec<Pos> = Vec::new();
        let mut current_box: Option<PeriodicBox> = None;
        let mut has_atoms = false;
        let mut line = String::new();

        loop {
            line.clear();
            let n = buf.read_line(&mut line)?;
            if n == 0 {
                break;
            }

            if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
                has_atoms = true;
                let x = parse_f32(&line, 30, 38).unwrap_or(0.0);
                let y = parse_f32(&line, 38, 46).unwrap_or(0.0);
                let z = parse_f32(&line, 46, 54).unwrap_or(0.0);
                coords.push(Pos::new(x * 0.1, y * 0.1, z * 0.1));

                let name = line.get(12..16).unwrap_or("").trim();
                let resname = line.get(17..20).unwrap_or("").trim();
                let chain = line
                    .get(21..22)
                    .and_then(|s| s.chars().next())
                    .unwrap_or(' ');
                let resid = parse_i32_opt(&line, 22, 26).unwrap_or(0);
                let occupancy = parse_f32_opt(&line, 54, 60).unwrap_or(1.0);
                let bfactor = parse_f32_opt(&line, 60, 66).unwrap_or(0.0);

                atoms.push(
                    Atom::new()
                        .with_name(name)
                        .with_resname(resname)
                        .with_resid(resid)
                        .with_chain(chain)
                        .with_occupancy(occupancy)
                        .with_bfactor(bfactor)
                        .guess()
                );
            } else if line.starts_with("CRYST1") {
                let a = parse_f32(&line, 6, 15).unwrap_or(0.0);
                let b = parse_f32(&line, 15, 24).unwrap_or(0.0);
                let c = parse_f32(&line, 24, 33).unwrap_or(0.0);
                let alpha = parse_f32(&line, 33, 40).unwrap_or(90.0);
                let beta = parse_f32(&line, 40, 47).unwrap_or(90.0);
                let gamma = parse_f32(&line, 47, 54).unwrap_or(90.0);
                current_box = PeriodicBox::from_vectors_angles(
                    a * 0.1,
                    b * 0.1,
                    c * 0.1,
                    alpha,
                    beta,
                    gamma,
                )
                .ok();
            } else if line.starts_with("MODEL") {
                if has_atoms {
                    // New MODEL starts — previous model is complete (ENDMDL absent)
                    break;
                }
            } else if line.starts_with("ENDMDL") {
                break;
            } else if line.starts_with("CONECT") {
                if let Some(a) = parse_serial_opt(&line, 6, 11) {
                    for (s, e) in [(11usize, 16usize), (16, 21), (21, 26), (26, 31)] {
                        if let Some(b) = parse_serial_opt(&line, s, e) {
                            if a != b {
                                let mut pair = [a, b];
                                pair.sort();
                                self.bonds.push(pair);
                            }
                        }
                    }
                }
            } else if line.starts_with("END") && !line.starts_with("ENDMDL") {
                break;
            }
        }

        if !has_atoms {
            return if self.at_least_one_state_read {
                Err(PdbHandlerError::Eof)?
            } else {
                Err(PdbHandlerError::Empty)?
            };
        }

        self.bonds.sort();
        self.bonds.dedup();

        let mut top = Topology::default();
        top.atoms = atoms;
        top.bonds = self.bonds.clone();
        top.assign_resindex();

        self.at_least_one_state_read = true;

        Ok((top, State { coords, time: 0.0, pbox: current_box, ..Default::default() }))
    }

    fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        if self.stored_topology.is_some() {
            Ok(self.stored_topology.take().unwrap())
        } else {
            let (top, st) = self.read()?;
            self.stored_state.get_or_insert(st);
            Ok(top)
        }
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        if self.stored_state.is_some() {
            Ok(self.stored_state.take().unwrap())
        } else {
            let (top, st) = self.read()?;
            self.stored_topology.get_or_insert(top);
            Ok(st)
        }
    }

    fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileFormatError> {
        let w = self.writer.as_mut().ok_or(FileFormatError::NotWritable)?;
        write_cryst1(w, data.get_box())?;
        for (i, (at, pos)) in data.iter_atoms_dyn().zip(data.iter_pos_dyn()).enumerate() {
            write_atom_record(w, (i % 99999) + 1, at, pos)?;
        }
        writeln!(w, "END")?;
        Ok(())
    }
}

// ── write helpers ────────────────────────────────────────────────────────────

fn write_cryst1(
    w: &mut BufWriter<File>,
    pbox: Option<&PeriodicBox>,
) -> Result<(), FileFormatError> {
    if let Some(b) = pbox {
        let (lengths, angles) = b.to_vectors_angles();
        writeln!(
            w,
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} P 1           1",
            lengths[0] * 10.0,
            lengths[1] * 10.0,
            lengths[2] * 10.0,
            angles[0],
            angles[1],
            angles[2]
        )?;
    }
    Ok(())
}

fn write_atom_record(
    w: &mut BufWriter<File>,
    serial: usize,
    at: &Atom,
    pos: &Pos,
) -> Result<(), FileFormatError> {
    let name = format_atom_name(at.name.as_str());
    let elem = element_symbol(at.atomic_number);
    writeln!(
        w,
        "ATOM  {:>5} {:<4} {:<3} {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
        serial,
        name,
        at.resname,
        at.chain,
        at.resid % 9999,
        pos.x * 10.0,
        pos.y * 10.0,
        pos.z * 10.0,
        at.occupancy,
        at.bfactor,
        elem
    )?;
    Ok(())
}
