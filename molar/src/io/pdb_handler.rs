use crate::atom::{element_symbol, AtomStr, ATOM_NAME_EXPECT, ATOM_RESNAME_EXPECT};
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
    bonds: Vec<[usize; 2]>,   // 0-indexed, collected inline when CONECT is encountered
    topology: Option<Topology>, // permanently stored after first parse; only cloned, not consumed
    stored_state: Option<State>, // cached state; consumed by read_state()
    topo_built: bool,
    exhausted: bool,
    // Writing
    writer: Option<BufWriter<File>>,
    serial: usize,              // running atom serial counter for write
    write_topo: Option<Vec<Atom>>, // buffered by write_topology() for use in write_state()
}

#[derive(Debug, Error)]
pub enum PdbHandlerError {
    #[error("can't open pdb file for reading")]
    OpenRead(#[source] std::io::Error),

    #[error("can't open pdb file for writing")]
    OpenWrite(#[source] std::io::Error),

    #[error("pdb file is empty or contains no ATOM/HETATM records")]
    Empty,

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

// ── core handler ────────────────────────────────────────────────────────────

impl PdbFileHandler {
    /// Stream-parse one MODEL block (or the whole file if no MODEL records).
    ///
    /// - On the first call (`!self.topo_built`): also builds `self.topology`.
    /// - Subsequent calls: only parse coordinates.
    /// - CONECT lines push bonds to `self.bonds`; they are applied to topology during the first call.
    /// - Returns `Ok(None)` when no ATOM/HETATM records were found (EOF or empty model).
    fn parse_next_model(&mut self) -> Result<Option<State>, PdbHandlerError> {
        let reader = match self.reader.as_mut() {
            Some(r) => r,
            None => return Ok(None),
        };

        let building_topo = !self.topo_built;
        let mut pending_atoms: Vec<Atom> = Vec::new();
        let mut coords: Vec<Pos> = Vec::new();
        let mut current_box: Option<PeriodicBox> = None;
        let mut has_atoms = false;
        let mut line = String::new();

        loop {
            line.clear();
            let n = reader.read_line(&mut line)?;
            if n == 0 {
                // True EOF
                self.exhausted = true;
                break;
            }

            if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
                has_atoms = true;
                let x = parse_f32(&line, 30, 38).unwrap_or(0.0);
                let y = parse_f32(&line, 38, 46).unwrap_or(0.0);
                let z = parse_f32(&line, 46, 54).unwrap_or(0.0);
                coords.push(Pos::new(x * 0.1, y * 0.1, z * 0.1));

                if building_topo {
                    let name = line.get(12..16).unwrap_or("").trim();
                    let resname = line.get(17..20).unwrap_or("").trim();
                    let chain = line
                        .get(21..22)
                        .and_then(|s| s.chars().next())
                        .unwrap_or(' ');
                    let resid = parse_i32_opt(&line, 22, 26).unwrap_or(0);
                    let occupancy = parse_f32_opt(&line, 54, 60).unwrap_or(1.0);
                    let bfactor = parse_f32_opt(&line, 60, 66).unwrap_or(0.0);

                    let mut at = Atom {
                        name: AtomStr::from_bytes(name.as_bytes()).expect(ATOM_NAME_EXPECT),
                        resname: AtomStr::from_bytes(resname.as_bytes())
                            .expect(ATOM_RESNAME_EXPECT),
                        resid,
                        chain,
                        occupancy,
                        bfactor,
                        type_name: AtomStr::from_bytes(b"").unwrap(),
                        ..Default::default()
                    };
                    at.guess_element_and_mass_from_name();
                    pending_atoms.push(at);
                }
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
                    // New model starts — previous model is complete (ENDMDL absent)
                    break;
                }
                // else: opening MODEL record before any atoms; keep going
            } else if line.starts_with("ENDMDL") {
                break;
            } else if line.starts_with("CONECT") {
                // Parse bond pairs: atom at [6..11] bonded to partners at [11..16],[16..21],[21..26],[26..31]
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
                self.exhausted = true;
                break;
            }
        }

        if !has_atoms {
            return Ok(None);
        }

        if building_topo {
            let mut top = Topology::default();
            top.atoms = pending_atoms;
            // Deduplicate and apply CONECT bonds
            self.bonds.sort();
            self.bonds.dedup();
            top.bonds = self.bonds.clone();
            top.assign_resindex();
            self.topology = Some(top);
            self.topo_built = true;
        }

        Ok(Some(State {
            coords,
            time: 0.0,
            pbox: current_box,
        }))
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
            topology: None,
            stored_state: None,
            topo_built: false,
            exhausted: false,
            serial: 0,
            write_topo: None,
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
            topology: None,
            stored_state: None,
            topo_built: false,
            exhausted: false,
            serial: 0,
            write_topo: None,
        })
    }

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        let state = self
            .parse_next_model()?
            .ok_or(PdbHandlerError::Empty)?;
        let top = self.topology.clone().unwrap();
        Ok((top, state))
    }

    fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        if !self.topo_built {
            let state = self
                .parse_next_model()?
                .ok_or(PdbHandlerError::Empty)?;
            self.stored_state.get_or_insert(state);
        }
        Ok(self.topology.clone().unwrap())
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        if let Some(s) = self.stored_state.take() {
            return Ok(s);
        }
        if self.exhausted {
            return Err(FileFormatError::Eof);
        }
        self.parse_next_model()?.ok_or(FileFormatError::Eof)
    }

    fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileFormatError> {
        let w = self.writer.as_mut().ok_or(FileFormatError::NotWritable)?;
        write_cryst1(w, data.get_box())?;
        let mut serial = self.serial;
        for (at, pos) in data.iter_atoms_dyn().zip(data.iter_pos_dyn()) {
            serial += 1;
            write_atom_record(w, serial % 99999, at, pos)?;
        }
        self.serial = serial;
        writeln!(w, "END")?;
        Ok(())
    }

    fn write_topology(&mut self, data: &dyn SaveTopology) -> Result<(), FileFormatError> {
        self.write_topo = Some(data.iter_atoms_dyn().cloned().collect());
        Ok(())
    }

    fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileFormatError> {
        let w = self.writer.as_mut().ok_or(FileFormatError::NotWritable)?;
        write_cryst1(w, data.get_box())?;
        let mut serial = self.serial;
        if let Some(ref atoms) = self.write_topo {
            for (at, pos) in atoms.iter().zip(data.iter_pos_dyn()) {
                serial += 1;
                write_atom_record(w, serial % 99999, at, pos)?;
            }
        } else {
            // No topology buffered — write minimal records with empty fields
            let empty = AtomStr::from_bytes(b"").unwrap();
            let dummy = Atom {
                name: empty,
                resname: empty,
                type_name: empty,
                ..Default::default()
            };
            for pos in data.iter_pos_dyn() {
                serial += 1;
                write_atom_record(w, serial % 99999, &dummy, pos)?;
            }
        }
        self.serial = serial;
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
