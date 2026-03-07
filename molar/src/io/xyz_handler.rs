use crate::atom::{element_symbol, AtomStr, ATOM_NAME_EXPECT};
use crate::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
    path::Path,
};
use thiserror::Error;

pub struct XyzFileHandler {
    // Reading
    reader: Option<BufReader<File>>,
    stored_topology: Option<Topology>, // cached topology; consumed by read_topology()
    stored_state: Option<State>,       // cached state; consumed by read_state()
    exhausted: bool,
    // Writing
    writer: Option<BufWriter<File>>,
    write_topo: Option<Vec<Atom>>, // buffered by write_topology() for write_state()
}

#[derive(Debug, Error)]
pub enum XyzHandlerError {
    #[error("can't open xyz file for reading")]
    OpenRead(#[source] std::io::Error),

    #[error("can't open xyz file for writing")]
    OpenWrite(#[source] std::io::Error),

    #[error("xyz file is empty or contains no atoms")]
    Empty,

    #[error("malformed atom count line")]
    BadAtomCount(#[from] ParseIntError),

    #[error("malformed coordinate on atom line {0}")]
    BadCoord(usize, #[source] ParseFloatError),

    #[error("truncated atom line {0}")]
    TruncatedLine(usize),

    #[error("io error")]
    Io(#[from] std::io::Error),
}

// ── core handler ────────────────────────────────────────────────────────────

impl XyzFileHandler {
    /// Stream-parse one XYZ frame (natoms / comment / atom lines).
    /// Topology is built and stored into `self.stored_topology` only when it is `None`.
    /// Returns `Ok(None)` at EOF.
    fn parse_next_frame(&mut self) -> Result<Option<State>, XyzHandlerError> {
        let reader = match self.reader.as_mut() {
            Some(r) => r,
            None => return Ok(None),
        };

        // Read atom count line; empty → EOF
        let mut line = String::new();
        let n = reader.read_line(&mut line)?;
        if n == 0 || line.trim().is_empty() {
            self.exhausted = true;
            return Ok(None);
        }
        let natoms: usize = line.trim().parse().map_err(XyzHandlerError::BadAtomCount)?;

        // Skip comment/title line
        line.clear();
        reader.read_line(&mut line)?;

        let build_topo = self.stored_topology.is_none();
        let mut atoms: Vec<Atom> = Vec::with_capacity(if build_topo { natoms } else { 0 });
        let mut coords: Vec<Pos> = Vec::with_capacity(natoms);

        for i in 0..natoms {
            line.clear();
            reader.read_line(&mut line)?;
            let mut tokens = line.split_whitespace();

            let elem = tokens.next().ok_or(XyzHandlerError::TruncatedLine(i))?;
            let x: f32 = tokens
                .next()
                .ok_or(XyzHandlerError::TruncatedLine(i))?
                .parse()
                .map_err(|e| XyzHandlerError::BadCoord(i, e))?;
            let y: f32 = tokens
                .next()
                .ok_or(XyzHandlerError::TruncatedLine(i))?
                .parse()
                .map_err(|e| XyzHandlerError::BadCoord(i, e))?;
            let z: f32 = tokens
                .next()
                .ok_or(XyzHandlerError::TruncatedLine(i))?
                .parse()
                .map_err(|e| XyzHandlerError::BadCoord(i, e))?;

            coords.push(Pos::new(x * 0.1, y * 0.1, z * 0.1));

            if build_topo {
                let mut at = Atom {
                    name: AtomStr::try_from_str(elem).expect(ATOM_NAME_EXPECT),
                    resname: AtomStr::try_from_str("MOL").unwrap(),
                    resid: 1,
                    chain: 'A',
                    occupancy: 1.0,
                    type_name: AtomStr::try_from_str("").unwrap(),
                    ..Default::default()
                };
                at.guess_element_and_mass_from_name();
                atoms.push(at);
            }
        }

        if build_topo {
            let mut top = Topology::default();
            top.atoms = atoms;
            top.assign_resindex();
            self.stored_topology = Some(top);
        }

        Ok(Some(State {
            coords,
            time: 0.0,
            pbox: None,
        }))
    }
}

// ── FileFormatHandler ────────────────────────────────────────────────────────

impl FileFormatHandler for XyzFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Ok(Self {
            reader: Some(BufReader::new(
                File::open(fname).map_err(XyzHandlerError::OpenRead)?,
            )),
            writer: None,
            stored_topology: None,
            stored_state: None,
            exhausted: false,
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
                File::create(fname).map_err(XyzHandlerError::OpenWrite)?,
            )),
            stored_topology: None,
            stored_state: None,
            exhausted: false,
            write_topo: None,
        })
    }

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        self.stored_topology = None; // force topology rebuild
        let st = self.parse_next_frame()?.ok_or(XyzHandlerError::Empty)?;
        let top = self.stored_topology.take().unwrap();
        Ok((top, st))
    }

    fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        if self.stored_topology.is_some() {
            return Ok(self.stored_topology.take().unwrap());
        }
        let st = self.parse_next_frame()?.ok_or(XyzHandlerError::Empty)?;
        self.stored_state.get_or_insert(st);
        Ok(self.stored_topology.take().unwrap())
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        if self.stored_state.is_some() {
            return Ok(self.stored_state.take().unwrap());
        }
        let st = self.parse_next_frame()?.ok_or(FileFormatError::Eof)?;
        Ok(st)
    }

    fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileFormatError> {
        let w = self.writer.as_mut().ok_or(FileFormatError::NotWritable)?;
        let n = data.len();
        writeln!(w, "{n}")?;
        writeln!(w)?; // blank comment line
        for (at, pos) in data.iter_atoms_dyn().zip(data.iter_pos_dyn()) {
            let elem = element_symbol(at.atomic_number);
            let sym = if elem.is_empty() { at.name.as_str() } else { elem };
            writeln!(
                w,
                "{} {:>12.6} {:>12.6} {:>12.6}",
                sym,
                pos.x * 10.0,
                pos.y * 10.0,
                pos.z * 10.0
            )?;
        }
        Ok(())
    }

    fn write_topology(&mut self, data: &dyn SaveTopology) -> Result<(), FileFormatError> {
        self.write_topo = Some(data.iter_atoms_dyn().cloned().collect());
        Ok(())
    }

    fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileFormatError> {
        let w = self.writer.as_mut().ok_or(FileFormatError::NotWritable)?;
        let n = data.len();
        writeln!(w, "{n}")?;
        writeln!(w)?;
        if let Some(ref atoms) = self.write_topo {
            for (at, pos) in atoms.iter().zip(data.iter_pos_dyn()) {
                let elem = element_symbol(at.atomic_number);
                let sym = if elem.is_empty() { at.name.as_str() } else { elem };
                writeln!(
                    w,
                    "{} {:>12.6} {:>12.6} {:>12.6}",
                    sym,
                    pos.x * 10.0,
                    pos.y * 10.0,
                    pos.z * 10.0
                )?;
            }
        } else {
            // No topology — write placeholder element
            for pos in data.iter_pos_dyn() {
                writeln!(
                    w,
                    "X {:>12.6} {:>12.6} {:>12.6}",
                    pos.x * 10.0,
                    pos.y * 10.0,
                    pos.z * 10.0
                )?;
            }
        }
        Ok(())
    }
}
