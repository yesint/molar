use crate::atom::element_symbol;
use crate::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
    path::Path,
};
use thiserror::Error;

pub struct XyzFileHandler {
    reader: Option<BufReader<File>>,
    writer: Option<BufWriter<File>>,

    // For separating IO of state and topology
    stored_topology: Option<Topology>,
    stored_state: Option<State>,

    at_least_one_state_read: bool,
}

#[derive(Debug, Error)]
pub enum XyzHandlerError {
    #[error("can't open xyz file for reading")]
    OpenRead(#[source] std::io::Error),

    #[error("can't open xyz file for writing")]
    OpenWrite(#[source] std::io::Error),

    #[error("xyz file is empty or contains no atoms")]
    Empty,

    #[error("unexpected end of file reached")]
    Eof,

    #[error("malformed atom count line")]
    BadAtomCount(#[from] ParseIntError),

    #[error("malformed coordinate on atom line {0}")]
    BadCoord(usize, #[source] ParseFloatError),

    #[error("truncated atom line {0}")]
    TruncatedLine(usize),

    #[error("io error")]
    Io(#[from] std::io::Error),
}

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
                File::create(fname).map_err(XyzHandlerError::OpenWrite)?,
            )),
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
                Err(XyzHandlerError::Eof)?
            } else {
                Err(XyzHandlerError::Empty)?
            };
        }

        // Read atom count line; empty line also signals end of trajectory
        let mut line = String::new();
        let n = buf.read_line(&mut line)?;
        if n == 0 || line.trim().is_empty() {
            return if self.at_least_one_state_read {
                Err(XyzHandlerError::Eof)?
            } else {
                Err(XyzHandlerError::Empty)?
            };
        }
        let natoms: usize = line.trim().parse().map_err(XyzHandlerError::BadAtomCount)?;

        // Skip comment/title line
        line.clear();
        buf.read_line(&mut line)?;

        let mut atoms: Vec<Atom> = Vec::with_capacity(natoms);
        let mut coords: Vec<Pos> = Vec::with_capacity(natoms);

        for i in 0..natoms {
            line.clear();
            buf.read_line(&mut line)?;
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

            atoms.push(
                Atom::new()
                    .with_name(elem)
                    .with_resname("MOL")
                    .with_resid(1)
                    .with_chain('A')
                    .guess()
            );
        }

        let mut top = Topology::default();
        top.atoms = atoms;
        top.assign_resindex();

        self.at_least_one_state_read = true;

        Ok((top, State { coords, time: 0.0, pbox: None }))
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
}
