use crate::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError}, path::Path,
};
use thiserror::Error;

pub struct GroFileHandler {
    reader: Option<BufReader<File>>,
    writer: Option<BufWriter<File>>,

    // For separating IO of state and topology
    stored_topology: Option<Topology>,
    stored_state: Option<State>,

    at_least_one_state_read: bool,
}

#[derive(Debug, Error)]
pub enum GroHandlerError {
    #[error("unexpected io error")]
    Io(#[from] std::io::Error),

    #[error("can't open gro file for reading")]
    OpenRead(#[source] std::io::Error),

    #[error("can't open gro file for writing")]
    OpenWrite(#[source] std::io::Error),

    #[error(transparent)]
    ParseInt(#[from] ParseIntError),

    #[error(transparent)]
    ParseFloat(#[from] ParseFloatError),

    #[error(transparent)]
    Pbc(#[from] PeriodicBoxError),

    #[error("atom {0} has incomplete {1} entry")]
    AtomEntry(usize, String),

    #[error("gro file is empty")]
    EmptyFile,

    #[error("end of file reached")]
    Eof,
}

impl GroFileHandler {
    pub fn open(fname: impl AsRef<Path>) -> Result<Self, GroHandlerError> {
        Ok(Self {
            reader: BufReader::new(File::open(fname)
                .map_err(|e| GroHandlerError::OpenRead(e))?).into(),
            writer: None,
            stored_state: None,
            stored_topology: None,
            at_least_one_state_read: false,
        })
    }

    pub fn create(fname: impl AsRef<Path>) -> Result<Self, GroHandlerError> {
        Ok(Self {
            writer: BufWriter::new(File::create(fname)
                .map_err(|e| GroHandlerError::OpenWrite(e))?).into(),
            reader: None,
            stored_state: None,
            stored_topology: None,
            at_least_one_state_read: false,
        })
    }

    pub fn write(
        &mut self,
        data: &(impl TopologyIoProvider + StateIoProvider),
    ) -> Result<(), GroHandlerError> {
        let natoms = data.len();
        let buf = self.writer.as_mut().unwrap();

        // Print title
        writeln!(buf, "Created by Molar, t= {:.3}", data.get_time())?;
        // Write number of atoms
        writeln!(buf, "{natoms}")?;
        // Write atom lines
        for (i, (at, pos)) in std::iter::zip(data.iter_atoms(), data.iter_pos()).enumerate() {
            let ind = (i % 99999) + 1; // Prevents overflow of index field. It's not used anyway.
            let resid = at.resid % 99999; // Prevents overflow of resid field.

            writeln!(
                buf,
                "{:>5.5}{:<5.5}{:>5.5}{:>5.5}{:>8.3}{:>8.3}{:>8.3}",
                resid, at.resname, at.name, ind, pos.x, pos.y, pos.z
            )?;
        }

        // Write periodic box
        if let Some(b) = data.get_box() {
            let m = b.get_matrix();
            // Diagonal elements
            // Use same format as Gromacs for consistency, but this is free format
            write!(
                buf,
                "{:>10.4} {:>10.4} {:>10.4}",
                m[(0, 0)],
                m[(1, 1)],
                m[(2, 2)]
            )?;

            // Write off-diagonal only for triclinic boxes
            if b.is_triclinic() {
                // note leading space added after diagonal
                write!(
                    buf,
                    " {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
                    m[(1, 0)],
                    m[(2, 0)],
                    m[(0, 1)],
                    m[(2, 1)],
                    m[(0, 2)],
                    m[(1, 2)]
                )?;
            }
            // Add training newline
            writeln!(buf)?;
        } else {
            // No box, write zero diagonal
            writeln!(buf, "0.0 0.0 0.0")?;
        }

        Ok(())
    }

    pub fn read(&mut self) -> Result<Option<(Topology, State)>, GroHandlerError> {
        let mut top = TopologyStorage::default();
        let mut state = StateStorage::default();

        let buf = self.reader.as_mut().unwrap();
        let mut line = String::new();
        
        // Check if we are at EOF
        if buf.fill_buf()?.is_empty() {
            // EOF reached. If we alread read some frames than this is end of trajectory
            // Otherwise this is an empty file.
            if self.at_least_one_state_read {
                return Ok(None);
            } else {
                return Err(GroHandlerError::EmptyFile);
            }
        }

        // Read the title 
        buf.read_line(&mut line).unwrap();

        // Try to extract time from trailing "t=..."
        state.time = if let Some(i) = line.rfind("t=") {
            line[i+2..].trim().parse::<f32>().unwrap_or(0.0)
        } else {
            0.0
        };
        
        // Read number of atoms
        line.clear();
        buf.read_line(&mut line).unwrap();
        let natoms = line.trim().parse::<usize>()?;

        // Go over atoms line by line
        for i in 0..natoms {
            line.clear();
            buf.read_line(&mut line).unwrap();
        
            let mut at = Atom {
                resid: line
                    .get(0..5)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "resid".into()))?
                    .trim()
                    .parse::<i32>()?,
                resname: line
                    .get(5..10)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "resname".into()))?
                    .trim()
                    .to_owned(),
                name: line
                    .get(10..15)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "name".into()))?
                    .trim()
                    .to_owned(),
                chain: ' ',
                type_name: "".into(),
                ..Default::default()
            };

            // We don't have element number and mass, so guess them
            at.guess_element_and_mass_from_name();
        
            // Add atom to topology
            top.atoms.push(at);

            // Read coordinates
            let v = Pos::new(
                line.get(20..28)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "x".into()))?
                    .trim()
                    .parse::<f32>()?,
                line.get(28..36)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "y".into()))?
                    .trim()
                    .parse::<f32>()?,
                line.get(36..44)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "z".into()))?
                    .trim()
                    .parse::<f32>()?,
            );
            state.coords.push(v);
        }
        /* Read the box
        Format: (https://manual.gromacs.org/archive/5.0.3/online/gro.html)
        v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)

        Our arrangement:
        v1(x) v2(x) v3(x)
        v1(y) v2(y) v3(y)
        v1(z) v2(z) v3(z)

        So, the sequence of reads is:
        (0,0) (1,1) (2,2) (1,0) (2,0) (0,1) (2,1) (0,2) (1,2)
        */
        line.clear();
        buf.read_line(&mut line).unwrap();
        let l = line.split_whitespace().map(|s| {
                Ok(s.parse::<f32>()?)
            })
            .collect::<Result<Vec<f32>, GroHandlerError>>()?;

        let mut m = Matrix3f::zeros();
        m[(0, 0)] = l[0];
        m[(1, 1)] = l[1];
        m[(2, 2)] = l[2];
        if l.len() == 9 {
            m[(1, 0)] = l[3];
            m[(2, 0)] = l[4];
            m[(0, 1)] = l[5];
            m[(2, 1)] = l[6];
            m[(0, 2)] = l[7];
            m[(1, 2)] = l[8];
        }
        state.pbox = Some(PeriodicBox::from_matrix(m).map_err(|e| GroHandlerError::Pbc(e))?);

        let state: State = state.into();
        let top: Topology = top.into();
        // Assign resindex
        top.assign_resindex();

        self.at_least_one_state_read = true;

        Ok(Some((top, state)))
    }

    pub fn read_topology(&mut self) -> Result<Topology, GroHandlerError> {
        if self.stored_topology.is_some() {
            Ok(self.stored_topology.take().unwrap())
        } else {
            let (top,st) = self.read()?.ok_or(GroHandlerError::Eof)?;
            self.stored_state.get_or_insert(st);
            Ok(top)
        }
    }

    pub fn read_state(&mut self) -> Result<Option<State>, GroHandlerError> {
        if self.stored_state.is_some() {
            Ok(Some(self.stored_state.take().unwrap()))
        } else {
            if let Some((top,st)) = self.read()? {
                self.stored_topology.get_or_insert(top);
                Ok(Some(st))
            } else {
                Ok(None)
            }
            
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::io::FileHandler;

    #[test]
    #[should_panic]
    fn invalid_file() {
        let mut fh = FileHandler::open("nonexisting.gro").unwrap();
        fh.read().unwrap();
    }
}
