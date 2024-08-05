use super::{FileFormatHandler, FileHandlerError, PeriodicBoxError, ReadTopAndState, State, StateProvider, Topology, TopologyProvider};
use crate::core::{Atom, Matrix3f, PeriodicBox, Pos, StateStorage, TopologyStorage};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
};
use thiserror::Error;

pub struct GroFileHandler {
    file_name: String,
}

#[derive(Debug, Error)]
pub enum GroHandlerError {
    #[error("unxpected io error")]
    Io(#[from] std::io::Error),

    #[error("unexpected end of GRO file at: {0}")]
    Eof(String),
    
    #[error("not an integer at {1}: {0}")]
    ParseInt(ParseIntError, String),
    
    #[error("not a float at {1}: {0}")]
    ParseFloat(ParseFloatError, String),
    
    #[error("atom entry {0} truncated at {1}")]
    AtomEntry(usize, String),

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),
}

impl FileFormatHandler for GroFileHandler {
    fn format_ext() -> String {
        "gro".into()
    }

    fn open(fname: &str) -> Result<Self, FileIoError> {
        Ok(Self {file_name: fname.to_owned()})
    }

    fn create(fname: &str) -> Result<Self, FileIoError> {
        Ok(Self {file_name: fname.to_owned()})
    }

    fn write(
        &mut self,
        data: &(impl TopologyProvider + StateProvider),
    ) -> Result<(), FileIoError> {
        // Open file for writing
        let mut buf = BufWriter::new(File::create(self.file_name.to_owned())?);
        let natoms = data.num_atoms();

        // Print title
        writeln!(buf, "Created by Molar")?;
        // Write number of atoms
        writeln!(buf, "{natoms}")?;
        // Write atom lines
        for (i, (at, pos)) in std::iter::zip(data.iter_atoms(), data.iter_pos()).enumerate() {
            let ind = (i % 99999) + 1; // Prevents overflow of index field. It's not used anyway.
            let resid = at.resid % 99999; // Prevents overflow of resid field.

            writeln!(
                buf,
                "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
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
    
    fn get_file_name(&self) -> &str {
        &self.file_name
    }

    fn read(&mut self) -> Result<(Topology, State), FileIoError> {
        let mut top = TopologyStorage::default();
        let mut state = StateStorage::default();

        let mut lines =
            BufReader::new(File::open(self.file_name.to_owned())
                .map_err(|e| GroHandlerError::Io(e))?)
                .lines();
        // Skip the title
        let _ = lines
            .next()
            .ok_or_else(|| GroHandlerError::Eof("title".into()))?;

        // Read number of atoms
        let natoms = lines
            .next()
            .ok_or_else(|| GroHandlerError::Eof("natoms".into()))?
            .map_err(|e| GroHandlerError::Io(e))?
            .parse::<usize>()
            .map_err(|e| GroHandlerError::ParseInt(e, "natoms".into()))?;

        // Go over atoms line by line
        for i in 0..natoms {
            let line = lines
                .next()
                .ok_or_else(|| GroHandlerError::Eof(format!("atom #{i}")))?
                .map_err(|e| GroHandlerError::Io(e))?;

            let at = Atom {
                resid: line
                    .get(0..5)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "resid".into()))?
                    .parse()
                    .map_err(|e| GroHandlerError::ParseInt(e, format!("atom #{i} resid")))?,
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
                type_id: 0,
                occupancy: 0.0,
                bfactor: 0.0,
                charge: 0.0,
                resindex: 0,
                ..Default::default()
            };
            // Add atom to topology
            top.atoms.push(at);

            // Read coordinates
            let v = Pos::new(
                line.get(20..28)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "x".into()))?
                    .parse()
                    .map_err(|e| GroHandlerError::ParseFloat(e, format!("atom #{i} y")))?,
                line.get(28..36)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "y".into()))?
                    .parse()
                    .map_err(|e| GroHandlerError::ParseFloat(e, format!("atom #{i} y")))?,
                line.get(36..44)
                    .ok_or_else(|| GroHandlerError::AtomEntry(i, "z".into()))?
                    .parse()
                    .map_err(|e| GroHandlerError::ParseFloat(e, format!("atom #{i} z")))?,
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
        let l = lines
            .next()
            .ok_or_else(|| GroHandlerError::Eof("pbc".into()))?
            .map_err(|e| GroHandlerError::Io(e))?
            .split(" ")
            .map(|s| {
                Ok(s.parse::<f32>()
                    .map_err(|e| GroHandlerError::ParseFloat(e, "pbc".into()))?)
            })
            .collect::<Result<Vec<f32>,GroHandlerError>>()?;

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
        state.pbox = Some(
            PeriodicBox::from_matrix(m)
            .map_err(|e| GroHandlerError::Pbc(e))?
        );

        let state: State = state.into();

        let mut top: Topology = top.into();
        // Assign resindex
        top.assign_resindex();

        Ok((top, state.into()))
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