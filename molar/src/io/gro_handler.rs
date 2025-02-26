use crate::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
};
use thiserror::Error;

pub struct GroFileHandler {
    file: File,
}

#[derive(Debug, Error)]
pub enum GroHandlerError {
    #[error("unxpected io error")]
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
}

impl GroFileHandler {
    pub fn open(fname: &str) -> Result<Self, GroHandlerError> {
        Ok(Self {
            file: File::open(fname.to_owned())
                .map_err(|e| GroHandlerError::OpenRead(e))?,
        })
    }

    pub fn create(fname: &str) -> Result<Self, GroHandlerError> {
        Ok(Self {
            file: File::create(fname.to_owned())
                .map_err(|e| GroHandlerError::OpenWrite(e))?,
        })
    }

    pub fn write(
        &mut self,
        data: &(impl TopologyProvider + StateProvider),
    ) -> Result<(), GroHandlerError> {
        // Open file for writing
        let mut buf = BufWriter::new(&self.file);
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

    // pub fn get_file_name(&self) -> &str {
    //     &self.file_name
    // }

    pub fn read(&mut self) -> Result<(Topology, State), GroHandlerError> {
        let mut top = TopologyStorage::default();
        let mut state = StateStorage::default();

        let mut reader = BufReader::new(&self.file);
        let mut line = String::new();
        
        // Skip the title
        reader.read_line(&mut line).unwrap();
        
        // Read number of atoms
        line.clear();
        reader.read_line(&mut line).unwrap();
        let natoms = line.trim().parse::<usize>()?;

        // Go over atoms line by line
        for i in 0..natoms {
            line.clear();
            reader.read_line(&mut line).unwrap();
        
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
        reader.read_line(&mut line).unwrap();
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
