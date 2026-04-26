use super::FileFormatError;
use crate::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
    path::Path,
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

    #[error("atom {0} has corrupted {1} entry")]
    AtomEntry(usize, String),

    #[error("gro file is empty")]
    EmptyFile,

    #[error("unexpcted end of file reached")]
    Eof,
}

impl GroFileHandler {
    fn read_inner(
        &mut self,
        _read_coords: bool,
        read_vels: bool,
        _read_forces: bool,
    ) -> Result<(Topology, State), FileFormatError> {
        let mut top = Topology::default();
        let mut state = State::default();

        let buf = self.reader.as_mut().unwrap();
        let mut line = String::new();

        // Check if we are at EOF
        if buf.fill_buf()?.is_empty() {
            if self.at_least_one_state_read {
                return Err(GroHandlerError::Eof)?;
            } else {
                return Err(GroHandlerError::EmptyFile)?;
            }
        }

        // Read the title
        buf.read_line(&mut line).unwrap();

        // Try to extract time from trailing "t=..."
        state.time = if let Some(i) = line.rfind("t=") {
            line[i + 2..].trim().parse::<Float>().unwrap_or(0.0)
        } else {
            0.0
        };

        // Read number of atoms
        line.clear();
        buf.read_line(&mut line).unwrap();
        let natoms = line
            .trim()
            .parse::<usize>()
            .map_err(GroHandlerError::ParseInt)?;

        state.coords.reserve(natoms);

        // Peek at the first atom line to detect whether the file has velocity columns.
        // This single check before the loop avoids per-atom branching.
        let file_has_vels = if read_vels && natoms > 0 {
            buf.fill_buf()?;
            // fill_buf gives us the buffer without consuming it; check the line length.
            // A GRO line with velocities has at least 68 bytes (cols 0-67).
            let peek = buf.fill_buf()?;
            let line_len = peek.iter().position(|&b| b == b'\n').map(|p| p + 1).unwrap_or(peek.len());
            line_len >= 68
        } else {
            false
        };

        if file_has_vels {
            let mut vels: Vec<Vel> = Vec::with_capacity(natoms);
            for i in 0..natoms {
                line.clear();
                buf.read_line(&mut line).unwrap();

                let resid = line.get(0..5).ok_or_else(|| GroHandlerError::AtomEntry(i, "resid".into()))?.trim().parse::<i32>().map_err(GroHandlerError::ParseInt)?;
                let resname = line.get(5..10).ok_or_else(|| GroHandlerError::AtomEntry(i, "resname".into()))?.trim();
                let name = line.get(10..15).ok_or_else(|| GroHandlerError::AtomEntry(i, "name".into()))?.trim();
                top.atoms.push(Atom::new().with_name(name).with_resname(resname).with_resid(resid).guess());

                state.coords.push(Pos::new(
                    line.get(20..28).ok_or_else(|| GroHandlerError::AtomEntry(i, "x".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                    line.get(28..36).ok_or_else(|| GroHandlerError::AtomEntry(i, "y".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                    line.get(36..44).ok_or_else(|| GroHandlerError::AtomEntry(i, "z".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                ));

                vels.push(Vel::new(
                    line.get(44..52).ok_or_else(|| GroHandlerError::AtomEntry(i, "vx".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                    line.get(52..60).ok_or_else(|| GroHandlerError::AtomEntry(i, "vy".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                    line.get(60..68).ok_or_else(|| GroHandlerError::AtomEntry(i, "vz".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                ));
            }
            state.velocities = vels;
        } else {
            for i in 0..natoms {
                line.clear();
                buf.read_line(&mut line).unwrap();

                let resid = line.get(0..5).ok_or_else(|| GroHandlerError::AtomEntry(i, "resid".into()))?.trim().parse::<i32>().map_err(GroHandlerError::ParseInt)?;
                let resname = line.get(5..10).ok_or_else(|| GroHandlerError::AtomEntry(i, "resname".into()))?.trim();
                let name = line.get(10..15).ok_or_else(|| GroHandlerError::AtomEntry(i, "name".into()))?.trim();
                top.atoms.push(Atom::new().with_name(name).with_resname(resname).with_resid(resid).guess());

                state.coords.push(Pos::new(
                    line.get(20..28).ok_or_else(|| GroHandlerError::AtomEntry(i, "x".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                    line.get(28..36).ok_or_else(|| GroHandlerError::AtomEntry(i, "y".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                    line.get(36..44).ok_or_else(|| GroHandlerError::AtomEntry(i, "z".into()))?.trim().parse::<Float>().map_err(GroHandlerError::ParseFloat)?,
                ));
            }
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
        let l = line
            .split_whitespace()
            .map(|s| Ok(s.parse::<Float>()?))
            .collect::<Result<Vec<Float>, GroHandlerError>>()?;

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
        state.pbox = Some(PeriodicBox::from_matrix(m).map_err(GroHandlerError::Pbc)?);

        // Assign resindex
        top.assign_resindex();

        self.at_least_one_state_read = true;

        Ok((top, state))
    }
}

impl FileFormatHandler for GroFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Ok(Self {
            reader: BufReader::new(File::open(fname).map_err(|e| GroHandlerError::OpenRead(e))?)
                .into(),
            writer: None,
            stored_state: None,
            stored_topology: None,
            at_least_one_state_read: false,
        })
    }

    fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Ok(Self {
            writer: BufWriter::new(File::create(fname).map_err(|e| GroHandlerError::OpenWrite(e))?)
                .into(),
            reader: None,
            stored_state: None,
            stored_topology: None,
            at_least_one_state_read: false,
        })
    }

    fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileFormatError> {
        let natoms = data.len();
        let buf = self.writer.as_mut().unwrap();

        // Print title
        writeln!(buf, "Created by Molar, t= {:.3}", data.get_time())?;
        // Write number of atoms
        writeln!(buf, "{natoms}")?;

        // Write atom lines — two separate loops to avoid per-atom velocity branch.
        let at_it = data.iter_atoms_dyn();
        let pos_it = data.iter_pos_dyn();
        let vel_it = data.iter_vel_dyn();

        if vel_it.len() > 0 {
            for (i, ((at, pos), v)) in at_it.zip(pos_it).zip(vel_it).enumerate() {
                let ind = (i % 99999) + 1;
                let resid = at.resid % 99999;
                write!(buf, "{:>5.5}{:<5.5}{:>5.5}{:>5.5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}",
                    resid, at.resname, at.name, ind, pos.x, pos.y, pos.z, v.x, v.y, v.z)?;
                writeln!(buf)?;
            }
        } else {
            for (i, (at, pos)) in at_it.zip(pos_it).enumerate() {
                let ind = (i % 99999) + 1;
                let resid = at.resid % 99999;
                write!(buf, "{:>5.5}{:<5.5}{:>5.5}{:>5.5}{:>8.3}{:>8.3}{:>8.3}",
                    resid, at.resname, at.name, ind, pos.x, pos.y, pos.z)?;
                writeln!(buf)?;
            }
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

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        self.read_inner(true, true, true)
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

    fn read_state_pick(
        &mut self,
        read_coords: bool,
        read_vels: bool,
        read_forces: bool,
    ) -> Result<State, FileFormatError> {
        if read_forces {
            return Err(FileFormatError::NoForces);
        }
        if self.stored_state.is_some() {
            // Stored state was produced by read_inner(all); trim unwanted fields.
            let mut st = self.stored_state.take().unwrap();
            if read_vels && st.velocities.is_empty() {
                return Err(FileFormatError::NoVelocities);
            }
            if !read_vels {
                st.velocities.clear();
            }
            if !read_coords {
                st.coords.clear();
            }
            Ok(st)
        } else {
            let (top, st) = self.read_inner(true, read_vels, false)?;
            self.stored_topology.get_or_insert(top);
            Ok(st)
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
