use super::{
    IoIndexProvider, IoTopologyProvider,
};
use crate::core::{Atom, Matrix3f, PeriodicBox, Pos, State, Topology};
use anyhow::{Result, anyhow};
use ascii::{AsciiChar, AsciiString};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
};

pub struct GroFileHandler {
    file_name: String,
    top: Option<Topology>,
    state: Option<State>,
    io_done: bool,
}

impl GroFileHandler {

    pub fn open(fname: &str) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self {
            file_name: fname.to_owned(),
            top: None,
            state: None,
            io_done: false,
        })
    }

    pub fn create(fname: &str) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self {
            file_name: fname.to_owned(),
            top: None,
            state: None,
            io_done: false,
        })
    }

    pub fn read(&mut self) -> Result<(Topology,State)> {
        self.io_done = true;
        self.do_read()
    }

    pub fn write(&mut self, data: &(impl IoIndexProvider + IoTopologyProvider + super::IoStateProvider)) -> Result<()> {
        self.io_done = true;
        self.do_write(data)
    }

    pub fn read_topology(&mut self) -> Result<Topology> {
        self.try_read()?;
        self.top.take().ok_or_else(|| anyhow!("Topology already read!"))
    }

    pub fn read_state(&mut self) -> Result<Option<State>> {
        self.try_read()?;
        Ok(self.state.take())
    }

    // Reads the file once
    fn try_read(&mut self) -> Result<()> {
        if !self.io_done {
            let (top,st) = self.do_read()?;
            self.top = Some(top);
            self.state = Some(st);
            self.io_done = true;
        }
        Ok(())
    }

    /*
    // Writes once
    fn try_write(
        &mut self,
        data: &(impl IoIndexProvider + IoTopologyProvider + super::IoStateProvider),
    ) -> Result<()> 
    {
        if !self.io_done {
            self.do_write(data)
            self.top = Some(top);
            self.state = Some(st);
            self.io_done = true;
        }
        Ok(())
    }
    */

    fn do_read(&mut self) -> Result<(Topology, State)> {
        let mut top = Topology::new();
        let mut state = State::new();

        let mut lines = BufReader::new(File::open(self.file_name.to_owned())?).lines();
        // Skip the title
        let _ = lines.next().unwrap();
        // Read number of atoms
        let natoms = lines.next().unwrap()?.parse::<usize>()?;
        // Go over atoms line by line
        for _ in 0..natoms {
            let line = lines.next().unwrap()?;

            let at = Atom {
                resid: line[0..5].parse()?,
                resname: AsciiString::from_ascii(line[5..10].trim().to_owned())?,
                name: AsciiString::from_ascii(line[10..15].trim().to_owned())?,
                chain: AsciiChar::Space,
                type_name: AsciiString::from_ascii("")?,
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
                line[20..28].parse()?,
                line[28..36].parse()?,
                line[36..44].parse()?,
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
            .unwrap()?
            .split(" ")
            .map(|s| s.parse().unwrap())
            .collect::<Vec<f32>>();
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
        state.box_ = Some(PeriodicBox::from_matrix(m)?);

        // Assign resindex
        top.assign_resindex();

        Ok((top, state))
    }

    fn do_write(
        &mut self,
        data: &(impl IoIndexProvider + IoTopologyProvider + super::IoStateProvider),
    ) -> Result<()> {
        // Open file for writing
        let mut buf = BufWriter::new(File::create(self.file_name.to_owned())?);
        let index = data.get_index();
        let natoms = index.len();
        let top = data.get_topology();
        let state = data.get_state();
        // Print title
        writeln!(buf, "Created by Molar")?;
        // Write number of atoms
        writeln!(buf, "{natoms}")?;
        // Write atom lines
        for i in 0..natoms {
            let ind = (i % 99999) + 1; // Prevents overflow of index field. It's not used anyway.
            let resid = top.atoms[i].resid % 99999; // Prevents overflow of resid field.

            writeln!(
                buf,
                "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                resid,
                top.atoms[i].resname,
                top.atoms[i].name,
                ind,
                state.coords[i].x,
                state.coords[i].y,
                state.coords[i].z
            )?;
        }

        // Write periodic box
        if let Some(b) = state.box_.as_ref() {
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
}