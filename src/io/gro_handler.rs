use super::{IoReader, IoTopologyReader, IoStateReader, IoTopologyWriter, IoWriter};
use crate::core::{Topology, Atom, Pos, State, Matrix3f, PeriodicBox};
use anyhow::{Result, anyhow};
use std::{
    fs::File,
    io::{BufRead, BufReader, Lines},
};
use ascii::{AsciiString, AsciiChar};

pub struct GroFileHandler {
    lines: Lines<BufReader<File>>,
    top: Option<Topology>,
    state: Option<State>,
    is_read: bool,
}

impl IoReader for GroFileHandler {
    fn new_reader(fname: &str) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self {
            lines: BufReader::new(File::open(fname)?).lines(),
            top: None,
            state: None,
            is_read: false,
        })
    }
}

impl IoWriter for GroFileHandler {
    fn new_writer(fname: &str) -> Result<Self>
        where
            Self: Sized 
    {
        todo!();
    }
}

impl IoTopologyReader for GroFileHandler {
    fn read_topology(&mut self) -> Result<Topology> {
        if !self.is_read {
            let (t,st) = self.read()?;
            self.top = Some(t);
            self.state = Some(st);
        }
        self.top.take().ok_or(anyhow!("Already read topology"))
    }
}

impl IoStateReader for GroFileHandler {
    fn read_state(&mut self) -> Result<Option<State>> {
        if !self.is_read {
            let (t,st) = self.read()?;
            self.top = Some(t);
            self.state = Some(st);
        }
        Ok(self.state.take())
    }
}

impl GroFileHandler {
    fn read(&mut self) -> Result<(Topology,State)> {
        let mut top = Topology::new();
        let mut state = State::new();

        // Skip the title
        let _ = self.lines.next().unwrap();
        // Read number of atoms
        let natoms = self.lines.next().unwrap()?.parse::<usize>()?;
        // Go over atoms line by line
        for i in 0..natoms {
            let line = self.lines.next().unwrap()?;

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
                line[36..44].parse()?
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
        let l = self.lines.next().unwrap()?.split(" ").map(|s| s.parse().unwrap()).collect::<Vec<f32>>();
        let mut m = Matrix3f::zeros();
        m[(0,0)] = l[0];
        m[(1,1)] = l[1];
        m[(2,2)] = l[2];
        if l.len()==9 {
            m[(1,0)] = l[3];
            m[(2,0)] = l[4];
            m[(0,1)] = l[5];
            m[(2,1)] = l[6];
            m[(0,2)] = l[7];
            m[(1,2)] = l[8];
        }
        state.box_ = Some(PeriodicBox::from_matrix(m)?);

        // Assign resindex
        top.assign_resindex();
        Ok((top,state))
    }
}

impl IoTopologyWriter for GroFileHandler {
    fn write_topology(&mut self, data: &impl super::IoIndexAndTopologyProvider) -> Result<()> {
        todo!()
    }
}