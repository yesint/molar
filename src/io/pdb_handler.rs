use super::{IoReader, IoTopologyReader};
use crate::core::{Topology, Atom, Pos, Vector3f, State};
use anyhow::Result;
use std::{
    fs::File,
    io::{BufRead, BufReader, Lines},
};
use ascii::{AsciiString, AsciiChar};

pub struct PdbFileHandler {
    lines: Lines<BufReader<File>>,
    file_name: String,
    natoms: usize,
    coord_buf: Vec<Pos>,
}

impl IoReader for PdbFileHandler {
    fn new_reader(fname: &str) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self {
            lines: BufReader::new(File::open(fname)?).lines(),
            file_name: fname.to_owned(),
            natoms: 0,
            coord_buf: vec![],
        })
    }
}

enum PdbRecord {
    Atom,
    Hetatm,
    Cryst,
    Conect,
    End,
    Ter,
    Other,
}

fn get_pdb_record(line: &str) -> PdbRecord {
    if line.starts_with("ATOM") {
        PdbRecord::Atom
    } else if line.starts_with("HETATM") {
        PdbRecord::Hetatm
    } else if line.starts_with("CRYST") {
        PdbRecord::Cryst
    } else if line.starts_with("END") {
        PdbRecord::End
    } else if line.starts_with("TER") {
        PdbRecord::Ter
    } else if line.starts_with("CONECT") {
        PdbRecord::Conect
    } else {
        PdbRecord::Other
    }
}

impl PdbFileHandler {
    fn read(&mut self) -> Result<(Topology,State)> {
        // Clear coordinate buffer
        self.coord_buf.clear();
        // Go line by line
        for line in self.lines {
            let line = line?;
            match get_pdb_record(&line) {
                PdbRecord::Atom | PdbRecord::Hetatm=> {
                    // Skip if altloc or insertion code is set
                    if line.as_bytes()[16] != b' ' || line.as_bytes()[26] != b' ' {
                        continue;
                    }

                    let at = Atom {
                        name: AsciiString::from_ascii(&line[12..=15])?,
                        resname: AsciiString::from_ascii(&line[17..=19])?,
                        chain: AsciiChar::from_ascii(line.as_bytes()[21])?,
                        resid: line[22..=25].parse()?,
                        occupancy: line[54..=59].parse()?,
                        bfactor: line[60..=65].parse()?,
                        atomic_number: 0, //line.as_bytes()[77],
                        ..Default::default()
                    };

                    // Save coordinates for later
                    let v = Pos::new(
                        line[30..=37].parse()?,
                        line[38..=45].parse()?,
                        line[46..=53].parse()?
                    );
                    self.coord_buf.push(v);
                    
                },
                _ => {}
            }
        }
        ()
    }
}
