use crate::core::{Atom, Pos, State, StateStorage, Topology, TopologyStorage, Vector3f};
use anyhow::Result;
use std::{
    fs::File,
    io::{BufRead, BufReader, Lines},
};

pub struct PdbFileHandler {
    lines: Lines<BufReader<File>>,
    file_name: String,
    natoms: usize,
    atoms_buf: Vec<Atom>,
    coord_buf: Vec<Pos>,
}

impl PdbFileHandler {
    fn open(fname: &str) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self {
            lines: BufReader::new(File::open(fname)?).lines(),
            file_name: fname.to_owned(),
            natoms: 0,
            atoms_buf: vec![],
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
        // Clear buffers
        self.atoms_buf.clear();
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
                        name: line[12..=15].into(),
                        resname: line[17..=19].into(),
                        chain: line.as_bytes()[21].into(),
                        resid: line[22..=25].parse()?,
                        occupancy: line[54..=59].parse()?,
                        bfactor: line[60..=65].parse()?,
                        atomic_number: 0, //line.as_bytes()[77],
                        ..Default::default()
                    };
                    self.atoms_buf.push(at);

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
        let top = TopologyStorage {
            atoms: self.atoms_buf,
            bonds: Default::default(),
            molecules: Default::default(),
        };

        let state = StateStorage {
            coords: self.coord_buf,
            pbox: None,
            time: 0.0,            
        };
        Ok((top.into(),state.into()))
    }
}
