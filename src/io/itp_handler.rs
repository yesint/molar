use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    num::{ParseFloatError, ParseIntError},
};
use thiserror::Error;

use super::{Topology, TopologyStorage};

pub struct ItpFileHandler {
    file: File,
}

impl ItpFileHandler {
    pub fn open(fname: &str) -> Result<Self, ItpHandlerError> {
        Ok(Self {
            file: File::open(fname.to_owned())
                .map_err(|e| ItpHandlerError::OpenRead(e))?,
        })
    }

    pub fn read_topology(&mut self) -> Result<Topology, ItpHandlerError> {
        let mut top = TopologyStorage::default();
        
        let mut reader = BufReader::new(&self.file);
        let mut line = String::new();
        
        //reader.read_line(&mut line).with_context(|| "a")?;
        // Go over atoms line by line
        // for i in 0..natoms {
        //     line.clear();
        //     reader.read_line(&mut line)?;
        // }
        todo!()
    }
}

#[derive(Debug, Error)]
pub enum ItpHandlerError {
    //#[error("unxpected io error")]
    //Io(#[from] std::io::Error),

    #[error("can't open itp file for reading")]
    OpenRead(#[source] std::io::Error),

    
}
