use std::{
    fs::File,
    io::{BufRead, BufReader},
    num::{ParseFloatError, ParseIntError},
};
use thiserror::Error;
use crate::core::Atom;

use super::{Topology, TopologyStorage};
use regex;

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
        
        // Loop over lines until we encounted moleculetype
        let moltype_re = regex::Regex::new(r"\[\s+moleculetype\s+\]").unwrap();
        loop {
            line.clear();
            let nread = reader.read_line(&mut line).unwrap();
            if nread == 0 {
                return Err(ItpHandlerError::NoMoleculetype);
            }
            if moltype_re.is_match(&line) {
                break;
            }
        }

        // We are now in moleculetype
        // Skip until we get atoms
        let atoms_re = regex::Regex::new(r"\[\s+atoms\s+\]").unwrap();
        loop {
            line.clear();
            let nread = reader.read_line(&mut line).unwrap();
            if nread == 0 {
                return Err(ItpHandlerError::NoAtoms);
            }
            if atoms_re.is_match(&line) {
                break;
            }
        }

        // We are in atoms. Read them
        loop {
            line.clear();
            let nread = reader.read_line(&mut line).unwrap();
            let l = line.trim();
            if nread == 0 {
                break;
            }
            // Skip comments and empty lines
            if l.trim_start().starts_with(";") || l.len() == 0 {
                continue
            }
            
            let fields = l.split_ascii_whitespace().collect::<Vec<_>>();
            if fields.len() < 8 {
                break
            }
            let mut at = Atom {
                name: fields[4].to_owned(),
                resname: fields[3].to_owned(),
                type_name: fields[1].to_owned(),
                resid: fields[2].parse()?,
                charge: fields[6].parse()?,
                mass: fields[7].parse()?,
                ..Default::default()
            };
            // We don't have element number, guess it
            at.guess_element_from_name();
            // Add atom to topology
            top.atoms.push(at);
        }

        Ok(top.into())
    }
}

#[derive(Debug, Error)]
pub enum ItpHandlerError {
    //#[error("unxpected io error")]
    //Io(#[from] std::io::Error),

    #[error("can't open itp file for reading")]
    OpenRead(#[source] std::io::Error),
    
    #[error("no moleculetype found")]
    NoMoleculetype,

    #[error("invalid atom entry")]
    InvalidAtomEntry(String),
    
    #[error("no atoms found")]
    NoAtoms,

    #[error("no resname found")]
    NoResname,

    #[error(transparent)]
    Regex(#[from] regex::Error),

    #[error(transparent)]
    ParseInt(#[from] ParseIntError),

    #[error(transparent)]
    ParseFloat(#[from] ParseFloatError),
}
