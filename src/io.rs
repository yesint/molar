use crate::core::{IndexIterator, State, Topology};
use anyhow::{anyhow, bail, Result};
use std::path::Path;

mod vmd_molfile_handler;
mod xtc_handler;

#[cfg(feature = "gromacs")]
mod tpr_handler;

// Reexports
#[cfg(feature = "gromacs")]
pub use tpr_handler::TprFileHandler;
pub use vmd_molfile_handler::VmdMolFileHandler;
pub use xtc_handler::XtcFileHandler;

//===============================
// Traits for file opening
//===============================
pub trait IoReader {
    fn new_reader(fname: &str) -> Result<Self>
    where
        Self: Sized;
}

pub trait IoWriter {
    fn new_writer(fname: &str) -> Result<Self>
    where
        Self: Sized;
}

//===============================
// Traits for writing
//===============================
pub trait IndexAndTopologyProvider {
    fn get_index_and_topology(&self) -> (impl IndexIterator, &Topology);
}

pub trait IndexAndStateProvider {
    fn get_index_and_state(&self) -> (impl IndexIterator, &State);
}

//===============================
// Traits for Topology IO
//===============================
pub trait IoTopologyReader: IoReader {
    fn read_topology(&mut self) -> Result<Topology>;
}

pub trait IoTopologyWriter: IoWriter {
    fn write_topology(
        &mut self,
        data: &impl IndexAndTopologyProvider
    ) -> Result<()>;
}

//===============================
// Traits for State IO
//===============================
pub trait IoStateReader: IoReader {
    fn read_next_state(&mut self) -> Result<Option<State>>;

    fn into_states_iter(self) -> IoStateIterator<Self>
    where
        Self: Sized,
    {
        IoStateIterator { reader: self }
    }
}

pub trait IoStateWriter: IoWriter {
    fn write_next_state(
        &mut self,
        data: &impl IndexAndStateProvider,
    ) -> Result<()>;
}

//==============================
// Random access trait
//==============================
pub trait IoRandomAccess: IoStateReader {
    fn seek_frame(&mut self, fr: usize) -> Result<()>;
    fn seek_time(&mut self, t: f32) -> Result<()>;
    fn tell_first(&self) -> Result<(usize, f32)>;
    fn tell_current(&self) -> Result<(usize, f32)>;
    fn tell_last(&self) -> Result<(usize, f32)>;
}

//==================================================================
// Iterator over the frames for any type implementing IoStateReader
//==================================================================
pub struct IoStateIterator<T>
where
    T: IoStateReader,
{
    reader: T,
}

impl<T> Iterator for IoStateIterator<T>
where
    T: IoStateReader,
{
    type Item = State;
    fn next(&mut self) -> Option<Self::Item> {
        self.reader.read_next_state().expect("Error reading state")
    }
}

//================================
// General type for file handlers
//================================

pub enum FileHandler<'a> {
    Pdb(VmdMolFileHandler<'a>),
    Dcd(VmdMolFileHandler<'a>),
    Xyz(VmdMolFileHandler<'a>),
    Xtc(XtcFileHandler),
    #[cfg(feature = "gromacs")]
    Tpr(TprFileHandler),
}

pub fn get_ext(fname: &str) -> Result<&str> {
    // Get extention
    Ok(Path::new(fname)
        .extension()
        .ok_or(anyhow!("File with extension expected, given {fname}"))?
        .to_str()
        .ok_or(anyhow!("Failed getting file extension from {fname}"))?)
}

impl<'a> IoReader for FileHandler<'a> {
    fn new_reader(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::new_reader(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::new_reader(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::new_reader(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::new_reader(fname)?)),
            #[cfg(feature = "gromacs")]
            "tpr" => Ok(Self::Tpr(TprFileHandler::new_reader(fname)?)),
            _ => bail!("Unrecognized extension {ext}"),
        }
    }
}

impl<'a> IoWriter for FileHandler<'a> {
    fn new_writer(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::new_writer(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::new_writer(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::new_writer(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::new_writer(fname)?)),
            _ => bail!("Unrecognized extension {ext}"),
        }
    }
}

impl<'a> IoTopologyReader for FileHandler<'a> {
    fn read_topology(&mut self) -> Result<Topology> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => h.read_topology(),
            #[cfg(feature = "gromacs")]
            Self::Tpr(ref mut h) => h.read_topology(),
            _ => bail!("Unable to read topology"),
        }
    }
}

impl<'a> IoTopologyWriter for FileHandler<'a> {
    fn write_topology(
        &mut self,
        data: &impl IndexAndTopologyProvider,
    ) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => {
                h.write_topology(data)
            }
            _ => bail!("Unable to write topology"),
        }
    }
}

impl<'a> IoStateReader for FileHandler<'a> {
    fn read_next_state(&mut self) -> Result<Option<State>> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.read_next_state()
            }
            Self::Xtc(ref mut h) => h.read_next_state(),
            #[cfg(feature = "gromacs")]
            Self::Tpr(ref mut h) => h.read_next_state(),
        }
    }
}

impl<'a> IoStateWriter for FileHandler<'a> {
    fn write_next_state(
        &mut self,
        data: &impl IndexAndStateProvider,
    ) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.write_next_state(data)
            }
            Self::Xtc(ref mut h) => h.write_next_state(data),
            _ => bail!("Unable to write state"),
        }
    }
}

impl<'a> IoRandomAccess for FileHandler<'a> {
    fn seek_frame(&mut self, fr: usize) -> Result<()> {
        match self {
            Self::Xtc(ref mut h) => h.seek_frame(fr),
            _ => bail!("Not a random access format!"),
        }
    }

    fn seek_time(&mut self, t: f32) -> Result<()> {
        match self {
            Self::Xtc(ref mut h) => h.seek_time(t),
            _ => bail!("Not a random access format!"),
        }
    }

    fn tell_first(&self) -> Result<(usize, f32)> {
        match self {
            Self::Xtc(ref h) => h.tell_first(),
            _ => bail!("Not a random access format!"),
        }
    }

    fn tell_current(&self) -> Result<(usize, f32)> {
        match self {
            Self::Xtc(ref h) => h.tell_current(),
            _ => bail!("Not a random access format!"),
        }
    }

    fn tell_last(&self) -> Result<(usize, f32)> {
        match self {
            Self::Xtc(ref h) => h.tell_last(),
            _ => bail!("Not a random access format!"),
        }
    }
}

#[test]
fn test_read() -> Result<()>{
    use super::io::*;

    let mut r = FileHandler::new_reader("tests/topol.tpr")?;
    let mut w = FileHandler::new_writer(concat!(env!("OUT_DIR"), "/1.pdb"))?;

    let st = r.read_topology()?;
    println!("{:?}", st.atoms);

    for fr in r.into_states_iter() {
        //println!("{:?}",fr);
        w.write_topology(&st)?;
        w.write_next_state(&fr)?;
        //w.write_structure(&st).unwrap();
        //w.write_next_state_subset(&fr,0..10).unwrap();
    }

    Ok(())
}

#[test]
fn test_traj() -> Result<()> {
    use super::io::*;

    let mut r = FileHandler::new_reader("tests/no_ATP.xtc")?;
    let (max_fr,max_t) = r.tell_last()?;
    println!("max: {max_fr}:{max_t}");
    
    let (cur_fr,cur_t) = r.tell_current()?;
    println!("cur: {cur_fr}:{cur_t}");
    
    r.seek_frame(2000)?;
    let (cur_fr,cur_t) = r.tell_current()?;
    println!("cur after seek to fr 2000: {cur_fr}:{cur_t}");

    //r.seek_time(250000.0)?;
    //let (cur_fr,cur_t) = r.tell_current()?;
    //println!("cur after seek to t 250k: {cur_fr}:{cur_t}");

    Ok(())
}