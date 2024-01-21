use crate::core::{IndexIterator, State, Topology};
use anyhow::{anyhow, bail, Result};
use std::{path::Path, ops::Deref};

mod vmd_molfile_handler;
mod xtc_handler;
mod gro_handler;
#[cfg(feature = "gromacs")]
mod tpr_handler;

// Reexports
#[cfg(feature = "gromacs")]
pub use tpr_handler::TprFileHandler;
pub use vmd_molfile_handler::VmdMolFileHandler;
pub use xtc_handler::XtcFileHandler;
pub use gro_handler::GroFileHandler;

//===============================
// Traits for file opening
//===============================
pub trait IoReader {
    fn open(fname: &str) -> Result<Self>
    where
        Self: Sized;
}

pub trait IoWriter {
    fn create(fname: &str) -> Result<Self>
    where
        Self: Sized;
}

// There are the following types of data file handlers:
// (1)  All-in-once files (GRO, TPR)
//      Topology+State is read and written at once
// (2)  State only multiple times (XTC, TRR)
//      These are classical trajectories
// (3)  Topology once + multiple states (PDB, TNG)

//===============================
// Traits for writing
//===============================

pub trait IoIndexProvider {
    fn get_index(&self) -> impl IndexIterator;
}

pub trait IoTopologyProvider {
    fn get_topology(&self) -> impl Deref<Target = Topology>;
}

pub trait IoStateProvider {
    fn get_state(&self) -> impl Deref<Target = State>;
}

//===================================================
// Traits for Topology IO (separate from trajectory)
//===================================================
pub trait IoTopologyReader: IoReader {
    fn read_topology(&mut self) -> Result<Topology>;
}

pub trait IoTopologyWriter: IoWriter {
    fn write_topology(&mut self, data: &(impl IoIndexProvider + IoTopologyProvider)) -> Result<()>;
}

//===============================
// Traits for trajectory IO
//===============================
pub trait IoTrajectoryReader: IoReader {
    fn read_state(&mut self) -> Result<Option<State>>;

    fn iter_states<'a>(&'a mut self) -> IoStateIterator<'a, Self>
    where
        Self: Sized,
    {
        IoStateIterator { reader: self }
    }
}

pub trait IoTrajectoryWriter: IoWriter {
    fn write_state(&mut self, data: &(impl IoIndexProvider+IoStateProvider)) -> Result<()>;
}

//===============================
// Traits for all-in-one formats
//===============================
pub trait IoOnceReader: IoReader {
    fn read(&mut self) -> Result<(Topology,State)>;
}

pub trait IoOnceWriter: IoWriter {
    fn write(&mut self, data: &(impl IoIndexProvider + IoTopologyProvider + IoStateProvider)) -> Result<()>;
}


//==============================
// Random access trait
//==============================
pub trait IoRandomAccess: IoTrajectoryReader {
    fn seek_frame(&mut self, fr: usize) -> Result<()>;
    fn seek_time(&mut self, t: f32) -> Result<()>;
    fn tell_first(&self) -> Result<(usize, f32)>;
    fn tell_current(&self) -> Result<(usize, f32)>;
    fn tell_last(&self) -> Result<(usize, f32)>;
}

//=======================================================================
// Iterator over the frames for any type implementing IoTrajectoryReader
//=======================================================================
pub struct IoStateIterator<'a, T>
where
    T: IoTrajectoryReader,
{
    reader: &'a mut T,
}

impl<'a, T> Iterator for IoStateIterator<'a, T>
where
    T: IoTrajectoryReader,
{
    type Item = State;
    fn next(&mut self) -> Option<Self::Item> {
        self.reader.read_state().expect("Error reading state")
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
    Gro(GroFileHandler),
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
    fn open(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::open(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::open(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::open(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::open(fname)?)),
            #[cfg(feature = "gromacs")]
            "tpr" => Ok(Self::Tpr(TprFileHandler::open(fname)?)),
            _ => bail!("Unrecognized extension for reading {ext}"),
        }
    }
}

impl<'a> IoWriter for FileHandler<'a> {
    fn create(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::create(fname)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::create(fname)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::create(fname)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::create(fname)?)),
            _ => bail!("Unrecognized extension for writing {ext}"),
        }
    }
}

impl IoOnceReader for FileHandler<'_> {
    fn read(&mut self) -> Result<(Topology,State)> {
        match self {
            #[cfg(feature = "gromacs")]
            Self::Tpr(ref mut h) => h.read(),
            Self::Gro(ref mut h) => h.read(),
            _ => bail!("Not a once-read format"),
        }
    }
}

impl IoOnceWriter for FileHandler<'_> {
    fn write(&mut self, data: &(impl IoIndexProvider + IoTopologyProvider + IoStateProvider)) -> Result<()> {
        match self {
            Self::Gro(ref mut h) => h.write(data),
            _ => bail!("Not a once-write format"),
        }
    }
}


impl<'a> IoTopologyReader for FileHandler<'a> {
    fn read_topology(&mut self) -> Result<Topology> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => h.read_topology(),
            _ => bail!("Unable to read topology"),
        }
    }
}

impl<'a> IoTopologyWriter for FileHandler<'a> {
    fn write_topology(&mut self, data: &(impl IoIndexProvider+IoTopologyProvider)) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => h.write_topology(data),
            _ => bail!("Unable to write topology"),
        }
    }
}

impl<'a> IoTrajectoryReader for FileHandler<'a> {
    fn read_state(&mut self) -> Result<Option<State>> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.read_state()
            }
            Self::Xtc(ref mut h) => h.read_state(),
            _ => bail!("Not a trajectory reader format!"),
        }
    }
}

impl<'a> IoTrajectoryWriter for FileHandler<'a> {
    fn write_state(&mut self, data: &(impl IoIndexProvider+IoStateProvider)) -> Result<()> {
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.write_state(data)
            }
            Self::Xtc(ref mut h) => h.write_state(data),
            _ => bail!("Not a trajectory writer format!"),
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

#[cfg(test)]
mod tests {
    use super::{FileHandler, IoReader, IoWriter, IoTrajectoryReader};
    use crate::{io::*, core::{SelectionAll, Select, Vector3f, ModifyPos}};
    use anyhow::Result;

    #[test]
    fn test_read() -> Result<()> {
        let mut r = FileHandler::open("tests/no_ATP.xtc")?;
        let mut w = FileHandler::create(concat!(env!("OUT_DIR"), "/1.xtc"))?;

        //let st = r.read_topology()?;
        //println!("{:?}", st.atoms);

        for fr in r.iter_states() {
            println!("{}", fr.time);
            //w.write_topology(&st)?;
            w.write_state(&fr)?;
            //w.write_structure(&st).unwrap();
            //w.write_next_state_subset(&fr,0..10).unwrap();
        }

        Ok(())
    }

    #[test]
    fn test_traj() -> Result<()> {
        let mut r = FileHandler::open("tests/no_ATP.xtc")?;
        let (max_fr, max_t) = r.tell_last()?;
        println!("max: {max_fr}:{max_t}");

        let (cur_fr, cur_t) = r.tell_current()?;
        println!("cur: {cur_fr}:{cur_t}");

        r.seek_frame(2000)?;
        let (cur_fr, cur_t) = r.tell_current()?;
        println!("cur after seek to fr 2000: {cur_fr}:{cur_t}");

        //r.seek_time(250000.0)?;
        //let (cur_fr,cur_t) = r.tell_current()?;
        //println!("cur after seek to t 250k: {cur_fr}:{cur_t}");

        Ok(())
    }
    
    #[test]
    fn test_pdb() -> Result<()> {
        let mut r = FileHandler::open("tests/no_ATP.pdb")?;
        let top1 = r.read_topology()?.to_rc();
        let st1 = r.read_state()?.unwrap().to_rc();
        let st2 = (*st1).borrow().clone().to_rc();
        println!("#1: {}",(*top1).borrow().atoms.len());

        let sel = SelectionAll::new().select(&top1,&st2)?;
        sel.modify().rotate(&Vector3f::x_axis(), 45.0_f32.to_radians());
        
        let outname = concat!(env!("OUT_DIR"), "/2.pdb");
        println!("{outname}");
        let mut w = FileHandler::create(outname)?;
        w.write_topology(&sel.query())?;
        w.write_state(&st1)?;
        w.write_state(&st2)?;

        //let top2 = r.read_topology()?;
        //let st2 = r.read_next_state()?.unwrap();
        //println!("#2: {}",top2.atoms.len());
        //for fr in r.iter_states() {
        //    println!("fr {}", fr.time);
        //}
        Ok(())
    }
}
