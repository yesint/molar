use crate::core::{Atom, GuardedQuery, IndexIterator, PeriodicBox, Pos, State, StateRc, Topology, TopologyRc};
use anyhow::{anyhow, bail, Result};
use std::path::Path;

mod vmd_molfile_handler;
mod xtc_handler;
mod gro_handler;
#[cfg(feature = "gromacs")]
mod tpr_handler;
mod file_content;
mod io_splitter;

// Reexports
#[cfg(feature = "gromacs")]
pub use tpr_handler::TprFileHandler;
pub use vmd_molfile_handler::{VmdMolFileHandler, VmdMolFileType};
pub use xtc_handler::XtcFileHandler;
pub use gro_handler::GroFileHandler;

pub use file_content::FileContent;
use self::io_splitter::IoSplitter;

//===============================
// Traits for file opening
//===============================

// There are the following types of data file handlers:
// (1)  All-in-once files (GRO, TPR)
//      Topology+State is read and written at once
// (2)  State only multiple times (XTC, TRR)
//      These are classical trajectories
// (3)  Topology once + multiple states (PDB, TNG)

//===============================
// Traits for writing
//===============================
pub trait IndexProvider {    
    fn iter_index(&self) -> impl Iterator<Item = usize>;    
}

pub trait TopologyProvider {    
    fn iter_atoms(&self) -> impl Iterator<Item = &Atom>;
    fn num_atoms(&self) -> usize;
}

pub trait StateProvider {
    fn iter_coords(&self) -> impl Iterator<Item = &Pos>;
    fn get_box(&self) -> Option<&PeriodicBox>;
    fn get_time(&self) -> f32;
    fn num_coords(&self) -> usize;
}

//=======================================================================
// Iterator over the frames for any type implementing IoTrajectoryReader
//=======================================================================
pub struct IoStateIterator<'a> {
    reader: FileHandler<'a>,
}

impl Iterator for IoStateIterator<'_> {
    type Item = StateRc;
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
    Tpr(IoSplitter<TprFileHandler>),
    Gro(IoSplitter<GroFileHandler>),
}

pub fn get_ext(fname: &str) -> Result<&str> {
    // Get extention
    Ok(Path::new(fname)
        .extension()
        .ok_or(anyhow!("File with extension expected, given {fname}"))?
        .to_str()
        .ok_or(anyhow!("Failed getting file extension from {fname}"))?
    )
}

impl FileHandler<'_> {
    pub fn open(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::open(fname,VmdMolFileType::Pdb)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::open(fname,VmdMolFileType::Dcd)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::open(fname,VmdMolFileType::Xyz)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::open(fname)?)),
            "gro" => Ok(Self::Gro(IoSplitter::new(GroFileHandler::open(fname)?))),
            #[cfg(feature = "gromacs")]
            "tpr" => Ok(Self::Tpr(IoSplitter::new(TprFileHandler::open(fname)?))),
            _ => bail!("Unrecognized extension for reading {ext}"),
        }
    }

    pub fn create(fname: &str) -> Result<Self> {
        let ext = get_ext(fname)?;
        match ext {
            "pdb" => Ok(Self::Pdb(VmdMolFileHandler::create(fname,VmdMolFileType::Pdb)?)),
            "dcd" => Ok(Self::Dcd(VmdMolFileHandler::create(fname,VmdMolFileType::Dcd)?)),
            "xyz" => Ok(Self::Xyz(VmdMolFileHandler::create(fname,VmdMolFileType::Xyz)?)),
            "xtc" => Ok(Self::Xtc(XtcFileHandler::create(fname)?)),
            "gro" => Ok(Self::Gro(IoSplitter::new(GroFileHandler::create(fname)?))),
            _ => bail!("Unrecognized extension for writing {ext}"),
        }
    }

    pub fn read_raw(&mut self) -> Result<(Topology,State)> {
        let (top,st) = match self {
            #[cfg(feature = "gromacs")]
            Self::Tpr(ref mut h) => h.read()?,
            Self::Gro(ref mut h) => h.read()?,
            Self::Pdb(ref mut h) => {
                let top = h.read_topology()?;
                let st = h.read_state()?.ok_or_else(|| anyhow!("Can't read first state!"))?;
                (top,st)
            },
            _ => bail!("Not a once-read format"),
        };
        Ok((top,st))
    }

    pub fn read(&mut self) -> Result<(TopologyRc,StateRc)> {
        let (top,st) = self.read_raw()?;
        Ok((top.to_rc(),st.to_rc()))
    }

    pub fn write<'a,T>(&mut self, data: &'a T) -> Result<()> 
    where T: GuardedQuery + 'a, T::Guard<'a>: TopologyProvider+StateProvider
    {
        let dp = data.guard();
        match self {
            Self::Gro(ref mut h) => h.handler.write(&dp),
            Self::Pdb(ref mut h) => {
                h.write_topology(&dp)?;
                h.write_state(&dp)?;
                Ok(())
            }
            _ => bail!("Not a once-write format"),
        }
    }

    pub fn read_topology_raw(&mut self) -> Result<Topology> {
        let top = match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => h.read_topology()?,
            _ => bail!("Unable to read topology"),
        };
        Ok(top)
    }

    pub fn read_topology(&mut self) -> Result<TopologyRc> {
        Ok(self.read_topology_raw()?.to_rc())
    }

    pub fn write_topology<'a,T>(&mut self, data: &'a T) -> Result<()> 
    where T: GuardedQuery + 'a, T::Guard<'a>: TopologyProvider
    {
        let dp = data.guard();
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) => h.write_topology(&dp),
            _ => bail!("Unable to write topology"),
        }
    }

    pub fn read_state_raw(&mut self) -> Result<Option<State>> {
        let st = match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.read_state()?
            },
            Self::Xtc(ref mut h) => {
                h.read_state()?
            },
            _ => bail!("Not a trajectory reader format!"),
        };
        Ok(st)
    }

    pub fn read_state(&mut self) -> Result<Option<StateRc>> {
        self.read_state_raw()?.map_or(Ok(None), |v| Ok(Some(v.to_rc())))
    }

    pub fn write_state<'a,T>(&mut self, data: &'a T) -> Result<()> 
    where T: GuardedQuery + 'a, T::Guard<'a>: StateProvider
    {
        let dp = data.guard();
        match self {
            Self::Pdb(ref mut h) | Self::Xyz(ref mut h) | Self::Dcd(ref mut h) => {
                h.write_state(&dp)
            }
            Self::Xtc(ref mut h) => h.write_state(&dp),
            _ => bail!("Not a trajectory writer format!"),
        }
    }

    pub fn seek_frame(&mut self, fr: usize) -> Result<()> {
        match self {
            Self::Xtc(ref mut h) => h.seek_frame(fr),
            _ => bail!("Not a random access format!"),
        }
    }

    pub fn seek_time(&mut self, t: f32) -> Result<()> {
        match self {
            Self::Xtc(ref mut h) => h.seek_time(t),
            _ => bail!("Not a random access format!"),
        }
    }

    pub fn tell_first(&self) -> Result<(usize, f32)> {
        match self {
            Self::Xtc(ref h) => h.tell_first(),
            _ => bail!("Not a random access format!"),
        }
    }

    pub fn tell_current(&self) -> Result<(usize, f32)> {
        match self {
            Self::Xtc(ref h) => h.tell_current(),
            _ => bail!("Not a random access format!"),
        }
    }

    pub fn tell_last(&self) -> Result<(usize, f32)> {
        match self {
            Self::Xtc(ref h) => h.tell_last(),
            _ => bail!("Not a random access format!"),
        }
    }
}

impl<'a> IntoIterator for FileHandler<'a> {
    type Item = StateRc;
    type IntoIter = IoStateIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        IoStateIterator{reader: self}
    }
}

#[cfg(test)]
mod tests {
    use super::FileHandler;
    use crate::core::{ModifyPos, Select, SelectionAll, Vector3f};
    use anyhow::Result;

    #[test]
    fn test_read() -> Result<()> {
        let r = FileHandler::open("tests/no_ATP.xtc")?;
        let mut w = FileHandler::create(concat!(env!("OUT_DIR"), "/1.xtc"))?;

        //let st = r.read_topology()?;
        //println!("{:?}", st.atoms);

        for fr in r {
            w.write_state(&fr)?;
            //let f = fr.into_inner();
            //println!("{}", f.time);
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
        let top1 = r.read_topology()?;
        let st1 = r.read_state()?.unwrap();
        let st2 = (*st1).borrow().clone().to_rc();
        println!("#1: {}",(*top1).borrow().atoms.len());

        let sel = SelectionAll::new().select(&top1,&st2)?;
        sel.rotate(&Vector3f::x_axis(), 45.0_f32.to_radians());
        
        let outname = concat!(env!("OUT_DIR"), "/2.pdb");
        println!("{outname}");
        let mut w = FileHandler::create(outname)?;
        w.write_topology(&sel)?;
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
