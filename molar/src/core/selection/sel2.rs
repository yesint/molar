use crate::prelude::*;
use sorted_vec::SortedSet;
use std::{default, error::Error, marker::PhantomData};
use triomphe::Arc;
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};

trait SelectionInterface {
    fn get_index(&self) -> &[usize];
    fn get_topology(&self) -> &Topology;
    fn get_state(&self) -> &State;
}

trait SerialSelectionInterface: SelectionInterface {}

//-----------------------------------------

pub struct SelSerial {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
    _phantom: PhantomData<*const ()>,
}

impl SelectionInterface for SelSerial {
    fn get_index(&self) -> &[usize] {
        &self.index_storage
    }

    fn get_state(&self) -> &State {
        self.state.as_ref()
    }

    fn get_topology(&self) -> &Topology {
        self.topology.as_ref()
    }
}

impl SelSerial {
    pub fn new(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        let n = top.len();
        Ok(Self {
            topology: Arc::new(top),
            state: Arc::new(st),
            index_storage: Arc::new(unsafe { SVec::from_sorted((0..n).collect::<Vec<_>>()) }),
            _phantom: Default::default(),
        })
    }

    pub fn select(&self, def: impl SelectionDef) -> Result<Self, SelectionError> {
        let ind = def.into_sel_index(
            &self.topology,
            &self.state,
            self.index_storage.as_slice().into(),
        )?;
        Ok(Self {
            topology: Arc::clone(&self.topology),
            state: Arc::clone(&self.state),
            index_storage: Arc::new(ind),
            _phantom: Default::default(),
        })
    }

    pub fn split_parallel<OP, R>(&self, op: OP) -> Result<ParSplit, SelectionError>
    where
        OP: Fn(&Particle) -> Option<R>,
        R: PartialOrd,
    {
        let mut it = self.iter_particle().filter_map(|p| {
            op(&p).map(|r| (p.id,r))
        });
        // Get first valid result
        let (mut cur_i, mut cur_val)= it.next().ok_or_else(|| SelectionError::EmptySplit)?;
        
        let mut selections = vec![];
        
        // Iterate until other value found
        let mut last_i = cur_i;
        while let Some((i,val)) = it.next() {
            last_i = i;
            if val != cur_val {
                // Create new selection
                let sel = self.select(cur_i..i)?;
                selections.push( unsafe{std::mem::transmute::<SelSerial,SelPar>(sel)} );
                // Reset cur
                cur_i = i;
                cur_val = val;
            }
        }

        // Add last selection
        let sel = self.select(cur_i..=last_i)?;
        selections.push( unsafe{std::mem::transmute::<SelSerial,SelPar>(sel)} );
        
        Ok(ParSplit {
            topology: Arc::clone(&self.topology),
            state: Arc::clone(&self.state),
            selections,
            _phantom: Default::default(),
        })
    }
}

impl LenProvider for SelSerial {
    fn len(&self) -> usize {
        self.index_storage.len()
    }
}

impl RandomParticleProvider for SelSerial {
    unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        Particle {
            atom: self.topology.nth_atom_unchecked(i),
            pos: self.state.nth_pos_unchecked(i),
            id: i,
        }
    }
}

impl RandomAtomProvider for SelSerial {
    fn num_atoms(&self) -> usize {
        self.topology.len()
    }

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        self.topology.nth_atom_unchecked(self.index_storage[i])
    }
}

impl RandomPosProvider for SelSerial {
    fn num_pos(&self) -> usize {
        self.state.len()
    }

    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        self.state.nth_pos_unchecked(self.index_storage[i])
    }
}

impl IndexProvider for SelSerial {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.index_storage.iter().cloned()
    }
}

impl ParticleIterProvider for SelSerial {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
        self.iter_index().map(|i| unsafe {
            Particle {
                atom: self.topology.nth_atom_unchecked(i),
                pos: self.state.nth_pos_unchecked(i),
                id: i,
            }
        })
    }
}
//--------------------------------------------
pub struct SelBuilder {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
    _phantom: PhantomData<*const ()>,
}

impl SelectionInterface for SelBuilder {
    fn get_index(&self) -> &[usize] {
        if self.index_storage.len() > self.topology.len()
            || self.index_storage.len() > self.state.len()
        {
            panic!("selection index is invalidated",);
        }
        &self.index_storage
    }

    fn get_state(&self) -> &State {
        self.state.as_ref()
    }

    fn get_topology(&self) -> &Topology {
        self.topology.as_ref()
    }
}

impl SelBuilder {
    pub fn select(&self, def: impl SelectionDef) -> Result<Self, SelectionError> {
        let ind = def.into_sel_index(&self.topology, &self.state, self.get_index().into())?;
        Ok(Self {
            topology: Arc::clone(&self.topology),
            state: Arc::clone(&self.state),
            index_storage: Arc::new(ind),
            _phantom: Default::default(),
        })
    }
}

pub struct Builder {
    topology: Arc<Topology>,
    state: Arc<State>,
    _phantom: PhantomData<*const ()>,
}

impl Builder {
    pub fn new(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        Ok(Self {
            topology: Arc::new(top),
            state: Arc::new(st),
            _phantom: Default::default(),
        })
    }

    pub fn select(&self, def: impl SelectionDef) -> Result<SelBuilder, SelectionError> {
        let ind = def.into_sel_index(&self.topology, &self.state, None)?;
        Ok(SelBuilder {
            topology: Arc::clone(&self.topology),
            state: Arc::clone(&self.state),
            index_storage: Arc::new(ind),
            _phantom: Default::default(),
        })
    }
}

//---------------------------------------------------

pub struct SelPar {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
}

impl SelectionInterface for SelPar {
    fn get_index(&self) -> &[usize] {
        &self.index_storage
    }

    fn get_state(&self) -> &State {
        self.state.as_ref()
    }

    fn get_topology(&self) -> &Topology {
        self.topology.as_ref()
    }
}

impl SelPar {
    // Subselection of SelPar is serial
    pub fn select(&self, def: impl SelectionDef) -> Result<SelSerial, SelectionError> {
        let ind = def.into_sel_index(&self.topology, &self.state, self.index_storage.as_slice().into())?;
        Ok(SelSerial {
            topology: Arc::clone(&self.topology),
            state: Arc::clone(&self.state),
            index_storage: Arc::new(ind),
            _phantom: Default::default(),
        })
    }
}

impl LenProvider for SelPar {
    fn len(&self) -> usize {
        self.index_storage.len()
    }
}

impl RandomAtomProvider for SelPar {
    fn num_atoms(&self) -> usize {
        self.index_storage.len()
    }

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        self.topology.nth_atom_unchecked(self.index_storage[i])
    }
}

pub struct ParSplit {
    topology: Arc<Topology>,
    state: Arc<State>,
    selections: Vec<SelPar>,
    _phantom: PhantomData<*const ()>,
}

impl ParSplit {
    /// Returns parallel iterator over stored parallel selections.
    pub fn par_iter(&mut self) -> rayon::slice::Iter<'_, SelPar> {
        self.selections.par_iter()
    }

    /// Returns parallel mutable iterator over stored parallel selections.
    pub fn par_iter_mut(&mut self) -> rayon::slice::IterMut<'_, SelPar> {
        self.selections.par_iter_mut()
    }
}

mod tests {
    use rayon::iter::ParallelIterator;

    use crate::{core::{selection::sel2::SelSerial, LenProvider, RandomAtomProvider, RandomPosProvider, SelectionError}, io::FileHandler};


    #[test]
    fn test1() -> anyhow::Result<()> {
        let (top,st) = FileHandler::open("tests/albumin.pdb")?.read()?;
        let sel = SelSerial::new(top, st)?;
        
        let mut par = sel.split_parallel(|p|{
            if p.atom.resname != "SOL" {
                Some(p.atom.resindex)
            } else {
                None
            }
        })?;

        par.par_iter().try_for_each(|sel|{
            println!("{} {}",sel.len(), sel.first_atom().resname);
            if sel.first_atom().resname == "ALA" {
                let subsel = sel.select("name CA")?;
                let pos = subsel.first_pos();
                println!("{}",pos);
            }
            Ok::<_,SelectionError>(())
        })?;

        // Add serial selection
        let ca = sel.select("name CA")?;
        println!("#ca: {}",ca.len());

        Ok(())
    }
}