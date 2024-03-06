use super::{
    measure::{GuardedQuery, MeasureMasses, MeasurePeriodic, MeasurePos},
    modify::{GuardedModify, ModifyPos, ModifyRandomAccess},
    providers::{
        AtomsProvider, BoxProvider, MassesProvider, PosMutProvider, PosProvider,
        RandomPosMutProvider,
    },
    Atom, AtomIterator, PeriodicBox, Pos, PosIterator, PosMutIterator, SelectionSplitIterator,
    State, StateRc, Topology, TopologyRc,
};
use crate::io::{FileHandler, IndexProvider, StateProvider, TopologyProvider};
use anyhow::{bail, Result};
use itertools::Itertools;
use std::{
    cell::{Ref, RefMut},
    collections::HashMap,
    rc::Rc,
};

pub use super::selection_parser::SelectionExpr;
//-----------------------------------------------------------------

/// Trait whic provides select method operating with generic smart pointers
pub trait Select {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection>;
}

//-----------------------------------------------------------------
fn check_sizes(topology: &Topology, state: &State) -> Result<()> {
    let n1 = topology.atoms.len();
    let n2 = state.coords.len();
    match n1 == n2 {
        true => Ok(()),
        false => bail!(
            "Structure and State are incompatible (sizes {} and {})",
            n1,
            n2
        ),
    }
}

pub struct SelectionAll {}

impl SelectionAll {
    pub fn new() -> Self {
        Self {}
    }
}

impl Select for SelectionAll {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        let st = state.borrow();
        check_sizes(&topology.borrow(), &st)?;
        let index: Vec<usize> = (0..st.coords.len()).collect();
        if index.len() > 0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for SelectionExpr {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        check_sizes(&topology.borrow(), &state.borrow())?;
        let index = self
            .apply_whole(&topology.borrow(), &state.borrow())
            .unwrap();
        if index.len() > 0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for str {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        let top = topology.borrow();
        let st = state.borrow();
        check_sizes(&top, &st)?;
        let index = SelectionExpr::try_from(self)
            .unwrap()
            .apply_whole(&top, &st)
            .unwrap();
        if index.len() > 0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for std::ops::Range<usize> {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        let top = topology.borrow();
        let st = state.borrow();
        check_sizes(&top, &st)?;
        let n = top.atoms.len();
        if self.start > n - 1 || self.end > n - 1 {
            bail!(
                "Index range {}:{} is invalid, 0:{} is allowed.",
                self.start,
                self.end,
                top.atoms.len()
            );
        }
        let index: Vec<usize> = self.clone().collect();
        if index.len() > 0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for Vec<usize> {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        let index: Vec<usize> = self.iter().cloned().sorted().dedup().collect();
        if index.len() > 0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

//---------------------------------------

/*
If Topology and State are wrapped into Rc<> then they
are guaranteed to be accessed from the single thread only and the only possible
case of memory unsafety is invalidation of references when adding or removing
elements from atoms and coords arrays.

If number of atoms is never modified (creation requires a dedicated builder
and modification consumes the State or Topology and returns a new one), then
it is impossible to get an undefined bahavior even with alising references
to *individual* atoms.

In principle even aliasing access from multiple threads to *individual* atoms is safe
in a sence that no UB could happen. The data races may result in incorrect values,
but no unsafe memory access is possible.

So Rc<UnsafeCell<State>> or Arc<UnsafeCell<State>> should be fine if API 
makes changing the number of atoms impossible.

- Should pbox be always immutable?
- Should we wrap into UnsafeCell<> only atoms and coords fields?
- What about bonds? Do we ever need them mutable?
*/

pub struct Selection {
    topology: TopologyRc,
    state: StateRc,
    index: Vec<usize>,
}

impl Selection {
    //===================
    // Subselections
    //===================

    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Selection> {
        let index = expr.apply_subset(
            &self.topology.borrow(),
            &self.state.borrow(),
            self.index.iter().cloned(),
        )?;
        if index.len() > 0 {
            Ok(Selection {
                topology: self.topology.clone(),
                state: self.state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }

    /// Subselection from string
    pub fn subsel_from_str(&self, sel_str: &str) -> Result<Selection> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_from_local_range(&self, range: std::ops::Range<usize>) -> Result<Selection> {
        if range.end >= self.index.len() {
            bail!(
                "Invalid local sub-range: {}:{}, valid range: 0:{}",
                range.start,
                range.end,
                self.index.len() - 1
            );
        }

        // Translate range of local indexes to global indexes
        let index: Vec<usize> = self
            .index
            .iter()
            .cloned()
            .skip(range.start)
            .take(range.len())
            .collect();

        Ok(Selection {
            topology: self.topology.clone(),
            state: self.state.clone(),
            index,
        })
    }

    /// Subselection from iterator that provides local selection indexes
    pub fn subsel_from_iter(
        &self,
        iter: impl ExactSizeIterator<Item = usize>,
    ) -> Result<Selection> {
        // Remove duplicates if any
        let index: Vec<usize> = iter.sorted().dedup().map(|i| self.index[i]).collect();
        // Now it's safe to call
        unsafe { self.subsel_from_vec_unchecked(index) }
    }

    // This method doesn't check if the vector has duplicates and thus unsafe
    pub unsafe fn subsel_from_vec_unchecked(&self, index: Vec<usize>) -> Result<Selection> {
        if index.len() > 0 {
            Ok(Selection {
                topology: self.topology.clone(),
                state: self.state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computer lazily when iterating.
    pub fn split_contig<T, F>(&self, func: F) -> SelectionSplitIterator<'_, T, F>
    where
        T: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> T,
    {
        SelectionSplitIterator::new(self, func)
    }

    /// Splits selection into pieces that could be disjoint
    /// according to the value of function.
    /// The number of selection correspond to the distinct values returned by `func`.
    /// Selections are stored in a provided container.
    pub fn split<T, F, C>(&self, func: F) -> C
    where
        T: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(usize, &Atom, &Pos) -> T,
        C: FromIterator<Selection> + Default,
    {
        let guard = self.guard();
        let mut ids = HashMap::<T, Vec<usize>>::default();

        for i in 0..guard.len() {
            let (ind, at, pos) = unsafe { guard.nth_unchecked(i) };
            let id = func(ind, at, pos);
            if let Some(el) = ids.get_mut(&id) {
                el.push(ind);
            } else {
                ids.insert(id, vec![ind]);
            }
        }

        C::from_iter(
            ids.into_values()
                .map(|ind| unsafe {
                    self.subsel_from_vec_unchecked(ind).unwrap() 
                }),
        )
    }

    //======================
    // Combining selections
    //======================
    pub fn union(&mut self, other: &Selection) -> Result<()> {
        if !Rc::ptr_eq(&self.topology, &other.topology) || !Rc::ptr_eq(&self.state, &other.state) {
            bail!("Can't combine selection pointing to different topologies or states!")
        };
        self.index.extend(other.index.iter());
        self.index = self.index.iter().cloned().sorted().dedup().collect();
        Ok(())
    }

    pub fn get_index_vec(&self) -> Vec<usize> {
        self.index.clone()
    }

    pub fn save(&self, fname: &str) -> Result<()> {
        let mut h = FileHandler::create(fname)?;
        h.write(self)
    }
}
//----------------------------------------------------

/// Scoped guard giving read-only access to selection
pub struct SelectionQueryGuard<'a> {
    topology_ref: Ref<'a, Topology>,
    state_ref: Ref<'a, State>,
    index: &'a Vec<usize>,
}

impl SelectionQueryGuard<'_> {
    pub fn len(&self) -> usize {
        self.index.len()
    }

    pub fn nth(&self, i: usize) -> (usize, &Atom, &Pos) {
        let ind = self.index[i];
        assert!(ind < self.len());
        (
            ind,
            &self.topology_ref.atoms[ind],
            &self.state_ref.coords[ind],
        )
    }

    pub unsafe fn nth_unchecked(&self, i: usize) -> (usize, &Atom, &Pos) {
        let ind = *self.index.get_unchecked(i);
        (
            ind,
            &self.topology_ref.atoms.get_unchecked(ind),
            &self.state_ref.coords.get_unchecked(ind),
        )
    }
}

/// Scoped guard giving read-write access to selection
#[allow(dead_code)]
pub struct SelectionModifyGuard<'a> {
    topology_ref: RefMut<'a, Topology>,
    state_ref: RefMut<'a, State>,
    index: &'a Vec<usize>,
}

//---------------------------------------------
// Implement traits for IO

impl<'a> IndexProvider for SelectionQueryGuard<'a> {
    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.index.iter().cloned()
    }
}

impl<'a> TopologyProvider for SelectionQueryGuard<'a> {
    fn iter_atoms(&self) -> impl Iterator<Item = &super::Atom> {
        self.index.iter().map(|i| &self.topology_ref.atoms[*i])
    }

    fn num_atoms(&self) -> usize {
        self.index.len()
    }
}

impl<'a> StateProvider for SelectionQueryGuard<'a> {
    fn iter_coords(&self) -> impl Iterator<Item = &Pos> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }

    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state_ref.pbox.as_ref()
    }

    fn get_time(&self) -> f32 {
        self.state_ref.time
    }

    fn num_coords(&self) -> usize {
        self.index.len()
    }
}

//==================================================================
// Implement analysis traits

impl GuardedQuery for Selection {
    type Guard<'a> = SelectionQueryGuard<'a>;
    fn guard<'a>(&'a self) -> Self::Guard<'a> {
        SelectionQueryGuard {
            topology_ref: self.topology.borrow(),
            state_ref: self.state.borrow(),
            index: &self.index,
        }
    }
}

impl BoxProvider for SelectionQueryGuard<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state_ref.pbox.as_ref()
    }
}

impl MeasurePeriodic for Selection {}

impl PosProvider for SelectionQueryGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }
}

impl MeasurePos for Selection {}

impl AtomsProvider for SelectionQueryGuard<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.index.iter().map(|i| &self.topology_ref.atoms[*i])
    }
}

impl MassesProvider for SelectionQueryGuard<'_> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.index.iter().map(|i| self.topology_ref.atoms[*i].mass)
    }
}

impl MeasureMasses for Selection {}

//-------------------------------------------------------
impl GuardedModify for Selection {
    type GuardMut<'a> = SelectionModifyGuard<'a>;
    fn guard_mut<'a>(&'a self) -> Self::GuardMut<'a> {
        SelectionModifyGuard {
            topology_ref: self.topology.borrow_mut(),
            state_ref: self.state.borrow_mut(),
            index: &self.index,
        }
    }
}

impl BoxProvider for SelectionModifyGuard<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state_ref.get_box()
    }
}

impl PosMutProvider for SelectionModifyGuard<'_> {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_> {
        self.index
            .iter()
            .map(|i| unsafe { &mut *self.state_ref.coords.as_mut_ptr().add(*i) })
    }
}

impl PosProvider for SelectionModifyGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }
}

impl RandomPosMutProvider for SelectionModifyGuard<'_> {
    #[inline(always)]
    fn nth_pos_mut(&mut self, i: usize) -> &mut Pos {
        &mut self.state_ref.coords[self.index[i]]
    }
}

impl ModifyPos for Selection {}
impl ModifyRandomAccess for Selection {}

//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {
    use super::{Selection, SelectionAll};
    use crate::{
        core::{
            providers::PosMutProvider, selection::Select, GuardedModify, GuardedQuery,
            MeasureMasses, MeasurePos, ModifyPos, ModifyRandomAccess, State, Topology, Vector3f,
            PBC_FULL,
        },
        io::*,
    };
    use lazy_static::lazy_static;

    fn read_test_pdb() -> (Topology, State) {
        let mut h = FileHandler::open("tests/no_ATP.pdb").unwrap();
        //let top = h.read_topology().unwrap();
        //let state = h.read_state().unwrap().unwrap();
        //(top, state)
        let (top, st) = h.read_raw().unwrap();
        (top, st)
    }

    // Read the test PDB file once and provide the content for tests
    lazy_static! {
        static ref SS: (Topology, State) = read_test_pdb();
    }

    fn make_sel() -> anyhow::Result<Selection> {
        let t = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = SelectionAll {}.select(&t, &s)?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<Selection> {
        let t = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = "not resname TIP3 POT CLA".select(&t, &s)?;
        Ok(sel)
    }

    #[test]
    fn test_make_sel() -> anyhow::Result<()> {
        make_sel().map(|_| ())
    }

    #[test]
    fn test_measure() -> anyhow::Result<()> {
        let sel = make_sel()?;
        println!("before {}", sel.guard().iter_coords().next().unwrap());

        let (minv, maxv) = sel.min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}", sel.guard().iter_coords().next().unwrap());
        println!("{:?}", sel.min_max());
        Ok(())
    }

    #[test]
    fn test_measure_pbc() -> anyhow::Result<()> {
        let sel = make_sel()?;

        let cm = sel.center_of_mass()?;
        println!("{cm}");
        Ok(())
    }

    #[test]
    fn test_translate() -> anyhow::Result<()> {
        let sel = make_sel()?;

        println!("before {}", sel.guard_mut().iter_pos_mut().next().unwrap());
        sel.translate(Vector3f::new(10.0, 10.0, 10.0));
        println!("after {}", sel.guard_mut().iter_pos_mut().next().unwrap());
        Ok(())
    }

    #[test]
    fn test_write_to_file() -> anyhow::Result<()> {
        let sel = make_sel()?;

        let mut h = FileHandler::create("f.pdb")?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.unwrap_connectivity_dim(0.2, &PBC_FULL)?;

        let mut h = FileHandler::create("unwrapped.pdb")?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn eigen_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let sel2 = make_sel_prot()?;

        sel2.rotate(&Vector3f::x_axis(), 80.0_f32.to_radians());

        sel2.save("tests/sel2.pdb")?;
        sel1.save("tests/sel1_before.pdb")?;
        println!("Initial RMSD:{}", Selection::rmsd_mw(&sel1, &sel2)?);

        let m = Selection::fit_transform(&sel1, &sel2)?;
        println!("{m}");

        sel1.apply_transform(&m);

        sel1.save("tests/sel1_after.pdb")?;
        println!("Final RMSD:{}", Selection::rmsd_mw(&sel1, &sel2)?);
        Ok(())
    }
}
