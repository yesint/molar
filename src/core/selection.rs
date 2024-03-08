use super::{
    measure::{MeasureMasses, MeasurePeriodic, MeasurePos},
    modify::{ModifyPos, ModifyRandomAccess},
    providers::*,
    Atom, AtomIterator, PeriodicBox, Pos, PosIterator, PosMutIterator, State, Topology,
};
use crate::io::{FileHandler, IndexProvider, StateProvider, TopologyProvider};
use anyhow::{anyhow, bail, Result};
use itertools::Itertools;
use std::{collections::HashMap, rc::Rc};

pub use super::selection_parser::SelectionExpr;
//-----------------------------------------------------------------

/// Trait which provides select method operating with generic smart pointers
pub trait Select {
    fn select(&self, topology: &Rc<Topology>, state: &Rc<State>) -> Result<Selection>;
}

//-----------------------------------------------------------------
fn check_sizes(topology: &Topology, state: &State) -> Result<()> {
    let n1 = topology.num_atoms();
    let n2 = state.num_coords();
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
    fn select(&self, topology: &Rc<Topology>, state: &Rc<State>) -> Result<Selection> {
        check_sizes(&topology, &state)?;
        let index: Vec<usize> = (0..state.num_coords()).collect();
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
    fn select(&self, topology: &Rc<Topology>, state: &Rc<State>) -> Result<Selection> {
        check_sizes(&topology, &state)?;
        let index = self.apply_whole(&topology, &state).unwrap();
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
    fn select(&self, topology: &Rc<Topology>, state: &Rc<State>) -> Result<Selection> {
        check_sizes(&topology, &state)?;
        let index = SelectionExpr::try_from(self)
            .unwrap()
            .apply_whole(&topology, &state)
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
    fn select(&self, topology: &Rc<Topology>, state: &Rc<State>) -> Result<Selection> {
        check_sizes(&topology, &state)?;
        let n = topology.num_atoms();
        if self.start > n - 1 || self.end > n - 1 {
            bail!(
                "Index range {}:{} is invalid, 0:{} is allowed.",
                self.start,
                self.end,
                topology.num_atoms()
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
    fn select(&self, topology: &Rc<Topology>, state: &Rc<State>) -> Result<Selection> {
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
    topology: Rc<Topology>,
    state: Rc<State>,
    index: Vec<usize>,
}

impl Selection {
    pub fn len(&self) -> usize {
        self.num_atoms()
    }

    //===================
    // Subselections
    //===================

    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Selection> {
        let index = expr.apply_subset(&self.topology, &self.state, self.index.iter().cloned())?;
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
        let index = iter
            .sorted()
            .dedup()
            .map(|i| {
                self.index.get(i)
                .cloned()
                .ok_or_else(|| anyhow!("Index {} is out of allowed range [0:{}]",i,self.index.len()))
            })
            .collect::<Result<Vec<usize>>>()?;
        // Now it's safe to call
        unsafe { self.subsel_from_vec_unchecked(index) }
    }

    // This method doesn't check if the vector has duplicates and thus unsafe
    unsafe fn subsel_from_vec_unchecked(&self, index: Vec<usize>) -> Result<Selection> {
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

    pub unsafe fn nth_unchecked(&self, i: usize) -> (usize, &Atom, &Pos) {
        let ind = *self.index.get_unchecked(i);
        (
            ind,
            self.topology.nth_atom_unchecked(ind),
            self.state.nth_pos_unchecked(ind),
        )
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
        let mut ids = HashMap::<T, Vec<usize>>::default();

        for i in 0..self.index.len() {
            let (ind, at, pos) = unsafe { self.nth_unchecked(i) };
            let id = func(ind, at, pos);
            if let Some(el) = ids.get_mut(&id) {
                el.push(ind);
            } else {
                ids.insert(id, vec![ind]);
            }
        }

        C::from_iter(
            ids.into_values()
                .map(|ind| unsafe { self.subsel_from_vec_unchecked(ind).unwrap() }),
            // This should never fail because `ind` can't be empty
        )
    }

    // Pre-defined splitters
    pub fn split_contig_resid(&self) -> SelectionSplitIterator<'_,i32,fn(usize, &Atom, &Pos)->i32> {
        self.split_contig(|_,at,_| at.resid)
    }

    pub fn split_resid<C>(&self) -> C
    where        
        C: FromIterator<Selection> + Default,
    {
        self.split(|_,at,_| at.resid)
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

    pub fn first(&self) -> (usize, &Atom, &Pos) {
        unsafe { self.nth_unchecked(0) }
    }

    pub fn nth(&self, i: usize) -> Result<(usize, &Atom, &Pos)> {
        if i>self.len() {
            bail!("Index {} is beyond the allowed range [0:{}]",i,self.len())
        }
        Ok(unsafe {self.nth_unchecked(i)})
    }

    pub fn iter(&self) -> SelectionIterator {
        SelectionIterator {
            sel: self,
            cur: 0,
        }
    }
}

//---------------------------------------------
pub struct SelectionIterator<'a> {
    sel: &'a Selection,
    cur: usize,
}

impl<'a> Iterator for SelectionIterator<'a> {
    type Item = (usize, &'a Atom, &'a Pos);
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.sel.len() {
            let ret = unsafe{ self.sel.nth_unchecked(self.cur) };
            self.cur += 1;
            Some(ret)
        } else {
            None
        }
    }
}
//---------------------------------------------
// Implement traits for IO

impl IndexProvider for Selection {
    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.index.iter().cloned()
    }
}

impl TopologyProvider for Selection {
    fn num_atoms(&self) -> usize {
        self.index.len()
    }
}

impl StateProvider for Selection {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.index.len()
    }
}

//==================================================================
// Implement analysis traits

impl BoxProvider for Selection {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl MeasurePeriodic for Selection {}

impl PosProvider for Selection {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe { self.index.iter().map(|i| self.state.nth_pos_unchecked(*i)) }
    }
}

impl MeasurePos for Selection {}

impl AtomsProvider for Selection {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.topology.nth_atom_unchecked(*i))
        }
    }
}

impl MassesProvider for Selection {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.topology.nth_atom_unchecked(*i).mass)
        }
    }
}

impl MeasureMasses for Selection {}

//-------------------------------------------------------

impl PosMutProvider for Selection {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.state.nth_pos_unchecked_mut(*i))
        }
    }
}

impl RandomPosMutProvider for Selection {
    #[inline(always)]
    unsafe fn nth_pos_unchecked_mut(&self, i: usize) -> &mut Pos {
        let ind = *self.index.get_unchecked(i);
        self.state.nth_pos_unchecked_mut(ind)
    }
}

impl ModifyPos for Selection {}
impl ModifyRandomAccess for Selection {}

//-------------------------------------------------------
// Splitting iterator
//-------------------------------------------------------

pub struct SelectionSplitIterator<'a, T, F> {
    sel: &'a Selection,
    func: F,
    counter: usize,
    id: T,
}

impl SelectionSplitIterator<'_, (), ()> {
    pub fn new<T, F>(sel: &Selection, func: F) -> SelectionSplitIterator<'_, T, F>
    where
        T: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> T,
    {
        SelectionSplitIterator {
            sel,
            func,
            counter: 0,
            id: T::default(),
        }
    }
}

impl<T, F> Iterator for SelectionSplitIterator<'_, T, F>
where
    T: Default + std::cmp::PartialEq,
    F: Fn(usize, &Atom, &Pos) -> T,
{
    type Item = Selection;
    fn next(&mut self) -> Option<Self::Item> {
        let mut index = vec![];
        while self.counter < self.sel.len() {
            let (i, at, pos) = unsafe { self.sel.nth_unchecked(self.counter) };
            let id = (self.func)(i, at, pos);

            if id == self.id {
                // Current selection continues. Add current index
                index.push(i);
            } else if index.is_empty() {
                // The very first id is not default, this is Ok, add index
                // and update self.id
                self.id = id;
                index.push(i);
            } else {
                // The end of current selection
                self.id = id; // Update self.id for the next selection
                return unsafe { Some(self.sel.subsel_from_vec_unchecked(index).unwrap()) };
            }
            // Next element
            self.counter += 1;
        }

        // Return any remaining index as last selection
        if !index.is_empty() {
            return unsafe { Some(self.sel.subsel_from_vec_unchecked(index).unwrap()) };
        }

        // If we are here stop iterating
        None
    }
}

//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {
    use super::{Selection, SelectionAll};
    use crate::{
        core::{
            providers::PosProvider, selection::{AtomsProvider, Select}, MeasureMasses, MeasurePos, ModifyPos,
            ModifyRandomAccess, State, Topology, Vector3f, PBC_FULL,
        },
        io::*,
    };

    fn read_test_pdb() -> (Topology, State) {
        let mut h = FileHandler::open("tests/no_ATP.pdb").unwrap();
        //let top = h.read_topology().unwrap();
        //let state = h.read_state().unwrap().unwrap();
        //(top, state)
        let (top, st) = h.read_raw().unwrap();
        (top, st)
    }

    fn make_sel() -> anyhow::Result<Selection> {
        let topst = read_test_pdb();
        let t = topst.0.clone().to_rc();
        let s = topst.1.clone().to_rc();
        let sel = SelectionAll {}.select(&t, &s)?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<Selection> {
        let topst = read_test_pdb();
        let t = topst.0.clone().to_rc();
        let s = topst.1.clone().to_rc();
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
        println!("before {}", sel.iter_pos().next().unwrap());

        let (minv, maxv) = sel.min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}", sel.iter_pos().next().unwrap());
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

        println!("before {}", sel.iter_pos().next().unwrap());
        sel.translate(Vector3f::new(10.0, 10.0, 10.0));
        println!("after {}", sel.iter_pos().next().unwrap());
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

    #[test]
    fn split_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        for res in sel1.split_contig(|_,at,_| at.resid) {
            println!("Res: {}",res.iter_atoms().next().unwrap().resid)
        }
        Ok(())
    }
}
