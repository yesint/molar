use crate::core::{Particle, SelectionExpr};
use anyhow::{bail, Result, anyhow};
use itertools::Itertools;

// Internal interface of selections implementation details
pub(super) mod selection_impl {
    use crate::prelude::*;
        
    pub trait SelectionImpl {
        type Subsel: super::Selection;

        fn get_index_iter(&self) -> impl ExactSizeIterator<Item = usize>;    
        fn get_topology(&self) -> &Topology;
        fn get_state(&self) -> &State;
        fn subselect(&self,index: &Vec<usize>) -> impl super::Selection;
        
        // Types can overload this default implementation to get higher performance
        fn nth_index(&self,i: usize) -> Option<usize> {
            self.get_index_iter().nth(i)
        }

        fn nth_particle(&self, i: usize) -> Option<Particle> {
            self.nth_index(i).map(|ind| {
                Particle {
                    id: ind,
                    atom: unsafe { self.get_topology().nth_atom_unchecked(ind) },
                    pos: unsafe {self.get_state().nth_pos_unchecked(ind) },
                }
            })
        }

        // Helper splitting function generic over returned selections kind
        fn split_gen<RT, F, C, S>(&self, func: F) -> C
        where
            RT: Default + std::hash::Hash + std::cmp::Eq,
            F: Fn(Particle) -> RT,
            C: FromIterator<S> + Default,
        {
            let mut ids = std::collections::HashMap::<RT, Vec<usize>>::default();

            for ind in self.get_index_iter() {
                let p = self.nth_particle(ind).unwrap();
                let i = p.id;
                let id = func(p);
                if let Some(el) = ids.get_mut(&id) {
                    el.push(i);
                } else {
                    ids.insert(id, vec![i]);
                }
            }

            C::from_iter(
                ids.into_values()
                    .map(|ind| self.subselect(&ind) ),
                // This should never fail because `ind` can't be empty
            )
        }
    }
}

// Public interface of selections
pub trait Selection: selection_impl::SelectionImpl + Sized {
    fn len(&self) -> usize {
        self.get_index_iter().len()
    }
    
    //===================
    // Subselections
    //===================

    /// Subselection from expression
    fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Self::Subsel> {
        let index = expr.apply_subset(self.get_topology(), self.get_state(), self.get_index_iter())?;
        if index.len() > 0 {
            Ok(self.subselect(&index))
        } else {
            bail!("Selection is empty")
        }
    }

    /// Subselection from string
    fn subsel_from_str(&self, sel_str: &str) -> Result<Self::Subsel> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    fn subsel_from_local_range(
        &self,
        range: std::ops::Range<usize>,
    ) -> Result<Self::Subsel> {
        if range.end >= self.len() {
            bail!(
                "Invalid local sub-range: {}:{}, valid range: 0:{}",
                range.start,
                range.end,
                self.len() - 1
            );
        }

        // Translate range of local indexes to global indexes
        let index: Vec<usize> = self
            .get_index_iter()
            .skip(range.start)
            .take(range.len())
            .collect();

        Ok(self.subselect(&index))
    }

    /// Subselection from iterator that provides local selection indexes
    fn subsel_from_iter(
        &self,
        iter: impl ExactSizeIterator<Item = usize>,
    ) -> Result<Self::Subsel> {
        // Remove duplicates if any
        let index = iter
            .sorted()
            .dedup()
            .map(|i| {
                self.nth_index(i).ok_or_else(|| {
                    anyhow!("Index {} is out of allowed range [0:{}]",i,self.len())
                })
            })
            .collect::<Result<Vec<usize>>>()?;
        // Now it's safe to call
        Ok(self.subselect(&index))
    }

    fn iter_particles(&self) -> SelectionParticleIterator<Self> {
        SelectionParticleIterator { sel: self, cur: 0 }
    }

    //============================
    // Splitting
    //============================

    

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is kept intact.
    /// 
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as subselections.
    fn split<RT, F, C>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> RT,
        C: FromIterator<Self::Subsel> + Default,
    {
        split_gen(self,func)
    }

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is consumed.
    /// 
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as parent selection.
    fn split_into<RT, F, C>(self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> RT,
        C: FromIterator<Self> + Default,
    {
        split_gen(&self,func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is left alive.
    fn split_resid<C>(&self) -> C
    where
        C: FromIterator<Self::Subsel> + Default,
    {
        self.split_gen(|p| p.atom.resid)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is consumed.
    fn split_resid_into<C>(self) -> C
    where
        C: FromIterator<Self> + Default,
    {
        self.split_gen(|p| p.atom.resid)
    }


}

//============================
// Iterators
//============================

/// Iterator over the [Particle]s from selection.
pub struct SelectionParticleIterator<'a, Sel> {
    sel: &'a Sel,
    cur: usize,
}

impl<'a, Sel: Selection> Iterator for SelectionParticleIterator<'a, Sel> {
    type Item = Particle<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.sel.len() {
            let ret = self.sel.nth_particle(self.cur);
            self.cur += 1;
            ret
        } else {
            None
        }
    }
}

impl<'a, Sel: Selection> ExactSizeIterator for SelectionParticleIterator<'a, Sel> {
    fn len(&self) -> usize {
        self.sel.len()
    }
}
