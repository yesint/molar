use std::marker::PhantomData;

use super::{
    selection_parser::{apply_ast_whole, generate_ast, SelectionAst},
    Atom, AtomIterator, AtomIteratorMut, IdAtomIterator, IdAtomIteratorMut, IdAtomPosIterator,
    IdAtomPosIteratorMut, IdPosIterator, IdPosIteratorMut, IndexIterator, Pos, PosIterator,
    PosIteratorMut, State, Structure,
};
use anyhow::{bail, Result};
use nalgebra::Point3;

#[derive(Debug)]
struct Selection<'a, I: IndexIterator> {
    ast: Option<SelectionAst>,
    index_iter: I,
    phantom: PhantomData<&'a I>,
}

impl<'a, I: IndexIterator> Selection<'a, I> {
    /// Creates selection from any iterator over usize
    fn from_iter(iter: I) -> Self {
        Self {
            index_iter: iter,
            ast: None,
            phantom: PhantomData,
        }
    }

    fn iter_pos(&self, state: &'a State) -> impl PosIterator<'a> + '_ {
        self.index_iter.clone().map(|i| &state.coords[i])
    }

    fn iter_pos_mut(&self, state: &'a mut State) -> impl PosIteratorMut<'a> + '_ {
        SubscriptedMutIter {
            iter: state.coords.iter_mut(),
            index_iter: self.index_iter.clone(),
            cur: 0,
        }
    }

    fn iter_index_pos(&self, state: &'a State) -> impl IdPosIterator<'a> + '_ {
        self.index_iter.clone().map(|i| (i, &state.coords[i]))
    }

    fn iter_index_pos_mut(&self, state: &'a mut State) -> impl IdPosIteratorMut<'a> + '_ {
        let it = SubscriptedMutIter {
            iter: state.coords.iter_mut(),
            index_iter: self.index_iter.clone(),
            cur: 0,
        };
        core::iter::zip(self.index_iter.clone(), it)
    }

    fn iter_atoms(&self, structure: &'a Structure) -> impl AtomIterator<'a> + '_ {
        self.index_iter.clone().map(|i| &structure.atoms[i])
    }

    fn iter_atoms_mut(&self, structure: &'a mut Structure) -> impl AtomIteratorMut + '_ {
        SubscriptedMutIter {
            iter: structure.atoms.iter_mut(),
            index_iter: self.index_iter.clone(),
            cur: 0,
        }
    }

    fn iter_index_atoms(&self, structure: &'a Structure) -> impl IdAtomIterator<'a> + '_ {
        self.index_iter.clone().map(|i| (i, &structure.atoms[i]))
    }

    fn iter_index_atoms_mut(&self, structure: &'a mut Structure) -> impl IdAtomIteratorMut + '_ {
        let it = SubscriptedMutIter {
            iter: structure.atoms.iter_mut(),
            index_iter: self.index_iter.clone(),
            cur: 0,
        };
        core::iter::zip(self.index_iter.clone(), it)
    }

    fn iter_index_atoms_pos(
        &self,
        structure: &'a Structure,
        state: &'a State,
    ) -> impl IdAtomPosIterator<'a> + '_ {
        self.index_iter
            .clone()
            .map(|i| (i, &structure.atoms[i], &state.coords[i]))
    }

    fn iter_index_atoms_pos_mut(
        &self,
        structure: &'a mut Structure,
        state: &'a mut State,
    ) -> impl IdAtomPosIteratorMut + '_ {
        let at = SubscriptedMutIter {
            iter: structure.atoms.iter_mut(),
            index_iter: self.index_iter.clone(),
            cur: 0,
        };
        let pos = SubscriptedMutIter {
            iter: state.coords.iter_mut(),
            index_iter: self.index_iter.clone(),
            cur: 0,
        };
        itertools::izip!(self.index_iter.clone(), at, pos)
    }
}

impl Selection<'_, <Vec<usize> as IntoIterator>::IntoIter> {
    fn from_sel_str(sel_str: &str, structure: &Structure, state: &State) -> Result<Self> {
        let ast = generate_ast(sel_str)?;
        let it = apply_ast_whole(&ast, structure, state)?;
        Ok(Self {
            index_iter: it.into_iter(),
            ast: Some(ast),
            phantom: PhantomData,
        })
    }

    fn apply(&mut self, structure: &Structure, state: &State) -> Result<()> {
        if let Some(ast) = self.ast.as_ref() {
            self.index_iter = apply_ast_whole(&ast, structure, state)?.into_iter()
        }
        Ok(())
    }
}

// Helper struct for creating subscripted mutable iterator
// from container iterator and index iterator
// IMPORTANT! Only works for **sorted** container iterator!
struct SubscriptedMutIter<'a, CI, T, II>
where
    T: 'a,                          // Element type
    CI: Iterator<Item = &'a mut T>, // Container iterator
    II: IndexIterator,              // Index iterator
{
    iter: CI,
    index_iter: II,
    cur: usize,
}

impl<'a, CI, T, II> Iterator for SubscriptedMutIter<'a, CI, T, II>
where
    T: 'a,
    CI: Iterator<Item = &'a mut T>,
    II: IndexIterator,
{
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        match self.index_iter.next() {
            Some(i) => {
                // Advance pos_iter by offset and yield
                let el = self.iter.nth(i - self.cur);
                // Advance current position
                self.cur = i + 1;
                el
            }
            None => None,
        }
    }
}

impl<'a, CI, T, II> ExactSizeIterator for SubscriptedMutIter<'a, CI, T, II>
where
    T: 'a,
    CI: Iterator<Item = &'a mut T>,
    II: IndexIterator,
{
    fn len(&self) -> usize {
        self.index_iter.len()
    }
}
