use std::iter::zip;

use super::providers::*;
use super::Atom;
use super::Matrix3f;
use super::Selection;
use super::{PbcDims, Pos, Vector3f};
use crate::distance_search::search::{DistanceSearcherSingle, SearchConnectivity};
use anyhow::{bail, Result, anyhow};
use itertools::izip;
use nalgebra::Rotation3;
use nalgebra::Unit;
use nalgebra::SVD;
use num_traits::Bounded;
use num_traits::Zero;

//---------------------------------------------------
// Free functions for computing properties that
// acceps a needed type of data provider as argument
//---------------------------------------------------


//-------------------------------------------------------
// Splitting iterator
//-------------------------------------------------------

pub struct SelectionSplitIterator<'a,T,F> {
    sel: &'a Selection,
    func: F,
    counter: usize,
    id: T,
}

impl SelectionSplitIterator<'_,(),()>
{
    pub fn new<T,F>(sel: &Selection, func: F) -> SelectionSplitIterator<'_,T,F> 
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

impl<T,F> Iterator for SelectionSplitIterator<'_,T,F> 
where 
    T: Default + std::cmp::PartialEq,
    F: Fn(usize, &Atom, &Pos) -> T,
{
    type Item = Selection;
    fn next(&mut self) -> Option<Self::Item> {
        let mut index = vec![];
        while self.counter < self.sel.len() {
            let (i,at,pos) = unsafe{ self.sel.nth_unchecked(self.counter) };
            let id = (self.func)(i,at,pos);

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
                return unsafe{ Some(self.sel.subsel_from_vec_unchecked(index).unwrap()) };
            }
            // Next element
            self.counter+=1;
        };

        // Return any remaining index as last selection
        if !index.is_empty() {
            return unsafe{ Some(self.sel.subsel_from_vec_unchecked(index).unwrap()) };
        }
        
        // If we are here stop iterating
        None
    }
}

