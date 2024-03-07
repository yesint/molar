use std::{cell::UnsafeCell, rc::Rc};
use crate::io::StateProvider;
use super::{providers::{BoxProvider, PosProvider}, PeriodicBox, Pos};
//use super::handle::{SharedHandle, Handle};

#[derive(Debug, Default,Clone)]
pub struct StateStorage {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub pbox: Option<PeriodicBox>,
}

#[derive(Debug)]
pub struct State(UnsafeCell<StateStorage>);

impl Clone for State {
    fn clone(&self) -> Self {
        Self(UnsafeCell::new(self.get().clone()))
    }
}

impl From<StateStorage> for State {
    fn from(value: StateStorage) -> Self {
        Self(UnsafeCell::new(value))
    }
}

impl State {
    pub fn to_rc(self) -> Rc<Self> {
        Rc::new(self)
    }

    #[inline(always)]
    fn get(&self) -> &StateStorage {
        unsafe {&*self.0.get()}
    }

    #[inline(always)]
    fn get_mut(&self) -> &mut StateStorage {
        unsafe {&mut *self.0.get()}
    }

    #[inline(always)]
    pub unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        self.get().coords.get_unchecked(i)
    }

    #[inline(always)]
    pub unsafe fn nth_pos_unchecked_mut(&self, i: usize) -> &mut Pos {
        self.get().coords.get_unchecked_mut(i)
    }

    #[inline(always)]
    pub fn nth_pos(&self, i: usize) -> Option<&Pos> {
        unsafe { self.get().coords.get(i) }
    }

    #[inline(always)]
    pub fn get_box(&self) -> Option<&PeriodicBox> {
        unsafe { self.get().pbox.as_ref() }
    }

    #[inline(always)]
    pub fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        unsafe { self.get_mut().pbox.as_mut() }
    }

    
}

impl StateProvider for State {
    fn get_time(&self) -> f32 {
        unsafe { self.get().time }
    }

    fn num_coords(&self) -> usize {
        unsafe { self.get().coords.len() }
    }
}

impl PosProvider for State {
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        unsafe { self.get().coords.iter() }
    }
}

impl BoxProvider for State {
    fn get_box(&self) -> Option<&PeriodicBox> {
        unsafe { self.get().pbox.as_ref() }
    }
}