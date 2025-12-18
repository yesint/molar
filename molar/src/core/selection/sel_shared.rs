use std::{cell::{self, Ref, RefCell, RefMut, UnsafeCell}, marker::PhantomData, path::Path, rc::Rc};

use nalgebra::{Const, IsometryMatrix3, Unit};
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use thiserror::Error;

use crate::prelude::*;

pub type SystemRc = Rc<UnsafeCell<System>>;

pub struct SystemShared(SystemRc);

impl SystemShared {
    pub fn from_file(fname: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let inner = System::from_file(fname)?;
        Ok(Self(Rc::new(UnsafeCell::new(inner))))
    }

    pub fn select(&self, def: impl SelectionDef) -> Result<SelShared, SelectionError> {
        let sb = unsafe{&*self.0.get()};
        Ok(SelShared{
            sys: Rc::clone(&self.0),
            index: def.into_sel_index(&sb.top, &sb.st, None)?,
            phantom_: Default::default(),
        })
    }

    pub fn select_all(&self) -> SelShared {
        let sb = unsafe{&*self.0.get()};
        SelShared{
            sys: Rc::clone(&self.0),
            index: unsafe {SVec::from_sorted( (0..sb.top.len()).into_iter().collect())},
            phantom_: Default::default(),
        }
        
    }
}

pub struct SelShared {
    sys: SystemRc,
    index: SVec,
    phantom_: PhantomData<*const ()>, // no Send
}

//------------------------
pub struct SelRef<'a> {
    sys: *const System,
    index: &'a SVec,
}

impl LenProvider for SelRef<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelRef<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelRef<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        unsafe{&*self.sys}.atoms_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        unsafe{&*self.sys}.coords_ptr()
    }
}

impl MeasurePosGuarded for SelShared {}

//------------------
pub struct SelRefMut<'a> {
    sys: *mut System,
    index: &'a SVec,
}

impl LenProvider for SelRefMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelRefMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelRefMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        unsafe{&*self.sys}.atoms_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        unsafe{&*self.sys}.coords_ptr()
    }
}

impl AtomPosAnalysisMut for SelRefMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        unsafe{&mut *self.sys}.atoms_ptr_mut()
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        unsafe{&mut *self.sys}.coords_ptr_mut()
    }
}

//-------------------------------


impl SelShared {
    // Since SelRef is only created inside the with_ref() closure and
    // SelShared is not Send only one of them exists  at any given time. 
    // SelRef stores raw pointer to System, which is only
    // turned to references to atoms and coords at the time of execution
    // so multiple &mut System never exist and the code is sound.
    fn get_ref(&self) -> SelRef<'_> {
        SelRef {
            sys: self.sys.get(),
            index: &self.index,
        }
    }

    fn get_ref_mut(&mut self) -> SelRefMut<'_> {
        SelRefMut {
            sys: self.sys.get(),
            index: &self.index,
        }
    }

    /// Used to execute any methods that provide references
    /// to atoms and coordinates by means of SelRef that exists inside 
    /// the closure and can't escape by any means.
    pub fn with_ref<'a,T,R>(&'a self, op: T) -> R
    where 
        T: Fn(SelRef<'a>) -> R,
    {
        op(self.get_ref())
    }

    pub fn with_ref_mut<'a,T,R>(&'a mut self, mut op: T) -> R
    where 
        T: FnMut(SelRefMut<'a>) -> R,
    {
        op(self.get_ref_mut())
    }
}

fn guard() {

}
// #[macro_export]
// macro_rules! with {
//     ($($sel:ident),+ , $body:block) => {{
//         $(let $sel = $sel.get_ref();)+
//         $body
//     }}
// }

// #[macro_export]
// macro_rules! with_mut {
//     ($($sel:ident),+ , $body:block) => {{
//         $(let $sel = $sel.borrow_mut();)+
//         $body
//     }}
// }

impl Guarded for SelShared {
    type Guard<'a> = SelRef<'a> where Self: 'a;
    fn guard(&self) -> Self::Guard<'_> {
        self.get_ref()
    }
}

impl GuardedMut for SelShared {
    type GuardMut<'a> = SelRefMut<'a> where Self: 'a;

    fn guard_mut(&mut self) -> Self::GuardMut<'_> {
        self.get_ref_mut()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() -> anyhow::Result<()> {
        let sys = SystemShared::from_file("tests/albumin.pdb")?;
        let ca = sys.select("name CA")?;
        let mut cb = sys.select("name CB")?;
        
        ca.with_ref(|s|{
            println!("{}", s.center_of_geometry());
            for p in s.iter_pos(){
                println!("{}",p);
            }
        });

        cb.with_ref_mut(|mut s|{
            println!("{}", s.center_of_geometry());
            for p in s.iter_pos_mut(){
                *p *= 2.0;
                println!("{}",p);
            }
        });

        println!("{}", ca.center_of_geometry());

        Ok(())
    }
}