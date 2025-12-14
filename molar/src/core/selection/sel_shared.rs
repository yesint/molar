use std::{cell::{self, Ref, RefCell, RefMut}, path::Path, rc::Rc};

use nalgebra::{Const, IsometryMatrix3, Unit};
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use thiserror::Error;

use crate::prelude::*;

pub struct SystemInner {
    top: Topology,
    st: State,
}

pub type SystemRc = Rc<RefCell<SystemInner>>;

pub struct SystemShared(SystemRc);

impl SystemShared {
    pub fn from_file(fname: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        let inner = SystemInner {top,st};
        Ok(Self(Rc::new(RefCell::new(inner))))
    }

    pub fn select(&self, def: impl SelectionDef) -> Result<SelShared, SelectionError> {
        let sb = self.0.borrow();
        Ok(SelShared{
            sys: Rc::clone(&self.0),
            index: def.into_sel_index(&sb.top, &sb.st, None)?,
        })
    }

    pub fn select_all(&self) -> SelShared {
        let sb = self.0.borrow();
        SelShared{
            sys: Rc::clone(&self.0),
            index: unsafe {SVec::from_sorted( (0..sb.top.len()).into_iter().collect())},
        }
        
    }
}

pub struct SelShared {
    sys: SystemRc,
    index: SVec,
}

pub struct SelRef<'a> {
    sys: Ref<'a,SystemInner>,
    index: &'a SVec,
}

impl LenProvider for SelRef<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

pub struct SelRefMut<'a> {
    sys: RefMut<'a,SystemInner>,
    index: &'a SVec,
}

impl SelShared {
    pub fn borrow(&self) -> SelRef<'_> {
        SelRef {
            sys: self.sys.borrow(),
            index: &self.index,
        }
    }

    pub fn borrow_mut(&self) -> SelRefMut<'_> {
        SelRefMut {
            sys: self.sys.borrow_mut(),
            index: &self.index,
        }
    }

    pub fn with<'a,T,R>(&'a self, op: T) -> R 
    where 
        T: Fn(SelRef<'a>)->R,
    {
        op(self.borrow())
    }

    pub fn with_mut<'a,T>(&'a mut self, mut op: T) -> Result<(), SelectionError>
    where 
        T: FnMut(SelRefMut<'a>),
    {
        op(self.borrow_mut());
        Ok(())
    }
}

#[macro_export]
macro_rules! with {
    ($($sel:ident),+ , $body:block) => {{
        $(let $sel = $sel.borrow();)+
        $body
    }}
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() -> anyhow::Result<()> {
        let sys = SystemShared::from_file("tests/albumin.pdb")?;
        let ca = sys.select("name CA")?;
        let cb = sys.select("name CB")?;
        let n = ca.with(|s|{
            println!("{}", s.len());
            s.len()
        });

        let n1 = with!(ca,cb,{
            println!("{} {}", ca.len(), cb.len());
            ca.len()
        });

        Ok(())
    }
}