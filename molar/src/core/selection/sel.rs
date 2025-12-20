use crate::prelude::*;

//===========================================================================
/// Selection index detached from any [System].
/// It is guaranteed to be non-empty.
//===========================================================================
pub struct Sel(pub(crate) SVec);

impl Sel {
    pub fn from_vec(index: Vec<usize>) -> Result<Self, SelectionError> {
        if index.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {
            Ok(Self(SVec::from_unsorted(index)))
        }
    }

    pub fn into_svec(self) -> SVec {
        self.0
    }

    pub fn from_svec(index: SVec) -> Result<Self, SelectionError> {
        if index.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {
            Ok(Self(index))
        }
    }

    pub(super) fn from_iter(iter: impl Iterator<Item = usize>) -> Result<Self, SelectionError> {
        let index: Vec<_> = iter.collect();
        Ok(Self::from_vec(index)?)
    }
}

impl LenProvider for Sel {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl IndexProvider for Sel {
    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.0.iter().cloned()
    }

    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.0.get_unchecked(i)
    }
}

//================================================
/// Read only subsystem that owns its index
/// Implements only read-only analysis traits
//================================================

#[derive(Clone, Debug)]
pub struct SelOwnBound<'a> {
    pub(super) sys: &'a System,
    pub(super) index: SVec,
}

impl SelOwnBound<'_> {
    /// Create new unbound sub-selection based on provided definition.
    pub fn select_as_index(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index.as_slice()),
        )?))
    }

    /// Create new bound sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Self, SelectionError> {
        Ok(Self {
            sys: self.sys,
            index: def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index.as_slice()))?,
        })
    }

    pub fn into_index(self) -> Sel {
        Sel(self.index)
    }

    pub fn clone_index(&self) -> Sel {
        Sel(self.index.clone())
    }

    pub fn split_bound<RT, F>(&self, func: F) -> impl Iterator<Item = SelOwnBound<'_>>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        self.split(func).map(|sel| SelOwnBound {
            index: sel.0,
            sys: &self.sys,
        })
    }

    pub fn split_resindex_bound(&self) -> impl Iterator<Item = SelOwnBound<'_>> {
        self.split_resindex().map(|sel| SelOwnBound {
            index: sel.0,
            sys: &self.sys,
        })
    }
}

impl LenProvider for SelOwnBound<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelOwnBound<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelOwnBound<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SelOwnBound<'_> {
    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }

    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }
}

impl TopologyWrite for SelOwnBound<'_> {}
impl StateWrite for SelOwnBound<'_> {}
impl TopologyStateWrite for SelOwnBound<'_> {}

//================================================
/// Read-write bound subsystem having access to
/// all fields of Topology and State
//================================================
pub struct SelOwnBoundMut<'a> {
    pub(super) sys: &'a mut System,
    pub(super) index: SVec,
}

impl SelOwnBoundMut<'_> {
    /// Create new sub-selection based on provided definition.
    pub fn unbind(self) -> Sel {
        Sel(self.index)
    }

    pub fn into_index(self) -> Sel {
        Sel(self.index)
    }

    pub fn clone_index(&self) -> Sel {
        Sel(self.index.clone())
    }
}

impl LenProvider for SelOwnBoundMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelOwnBoundMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelOwnBoundMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SelOwnBoundMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.sys.top.atoms.as_mut_ptr()
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.sys.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SelOwnBoundMut<'_> {
    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }

    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }
}

impl NonAtomPosAnalysisMut for SelOwnBoundMut<'_> {
    fn st_ptr_mut(&mut self) -> *mut State {
        &mut self.sys.st
    }

    fn top_ptr_mut(&mut self) -> *mut Topology {
        &mut self.sys.top
    }
}

//================================================
/// Read only selection that borrows its index
/// Implements only read-only analysis traits
//================================================

#[derive(Clone, Debug)]
pub struct SelBound<'a> {
    pub(super) sys: &'a System,
    pub(super) index: &'a [usize],
}

impl SelBound<'_> {
    /// Create new owned sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<SelOwnBound<'_>, SelectionError> {
        let index = def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index))?;
        Ok(SelOwnBound {
            index,
            sys: &self.sys,
        })
    }

    pub fn select_as_index(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index))?))
    }

    pub fn clone_index(&self) -> Sel {
        Sel(SVec::from_iter(self.index.iter().cloned()))
    }
}

impl LenProvider for SelBound<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelBound<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelBound<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SelBound<'_> {
    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }

    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }
}

impl TopologyWrite for SelBound<'_> {}
impl StateWrite for SelBound<'_> {}
impl TopologyStateWrite for SelBound<'_> {}

//================================================
/// Read-write selection that borrows its index
//================================================
pub struct SelBoundMut<'a> {
    pub(super) sys: &'a mut System,
    pub(super) index: &'a [usize],
}

impl SelBoundMut<'_> {
    /// Create new sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<SelOwnBound<'_>, SelectionError> {
        let index = def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index))?;
        Ok(SelOwnBound {
            index,
            sys: &self.sys,
        })
    }

    pub fn clone_index(&self) -> Sel {
        Sel(SVec::from_iter(self.index.iter().cloned()))
    }
}

impl LenProvider for SelBoundMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelBoundMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelBoundMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SelBoundMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.sys.top.atoms.as_mut_ptr()
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.sys.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SelBoundMut<'_> {
    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }

    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }
}

impl NonAtomPosAnalysisMut for SelBoundMut<'_> {
    fn st_ptr_mut(&mut self) -> *mut State {
        &mut self.sys.st
    }

    fn top_ptr_mut(&mut self) -> *mut Topology {
        &mut self.sys.top
    }
}

#[macro_export]
macro_rules! bind {
    ($sys:expr, $($sel:ident),+ , $body:block) => {{
        $(let $sel = $sys.bind(&$sel);)+
        $body
    }}
}

#[macro_export]
macro_rules! bind_mut {
    ($sys:expr, $($sel:ident),+ , $body:block) => {{
        $(let $sel = $sys.bind_mut(&$sel);)+
        $body
    }}
}

//====================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    //use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};

    #[test]
    fn test1() -> anyhow::Result<()> {
        let mut h = FileHandler::open("tests/protein.pdb").unwrap();
        let (top, st) = h.read().unwrap();
        let mut sys = System::new(top, st)?;
        let sel1 = sys.select_bound(vec![1, 2, 6, 7])?;

        for at in sel1.iter_atoms() {
            println!("{} {}", at.name, at.resname);
        }

        let ind = sel1.into_index();
        for at in sys.try_bind_mut(&ind)?.iter_atoms_mut() {
            at.bfactor += 1.0;
        }

        // let lock2 = sel1.bind(&sys);
        // for at in lock2.iter_atoms() {
        //     println!("{} ", at.bfactor);
        // }

        //drop(lock);

        Ok(())
    }

    #[test]
    fn test1_1() -> anyhow::Result<()> {
        let mut sys = System::from_file("tests/albumin.pdb")?;

        let par = sys.split_par(|p| {
            if p.atom.resname != "SOL" {
                Some(p.atom.resindex)
            } else {
                None
            }
        })?;

        sys.bind_par_mut(&par)?
            .par_iter_mut()
            .try_for_each(|mut sel| {
                //println!("{} {}", sel.len(), sel.first_atom().resname);
                if sel.first_atom().resname == "ALA" {
                    sel.first_pos_mut().coords[1] = 555.5;
                    println!("{}", sel.first_pos());
                } else {
                    sel.first_pos_mut().coords[1] = 777.7;
                }
                Ok::<_, SelectionError>(())
            })?;

        // Add serial selection
        let ca = sys.select_bound("name CA")?;
        let cb = sys.select_bound("name CB")?;
        println!("#ca: {} {}", ca.len(), cb.center_of_mass()?);

        for a in cb.iter_atoms().take(10) {
            println!("{}", a.name);
        }

        //Iter serial views
        let mut b = sys.bind_par_mut(&par)?;
        let serials: Vec<_> = b.iter_mut().collect();
        println!("serial #5: {}", serials[5].first_particle().pos);

        Ok(())
    }
}
