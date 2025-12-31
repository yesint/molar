use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};

use crate::{prelude::*, selection::utils::{difference_sorted, intersection_sorted, union_sorted}};

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

impl IndexSliceProvider for Sel {
    fn get_index_slice(&self) -> &[usize] {
        self.0.as_slice()
    }
}

impl SelectionLogic for Sel {
    type DerivedSel = Sel;
    fn clone_with_index(&self, index: SVec) -> Self::DerivedSel {
        Sel(index)
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

impl Selectable for SelOwnBound<'_> {
/// Create new unbound sub-selection based on provided definition.
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index.as_slice()),
        )?))
    }
}

impl<'a> SelectionLogic for SelOwnBound<'a> {
    type DerivedSel = SelOwnBound<'a>;
    fn clone_with_index(&self, index: SVec) -> Self::DerivedSel {
        SelOwnBound {
            sys: &self.sys,
            index
        }
    }
}

impl SelOwnBound<'_> {
    /// Create new bound sub-selection based on provided definition.
    pub fn select_bound(&self, def: impl SelectionDef) -> Result<Self, SelectionError> {
        Ok(Self {
            sys: self.sys,
            index: def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index.as_slice()))?,
        })
    }

    pub fn into_unbound(self) -> Sel {
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

impl IndexSliceProvider for SelOwnBound<'_> {
    fn get_index_slice(&self) -> &[usize] {
        self.index.as_slice()
    }
}

impl SystemProvider for SelOwnBound<'_> {
    fn get_system_ptr(&self) -> *const System {
        self.sys
    }
}

impl SaveTopology for SelOwnBound<'_> {}
impl SaveState for SelOwnBound<'_> {}
impl SaveTopologyState for SelOwnBound<'_> {}

//================================================
/// Read-write bound subsystem having access to
/// all fields of Topology and State
//================================================
pub struct SelOwnBoundMut<'a> {
    pub(super) sys: &'a mut System,
    pub(super) index: SVec,
}

impl Selectable for SelOwnBoundMut<'_> {
/// Create new unbound sub-selection based on provided definition.
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index.as_slice()),
        )?))
    }
}

impl SelOwnBoundMut<'_> {
    /// Create new sub-selection based on provided definition.
    pub fn into_unbound(self) -> Sel {
        Sel(self.index)
    }

    pub fn into_index(self) -> Sel {
        Sel(self.index)
    }

    pub fn clone_index(&self) -> Sel {
        Sel(self.index.clone())
    }
}

impl IndexSliceProvider for SelOwnBoundMut<'_> {
    fn get_index_slice(&self) -> &[usize] {
        self.index.as_slice()
    }
}

impl SystemProvider for SelOwnBoundMut<'_> {
    fn get_system_ptr(&self) -> *const System {
        self.sys
    }
}

impl SystemMutProvider for SelOwnBoundMut<'_> {}

impl SelectionLogic for SelOwnBoundMut<'_> {
    type DerivedSel = Sel;
    fn clone_with_index(&self, index: SVec) -> Self::DerivedSel {
        Sel(index)
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

impl Selectable for SelBound<'_> {
/// Create new unbound sub-selection based on provided definition.
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index),
        )?))
    }
}

impl<'a> SelectionLogic for SelBound<'a> {
    type DerivedSel = SelOwnBound<'a>;
    fn clone_with_index(&self, index: SVec) -> Self::DerivedSel {
        SelOwnBound {
            sys: &self.sys,
            index
        }
    }
}

impl SelBound<'_> {
    /// Create new owned sub-selection based on provided definition.
    pub fn select_bound(&self, def: impl SelectionDef) -> Result<SelOwnBound<'_>, SelectionError> {
        let index = def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index))?;
        Ok(SelOwnBound {
            index,
            sys: &self.sys,
        })
    }

    pub fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index))?))
    }

    pub fn clone_index(&self) -> Sel {
        Sel(SVec::from_iter(self.index.iter().cloned()))
    }
}

impl IndexSliceProvider for SelBound<'_> {
    fn get_index_slice(&self) -> &[usize] {
        self.index
    }
}

impl SystemProvider for SelBound<'_> {
    fn get_system_ptr(&self) -> *const System {
        self.sys
    }
}

impl SaveTopology for SelBound<'_> {}
impl SaveState for SelBound<'_> {}
impl SaveTopologyState for SelBound<'_> {}

//================================================
/// Read-write selection that borrows its index
//================================================
pub struct SelBoundMut<'a> {
    pub(super) sys: &'a mut System,
    pub(super) index: &'a [usize],
}

impl Selectable for SelBoundMut<'_> {
    /// Create new unbound sub-selection based on provided definition.
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index),
        )?))
    }
}

impl SelBoundMut<'_> {
    pub fn clone_index(&self) -> Sel {
        Sel(SVec::from_iter(self.index.iter().cloned()))
    }
}

impl IndexSliceProvider for SelBoundMut<'_> {
    fn get_index_slice(&self) -> &[usize] {
        self.index
    }
}

impl SystemProvider for SelBoundMut<'_> {
    fn get_system_ptr(&self) -> *const System {
        self.sys
    }
}

impl SystemMutProvider for SelBoundMut<'_> {}

impl SelectionLogic for SelBoundMut<'_> {
    type DerivedSel = Sel;
    fn clone_with_index(&self, index: SVec) -> Self::DerivedSel {
        Sel(index)
    }
}

//-------------------------------------------------------------------
/// Convenience macro for binding several selections ar once
#[macro_export]
macro_rules! bind {
    ($sys:expr, $($sel:ident),+ , $body:block) => {{
        $(let $sel = $sys.bind(&$sel);)+
        $body
    }}
}

/// Convenience macro for binding several selections ar once mutably
#[macro_export]
macro_rules! bind_mut {
    ($sys:expr, $($sel:ident),+ , $body:block) => {{
        $(let $sel = $sys.bind_mut(&$sel);)+
        $body
    }}
}

// Operator overloads for selections
impl std::ops::BitOr for &Sel {
    type Output = Sel;
    fn bitor(self, rhs: Self) -> Self::Output {
        Sel(unsafe{union_sorted(&self.0, &rhs.0)})
    }
}

impl std::ops::BitAnd for &Sel {
    type Output = Sel;
    fn bitand(self, rhs: Self) -> Self::Output {
        let ind = unsafe{intersection_sorted(&self.0, &rhs.0)};
        if ind.is_empty() { panic!("empty selection intersection") }
        Sel(ind)
    }
}

impl std::ops::Sub for &Sel {
    type Output = Sel;
    fn sub(self, rhs: Self) -> Self::Output {
        let ind = unsafe{difference_sorted(&self.0, &rhs.0)};
        if ind.is_empty() { panic!("empty selection difference") }
        Sel(ind)
    }
}

//====================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::iter::ParallelIterator;

    #[test]
    fn test1() -> anyhow::Result<()> {
        let mut h = FileHandler::open("tests/protein.pdb").unwrap();
        let (top, st) = h.read().unwrap();
        let mut sys = System::new(top, st)?;
        let sel1 = sys.select_bound(vec![1, 2, 6, 7])?;

        for at in sel1.iter_atoms() {
            println!("{} {}", at.name, at.resname);
        }

        let ind = sel1.into_unbound();
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

        sys.iter_par_split_mut(&par)
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

        Ok(())
    }

    #[test]
    fn test_pos_par_iter() -> anyhow::Result<()> {
        let sys = System::from_file("tests/albumin.pdb")?;
        let sel = sys.select_bound(vec![0, 1, 2, 3, 4])?;

        // Collect from sequential iterator
        let seq_coords: Vec<_> = sel.iter_pos().collect();

        // Collect from parallel iterator
        let par_coords: Vec<_> = sel.par_iter_pos().collect();

        println!("{:?}", seq_coords);
        println!("{:?}", par_coords);

        // Verify they contain the same elements (order might differ due to parallelism)
        assert_eq!(seq_coords.len(), par_coords.len());
        assert_eq!(seq_coords.len(), 5);

        // Since parallel iterator may reorder, we compare sorted by first coordinate
        let mut seq_sorted = seq_coords.clone();
        let mut par_sorted = par_coords.clone();
        seq_sorted.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap_or(std::cmp::Ordering::Equal));
        par_sorted.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap_or(std::cmp::Ordering::Equal));

        for (s, p) in seq_sorted.iter().zip(par_sorted.iter()) {
            assert_eq!(s[0], p[0]);
            assert_eq!(s[1], p[1]);
            assert_eq!(s[2], p[2]);
        }

        Ok(())
    }

    #[test]
    fn test_atom_par_iter() -> anyhow::Result<()> {
        let sys = System::from_file("tests/albumin.pdb")?;
        let sel = sys.select_bound(vec![0, 1, 2, 3, 4])?;

        // Collect names from sequential iterator
        let seq_names: Vec<String> = sel.iter_atoms().map(|a| a.name.to_string()).collect();

        // Collect names from parallel iterator
        let par_names: Vec<String> = sel.par_iter_atoms().map(|a| a.name.to_string()).collect();

        // Verify they contain the same elements
        assert_eq!(seq_names.len(), par_names.len());
        assert_eq!(seq_names.len(), 5);

        // Sort and compare
        let mut seq_sorted = seq_names.clone();
        let mut par_sorted = par_names.clone();
        seq_sorted.sort();
        par_sorted.sort();

        assert_eq!(seq_sorted, par_sorted);

        Ok(())
    }

    #[test]
    fn test_particle_par_iter() -> anyhow::Result<()> {
        let sys = System::from_file("tests/albumin.pdb")?;
        let sel = sys.select_bound(vec![0, 1, 2, 3, 4])?;

        // Collect from sequential iterator
        let seq_names: Vec<(usize, String)> = sel
            .iter_particle()
            .map(|p| (p.id, p.atom.name.to_string()))
            .collect();

        // Collect from parallel iterator
        let par_names: Vec<(usize, String)> = sel
            .par_iter_particle()
            .map(|p| (p.id, p.atom.name.to_string()))
            .collect();

        // Verify they contain the same elements
        assert_eq!(seq_names.len(), par_names.len());
        assert_eq!(seq_names.len(), 5);

        // Sort and compare
        let mut seq_sorted = seq_names.clone();
        let mut par_sorted = par_names.clone();
        seq_sorted.sort_by_key(|x| x.0);
        par_sorted.sort_by_key(|x| x.0);

        assert_eq!(seq_sorted, par_sorted);

        Ok(())
    }
}
