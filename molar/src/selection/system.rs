use crate::par::*;

use super::utils::check_topology_state_sizes;
use crate::prelude::*;
use std::path::Path;

//================================================
/// System that stores Topology and State
//================================================
#[derive(Debug, Default)]
pub struct System {
    pub(super) top: Topology,
    pub(super) st: State,
}

impl System {
    pub fn new(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        Ok(Self { top, st })
    }

    pub fn from_file(fname: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Self::new(top, st)?)
    }

    /// Create unbound sub-selection
    pub fn sub_select(&self, ind: &impl IndexSliceProvider, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        unsafe {Ok(self.try_bind_sorted_slice(ind.get_index_slice())?.select(def)?)}
    }

    /// Create all detached
    pub fn select_all(&self) -> Sel {
        Sel(unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) })
    }

    pub fn select_all_bound(&self) -> SelOwnBound<'_> {
        SelOwnBound {
            sys: self,
            index: unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) },
        }
    }

    pub fn select_bound_mut(
        &mut self,
        def: impl SelectionDef,
    ) -> Result<SelOwnBoundMut<'_>, SelectionError> {
        let index = def.into_sel_index(self, None)?;
        Ok(SelOwnBoundMut { sys: self, index })
    }

    pub fn select_all_bound_mut(&mut self) -> SelOwnBoundMut<'_> {
        let index = unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) };
        SelOwnBoundMut { sys: self, index }
    }

    // Internal function. Slice is supposed to be sorted!
    unsafe fn try_bind_sorted_slice<'a>(&'a self, sel: &'a [usize]) -> Result<SelBound<'a>, SelectionError> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.get_unchecked(sel.len() - 1) };
        if last >= self.top.len() {
            Err(SelectionError::IndexValidation(
                *sel.first().unwrap(),
                last,
                self.top.len() - 1,
            ))
        } else {
            Ok(SelBound {
                sys: self,
                index: sel,
            })
        }
    }

    /// Binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn bind<'a>(&'a self, sel: &'a Sel) -> SelBound<'a> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            panic!("selection is out of bounds");
        } else {
            SelBound {
                sys: self,
                index: sel.0.as_slice(),
            }
        }
    }

    /// Borrow the current state (coordinates, box, time).
    pub fn state(&self) -> &State {
        &self.st
    }

    /// Borrow the topology.
    pub fn topology(&self) -> &Topology {
        &self.top
    }

    /// Replace the topology's bond table (see [`Topology::set_bonds`]).
    ///
    /// Bonds are the one part of a topology that can be swapped on a live system: they take
    /// no part in the topology/state size invariant, so no existing selection is
    /// invalidated. This is how connectivity that did **not** come from the structure file —
    /// distance-based perception, an interactive editor, a separately read topology — reaches
    /// everything that reads bonds: the `polh` / `apolh` selection keywords,
    /// [`perceive`](Self::perceive), and `molar_ff`'s force-field typing / charge assignment.
    pub fn set_bonds(&mut self, bonds: BondStorage) -> Result<(), SelectionError> {
        Ok(self.top.set_bonds(bonds)?)
    }

    /// Perceive connectivity from the current elements and coordinates.
    ///
    /// This replaces the complete bond table. New bonds have unspecified order.
    pub fn perceive_connectivity(
        &mut self,
        options: &ConnectivityOptions,
    ) -> Result<usize, BondPerceptionError> {
        let bonds = crate::perception::perceive_connectivity(&*self, options)?;
        let count = bonds.len();
        // The free function returns indices local to its input. `System` has an
        // identity index, so every endpoint was checked by construction.
        self.top.bonds = bonds;
        Ok(count)
    }

    /// Perceive bond orders and formal charges for the current topology, returning a validated
    /// [`BondAssignment`]. This does not change the system; apply the result with
    /// [`System::apply_bond_assignment`].
    pub fn assign_bond_orders(
        &self,
        options: &BondOrderOptions,
    ) -> Result<BondAssignment, BondPerceptionError> {
        crate::perception::assign_bond_orders(self.topology(), Some(&self.st.coords), options)
    }

    /// Apply a bond-order and formal-charge result to this system's topology.
    pub fn apply_bond_assignment(
        &mut self,
        assignment: &BondAssignment,
    ) -> Result<(), BondPerceptionError> {
        assignment.apply_to(&mut self.top)
    }

    /// Add the explicit hydrogens of a [`HydrogenAddition`] plan as one transaction: the new
    /// atoms, their coordinates, and their bonds are appended together. New atoms go at the end,
    /// so existing atom indices — and any selection built on them — stay valid. All validation
    /// happens before the first change, so a rejected plan leaves the system untouched.
    pub fn add_hydrogens(
        &mut self,
        plan: &HydrogenAddition,
    ) -> Result<(), BondPerceptionError> {
        let n = self.top.atoms.len();
        if n != plan.source_atom_count() {
            return Err(BondPerceptionError::AtomCountChanged {
                expected: plan.source_atom_count(),
                actual: n,
            });
        }
        for &parent in plan.parents() {
            if parent >= n {
                return Err(BondPerceptionError::NoValidAssignment { atom: parent });
            }
        }
        let has_vel = self.st.has_vel();
        let has_force = self.st.has_force();
        if !plan.zero_fill_dynamics() && (has_vel || has_force) {
            return Err(BondPerceptionError::DynamicsPresent);
        }

        // Validation passed; append atoms, coordinates, and bonds together.
        for (i, (&parent, &pos)) in plan.parents().iter().zip(plan.positions()).enumerate() {
            let src = self.top.atoms.get(parent).unwrap();
            let hydrogen = Atom::new()
                .with_atomic_number(1)
                .with_mass(1.008)
                .with_name("H")
                .with_resname(src.get_resname())
                .with_resid(src.get_resid() as i32)
                .with_resindex(src.get_resindex())
                .with_chain(src.get_chain());
            self.top.atoms.push(&hydrogen);
            self.st.coords.push(pos);
            if has_vel {
                self.st.velocities.push(Vel::zeros());
            }
            if has_force {
                self.st.forces.push(Force::zeros());
            }
            self.top.bonds.push(&Bond::with_order(parent, n + i, BondOrder::Single));
        }
        // A single molecule range grows to cover the new atoms; leave a multi-molecule table
        // (only TPR sets one) untouched.
        if let [range] = self.top.molecules.as_mut_slice() {
            range[1] = self.top.atoms.len();
        }
        self.top.bonds.invalidate_adjacency();
        Ok(())
    }

    /// Perceive rings + aromaticity, annotating this system's topology in place: sets
    /// `BondOrder::Aromatic` on aromatic-ring bonds and the in-ring/aromatic flag bits on
    /// the atoms (see [`crate::perception`]). Returns the [`Perception`] (SSSR rings + net
    /// charge). Coordinate-free — works on the connection table alone.
    pub fn perceive(&mut self) -> Perception {
        crate::perception::perceive(&mut self.top)
    }

    /// Bind `sel` using this system's **topology** but coordinates from an
    /// **external** `state` (e.g. a trajectory frame), without copying the state
    /// into the system. The disjoint counterpart of [`bind`](Self::bind), which
    /// uses the system's own state. The selection's indices must be valid for the
    /// topology (guaranteed when `sel` was made from this system); `state` must
    /// have the same atom count.
    pub fn bind_with_state<'a>(&'a self, sel: &'a Sel, state: &'a State) -> SelBoundParts<'a> {
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            panic!("selection is out of bounds");
        }
        SelBoundParts {
            top: &self.top,
            st: state,
            index: sel.0.as_slice(),
        }
    }

    /// Binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn try_bind<'a>(&'a self, sel: &'a Sel) -> Result<SelBound<'a>, SelectionError> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            Err(SelectionError::IndexValidation(
                *sel.0.first().unwrap(),
                last,
                self.top.len() - 1,
            ))
        } else {
            Ok(SelBound {
                sys: self,
                index: sel.0.as_slice(),
            })
        }
    }

    /// Mutably binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn bind_mut<'a>(&'a mut self, sel: &'a Sel) -> SelBoundMut<'a> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            panic!("selection is out of bounds");
        } else {
            SelBoundMut {
                sys: self,
                index: sel.0.as_slice(),
            }
        }
    }

    /// Mutably binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn try_bind_mut<'a>(&'a mut self, sel: &'a Sel) -> Result<SelBoundMut<'a>, SelectionError> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            Err(SelectionError::IndexValidation(
                *sel.0.first().unwrap(),
                last,
                self.top.len() - 1,
            ))
        } else {
            Ok(SelBoundMut {
                sys: self,
                index: sel.0.as_slice(),
            })
        }
    }

    /// Returns mutable parallel iterator over parallel selections.
    pub fn iter_par_split_mut<'a>(
        &'a mut self,
        par: &'a ParSplit,
    ) -> impl IndexedParallelIterator<Item = SelParMut<'a>> {
        par.check_bounds(self);
        let ptr = self as *mut System as usize;
        par.selections
            .par_iter()
            .map(move |sel| SelParMut::new(ptr as *mut System, &sel.0))
    }

    /// Returns parallel iterator over parallel selections.
    pub fn iter_par_split<'a>(
        &'a mut self,
        par: &'a ParSplit,
    ) -> impl IndexedParallelIterator<Item = SelPar<'a>> {
        par.check_bounds(self);
        par.selections
            .par_iter()
            .map(|sel| SelPar::new(self, &sel.0))
    }

    pub fn split_bound<'a, RT, F>(&'a self, func: F) -> impl Iterator<Item = SelOwnBound<'a>>
    where
        RT: Default + std::cmp::PartialEq + 'a,
        F: Fn(Particle) -> Option<RT> + 'a,
    {
        self.split(func).map(|sel| SelOwnBound {
            sys: self,
            index: sel.0,
        })
    }

    pub fn split_resindex_bound(&self) -> impl Iterator<Item = SelOwnBound<'_>> {
        self.split_bound(|p| Some(p.atom.get_resindex()))
    }

    pub fn set_state(&mut self, st: State) -> Result<State, SelectionError> {
        if self.len() != st.len() {
            Err(SelectionError::IncompatibleState)
        } else {
            Ok(std::mem::replace(&mut self.st, st))
        }
    }

    pub fn set_topology(&mut self, top: Topology) -> Result<Topology, SelectionError> {
        if self.len() != top.len() {
            Err(SelectionError::IncompatibleTopology)
        } else {
            Ok(std::mem::replace(&mut self.top, top))
        }
    }

    pub fn release(self) -> (Topology, State) {
        (self.top, self.st)
    }

    //===============
    // Modifying
    //===============

    /// Append selection derived from self
    pub fn append_from_self(&mut self, sel: &Sel) -> Result<Sel, SelectionError> {
        let old_last = self.len() - 1;
        let last_ind = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last_ind >= self.top.len() {
            return Err(SelectionError::IndexValidation(
                sel.0[0],
                last_ind,
                self.top.len(),
            ))?;
        }
        let pos: Vec<_> = sel.0.iter().map(|i| &self.st.coords[*i]).cloned().collect();
        let atoms: Vec<_> = sel.0.iter().map(|i| self.top.atoms.to_atom(*i)).collect();
        self.st.add_coords(pos.into_iter());
        self.top.add_atoms(atoms.into_iter());
        Ok(Sel::from_iter(old_last + 1..self.len())?)
    }

    pub fn append_atoms<'a>(
        &mut self,
        atoms: impl Iterator<Item = Atom>,
        coords: impl Iterator<Item = &'a Pos>,
    ) -> Result<Sel, SelectionError> {
        let old_last = self.len() - 1;
        self.st.add_coords(coords.cloned());
        self.top.add_atoms(atoms);
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(Sel::from_iter(old_last + 1..self.len())?)
    }

    pub fn append_atom(&mut self, atom: &Atom, pos: &Pos) -> Result<Sel, SelectionError> {
        self.append_atoms(std::iter::once(atom.clone()), std::iter::once(pos))
    }

    pub fn append(
        &mut self,
        data: &(impl AtomProvider + PosProvider),
    ) -> Result<Sel, SelectionError> {
        let old_last = self.len() - 1;
        self.st.add_coords(data.iter_pos().cloned());
        self.top.add_atoms(data.iter_atoms().map(|a| Atom::from(&a)));
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(Sel::from_iter(old_last + 1..self.len())?)
    }

    pub fn remove(
        &mut self,
        removed: impl Iterator<Item = usize> + Clone,
    ) -> Result<(), BuilderError> {
        self.st.remove_coords(removed.clone())?;
        self.top.remove_atoms(removed)?;
        Ok(())
    }

    pub fn set_box_from(&mut self, src: &impl BoxProvider) {
        self.st.pbox = src.get_box().cloned();
    }

    pub fn multiply_periodically(&mut self, nbox: [usize; 3]) -> Result<(), SelectionError> {
        if self.get_box().is_none() {
            return Err(PeriodicBoxError::NoPbc)?;
        }
        let m = self.require_box()?.get_matrix();
        let all = self.select_all();
        for x in 0..=nbox[0] {
            for y in 0..=nbox[1] {
                for z in 0..=nbox[2] {
                    if x == 0 && y == 0 && z == 0 {
                        continue;
                    }
                    let added = self.append_from_self(&all)?;
                    let shift =
                        m.column(0) * x as Float + m.column(1) * y as Float + m.column(2) * z as Float;
                    self.select_bound_mut(added)?.translate(&shift);
                }
            }
        }
        // Scale the box
        self.get_box_mut().unwrap().scale_vectors([
            nbox[0] as Float,
            nbox[1] as Float,
            nbox[2] as Float,
        ])?;

        // Re-assign resindex
        self.top.assign_resindex();
        Ok(())
    }

    pub fn assign_resindex(&mut self) {
        self.top.assign_resindex();
    }
}

impl Selectable for System {
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(self, None)?))
    }
}

impl SelectableBound for System {
    fn select_bound(&self, def: impl SelectionDef) -> Result<SelOwnBound<'_>, SelectionError> {
        Ok(SelOwnBound {
            index: def.into_sel_index(self, None)?,
            sys: unsafe{&*self.get_system_ptr()},
        })
    }
}

impl SaveTopology for System {
    fn iter_atoms_dyn(&self) -> Box<dyn Iterator<Item = AtomRef<'_>> + '_> {
        Box::new(self.iter_atoms())
    }
    fn iter_bonds_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = BondRef<'a>> + 'a> {
        Box::new(BondProvider::iter_bonds(self))
    }
    fn num_bonds(&self) -> usize {
        BondProvider::num_bonds(self)
    }
}

impl SaveState for System {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Pos> + 'a> {
        Box::new(self.st.coords.iter())
    }

    fn iter_vel_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Vel> + 'a> {
        Box::new(self.st.velocities.iter())
    }

    fn iter_force_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Force> + 'a> {
        Box::new(self.st.forces.iter())
    }
}

impl SaveTopologyState for System {}

impl SystemProvider for System {
    fn get_system_ptr(&self) -> *const System {
        self
    }
}

impl SystemMutProvider for System {}

impl LenProvider for System {
    fn len(&self) -> usize {
        self.top.atoms.len()
    }
}

impl IndexProvider for System {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        (0..self.len()).into_iter()
    }
}

impl IndexParProvider for System {
    fn par_iter_index(&self) -> impl IndexedParallelIterator<Item = usize> {
        (0..self.len()).into_par_iter()
    }
}

// Analog of `bind` with `&sel >> &system` syntax
impl<'a> std::ops::Shr<&'a System> for &'a Sel {
    type Output = SelBound<'a>;
    fn shr(self, rhs: &'a System) -> Self::Output {
        rhs.bind(self)
    }
}

// Analog of `bind_mut` with `&sel >> &mut system` syntax
impl<'a> std::ops::Shr<&'a mut System> for &'a Sel {
    type Output = SelBoundMut<'a>;
    fn shr(self, rhs: &'a mut System) -> Self::Output {
        rhs.bind_mut(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    /// `bind_with_state` must read coordinates from the supplied external state
    /// (atoms still from the system's topology), and must NOT touch the system.
    #[test]
    fn bind_with_state_reads_external_coords() -> anyhow::Result<()> {
        let sys = System::from_file("tests/2lao.pdb")?;
        let sel = sys.select("name CA")?;

        // Baseline coords via the normal (system-owned) bind.
        let base: Vec<Pos> = sys.bind(&sel).iter_pos().cloned().collect();
        assert!(!base.is_empty());

        // An external state = the system's state shifted by +1 nm in x.
        let mut frame = sys.state().clone();
        for p in frame.coords.iter_mut() {
            p.coords.x += 1.0;
        }

        // bind_with_state reads the shifted coords (not the system's own).
        let got: Vec<Pos> = sys.bind_with_state(&sel, &frame).iter_pos().cloned().collect();
        assert_eq!(base.len(), got.len());
        for (b, g) in base.iter().zip(&got) {
            assert!((g.coords.x - (b.coords.x + 1.0)).abs() < 1e-4);
            assert_eq!(b.coords.y, g.coords.y);
            assert_eq!(b.coords.z, g.coords.z);
        }

        // The system's own state is untouched.
        let after: Vec<Pos> = sys.bind(&sel).iter_pos().cloned().collect();
        assert_eq!(base, after);
        Ok(())
    }

    /// Water (O-H, O-H) then methane-ish (C-H), as a bond-less system: the atoms are there
    /// but nothing says what is attached to what.
    fn unbonded_h_system() -> System {
        let mut top = Topology::default();
        for name in ["O", "H1", "H2", "C", "H3"] {
            top.atoms.push(&Atom::new().with_name(name).with_resname("MOL").with_resid(1).guess());
        }
        top.assign_resindex();
        let st = State {
            coords: (0..5).map(|i| Pos::new(i as Float * 0.1, 0.0, 0.0)).collect(),
            ..Default::default()
        };
        System::new(top, st).unwrap()
    }

    /// Connectivity from outside the structure file reaches the bond-reading machinery:
    /// `polh` / `apolh` classify hydrogens by what they are bonded to, so they match
    /// nothing until the bonds are installed.
    #[test]
    fn set_bonds_feeds_the_bond_graph() -> anyhow::Result<()> {
        let mut sys = unbonded_h_system();
        assert!(sys.select("polh").is_err(), "no bonds -> no polar hydrogens");
        assert!(sys.select("apolh").is_err(), "no bonds -> no apolar hydrogens");

        let mut bonds = BondStorage::default();
        for pair in [[0, 1], [0, 2], [3, 4]] {
            bonds.push(&Bond::new(pair[0], pair[1]));
        }
        sys.set_bonds(bonds)?;

        // H1/H2 hang off the O; H3 off the C.
        assert_eq!(sys.select("polh")?.iter_index().collect::<Vec<_>>(), vec![1, 2]);
        assert_eq!(sys.select("apolh")?.iter_index().collect::<Vec<_>>(), vec![4]);
        Ok(())
    }

    /// An out-of-range or self-referencing pair is rejected, leaving the old table in place.
    #[test]
    fn set_bonds_validates_indices() -> anyhow::Result<()> {
        let mut sys = unbonded_h_system();
        let mut good = BondStorage::default();
        good.push(&Bond::new(0, 1));
        sys.set_bonds(good)?;

        let mut bad = BondStorage::default();
        bad.push(&Bond::new(0, 999_999));
        assert!(sys.set_bonds(bad).is_err(), "out-of-range endpoint must be rejected");

        let mut selfy = BondStorage::default();
        selfy.push(&Bond::new(3, 3));
        assert!(sys.set_bonds(selfy).is_err(), "self-bond must be rejected");

        assert_eq!(sys.topology().bonds.len(), 1, "a rejected table changes nothing");
        Ok(())
    }
}
