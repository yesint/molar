use super::algorithms::*;
use super::providers::*;
use super::Pos;
use anyhow::Result;
//==============================================================
// Traits for measuring (immutable access)
//==============================================================

// Main trait giving scoped read-only data provider
pub trait GuardedQuery {
    type Guard<'a>
    where
        Self: 'a;
    fn guard<'a>(&'a self) -> Self::Guard<'a>;
}

/// Trait for analysis requiring only positions
pub trait MeasurePos: GuardedQuery
where
    for<'a> Self::Guard<'a>: PosProvider,
{
    fn min_max(&self) -> (Pos, Pos) {
        min_max(&self.guard())
    }

    fn center_of_geometry(&self) -> Pos {
        center_of_geometry(&self.guard())
    }

    fn rmsd_(sel1: &Self, sel2: &Self) -> Result<f32> {
        rmsd(&sel1.guard(), &sel2.guard())
    }
}

/// Trait for analysis requiring positions and masses
pub trait MeasureMasses: GuardedQuery
where
    for<'a> Self::Guard<'a>: PosProvider + MassesProvider,
{
    fn center_of_mass(&self) -> Result<Pos> {
        center_of_mass(&self.guard())
    }

    fn fit_transform(sel1: &Self, sel2: &Self) -> Result<nalgebra::IsometryMatrix3<f32>> {
        fit_transform(&sel1.guard(), &sel2.guard())
    }

    fn fit_transform_at_origin(sel1: &Self, sel2: &Self) -> Result<nalgebra::IsometryMatrix3<f32>> {
        fit_transform_at_origin(&sel1.guard(), &sel2.guard())
    }

    fn rmsd_mw(sel1: &Self, sel2: &Self) -> Result<f32> {
        rmsd_mw(&sel1.guard(), &sel2.guard())
    }
}

/// Trait for analysis requiring positions, masses and pbc
pub trait MeasurePeriodic: GuardedQuery
where
    for<'a> Self::Guard<'a>: PosProvider + MassesProvider + BoxProvider,
{
    fn center_of_mass_pbc(&self) -> Result<Pos> {
        center_of_mass_pbc(&self.guard())
    }
}

/*
pub trait MeasureIndexed: GuardedQuery
where
    for<'a> Self::Guard<'a>: IoIndexProvider + AtomsProvider + PosProvider,
{
    fn split<'a>(
        &self,
        func: fn(usize, &Atom, &Pos) -> usize,
    ) -> SelectionSplitter<impl IoIndexProvider + AtomsProvider + PosProvider,impl Iterator<Item = (usize,usize)>> {
        SelectionSplitter::new(self.guard(), func)
    }
}
*/