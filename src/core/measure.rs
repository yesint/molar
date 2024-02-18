use super::providers::*;
use super::algorithms::*;
use super::Pos;
use anyhow::Result;
//==============================================================
// Traits for measuring (immutable access)
//==============================================================

// Main trait giving scoped read-only data provider
pub trait GuardedQuery {
    type Guard<'a> where Self: 'a;
    fn guard<'a>(&'a self) -> Self::Guard<'a>;
}

/// Trait for analysis requiring only positions
pub trait MeasurePos: GuardedQuery
    where for<'a> Self::Guard<'a>: PosProvider
{
    fn min_max(&self) -> (Pos, Pos) {
        min_max(&self.guard())
    }

    fn center_of_geometry(&self) -> Pos {
        center_of_geometry(&self.guard())
    }

    fn rmsd_(sel1: &Self, sel2: &Self) -> Result<f32> {
        let dp1 = sel1.guard();
        let dp2 = sel2.guard();
        rmsd(&dp1, &dp2)
    }
}

/// Trait for analysis requiring positions and masses
pub trait MeasureMasses: GuardedQuery
    where for<'a> Self::Guard<'a>: PosProvider + MassesProvider
{
    fn center_of_mass(&self) -> Result<Pos> {
        center_of_mass(&self.guard())
    }

    fn fit_transform(sel1: &Self, sel2: &Self) -> Result<nalgebra::IsometryMatrix3<f32>> {
        let dp1 = sel1.guard();
        let dp2 = sel2.guard();
        fit_transform(&dp1, &dp2)
    }

    fn fit_transform_at_origin(sel1: &Self, sel2: &Self) -> Result<nalgebra::IsometryMatrix3<f32>> {
        let dp1 = sel1.guard();
        let dp2 = sel2.guard();
        fit_transform_at_origin(&dp1, &dp2)
    }

    fn  rmsd_mw(sel1: &Self, sel2: &Self) -> Result<f32> {
        let dp1 = sel1.guard();
        let dp2 = sel2.guard();
        rmsd_mw(&dp1, &dp2)
    }
}

/// Trait for analysis requiring positions, masses and pbc
pub trait MeasurePeriodic: GuardedQuery
    where for<'a> Self::Guard<'a>: PosProvider + MassesProvider + BoxProvider
{
    fn center_of_mass_pbc(&self) -> Result<Pos> {
        center_of_mass_pbc(&self.guard())
    }
}
