use nalgebra::Unit;
use anyhow::Result;
use super::providers::*;
use super::algorithms::*;
use super::PbcDims;
use super::Vector3f;
use super::PBC_FULL;

//==============================================================
// Traits for modification (mutable access)
//==============================================================

pub trait GuardedModify {
    type GuardMut<'a> where Self: 'a;
    fn guard_mut<'a>(&'a self) -> Self::GuardMut<'a>;
}

/// Trait for modification requiring only positions
pub trait ModifyPos: GuardedModify 
    where for<'a> Self::GuardMut<'a>: PosMutProvider + PosProvider
{
    fn translate(&self, shift: Vector3f) {
        translate(&mut self.guard_mut(), shift)
    }

    fn rotate(&self, ax: &Unit<Vector3f>, ang: f32) {
        rotate(&mut self.guard_mut(), ax, ang)
    }

    fn apply_transform(&self, tr: &nalgebra::IsometryMatrix3<f32>) {
        apply_transform(&mut self.guard_mut(), tr)
    }
}

/// Trait for modification requiring positions and pbc
pub trait ModifyPeriodic: GuardedModify 
    where for<'a> Self::GuardMut<'a>: PosMutProvider + BoxProvider
{
    fn unwrap_simple_dim(&self, dims: PbcDims) -> Result<()> {
        unwrap_simple_dim(&mut self.guard_mut(), dims)
    }

    fn unwrap_simple(&self) -> Result<()> {
        self.unwrap_simple_dim(PBC_FULL)
    }
}

/// Trait for modification requiring random access positions and pbc
pub trait ModifyRandomAccess: GuardedModify
    where for<'a> Self::GuardMut<'a>: PosMutProvider + PosProvider + BoxProvider + RandomPosMutProvider
{
    fn unwrap_connectivity(&self, cutoff: f32) -> Result<()> {
        unwrap_connectivity_dim(&mut self.guard_mut(), cutoff, &PBC_FULL)
    }

    fn unwrap_connectivity_dim(&self, cutoff: f32, dims: &PbcDims) -> Result<()> {
        unwrap_connectivity_dim(&mut self.guard_mut(), cutoff, dims)
    }
}

/// Trait for modification requiring atoms
pub trait ModifyAtoms: GuardedModify
    where for<'a> Self::GuardMut<'a>: AtomsMutProvider
{
    fn assign_resindex(&self) {
        assign_resindex(&mut self.guard_mut())
    }
}

