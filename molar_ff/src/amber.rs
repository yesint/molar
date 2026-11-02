//! Amber mechanics, automatic GAFF/GAFF2 parameter completion, and minimization.
//!
//! Coordinates use nm; energies use kJ/mol; forces use kJ/mol/nm. Amber's
//! harmonic constants multiply the squared displacement without a factor 1/2.
//! [`AmberModel`] adds cutoff LJ/Coulomb and L-BFGS minimization to [`BondedModel`].
//! Ordinary LJ uses a potential switch; Coulomb uses a shifted force. There is
//! no PME, implicit solvent, periodic wrapping, or dispersion correction.
//!
//! ```no_run
//! use molar::prelude::*;
//! use molar_ff::{FFType, amber::BondedModel};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let sys = System::from_file("ligand.sdf")?;
//! let model = BondedModel::gaff(&sys, FFType::Gaff2)?;
//! let positions: Vec<[f64; 3]> = sys.state().coords.iter()
//!     .map(|p| [p.x as f64, p.y as f64, p.z as f64]).collect();
//! let result = model.evaluate(&positions)?;
//! println!("Bonded energy: {} kJ/mol", result.energy.total());
//! for (kind, atoms, resolution) in model.estimated_terms() {
//!     println!("{kind:?} {atoms:?}: {resolution:?}");
//! }
//! # Ok(())
//! # }
//! ```
//! Default cutoff mechanics with Espaloma charges, followed by minimization.
//!
//! ```no_run
//! use molar::prelude::*;
//! use molar_ff::{ApplyCharges, ChargeModel, FFType, amber::*};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut system = System::from_file("ligand.sdf")?;
//! system.apply_charges(ChargeModel::Espaloma)?;
//! let model = AmberModel::gaff(&system, FFType::Gaff2, NonbondedOptions::default())?;
//! let mut positions: Vec<[f64;3]> = system.state().coords.iter()
//!     .map(|p| [p.x as f64, p.y as f64, p.z as f64]).collect();
//! let movable = vec![true; positions.len()];
//! let report = model.minimize(&mut positions, &movable, MinimizeOptions::default())?;
//! println!("{:?}: {} kJ/mol", report.reason, report.final_energy.total());
//! # Ok(())
//! # }
//! ```

mod energy;
mod full;
mod minimize;
mod model;
mod nonbonded;
mod parameters;
mod resolve;

pub use energy::{BondedEnergy, BondedEvaluation, EvaluationError};
pub use model::{BondedError, BondedModel, BondedTerm, MissingParameter, ParameterPolicy};
pub use parameters::{
    AmberParameters, AngleParameter, BondParameter, FourierTerm, LennardJonesParameter,
    ParameterError, ParameterOrigin, Resolution, ResolutionMethod, TermKind,
};

pub use full::{AmberEnergy, AmberError, AmberEvaluation, AmberModel};
pub use minimize::{MinimizeError, MinimizeOptions, MinimizeResult, StopReason};
pub use nonbonded::{
    COULOMB_CONSTANT, NonbondedAtom, NonbondedEnergy, NonbondedError, NonbondedEvaluation,
    NonbondedException, NonbondedModel, NonbondedOptions,
};
