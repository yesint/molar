//! Amber bonded mechanics, with automatic GAFF/GAFF2 parameter completion.
//!
//! Coordinates use nm; energies use kJ/mol; forces use kJ/mol/nm. Amber's
//! harmonic constants multiply the squared displacement without a factor 1/2.
//! This module contains no nonbonded energy or minimizer.
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
mod energy;
mod model;
mod parameters;
mod resolve;

pub use energy::{BondedEnergy, BondedEvaluation, EvaluationError};
pub use model::{BondedError, BondedModel, BondedTerm, MissingParameter, ParameterPolicy};
pub use parameters::{
    AmberParameters, AngleParameter, BondParameter, FourierTerm, ParameterError, ParameterOrigin,
    Resolution, ResolutionMethod, TermKind,
};
