/// Information about the atom except its coordinates.
#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Atom {
    /// Atom name.
    pub name: String,
    /// Residue name.
    pub resname: String,
    /// Residue id (aka residue number). This could be negative!
    pub resid: i32, // Could be negative.
    /// Residue index. Assigned when reading the [topology](super::Topology).
    /// Unique for each contigous span of resid. Starts from zero.
    pub resindex: usize,
    /// Atomic number in the periodic table.
    pub atomic_number: u8,
    /// Mass in atomic units
    pub mass: f32,
    /// Charge in electroc charges.
    pub charge: f32,
    /// Name of the atom type.
    pub type_name: String,
    /// Unique id of the atom type.
    pub type_id: u32,
    // PDB chain identifier.
    pub chain: char,
    // PDB B-factor.
    pub bfactor: f32,
    /// PDB occupancy.
    pub occupancy: f32,
}

impl Atom {
    pub fn new() -> Self {
        Default::default()
    }

    /// Naive guessing of the mass and element from the atom name. 
    pub fn guess_element_and_mass_from_name(&mut self) {
        (self.atomic_number, self.mass) = match
            self.name.as_str()
            .trim_start_matches(char::is_numeric)
            .chars().next().unwrap()
        {
            'C' => (6, 12.0),
            'O' => (8, 16.0),
            'N' => (7, 14.0),
            'S' => (16, 32.0),
            'H' => (1, 1.0),
            'P' => (15, 31.0),
            'F' => (9, 19.0),
            'B' => (5, 11.0),
            _ => (0, 1.0), // Unknown atom
        }
    }

    /// Returns atom's Van der Waals radius based on its atomic number.
    /// If the element is not recognized returns a default value of 0.15 nm.
    pub fn vdw(&self) -> f32 {
        use super::periodic_table::*;

        if self.atomic_number<=0 {
            get_vdw_from_atom_name(&self.name)
        } else {
            ELEMENT_VDW[self.atomic_number as usize] * 0.1
        }
    }
}
