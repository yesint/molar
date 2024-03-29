
#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Atom {
    // Mandatory fields
    pub name: String,
    pub resname: String,
    pub resid: i32, // Could be negative
    pub resindex: usize,
    // Atom physical properties from topology
    pub atomic_number: u8,
    pub mass: f32,
    pub charge: f32,
    pub type_name: String,
    pub type_id: u32,
    // Specific PDB fields
    pub chain: char,
    pub bfactor: f32,
    pub occupancy: f32,
}

impl Atom {
    pub fn new() -> Self {
        Default::default()
    }

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

    pub fn vdw(&self) -> f32 {
        use super::periodic_table::*;

        if self.atomic_number<=0 {
            get_vdw_from_atom_name(&self.name)
        } else {
            ELEMENT_VDW[self.atomic_number as usize] * 0.1
        }
    }
}
