use super::{ELEMENT_MASS, ELEMENT_NAME, ELEMENT_NAME_UPPER, ELEMENT_VDW};

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

    pub fn guess_element_from_name(&mut self) {
        self.atomic_number = 0;
        // Index of the first letter in atom name
        if let Some(i) = self.name.find(|c: char| c.is_ascii_alphabetic()) {
            // Match special cases when atom name doesn't 
            // start with the element name at all
            match self.name.as_str() {
                "SOD" => self.atomic_number = 11, //Sodium
                "POT" => self.atomic_number = 19, //Potassium
                _ => (),
            } 

            // Find matching element name in periodic table
            // Attempt 2-letter matching if possible
            if self.atomic_number ==0 && self.name.len() >= 2 {
                let c2 = self.name[i..=i + 1].to_ascii_uppercase();
                for an in 1..ELEMENT_NAME_UPPER.len() {
                    let el = ELEMENT_NAME_UPPER[an];
                    if el.len() == 2 && el == c2 {
                        // If the first letters are C,N,O,H,P be extra careful
                        // and only match to two-letter elements if the resname is the 
                        // same as name (like in ions CA and CL).
                        // Otherwise skip to single-letter matching
                        match el.chars().next().unwrap() {
                            'C'|'N'|'O'|'H'|'P' => {
                                if self.name == self.resname {
                                    self.atomic_number = an as u8;
                                }
                            },
                            _ => {
                                self.atomic_number = an as u8;
                            },
                        }
                    }
                }
            }

            // If atomic_number is still 0 try 1-letter matching
            if self.atomic_number ==0 {
                for an in 1..ELEMENT_NAME.len() {
                    let el = ELEMENT_NAME[an];
                    if el.len() == 1 && el == &self.name[i..=i] {
                        self.atomic_number = an as u8;
                    }
                }
            }
        }
    }

    pub fn guess_element_and_mass_from_name(&mut self) {
        self.guess_element_from_name();
        // Fill mass field
        self.mass = ELEMENT_MASS[self.atomic_number as usize];
    }

    /// Naive guessing of the mass and element from the atom name.
    // pub fn guess_element_and_mass_from_name(&mut self) {
    //     (self.atomic_number, self.mass) = match self
    //         .name
    //         .as_str()
    //         .trim_start_matches(char::is_numeric)
    //         .chars()
    //         .next()
    //         .unwrap()
    //     {
    //         'C' => (6, ELEMENT_MASS[6]),
    //         'O' => (8, ELEMENT_MASS[8]),
    //         'N' => (7, ELEMENT_MASS[7]),
    //         'S' => (16, ELEMENT_MASS[16]),
    //         'H' => (1, ELEMENT_MASS[1]),
    //         'P' => (15, ELEMENT_MASS[15]),
    //         'F' => (9, ELEMENT_MASS[9]),
    //         'B' => (5, ELEMENT_MASS[5]),
    //         _ => (0, 1.0), // Unknown atom
    //     }
    // }

    /// Returns atom's Van der Waals radius based on its atomic number.
    /// If the element is not recognized returns a default value of 0.15 nm (atomnum '0').
    pub fn vdw(&self) -> f32 {
        ELEMENT_VDW[self.atomic_number as usize] * 0.1
    }
}
