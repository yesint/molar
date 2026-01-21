use super::periodic_table::{ELEMENT_MASS, ELEMENT_NAME, ELEMENT_NAME_UPPER, ELEMENT_VDW};
use compact_str::CompactString;

pub trait AtomLike {
    /// Atom name.
    fn get_name(&self) -> &str;
    fn set_name(&mut self, name: &str);

    /// Residue name.
    fn get_resname(&self) -> &str;
    fn set_resname(&mut self, resname: &str);

    /// Residue id (aka residue number). This could be negative!
    fn get_resid(&self) -> isize;
    fn set_resid(&mut self, resid: isize);

    /// Residue index. Assigned when reading the topology.
    /// Unique for each contiguous span of resid. Starts from zero.
    fn get_resindex(&self) -> usize;
    fn set_resindex(&mut self, resindex: usize);

    /// Atomic number in the periodic table.
    fn get_atomic_number(&self) -> u8;
    fn set_atomic_number(&mut self, atomic_number: u8);

    /// Mass in atomic units
    fn get_mass(&self) -> f32;
    fn set_mass(&mut self, mass: f32);

    /// Charge in electric charges.
    fn get_charge(&self) -> f32;
    fn set_charge(&mut self, charge: f32);

    /// Name of the atom type.
    fn get_type_name(&self) -> &str;
    fn set_type_name(&mut self, type_name: &str);

    /// Unique id of the atom type.
    fn get_type_id(&self) -> u32;
    fn set_type_id(&mut self, type_id: u32);

    /// PDB chain identifier.
    fn get_chain(&self) -> char;
    fn set_chain(&mut self, chain: char);

    /// PDB B-factor.
    fn get_bfactor(&self) -> f32;
    fn set_bfactor(&mut self, bfactor: f32);

    /// PDB occupancy.
    fn get_occupancy(&self) -> f32;
    fn set_occupancy(&mut self, occupancy: f32);
}

/// Information about the atom except its coordinates.
#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Atom {
    /// Atom name.
    pub name: CompactString,
    /// Residue name.
    pub resname: CompactString,
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
    pub type_name: CompactString,
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
            if self.atomic_number == 0 && self.name.len() >= 2 {
                let c2 = self.name[i..=i + 1].to_ascii_uppercase();
                for an in 1..ELEMENT_NAME_UPPER.len() {
                    let el = ELEMENT_NAME_UPPER[an];
                    if el.len() == 2 && el == c2 {
                        // If the first letters are C,N,O,H,P be extra careful
                        // and only match to two-letter elements if the resname is the
                        // same as name (like in ions CA and CL).
                        // Otherwise skip to single-letter matching
                        match el.chars().next().unwrap() {
                            'C' | 'N' | 'O' | 'H' | 'P' => {
                                if self.name == self.resname {
                                    self.atomic_number = an as u8;
                                }
                            }
                            _ => {
                                self.atomic_number = an as u8;
                            }
                        }
                    }
                }
            }

            // If atomic_number is still 0 try 1-letter matching
            if self.atomic_number == 0 {
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

impl AtomLike for Atom {
    // Atom name
    fn get_name(&self) -> &str {
        self.name.as_str()
    }
    fn set_name(&mut self, name: &str) {
        self.name = name.into();
    }

    // Residue name
    fn get_resname(&self) -> &str {
        self.resname.as_str()
    }
    fn set_resname(&mut self, resname: &str) {
        self.resname = resname.into();
    }

    // Residue id
    fn get_resid(&self) -> isize {
        self.resid as isize
    }
    fn set_resid(&mut self, resid: isize) {
        self.resid = resid as i32
    }

    // Residue index
    fn get_resindex(&self) -> usize {
        self.resindex
    }
    fn set_resindex(&mut self, resindex: usize) {
        self.resindex = resindex;
    }

    // Atomic number
    fn get_atomic_number(&self) -> u8 {
        self.atomic_number
    }
    fn set_atomic_number(&mut self, atomic_number: u8) {
        self.atomic_number = atomic_number;
    }

    // Mass
    fn get_mass(&self) -> f32 {
        self.mass
    }
    fn set_mass(&mut self, mass: f32) {
        self.mass = mass;
    }

    // Charge
    fn get_charge(&self) -> f32 {
        self.charge
    }
    fn set_charge(&mut self, charge: f32) {
        self.charge = charge;
    }

    // Type name
    fn get_type_name(&self) -> &str {
        self.type_name.as_str()
    }
    fn set_type_name(&mut self, type_name: &str) {
        self.type_name = type_name.into();
    }

    // Type id
    fn get_type_id(&self) -> u32 {
        self.type_id
    }
    fn set_type_id(&mut self, type_id: u32) {
        self.type_id = type_id;
    }

    // Chain
    fn get_chain(&self) -> char {
        self.chain
    }
    fn set_chain(&mut self, chain: char) {
        self.chain = chain;
    }

    // B-factor
    fn get_bfactor(&self) -> f32 {
        self.bfactor
    }
    fn set_bfactor(&mut self, bfactor: f32) {
        self.bfactor = bfactor;
    }
    // Occupancy
    fn get_occupancy(&self) -> f32 {
        self.occupancy
    }
    fn set_occupancy(&mut self, occupancy: f32) {
        self.occupancy = occupancy;
    }
}  

