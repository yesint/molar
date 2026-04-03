use super::periodic_table::{ELEMENT_MASS, ELEMENT_NAME, ELEMENT_NAME_UPPER, ELEMENT_VDW};
use tinystr::TinyAsciiStr;

/// Stack-allocated ASCII atom string (max 8 bytes, no heap allocation).
pub type AtomStr = TinyAsciiStr<8>;

pub(crate) const ATOM_NAME_EXPECT:      &str = "atom name fits in 8 bytes";
pub(crate) const ATOM_RESNAME_EXPECT:   &str = "residue name fits in 8 bytes";
pub(crate) const ATOM_TYPE_NAME_EXPECT: &str = "atom type name fits in 8 bytes";

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
#[derive(Debug, Clone)]
pub struct Atom {
    /// Atom name.
    pub name: AtomStr,
    /// Residue name.
    pub resname: AtomStr,
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
    pub type_name: AtomStr,
    /// Unique id of the atom type.
    pub type_id: u32,
    // PDB chain identifier.
    pub chain: char,
    // PDB B-factor.
    pub bfactor: f32,
    /// PDB occupancy.
    pub occupancy: f32,
}

impl Default for Atom {
    fn default() -> Self {
        let empty = AtomStr::try_from_str("").unwrap();
        Atom {
            name: empty,
            resname: empty,
            type_name: empty,
            resid: 0,
            resindex: 0,
            atomic_number: 0,
            mass: 0.0,
            charge: 0.0,
            type_id: 0,
            chain: ' ',
            bfactor: 0.0,
            occupancy: 1.0,
        }
    }
}

impl Atom {
    pub fn new() -> Self {
        Default::default()
    }

    // Fluent setters
    pub fn with_name(mut self, name: &str) -> Self {
        self.name = AtomStr::try_from_str(name).expect(ATOM_NAME_EXPECT);
        self
    }
    pub fn with_resname(mut self, resname: &str) -> Self {
        self.resname = AtomStr::try_from_str(resname).expect(ATOM_RESNAME_EXPECT);
        self
    }
    pub fn with_resid(mut self, resid: i32) -> Self { self.resid = resid; self }
    pub fn with_resindex(mut self, resindex: usize) -> Self { self.resindex = resindex; self }
    pub fn with_atomic_number(mut self, n: u8) -> Self { self.atomic_number = n; self }
    pub fn with_mass(mut self, mass: f32) -> Self { self.mass = mass; self }
    pub fn with_charge(mut self, charge: f32) -> Self { self.charge = charge; self }
    pub fn with_type_name(mut self, type_name: &str) -> Self {
        self.type_name = AtomStr::try_from_str(type_name).expect(ATOM_TYPE_NAME_EXPECT);
        self
    }
    pub fn with_type_id(mut self, type_id: u32) -> Self { self.type_id = type_id; self }
    pub fn with_chain(mut self, chain: char) -> Self { self.chain = chain; self }
    pub fn with_bfactor(mut self, bfactor: f32) -> Self { self.bfactor = bfactor; self }
    pub fn with_occupancy(mut self, occupancy: f32) -> Self { self.occupancy = occupancy; self }

    /// Chainable version of `guess_element_and_mass_from_name()`.
    pub fn guess(mut self) -> Self {
        self.guess_element_and_mass_from_name();
        self
    }

    // Element constructors
    pub fn hydrogen() -> Self { Atom::new().with_name("H").guess() }
    pub fn carbon() -> Self { Atom::new().with_name("C").guess() }
    pub fn nitrogen() -> Self { Atom::new().with_name("N").guess() }
    pub fn oxygen() -> Self { Atom::new().with_name("O").guess() }
    pub fn phosphorus() -> Self { Atom::new().with_name("P").guess() }
    pub fn sulfur() -> Self { Atom::new().with_name("S").guess() }

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

/// Returns the uppercase element symbol for the given atomic number (e.g. `"FE"`, `"C"`, `"HE"`).
/// Returns `""` for atomic number 0 (unknown).
pub(crate) fn element_symbol(atomic_number: u8) -> &'static str {
    ELEMENT_NAME_UPPER
        .get(atomic_number as usize)
        .copied()
        .unwrap_or("")
}

impl<T: AtomLike> From<&T> for Atom {
    fn from(a: &T) -> Self {
        Atom::new()
            .with_name(a.get_name())
            .with_resname(a.get_resname())
            .with_resid(a.get_resid() as i32)
            .with_resindex(a.get_resindex())
            .with_atomic_number(a.get_atomic_number())
            .with_mass(a.get_mass())
            .with_charge(a.get_charge())
            .with_type_name(a.get_type_name())
            .with_type_id(a.get_type_id())
            .with_chain(a.get_chain())
            .with_bfactor(a.get_bfactor())
            .with_occupancy(a.get_occupancy())
    }
}

impl AtomLike for Atom {
    // Atom name
    fn get_name(&self) -> &str {
        self.name.as_str()
    }
    fn set_name(&mut self, name: &str) {
        self.name = AtomStr::try_from_str(name).expect(ATOM_NAME_EXPECT);
    }

    // Residue name
    fn get_resname(&self) -> &str {
        self.resname.as_str()
    }
    fn set_resname(&mut self, resname: &str) {
        self.resname = AtomStr::try_from_str(resname).expect(ATOM_RESNAME_EXPECT);
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
        self.type_name = AtomStr::try_from_str(type_name).expect(ATOM_TYPE_NAME_EXPECT);
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

