use super::periodic_table::{ELEMENT_MASS, ELEMENT_NAME, ELEMENT_NAME_UPPER, ELEMENT_VDW};
use crate::aliases::Float;
use tinystr::TinyAsciiStr;

/// Stack-allocated ASCII atom string (max 8 bytes, no heap allocation).
pub type AtomStr = TinyAsciiStr<8>;

pub(crate) const ATOM_NAME_EXPECT:      &str = "atom name fits in 8 bytes";
pub(crate) const ATOM_RESNAME_EXPECT:   &str = "residue name fits in 8 bytes";
pub(crate) const ATOM_TYPE_NAME_EXPECT: &str = "atom type name fits in 8 bytes";

/// Per-atom perceived-chemistry flags (ring membership, aromaticity, …).
///
/// Stored in an optional column (see the SoA `AtomStorage`); this newtype replaces the
/// former bit-packing of ring/aromatic flags into the top two bits of `type_id`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AtomFlags(u8);

impl AtomFlags {
    const IN_RING: u8 = 1 << 0;
    const AROMATIC: u8 = 1 << 1;

    pub fn is_in_ring(&self) -> bool {
        self.0 & Self::IN_RING != 0
    }
    pub fn set_in_ring(&mut self, v: bool) {
        if v {
            self.0 |= Self::IN_RING;
        } else {
            self.0 &= !Self::IN_RING;
        }
    }
    pub fn is_aromatic(&self) -> bool {
        self.0 & Self::AROMATIC != 0
    }
    pub fn set_aromatic(&mut self, v: bool) {
        if v {
            self.0 |= Self::AROMATIC;
        } else {
            self.0 &= !Self::AROMATIC;
        }
    }
}

/// Read-only access to atom properties.
///
/// Implemented by the owned [`Atom`] and by the borrowed column proxies
/// [`AtomRef`](crate::AtomRef) / [`AtomRefMut`](crate::AtomRefMut). Optional (force-field /
/// chemistry) getters return `None` when the property was never assigned.
pub trait AtomLike {
    /// Atom name.
    fn get_name(&self) -> &str;
    /// Residue name.
    fn get_resname(&self) -> &str;
    /// Residue id (aka residue number). This could be negative!
    fn get_resid(&self) -> isize;
    /// Residue index. Unique for each contiguous span of resid. Starts from zero.
    fn get_resindex(&self) -> usize;
    /// Atomic number in the periodic table.
    fn get_atomic_number(&self) -> u8;
    /// Mass in atomic units.
    fn get_mass(&self) -> Float;
    /// Partial (working) charge in electron charges.
    fn get_charge(&self) -> Float;
    /// PDB chain identifier.
    fn get_chain(&self) -> char;
    /// PDB B-factor.
    fn get_bfactor(&self) -> Float;
    /// PDB occupancy.
    fn get_occupancy(&self) -> Float;

    // --- Optional (force-field / chemistry) properties: `None` if never assigned ---
    /// Name of the atom type (force field).
    fn get_type_name(&self) -> Option<&str>;
    /// Unique id of the atom type.
    fn get_type_id(&self) -> Option<u32>;
    /// Integer formal charge (e.g. from an SDF `M  CHG` record).
    fn get_formal_charge(&self) -> Option<i32>;
    /// Perceived-chemistry flags.
    fn get_flags(&self) -> Option<AtomFlags>;

    /// Whether the atom is a member of some ring (false when flags are unset).
    fn is_in_ring(&self) -> bool {
        self.get_flags().map_or(false, |f| f.is_in_ring())
    }
    /// Whether the atom belongs to an aromatic ring (false when flags are unset).
    fn is_aromatic(&self) -> bool {
        self.get_flags().map_or(false, |f| f.is_aromatic())
    }

    /// Van der Waals radius (nm) from the atomic number; 0.15 nm for unknown (Z=0).
    fn vdw(&self) -> Float {
        ELEMENT_VDW[self.get_atomic_number() as usize] as Float * 0.1
    }
}

/// Mutable access to atom properties. Setters for the optional (force-field / chemistry)
/// properties assign `Some(..)`.
pub trait AtomLikeMut: AtomLike {
    fn set_name(&mut self, name: &str);
    fn set_resname(&mut self, resname: &str);
    fn set_resid(&mut self, resid: isize);
    fn set_resindex(&mut self, resindex: usize);
    fn set_atomic_number(&mut self, atomic_number: u8);
    fn set_mass(&mut self, mass: Float);
    fn set_charge(&mut self, charge: Float);
    fn set_chain(&mut self, chain: char);
    fn set_bfactor(&mut self, bfactor: Float);
    fn set_occupancy(&mut self, occupancy: Float);

    fn set_type_name(&mut self, type_name: &str);
    fn set_type_id(&mut self, type_id: u32);
    fn set_formal_charge(&mut self, formal_charge: i32);
    fn set_flags(&mut self, flags: AtomFlags);

    fn set_in_ring(&mut self, v: bool) {
        let mut f = self.get_flags().unwrap_or_default();
        f.set_in_ring(v);
        self.set_flags(f);
    }
    fn set_aromatic(&mut self, v: bool) {
        let mut f = self.get_flags().unwrap_or_default();
        f.set_aromatic(v);
        self.set_flags(f);
    }
}

/// Information about the atom except its coordinates.
///
/// Retained as the detached, densely-packed construction/interchange row even though
/// [`Topology`](crate::Topology) stores atoms column-wise (Struct-of-Arrays). The four
/// optional force-field/chemistry fields are `Option` — `None` means "never assigned".
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
    pub mass: Float,
    /// Partial (working) charge in electron charges. The integer formal charge is stored
    /// separately in [`formal_charge`](Self::formal_charge).
    pub charge: Float,
    // PDB chain identifier.
    pub chain: char,
    // PDB B-factor.
    pub bfactor: Float,
    /// PDB occupancy.
    pub occupancy: Float,
    // --- optional force-field / chemistry properties ---
    /// Name of the atom type (force field).
    pub type_name: Option<AtomStr>,
    /// Unique id of the atom type.
    pub type_id: Option<u32>,
    /// Integer formal charge (e.g. from an SDF `M  CHG` record).
    pub formal_charge: Option<i32>,
    /// Perceived-chemistry flags (ring membership, aromaticity, …).
    pub flags: Option<AtomFlags>,
}

impl Default for Atom {
    fn default() -> Self {
        let empty = AtomStr::try_from_str("").unwrap();
        Atom {
            name: empty,
            resname: empty,
            resid: 0,
            resindex: 0,
            atomic_number: 0,
            mass: 0.0,
            charge: 0.0,
            chain: ' ',
            bfactor: 0.0,
            occupancy: 1.0,
            type_name: None,
            type_id: None,
            formal_charge: None,
            flags: None,
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
    pub fn with_mass(mut self, mass: Float) -> Self { self.mass = mass; self }
    pub fn with_charge(mut self, charge: Float) -> Self { self.charge = charge; self }
    pub fn with_chain(mut self, chain: char) -> Self { self.chain = chain; self }
    pub fn with_bfactor(mut self, bfactor: Float) -> Self { self.bfactor = bfactor; self }
    pub fn with_occupancy(mut self, occupancy: Float) -> Self { self.occupancy = occupancy; self }
    pub fn with_type_name(mut self, type_name: &str) -> Self {
        self.type_name = Some(AtomStr::try_from_str(type_name).expect(ATOM_TYPE_NAME_EXPECT));
        self
    }
    pub fn with_type_id(mut self, type_id: u32) -> Self { self.type_id = Some(type_id); self }
    pub fn with_formal_charge(mut self, formal_charge: i32) -> Self {
        self.formal_charge = Some(formal_charge);
        self
    }
    pub fn with_flags(mut self, flags: AtomFlags) -> Self { self.flags = Some(flags); self }

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
            // Attempt 2-letter matching if possible (only when >= 2 chars remain from i)
            if self.atomic_number == 0 && i + 1 < self.name.len() {
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
        // Fill mass field — periodic table is stored as Float for compactness; cast at lookup.
        self.mass = ELEMENT_MASS[self.atomic_number as usize] as Float;
    }

    // Naive guessing of the mass and element from the atom name.
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
        let mut at = Atom::new()
            .with_name(a.get_name())
            .with_resname(a.get_resname())
            .with_resid(a.get_resid() as i32)
            .with_resindex(a.get_resindex())
            .with_atomic_number(a.get_atomic_number())
            .with_mass(a.get_mass())
            .with_charge(a.get_charge())
            .with_chain(a.get_chain())
            .with_bfactor(a.get_bfactor())
            .with_occupancy(a.get_occupancy());
        if let Some(t) = a.get_type_name() {
            at = at.with_type_name(t);
        }
        if let Some(id) = a.get_type_id() {
            at = at.with_type_id(id);
        }
        if let Some(fc) = a.get_formal_charge() {
            at = at.with_formal_charge(fc);
        }
        if let Some(f) = a.get_flags() {
            at = at.with_flags(f);
        }
        at
    }
}

impl AtomLike for Atom {
    fn get_name(&self) -> &str {
        self.name.as_str()
    }
    fn get_resname(&self) -> &str {
        self.resname.as_str()
    }
    fn get_resid(&self) -> isize {
        self.resid as isize
    }
    fn get_resindex(&self) -> usize {
        self.resindex
    }
    fn get_atomic_number(&self) -> u8 {
        self.atomic_number
    }
    fn get_mass(&self) -> Float {
        self.mass
    }
    fn get_charge(&self) -> Float {
        self.charge
    }
    fn get_chain(&self) -> char {
        self.chain
    }
    fn get_bfactor(&self) -> Float {
        self.bfactor
    }
    fn get_occupancy(&self) -> Float {
        self.occupancy
    }
    fn get_type_name(&self) -> Option<&str> {
        self.type_name.as_ref().map(|s| s.as_str())
    }
    fn get_type_id(&self) -> Option<u32> {
        self.type_id
    }
    fn get_formal_charge(&self) -> Option<i32> {
        self.formal_charge
    }
    fn get_flags(&self) -> Option<AtomFlags> {
        self.flags
    }
}

impl AtomLikeMut for Atom {
    fn set_name(&mut self, name: &str) {
        self.name = AtomStr::try_from_str(name).expect(ATOM_NAME_EXPECT);
    }
    fn set_resname(&mut self, resname: &str) {
        self.resname = AtomStr::try_from_str(resname).expect(ATOM_RESNAME_EXPECT);
    }
    fn set_resid(&mut self, resid: isize) {
        self.resid = resid as i32
    }
    fn set_resindex(&mut self, resindex: usize) {
        self.resindex = resindex;
    }
    fn set_atomic_number(&mut self, atomic_number: u8) {
        self.atomic_number = atomic_number;
    }
    fn set_mass(&mut self, mass: Float) {
        self.mass = mass;
    }
    fn set_charge(&mut self, charge: Float) {
        self.charge = charge;
    }
    fn set_chain(&mut self, chain: char) {
        self.chain = chain;
    }
    fn set_bfactor(&mut self, bfactor: Float) {
        self.bfactor = bfactor;
    }
    fn set_occupancy(&mut self, occupancy: Float) {
        self.occupancy = occupancy;
    }
    fn set_type_name(&mut self, type_name: &str) {
        self.type_name = Some(AtomStr::try_from_str(type_name).expect(ATOM_TYPE_NAME_EXPECT));
    }
    fn set_type_id(&mut self, type_id: u32) {
        self.type_id = Some(type_id);
    }
    fn set_formal_charge(&mut self, formal_charge: i32) {
        self.formal_charge = Some(formal_charge);
    }
    fn set_flags(&mut self, flags: AtomFlags) {
        self.flags = Some(flags);
    }
}
