//! MDL **Molfile / SDF** (`.mol` / `.sdf`) reader & writer — a standalone handler,
//! no external cheminformatics dependency (the format is simple and molar's IO is
//! deliberately self-contained).
//!
//! Supports the **V2000** connection table: a counts line, an atom block
//! (`x y z element …`, coordinates in Ångströms), and a bond block
//! (`atom1 atom2 type …`, 1-based atom indices). The bond `type` maps to molar's
//! [`BondOrder`] (1→Single, 2→Double, 3→Triple, 4→Aromatic) — this is the format
//! that motivates bonds carrying an order. A `.sdf` may hold several records
//! separated by `$$$$`; successive [`read`](FileFormatHandler::read) calls return
//! successive records (each its own topology), like reading multi-MODEL PDB.
//!
//! **Out of scope (deliberate):**
//! - **V3000** connection tables (detected and rejected with a clear error).
//! - **Data fields** (the `> <name>` blocks after `M  END`) are skipped on read.
//!   Exposing them is future work and will get its own *non-trait* accessor rather
//!   than being forced through the generic `read` interface.

use crate::atom::element_symbol;
use crate::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
    path::Path,
};
use thiserror::Error;

/// Molfile coordinates are in Ångströms; molar works in nanometres.
const ANGSTROM_TO_NM: Float = 0.1;
const NM_TO_ANGSTROM: Float = 10.0;

pub struct SdfFileHandler {
    reader: Option<BufReader<DynSource>>,
    writer: Option<BufWriter<File>>,
    /// Write `$$$$` record terminators (true for `.sdf`, false for a bare `.mol`).
    sdf: bool,
    /// For the separate topology/state read path.
    stored_topology: Option<Topology>,
    stored_state: Option<State>,
    at_least_one_state_read: bool,
}

#[derive(Debug, Error)]
pub enum SdfHandlerError {
    #[error("can't open sdf/mol file for reading")]
    OpenRead(#[source] std::io::Error),

    #[error("can't open sdf/mol file for writing")]
    OpenWrite(#[source] std::io::Error),

    #[error("sdf/mol file is empty")]
    Empty,

    #[error("truncated molfile header (need title, program, comment, counts lines)")]
    TruncatedHeader,

    #[error("malformed counts line: {0:?}")]
    BadCounts(String),

    #[error("V3000 molfiles are not supported yet (only V2000)")]
    V3000Unsupported,

    #[error("truncated atom block at atom {0}")]
    TruncatedAtoms(usize),

    #[error("truncated bond block at bond {0}")]
    TruncatedBonds(usize),

    #[error("malformed number in atom {0}")]
    BadAtom(usize, #[source] ParseFloatError),

    #[error("malformed index/order in bond {0}")]
    BadBond(usize, #[source] ParseIntError),

    #[error("io error")]
    Io(#[from] std::io::Error),
}

/// Read one line; `Ok(None)` at end of file.
fn next_line(buf: &mut impl BufRead) -> std::io::Result<Option<String>> {
    let mut s = String::new();
    if buf.read_line(&mut s)? == 0 {
        Ok(None)
    } else {
        Ok(Some(s))
    }
}

/// Parse a fixed-width integer field `[start, start+len)` of a molfile line
/// (right-justified, possibly run together with its neighbour — so column slicing,
/// not whitespace splitting). A short line yields 0 (trailing fields are optional).
fn int_field(line: &str, start: usize, len: usize) -> Result<i64, ParseIntError> {
    let end = (start + len).min(line.len());
    let field = line.get(start..end).unwrap_or("").trim();
    if field.is_empty() {
        Ok(0)
    } else {
        field.parse()
    }
}

impl SdfFileHandler {
    pub(crate) fn from_source(src: DynSource) -> Result<Self, FileFormatError> {
        Ok(Self {
            reader: Some(BufReader::new(src)),
            writer: None,
            sdf: true,
            stored_topology: None,
            stored_state: None,
            at_least_one_state_read: false,
        })
    }
}

impl FileFormatHandler for SdfFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let file = File::open(fname).map_err(SdfHandlerError::OpenRead)?;
        Self::from_source(DynSource(Box::new(file)))
    }

    fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let fname = fname.as_ref();
        let sdf = fname
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("sdf"))
            .unwrap_or(false);
        Ok(Self {
            reader: None,
            writer: Some(BufWriter::new(
                File::create(fname).map_err(SdfHandlerError::OpenWrite)?,
            )),
            sdf,
            stored_topology: None,
            stored_state: None,
            at_least_one_state_read: false,
        })
    }

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        let buf = self.reader.as_mut().unwrap();

        // True EOF (no bytes left at all).
        if buf.fill_buf()?.is_empty() {
            return if self.at_least_one_state_read {
                Err(FileFormatError::Eof)?
            } else {
                Err(SdfHandlerError::Empty)?
            };
        }

        // The 4-line header: title (line 1, *may be blank* — so we don't skip blanks),
        // program/timestamp, comment, counts. If the file ends mid-header after only
        // blank lines, that's trailing whitespace past the last record → EOF.
        let mut header = [String::new(), String::new(), String::new(), String::new()];
        let mut saw_content = false;
        for slot in header.iter_mut() {
            match next_line(buf)? {
                Some(l) => {
                    saw_content |= !l.trim().is_empty();
                    *slot = l;
                }
                None if saw_content => return Err(SdfHandlerError::TruncatedHeader)?,
                None if self.at_least_one_state_read => return Err(FileFormatError::Eof)?,
                None => return Err(SdfHandlerError::Empty)?,
            }
        }
        let counts = &header[3];

        if counts.contains("V3000") {
            return Err(SdfHandlerError::V3000Unsupported)?;
        }
        let natoms = int_field(counts, 0, 3)
            .map_err(|_| SdfHandlerError::BadCounts(counts.clone()))? as usize;
        let nbonds = int_field(counts, 3, 3)
            .map_err(|_| SdfHandlerError::BadCounts(counts.clone()))? as usize;
        if natoms == 0 {
            return Err(SdfHandlerError::BadCounts(counts.clone()))?;
        }

        // Atom block: `x y z element …` (coords whitespace-separated in practice).
        let mut atoms: Vec<Atom> = Vec::with_capacity(natoms);
        let mut coords: Vec<Pos> = Vec::with_capacity(natoms);
        for i in 0..natoms {
            let line = next_line(buf)?.ok_or(SdfHandlerError::TruncatedAtoms(i))?;
            let mut t = line.split_whitespace();
            let mut coord = [0.0 as Float; 3];
            for c in coord.iter_mut() {
                *c = t
                    .next()
                    .ok_or(SdfHandlerError::TruncatedAtoms(i))?
                    .parse::<Float>()
                    .map_err(|e| SdfHandlerError::BadAtom(i, e))?
                    * ANGSTROM_TO_NM;
            }
            let elem = t.next().ok_or(SdfHandlerError::TruncatedAtoms(i))?;
            coords.push(Pos::new(coord[0], coord[1], coord[2]));
            atoms.push(
                Atom::new()
                    .with_name(elem)
                    .with_resname("MOL")
                    .with_resid(1)
                    .with_chain('A')
                    .guess(),
            );
        }

        // Bond block: `atom1 atom2 type …`, fixed 3-wide columns, 1-based indices.
        let mut bonds: Vec<Bond> = Vec::with_capacity(nbonds);
        for i in 0..nbonds {
            let line = next_line(buf)?.ok_or(SdfHandlerError::TruncatedBonds(i))?;
            let a1 = int_field(&line, 0, 3).map_err(|e| SdfHandlerError::BadBond(i, e))?;
            let a2 = int_field(&line, 3, 3).map_err(|e| SdfHandlerError::BadBond(i, e))?;
            let ty = int_field(&line, 6, 3).map_err(|e| SdfHandlerError::BadBond(i, e))?;
            if a1 < 1 || a2 < 1 || a1 as usize > natoms || a2 as usize > natoms {
                return Err(SdfHandlerError::TruncatedBonds(i))?;
            }
            let order = match ty {
                2 => BondOrder::Double,
                3 => BondOrder::Triple,
                4 => BondOrder::Aromatic,
                _ => BondOrder::Single, // 1, plus the "query" orders 5–8 → single
            };
            bonds.push(Bond::with_order(a1 as usize - 1, a2 as usize - 1, order));
        }

        // Skip the rest of the record: properties (`M  …`, `M  END`) and any data
        // fields, stopping after the `$$$$` separator (SDF) or at EOF (.mol).
        while let Some(l) = next_line(buf)? {
            if l.trim_end() == "$$$$" {
                break;
            }
        }

        let mut top = Topology::default();
        top.atoms = atoms;
        top.bonds = bonds;
        top.assign_resindex();
        self.at_least_one_state_read = true;
        Ok((top, State { coords, time: 0.0, pbox: None, ..Default::default() }))
    }

    fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        if let Some(top) = self.stored_topology.take() {
            Ok(top)
        } else {
            let (top, st) = self.read()?;
            self.stored_state.get_or_insert(st);
            Ok(top)
        }
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        if let Some(st) = self.stored_state.take() {
            Ok(st)
        } else {
            let (top, st) = self.read()?;
            self.stored_topology.get_or_insert(top);
            Ok(st)
        }
    }

    fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileFormatError> {
        let sdf = self.sdf;
        let w = self.writer.as_mut().ok_or(FileFormatError::NotWritable)?;
        let natoms = data.len();
        let nbonds = data.num_bonds();

        // Header: title / program / comment, then the V2000 counts line.
        writeln!(w)?; // title (blank — molar has no per-molecule name here)
        writeln!(w, "  molar")?; // program/timestamp line
        writeln!(w)?; // comment
        writeln!(
            w,
            "{natoms:>3}{nbonds:>3}  0  0  0  0  0  0  0  0999 V2000"
        )?;

        for (at, pos) in data.iter_atoms_dyn().zip(data.iter_pos_dyn()) {
            let elem = element_symbol(at.atomic_number);
            let sym = if elem.is_empty() { at.name.as_str() } else { elem };
            writeln!(
                w,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                pos.x * NM_TO_ANGSTROM,
                pos.y * NM_TO_ANGSTROM,
                pos.z * NM_TO_ANGSTROM,
                sym,
            )?;
        }

        for b in data.iter_bonds_dyn() {
            let ty = match b.order {
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                BondOrder::Single | BondOrder::Unspecified => 1,
            };
            // 1-based atom indices.
            writeln!(w, "{:>3}{:>3}{:>3}  0  0  0  0", b.i1 + 1, b.i2 + 1, ty)?;
        }

        writeln!(w, "M  END")?;
        if sdf {
            writeln!(w, "$$$$")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// A minimal ethene V2000 molfile (C=C + 4 H), coordinates in Å.
    const ETHENE: &str = "\
ethene
  test

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3300    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5600    0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5600   -0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8900    0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8900   -0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  2  5  1  0  0  0  0
  2  6  1  0  0  0  0
M  END
$$$$
";

    #[test]
    fn reads_atoms_bonds_orders() {
        let mut h = SdfFileHandler::from_source(DynSource(Box::new(Cursor::new(
            ETHENE.as_bytes().to_vec(),
        ))))
        .unwrap();
        let (top, st) = h.read().unwrap();
        assert_eq!(top.atoms.len(), 6);
        assert_eq!(st.coords.len(), 6);
        assert_eq!(top.bonds.len(), 5);
        // First two atoms are carbons; coordinates converted Å → nm.
        assert_eq!(top.atoms[0].atomic_number, 6);
        assert!((st.coords[1].x - 0.133).abs() < 1e-5, "C–C along x ≈ 0.133 nm");
        // The C=C bond carries a double order; the C–H bonds are single.
        assert_eq!(top.bonds[0], Bond::with_order(0, 1, BondOrder::Double));
        assert_eq!(top.bonds[1].order, BondOrder::Single);
        // A second read hits EOF (single record).
        assert!(matches!(h.read(), Err(FileFormatError::Eof)));
    }

    #[test]
    fn write_read_round_trip_preserves_order() {
        // Build a tiny molecule: C=O with two single C–H, then write + re-read.
        let mut top = Topology::default();
        for sym in ["C", "O", "H", "H"] {
            top.atoms.push(Atom::new().with_name(sym).with_resname("MOL").with_resid(1).guess());
        }
        top.bonds = vec![
            Bond::with_order(0, 1, BondOrder::Double),
            Bond::with_order(0, 2, BondOrder::Single),
            Bond::with_order(0, 3, BondOrder::Single),
        ];
        top.assign_resindex();
        let st = State {
            coords: vec![
                Pos::new(0.0, 0.0, 0.0),
                Pos::new(0.123, 0.0, 0.0),
                Pos::new(-0.05, 0.09, 0.0),
                Pos::new(-0.05, -0.09, 0.0),
            ],
            ..Default::default()
        };
        let sys = System::new(top, st).unwrap();

        // Write to an in-memory buffer via a temp file (the writer is File-backed),
        // then read it back.
        let dir = std::env::temp_dir();
        let path = dir.join("molar_sdf_roundtrip.sdf");
        {
            let mut wh = SdfFileHandler::create(&path).unwrap();
            wh.write(&sys).unwrap();
        }
        let mut rh = SdfFileHandler::open(&path).unwrap();
        let (top2, st2) = rh.read().unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(top2.atoms.len(), 4);
        assert_eq!(top2.bonds.len(), 3);
        assert_eq!(top2.bonds[0].order, BondOrder::Double);
        assert_eq!(top2.bonds[1].order, BondOrder::Single);
        assert_eq!(top2.atoms[1].atomic_number, 8); // O survived
        assert!((st2.coords[1].x - 0.123).abs() < 1e-3, "C–O distance preserved");
    }

    #[test]
    fn rejects_v3000() {
        let v3000 = "\
mol

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  END
$$$$
";
        let mut h = SdfFileHandler::from_source(DynSource(Box::new(Cursor::new(
            v3000.as_bytes().to_vec(),
        ))))
        .unwrap();
        assert!(h.read().is_err());
    }
}
