//! Partial LAMMPS data support for KG polymers. Coordinates are unwrapped and
//! translated to a zero box origin. Interaction parameters are not retained.
use crate::prelude::*;
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};
use thiserror::Error;

/// Conversion from LAMMPS input units to MolAR units. Writing uses the inverse.
/// LJ units have no fixed physical scale: the default assigns 1 nm and 1 atomic
/// mass unit to one input unit. The data-file title is not used to infer units.
#[derive(Debug, Clone, Copy)]
pub struct LammpsOptions {
    /// Nanometers per input distance unit.
    pub length_scale: Float,
    /// Atomic mass units per input mass unit.
    pub mass_scale: Float,
}
impl Default for LammpsOptions {
    fn default() -> Self {
        Self {
            length_scale: 1.0,
            mass_scale: 1.0,
        }
    }
}
impl LammpsOptions {
    pub(crate) fn validate(self) -> Result<Self, LammpsHandlerError> {
        if !self.length_scale.is_finite()
            || self.length_scale <= 0.0
            || !self.mass_scale.is_finite()
            || self.mass_scale <= 0.0
        {
            return Err(LammpsHandlerError::InvalidScales);
        }
        Ok(self)
    }
}

#[derive(Debug, Error)]
pub enum LammpsHandlerError {
    #[error("LAMMPS {section}, line {line}: {message}")]
    Parse {
        section: String,
        line: usize,
        message: String,
    },
    #[error("LAMMPS length and mass scales must be finite and positive")]
    InvalidScales,
    #[error("cannot write LAMMPS data: {0}")]
    Write(String),
    #[error("LAMMPS IO error")]
    Io(#[from] std::io::Error),
}
fn bad(section: &str, line: usize, message: impl Into<String>) -> LammpsHandlerError {
    LammpsHandlerError::Parse {
        section: section.into(),
        line,
        message: message.into(),
    }
}
fn number<T: std::str::FromStr>(
    s: &str,
    section: &str,
    line: usize,
) -> Result<T, LammpsHandlerError> {
    s.parse()
        .map_err(|_| bad(section, line, format!("invalid number {s:?}")))
}
fn real(s: &str, section: &str, line: usize) -> Result<Float, LammpsHandlerError> {
    let v = number::<Float>(s, section, line)?;
    if !v.is_finite() {
        return Err(bad(section, line, "number must be finite"));
    }
    Ok(v)
}
fn type_id(s: &str, count: usize, section: &str, line: usize) -> Result<u32, LammpsHandlerError> {
    let id = number::<u32>(s, section, line)?;
    if id == 0 || id as usize > count {
        return Err(bad(section, line, "type ID outside declared range"));
    }
    Ok(id)
}
fn section_count(name: &str) -> Option<&'static str> {
    Some(match name {
        "Atoms" | "Velocities" => "atoms",
        "Masses" | "Pair Coeffs" => "atom types",
        "Bonds" => "bonds",
        "Angles" => "angles",
        "Dihedrals" => "dihedrals",
        "Impropers" => "impropers",
        "Bond Coeffs" => "bond types",
        "Angle Coeffs" | "BondBond Coeffs" | "BondAngle Coeffs" => "angle types",
        "Dihedral Coeffs"
        | "MiddleBondTorsion Coeffs"
        | "EndBondTorsion Coeffs"
        | "AngleTorsion Coeffs"
        | "AngleAngleTorsion Coeffs"
        | "BondBond13 Coeffs" => "dihedral types",
        "Improper Coeffs" | "AngleAngle Coeffs" => "improper types",
        "Ellipsoids" => "ellipsoids",
        "Lines" => "lines",
        "Triangles" => "triangles",
        "Bodies" => "bodies",
        "PairIJ Coeffs" => "pairij",
        _ => return None,
    })
}
struct AtomRow {
    id: u64,
    mol: i32,
    ty: u32,
    pos: Pos,
    image: Vector3f,
    line: usize,
}

pub(crate) struct LammpsFileHandler {
    reader: Option<BufReader<DynSource>>,
    writer: Option<BufWriter<File>>,
    options: LammpsOptions,
    consumed: bool,
    written: bool,
    stored_topology: Option<Topology>,
    stored_state: Option<State>,
}
impl LammpsFileHandler {
    pub(crate) fn from_source(
        src: DynSource,
        options: LammpsOptions,
    ) -> Result<Self, FileFormatError> {
        Ok(Self {
            reader: Some(BufReader::new(src)),
            writer: None,
            options: options.validate()?,
            consumed: false,
            written: false,
            stored_topology: None,
            stored_state: None,
        })
    }
    pub(crate) fn open_with_options(
        path: &Path,
        options: LammpsOptions,
    ) -> Result<Self, FileFormatError> {
        let options = options.validate()?;
        Self::from_source(
            DynSource(Box::new(File::open(path).map_err(LammpsHandlerError::Io)?)),
            options,
        )
    }
    pub(crate) fn create_with_options(
        path: &Path,
        options: LammpsOptions,
    ) -> Result<Self, FileFormatError> {
        let options = options.validate()?;
        Ok(Self {
            reader: None,
            writer: Some(BufWriter::new(
                File::create(path).map_err(LammpsHandlerError::Io)?,
            )),
            options,
            consumed: false,
            written: false,
            stored_topology: None,
            stored_state: None,
        })
    }

    fn parse(&mut self) -> Result<(Topology, State), LammpsHandlerError> {
        let reader = self.reader.as_mut().expect("read path checked");
        let mut text = String::new();
        if reader.read_line(&mut text)? == 0 {
            return Err(bad("header", 1, "empty file"));
        }
        let mut counts = HashMap::<String, usize>::new();
        let mut bounds = [None; 3];
        let mut tilt = Vector3f::zeros();
        let mut section = String::new();
        let mut section_line = 1;
        let mut rows = 0usize;
        let mut expected = 0usize;
        let mut seen = HashSet::new();
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut atom_ids = HashSet::new();
        let mut bond_ids = HashSet::new();
        let mut masses = HashMap::new();
        let mut atom_width = None;
        let mut line = 1usize;
        loop {
            text.clear();
            if reader.read_line(&mut text)? == 0 {
                break;
            }
            line += 1;
            let (content, comment) = text.split_once('#').unwrap_or((&text, ""));
            let content = content.trim();
            if content.is_empty() {
                continue;
            }
            let fields: Vec<_> = content.split_whitespace().collect();
            if let Some(key) = section_count(content) {
                if !section.is_empty() && rows != expected {
                    return Err(bad(
                        &section,
                        line,
                        format!("expected {expected} rows, found {rows}"),
                    ));
                }
                if !seen.insert(content.to_owned()) {
                    return Err(bad(content, line, "duplicate section"));
                }
                if bounds.iter().any(Option::is_none) {
                    return Err(bad("header", line, "all three box bounds are required"));
                }
                if counts.get("atoms").copied().unwrap_or(0) == 0
                    || counts.get("atom types").copied().unwrap_or(0) == 0
                {
                    return Err(bad(
                        "header",
                        line,
                        "positive atom and atom type counts are required",
                    ));
                }
                if content == "Atoms" {
                    let style = comment.split_whitespace().next().unwrap_or("");
                    if !["", "molecular", "bond", "angle", "id", "atom-ID"].contains(&style) {
                        return Err(bad(
                            content,
                            line,
                            format!("unsupported atom style/comment {style:?}"),
                        ));
                    }
                }
                if content == "Bodies" {
                    return Err(bad(
                        content,
                        line,
                        "variable-length Bodies section is unsupported",
                    ));
                }
                expected = if key == "pairij" {
                    let n = counts.get("atom types").copied().unwrap_or(0);
                    n.checked_add(1)
                        .and_then(|v| n.checked_mul(v))
                        .map(|v| v / 2)
                        .ok_or_else(|| bad("header", line, "type count overflow"))?
                } else {
                    counts.get(key).copied().unwrap_or(0)
                };
                section = content.into();
                section_line = line;
                rows = 0;
                continue;
            }
            if section.is_empty() {
                if fields.iter().any(|s| ["avec", "bvec", "cvec"].contains(s))
                    || content.ends_with("abc origin")
                {
                    return Err(bad(
                        "header",
                        line,
                        "general triclinic boxes are unsupported",
                    ));
                }
                let mut handled = false;
                for (axis, suffix) in ["xlo xhi", "ylo yhi", "zlo zhi"].iter().enumerate() {
                    if content.ends_with(suffix) && fields.len() == 4 {
                        if bounds[axis].is_some() {
                            return Err(bad("header", line, "duplicate box bounds"));
                        }
                        let lo = real(fields[0], "header", line)?;
                        let hi = real(fields[1], "header", line)?;
                        if hi <= lo {
                            return Err(bad(
                                "header",
                                line,
                                "box upper bound must exceed lower bound",
                            ));
                        }
                        bounds[axis] = Some((lo, hi));
                        handled = true;
                    }
                }
                if handled {
                    continue;
                }
                if content.ends_with("xy xz yz") && fields.len() == 6 {
                    if !seen.insert("tilt".into()) {
                        return Err(bad("header", line, "duplicate tilt factors"));
                    }
                    tilt = Vector3f::new(
                        real(fields[0], "header", line)?,
                        real(fields[1], "header", line)?,
                        real(fields[2], "header", line)?,
                    );
                    continue;
                }
                let key = fields[1..].join(" ");
                if [
                    "atoms",
                    "bonds",
                    "angles",
                    "dihedrals",
                    "impropers",
                    "atom types",
                    "bond types",
                    "angle types",
                    "dihedral types",
                    "improper types",
                    "ellipsoids",
                    "lines",
                    "triangles",
                    "bodies",
                ]
                .contains(&key.as_str())
                {
                    let count = number(fields[0], "header", line)?;
                    if counts.insert(key, count).is_some() {
                        return Err(bad("header", line, "duplicate count"));
                    }
                } else if key.starts_with("extra ") && key.ends_with("per atom") {
                    let _: usize = number(fields[0], "header", line)?;
                } else {
                    return Err(bad(
                        "header",
                        line,
                        format!("unsupported header or section {content:?}"),
                    ));
                }
                continue;
            }
            if rows >= expected {
                return Err(bad(
                    &section,
                    line,
                    format!("extra row or unsupported section {content:?}"),
                ));
            }
            match section.as_str() {
                "Atoms" => {
                    if ![6, 9].contains(&fields.len())
                        || atom_width.is_some_and(|w| w != fields.len())
                    {
                        return Err(bad(
                            &section,
                            line,
                            "expected consistently six or nine columns",
                        ));
                    }
                    atom_width = Some(fields.len());
                    let id = number::<u64>(fields[0], &section, line)?;
                    if id == 0 || !atom_ids.insert(id) {
                        return Err(bad(&section, line, "zero or duplicate atom ID"));
                    }
                    let mol = number::<i32>(fields[1], &section, line)?;
                    if mol < 0 {
                        return Err(bad(
                            &section,
                            line,
                            "molecule ID must be nonnegative and fit i32",
                        ));
                    }
                    let ty = type_id(fields[2], counts["atom types"], &section, line)?;
                    let pos = Pos::new(
                        real(fields[3], &section, line)?,
                        real(fields[4], &section, line)?,
                        real(fields[5], &section, line)?,
                    );
                    let mut image = Vector3f::zeros();
                    if fields.len() == 9 {
                        for i in 0..3 {
                            image[i] = number::<i32>(fields[6 + i], &section, line)? as Float;
                        }
                    }
                    atoms.push(AtomRow {
                        id,
                        mol,
                        ty,
                        pos,
                        image,
                        line,
                    });
                }
                "Masses" => {
                    if fields.len() != 2 {
                        return Err(bad(&section, line, "expected type ID and mass"));
                    }
                    let ty = type_id(fields[0], counts["atom types"], &section, line)?;
                    let mass = real(fields[1], &section, line)? * self.options.mass_scale;
                    if !mass.is_finite() || mass <= 0.0 || masses.insert(ty, mass).is_some() {
                        return Err(bad(
                            &section,
                            line,
                            "mass must be positive, finite, and unique per type",
                        ));
                    }
                }
                "Bonds" => {
                    if fields.len() != 4 {
                        return Err(bad(
                            &section,
                            line,
                            "expected bond ID, type ID, and two atom IDs",
                        ));
                    }
                    let id = number::<u64>(fields[0], &section, line)?;
                    if id == 0 || !bond_ids.insert(id) {
                        return Err(bad(&section, line, "zero or duplicate bond ID"));
                    }
                    type_id(
                        fields[1],
                        counts.get("bond types").copied().unwrap_or(0),
                        &section,
                        line,
                    )?;
                    bonds.push((
                        number::<u64>(fields[2], &section, line)?,
                        number::<u64>(fields[3], &section, line)?,
                        line,
                    ));
                }
                _ => {} // Known sections outside the supported scope are discarded.
            }
            rows += 1;
        }
        if !section.is_empty() && rows != expected {
            return Err(bad(
                &section,
                line,
                format!("expected {expected} rows, found {rows}"),
            ));
        }
        if atoms.is_empty() {
            return Err(bad("Atoms", line, "missing or empty section"));
        }
        if bonds.len() != counts.get("bonds").copied().unwrap_or(0) {
            return Err(bad("Bonds", line, "missing section or incorrect count"));
        }
        let mut matrix = Matrix3f::zeros();
        let mut origin = Vector3f::zeros();
        for i in 0..3 {
            let (lo, hi) = bounds[i].unwrap();
            matrix[(i, i)] = hi - lo;
            origin[i] = lo;
        }
        matrix[(0, 1)] = tilt[0];
        matrix[(0, 2)] = tilt[1];
        matrix[(1, 2)] = tilt[2];
        let scaled = matrix * self.options.length_scale;
        let pbox = PeriodicBox::from_matrix(scaled)
            .map_err(|e| bad("header", section_line, e.to_string()))?;
        atoms.sort_unstable_by_key(|a| (a.mol, a.id));
        let mut top = Topology::default();
        let mut coords = Vec::with_capacity(atoms.len());
        let mut id_map = HashMap::new();
        for (i, row) in atoms.iter().enumerate() {
            let mass = if seen.contains("Masses") {
                *masses
                    .get(&row.ty)
                    .ok_or_else(|| bad("Masses", row.line, "missing mass for atom type"))?
            } else {
                self.options.mass_scale
            };
            let mut atom = Atom::new()
                .with_name("B")
                .with_resname("MOL")
                .with_type_id(row.ty)
                .with_mass(mass);
            atom.resid = row.mol;
            top.atoms.push(&atom);
            let pos = (row.pos.coords - origin + matrix * row.image) * self.options.length_scale;
            if !pos.iter().all(|v| v.is_finite()) {
                return Err(bad("Atoms", row.line, "scaled position is not finite"));
            }
            coords.push(Pos::from(pos));
            id_map.insert(row.id, i);
            if row.mol > 0 {
                if i > 0 && atoms[i - 1].mol == row.mol {
                    top.molecules.last_mut().unwrap()[1] = i;
                } else {
                    top.molecules.push([i, i]);
                }
            }
        }
        for (a, b, line) in bonds {
            let i = *id_map
                .get(&a)
                .ok_or_else(|| bad("Bonds", line, "unknown atom ID"))?;
            let j = *id_map
                .get(&b)
                .ok_or_else(|| bad("Bonds", line, "unknown atom ID"))?;
            if i == j {
                return Err(bad("Bonds", line, "self bond"));
            }
            top.bonds.push(&Bond::new(i, j));
        }
        top.assign_resindex();
        Ok((
            top,
            State {
                coords,
                pbox: Some(pbox),
                ..Default::default()
            },
        ))
    }
}
impl FileFormatHandler for LammpsFileHandler {
    fn open(path: impl AsRef<Path>) -> Result<Self, FileFormatError> {
        Self::open_with_options(path.as_ref(), LammpsOptions::default())
    }
    fn create(path: impl AsRef<Path>) -> Result<Self, FileFormatError> {
        Self::create_with_options(path.as_ref(), LammpsOptions::default())
    }
    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        if self.reader.is_none() {
            return Err(FileFormatError::NotReadable);
        }
        if self.consumed {
            return Err(FileFormatError::Eof);
        }
        self.consumed = true;
        Ok(self.parse()?)
    }
    fn read_topology(&mut self) -> Result<Topology, FileFormatError> {
        if let Some(top) = self.stored_topology.take() {
            return Ok(top);
        }
        let (top, state) = self.read()?;
        self.stored_state = Some(state);
        Ok(top)
    }
    fn read_state(&mut self) -> Result<State, FileFormatError> {
        if let Some(state) = self.stored_state.take() {
            return Ok(state);
        }
        let (top, state) = self.read()?;
        self.stored_topology = Some(top);
        Ok(state)
    }
    fn write_topology(&mut self, data: &dyn SaveTopology) -> Result<(), FileFormatError> {
        if self.writer.is_none() {
            return Err(FileFormatError::NotWritable);
        }
        if self.written || self.stored_topology.is_some() {
            return Err(LammpsHandlerError::Write("topology already supplied".into()).into());
        }
        let atoms: Vec<_> = data.iter_atoms_dyn().collect();
        let top = Topology {
            bonds: data.bonds_for_write()?.into_iter().collect(),
            atoms: atoms.iter().map(Atom::from).collect(),
            ..Default::default()
        };
        self.stored_topology = Some(top);
        Ok(())
    }
    fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileFormatError> {
        if self.writer.is_none() {
            return Err(FileFormatError::NotWritable);
        }
        let top = self.stored_topology.take().ok_or_else(|| {
            LammpsHandlerError::Write("write_topology must precede write_state".into())
        })?;
        let state = State {
            coords: data.iter_pos_dyn().copied().collect(),
            pbox: data.get_box().cloned(),
            ..Default::default()
        };
        let system =
            System::new(top, state).map_err(|e| LammpsHandlerError::Write(e.to_string()))?;
        self.write(&system)
    }
    fn write(&mut self, data: &dyn SaveTopologyState) -> Result<(), FileFormatError> {
        if self.writer.is_none() {
            return Err(FileFormatError::NotWritable);
        }
        if self.written {
            return Err(LammpsHandlerError::Write(
                "only one structure per file is supported".into(),
            )
            .into());
        }
        let fail = |s: &str| LammpsHandlerError::Write(s.into());
        let atoms: Vec<_> = data.iter_atoms_dyn().collect();
        let positions: Vec<_> = data.iter_pos_dyn().collect();
        if atoms.is_empty() || atoms.len() != positions.len() || atoms.len() != data.len() {
            return Err(fail("atom and coordinate counts must agree and be positive").into());
        }
        let matrix = data
            .get_box()
            .ok_or_else(|| fail("a periodic box is required"))?
            .get_matrix();
        if !matrix.iter().all(|v| v.is_finite())
            || (0..3).any(|i| matrix[(i, i)] <= 0.0)
            || matrix[(1, 0)] != 0.0
            || matrix[(2, 0)] != 0.0
            || matrix[(2, 1)] != 0.0
        {
            return Err(fail("box must be orthogonal or restricted triclinic").into());
        }
        let inverse = matrix.try_inverse().ok_or_else(|| fail("invalid box"))?;
        let output_box = matrix / self.options.length_scale;
        if !output_box.iter().all(|v| v.is_finite()) || (0..3).any(|i| output_box[(i, i)] <= 0.0) {
            return Err(fail("scaled box is invalid").into());
        }
        let mut masses = BTreeMap::new();
        let mut output_atoms = Vec::new();
        for (atom, pos) in atoms.iter().zip(&positions) {
            let ty = atom.get_type_id().unwrap_or(1);
            let mass = atom.get_mass() / self.options.mass_scale;
            if ty == 0 || !mass.is_finite() || mass <= 0.0 {
                return Err(
                    fail("type IDs and masses must be positive; masses must be finite").into(),
                );
            }
            if let Some(previous) = masses.insert(ty, mass) {
                if previous != mass {
                    return Err(fail("atoms of one type have different masses").into());
                }
            }
            if atom.get_resid() < 0 {
                return Err(fail("molecule IDs must be nonnegative").into());
            }
            if !pos.coords.iter().all(|v| v.is_finite()) {
                return Err(fail("coordinates must be finite").into());
            }
            let image = (inverse * pos.coords).map(Float::floor);
            if !image
                .iter()
                .all(|v| *v >= i32::MIN as Float && (*v as f64) <= i32::MAX as f64)
            {
                return Err(fail("image flags must fit i32").into());
            }
            let wrapped = (pos.coords - matrix * image) / self.options.length_scale;
            if !wrapped.iter().all(|v| v.is_finite()) {
                return Err(fail("scaled position is not finite").into());
            }
            output_atoms.push((ty, atom.get_resid(), wrapped, image));
        }
        let bonds = data.bonds_for_write()?;
        let type_count = *masses.last_key_value().unwrap().0;
        let w = self.writer.as_mut().unwrap();
        self.written = true; // Never append another structure, including after an IO failure.
        writeln!(
            w,
            "MolAR LAMMPS data; length scale = {}, mass scale = {}\n",
            self.options.length_scale, self.options.mass_scale
        )?;
        writeln!(
            w,
            "{} atoms\n{} bonds\n\n{} atom types\n{} bond types\n",
            atoms.len(),
            bonds.len(),
            type_count,
            usize::from(!bonds.is_empty())
        )?;
        for (i, key) in ["xlo xhi", "ylo yhi", "zlo zhi"].iter().enumerate() {
            writeln!(w, "0 {} {key}", output_box[(i, i)])?;
        }
        if output_box[(0, 1)] != 0.0 || output_box[(0, 2)] != 0.0 || output_box[(1, 2)] != 0.0 {
            writeln!(
                w,
                "{} {} {} xy xz yz",
                output_box[(0, 1)],
                output_box[(0, 2)],
                output_box[(1, 2)]
            )?;
        }
        writeln!(w, "\nMasses\n")?;
        for ty in 1..=type_count {
            writeln!(w, "{ty} {}", masses.get(&ty).copied().unwrap_or(1.0))?;
        }
        writeln!(w, "\nAtoms # molecular\n")?;
        for (i, (ty, mol, pos, image)) in output_atoms.iter().enumerate() {
            writeln!(
                w,
                "{} {mol} {ty} {} {} {} {} {} {}",
                i + 1,
                pos[0],
                pos[1],
                pos[2],
                image[0] as i32,
                image[1] as i32,
                image[2] as i32
            )?;
        }
        if !bonds.is_empty() {
            writeln!(w, "\nBonds\n")?;
            for (i, bond) in bonds.iter().enumerate() {
                writeln!(w, "{} 1 {} {}", i + 1, bond.i1 + 1, bond.i2 + 1)?;
            }
        }
        w.flush()?;
        Ok(())
    }
}
