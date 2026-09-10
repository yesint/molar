//! PDBx/mmCIF text reader and writer.
//!
//! This handler intentionally reads the structure-sized subset that MolAR can represent:
//! `_atom_site`, `_cell`, and covalent `_struct_conn` records.  The lexer still implements
//! CIF 1.1 quoting and multiline text rules, so unrelated categories can be skipped safely.

use crate::atom::{atomic_number_from_symbol, element_symbol};
use crate::periodic_table::ELEMENT_MASS;
use crate::prelude::*;
use std::{
    collections::{HashMap, VecDeque},
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    path::Path,
};
use thiserror::Error;

const ANGSTROM_TO_NM: Float = 0.1;
const NM_TO_ANGSTROM: Float = 10.0;

pub struct CifFileHandler {
    reader: Option<BufReader<DynSource>>,
    writer: Option<BufWriter<File>>,
    frames: Option<VecDeque<(Topology, State)>>,
    stored_topology: Option<Topology>,
    stored_state: Option<State>,
    at_least_one_state_read: bool,
}

#[derive(Debug, Error)]
pub enum CifHandlerError {
    #[error("can't open cif file for reading")]
    OpenRead(#[source] std::io::Error),
    #[error("can't open cif file for writing")]
    OpenWrite(#[source] std::io::Error),
    #[error("cif file is empty or has no _atom_site records")]
    Empty,
    #[error("invalid UTF-8 in cif file")]
    Utf8(#[from] std::string::FromUtf8Error),
    #[error("cif syntax error at line {line}, column {column}: {message}")]
    Syntax {
        line: usize,
        column: usize,
        message: String,
    },
    #[error("missing required cif field {0}")]
    MissingField(&'static str),
    #[error("invalid value {value:?} for cif field {field}")]
    InvalidValue { field: &'static str, value: String },
    #[error("{field} value {value:?} does not fit MolAR's 8-byte atom string")]
    StringTooLong { field: &'static str, value: String },
    #[error("cif chain identifier {0:?} cannot be represented by MolAR's one-character chain")]
    UnsupportedChain(String),
    #[error("cif chain identifiers {first:?} and {second:?} both map to {mapped:?}")]
    ChainCollision {
        first: String,
        second: String,
        mapped: char,
    },
    #[error("model {model} does not have the same atom identities and order as the first model")]
    ModelMismatch { model: i32 },
    #[error("unsupported cif bond order {0:?}")]
    UnsupportedBondOrder(String),
    #[error("conflicting orders for cif bond {0}-{1}")]
    ConflictingBondOrder(usize, usize),
    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),
    #[error("io error")]
    Io(#[from] std::io::Error),
}

#[derive(Clone, Debug)]
struct Token {
    text: String,
    line: usize,
    column: usize,
}

struct Lexer {
    input: Vec<u8>,
    pos: usize,
    line: usize,
    column: usize,
}

impl Lexer {
    fn new(input: String) -> Self {
        Self {
            input: input.into_bytes(),
            pos: 0,
            line: 1,
            column: 1,
        }
    }

    fn bump(&mut self) -> Option<u8> {
        let b = *self.input.get(self.pos)?;
        self.pos += 1;
        if b == b'\n' {
            self.line += 1;
            self.column = 1;
        } else {
            self.column += 1;
        }
        Some(b)
    }

    fn syntax<T>(
        &self,
        line: usize,
        column: usize,
        message: impl Into<String>,
    ) -> Result<T, CifHandlerError> {
        Err(CifHandlerError::Syntax {
            line,
            column,
            message: message.into(),
        })
    }

    fn next_token(&mut self) -> Result<Option<Token>, CifHandlerError> {
        loop {
            while self
                .input
                .get(self.pos)
                .is_some_and(|b| b.is_ascii_whitespace())
            {
                self.bump();
            }
            if self.input.get(self.pos) == Some(&b'#') {
                while let Some(b) = self.bump() {
                    if b == b'\n' {
                        break;
                    }
                }
                continue;
            }
            break;
        }
        if self.pos == self.input.len() {
            return Ok(None);
        }

        let line = self.line;
        let column = self.column;
        let start = self.pos;
        let first = self.bump().unwrap();

        // A semicolon in column one starts a multiline text field.  Its closing semicolon
        // must also be the first byte on a line.
        if first == b';' && column == 1 {
            let value_start = self.pos;
            loop {
                if self.pos >= self.input.len() {
                    return self.syntax(line, column, "unterminated semicolon text field");
                }
                if self.column == 1 && self.input[self.pos] == b';' {
                    let value_end = self.pos;
                    self.bump();
                    while self
                        .input
                        .get(self.pos)
                        .is_some_and(|b| *b != b'\n' && *b != b'\r')
                    {
                        if !self.input[self.pos].is_ascii_whitespace() {
                            return self.syntax(
                                self.line,
                                self.column,
                                "text-field delimiter has trailing data",
                            );
                        }
                        self.bump();
                    }
                    let text = String::from_utf8(self.input[value_start..value_end].to_vec())?;
                    return Ok(Some(Token { text, line, column }));
                }
                self.bump();
            }
        }

        if first == b'\'' || first == b'"' {
            let quote = first;
            let value_start = self.pos;
            loop {
                let Some(&b) = self.input.get(self.pos) else {
                    return self.syntax(line, column, "unterminated quoted value");
                };
                if b == quote
                    && self
                        .input
                        .get(self.pos + 1)
                        .is_none_or(|next| next.is_ascii_whitespace())
                {
                    let value_end = self.pos;
                    self.bump();
                    let text = String::from_utf8(self.input[value_start..value_end].to_vec())?;
                    return Ok(Some(Token { text, line, column }));
                }
                self.bump();
            }
        }

        while self
            .input
            .get(self.pos)
            .is_some_and(|b| !b.is_ascii_whitespace())
        {
            self.bump();
        }
        let text = String::from_utf8(self.input[start..self.pos].to_vec())?;
        Ok(Some(Token { text, line, column }))
    }
}

#[derive(Default)]
struct CifDocument {
    scalars: HashMap<String, String>,
    loops: Vec<CifLoop>,
}

struct CifLoop {
    tags: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn is_control(text: &str) -> bool {
    let lower = text.to_ascii_lowercase();
    text.starts_with('_')
        || lower == "loop_"
        || lower == "stop_"
        || lower == "global_"
        || lower.starts_with("data_")
        || lower.starts_with("save_")
}

fn wanted_tag(tag: &str) -> bool {
    tag.starts_with("_atom_site.") || tag.starts_with("_cell.") || tag.starts_with("_struct_conn.")
}

fn parse_document(input: String) -> Result<CifDocument, CifHandlerError> {
    let mut lexer = Lexer::new(input);
    let mut pending: Option<Token> = None;
    let mut doc = CifDocument::default();
    loop {
        let token = match pending.take() {
            Some(t) => t,
            None => match lexer.next_token()? {
                Some(t) => t,
                None => break,
            },
        };
        let lower = token.text.to_ascii_lowercase();
        if lower == "loop_" {
            let mut tags = Vec::new();
            while let Some(t) = lexer.next_token()? {
                if t.text.starts_with('_') {
                    tags.push(t.text.to_ascii_lowercase());
                } else {
                    pending = Some(t);
                    break;
                }
            }
            if tags.is_empty() {
                return Err(CifHandlerError::Syntax {
                    line: token.line,
                    column: token.column,
                    message: "loop_ has no data names".into(),
                });
            }
            let keep = tags.iter().any(|t| wanted_tag(t));
            let mut values = Vec::new();
            loop {
                let t = match pending.take() {
                    Some(t) => t,
                    None => match lexer.next_token()? {
                        Some(t) => t,
                        None => break,
                    },
                };
                if is_control(&t.text) {
                    pending = Some(t);
                    break;
                }
                if keep {
                    values.push(t.text);
                }
            }
            if keep {
                if values.len() % tags.len() != 0 {
                    return Err(CifHandlerError::Syntax {
                        line: token.line,
                        column: token.column,
                        message: format!(
                            "loop has {} values for {} columns",
                            values.len(),
                            tags.len()
                        ),
                    });
                }
                let rows = values.chunks(tags.len()).map(<[String]>::to_vec).collect();
                doc.loops.push(CifLoop { tags, rows });
            }
        } else if token.text.starts_with('_') {
            let value = lexer.next_token()?.ok_or_else(|| CifHandlerError::Syntax {
                line: token.line,
                column: token.column,
                message: "data name has no value".into(),
            })?;
            if is_control(&value.text) {
                return Err(CifHandlerError::Syntax {
                    line: value.line,
                    column: value.column,
                    message: "data name has no value".into(),
                });
            }
            if wanted_tag(&lower) {
                doc.scalars.insert(lower, value.text);
            }
        }
    }
    Ok(doc)
}

fn missing(v: &str) -> bool {
    v == "." || v == "?"
}

fn column(tags: &[String], name: &str) -> Option<usize> {
    tags.iter().position(|t| t == name)
}

fn value(row: &[String], col: Option<usize>) -> Option<&str> {
    col.and_then(|i| row.get(i))
        .map(String::as_str)
        .filter(|v| !missing(v))
}

fn required<'a>(
    row: &'a [String],
    col: Option<usize>,
    field: &'static str,
) -> Result<&'a str, CifHandlerError> {
    value(row, col).ok_or(CifHandlerError::MissingField(field))
}

fn parse_value<T: std::str::FromStr>(s: &str, field: &'static str) -> Result<T, CifHandlerError> {
    s.parse().map_err(|_| CifHandlerError::InvalidValue {
        field,
        value: s.into(),
    })
}

fn check_atom_str(s: &str, field: &'static str) -> Result<(), CifHandlerError> {
    if s.is_ascii() && s.len() <= 8 {
        Ok(())
    } else {
        Err(CifHandlerError::StringTooLong {
            field,
            value: s.into(),
        })
    }
}

fn chain_char(id: &str) -> Result<char, CifHandlerError> {
    let mut chars = id.chars();
    if let (Some(c), None) = (chars.next(), chars.next()) {
        return Ok(c);
    }
    if let Some(last) = id.rsplit('.').next() {
        let mut chars = last.chars();
        if let (Some(c), None) = (chars.next(), chars.next()) {
            return Ok(c);
        }
    }
    Err(CifHandlerError::UnsupportedChain(id.into()))
}

#[derive(Clone, PartialEq, Eq)]
struct AtomIdentity {
    label_asym: String,
    label_seq: String,
    label_comp: String,
    label_atom: String,
    alt: String,
}

fn identity_key(id: &AtomIdentity) -> String {
    format!(
        "{}\x1f{}\x1f{}\x1f{}\x1f{}",
        id.label_asym, id.label_seq, id.label_comp, id.label_atom, id.alt
    )
}

fn bond_order(s: Option<&str>) -> Result<BondOrder, CifHandlerError> {
    match s.map(str::to_ascii_lowercase).as_deref() {
        None => Ok(BondOrder::Unspecified),
        Some("sing") => Ok(BondOrder::Single),
        Some("doub") => Ok(BondOrder::Double),
        Some("trip") => Ok(BondOrder::Triple),
        Some("arom") => Ok(BondOrder::Aromatic),
        Some(v) => Err(CifHandlerError::UnsupportedBondOrder(v.into())),
    }
}

fn parse_cell(doc: &CifDocument) -> Option<PeriodicBox> {
    let get = |name: &str| {
        doc.scalars
            .get(name)
            .filter(|v| !missing(v))?
            .parse::<Float>()
            .ok()
    };
    PeriodicBox::from_vectors_angles(
        get("_cell.length_a")? * ANGSTROM_TO_NM,
        get("_cell.length_b")? * ANGSTROM_TO_NM,
        get("_cell.length_c")? * ANGSTROM_TO_NM,
        get("_cell.angle_alpha")?,
        get("_cell.angle_beta")?,
        get("_cell.angle_gamma")?,
    )
    .ok()
}

fn parse_frames(doc: &CifDocument) -> Result<VecDeque<(Topology, State)>, CifHandlerError> {
    let atom_loop = doc
        .loops
        .iter()
        .find(|l| l.tags.iter().any(|t| t.starts_with("_atom_site.")))
        .ok_or(CifHandlerError::Empty)?;
    if atom_loop.rows.is_empty() {
        return Err(CifHandlerError::Empty);
    }

    macro_rules! col {
        ($name:literal) => {
            column(&atom_loop.tags, concat!("_atom_site.", $name))
        };
    }
    let group = col!("group_pdb");
    let label_atom = col!("label_atom_id");
    let label_comp = col!("label_comp_id");
    let label_asym = col!("label_asym_id");
    let label_seq = col!("label_seq_id");
    let alt = col!("label_alt_id");
    let ins = col!("pdbx_pdb_ins_code");
    let auth_atom = col!("auth_atom_id");
    let auth_comp = col!("auth_comp_id");
    let auth_asym = col!("auth_asym_id");
    let auth_seq = col!("auth_seq_id");
    let elem = col!("type_symbol");
    let x = col!("cartn_x");
    let y = col!("cartn_y");
    let z = col!("cartn_z");
    let occupancy = col!("occupancy");
    let bfactor = col!("b_iso_or_equiv");
    let formal_charge = col!("pdbx_formal_charge");
    let model = col!("pdbx_pdb_model_num");

    let mut model_order = Vec::<i32>::new();
    let mut grouped: HashMap<i32, Vec<&[String]>> = HashMap::new();
    for row in &atom_loop.rows {
        if let Some(g) = value(row, group)
            && !g.eq_ignore_ascii_case("ATOM")
            && !g.eq_ignore_ascii_case("HETATM")
        {
            continue;
        }
        let m = value(row, model)
            .map(|s| parse_value(s, "_atom_site.pdbx_PDB_model_num"))
            .transpose()?
            .unwrap_or(1);
        match grouped.entry(m) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                model_order.push(m);
                entry.insert(vec![row]);
            }
            std::collections::hash_map::Entry::Occupied(mut entry) => {
                entry.get_mut().push(row);
            }
        }
    }
    if model_order.is_empty() {
        return Err(CifHandlerError::Empty);
    }

    let pbox = parse_cell(doc);
    let mut chain_sources: HashMap<char, String> = HashMap::new();
    let mut first_ids: Option<Vec<AtomIdentity>> = None;
    let mut frames = VecDeque::new();
    for m in model_order {
        let rows = &grouped[&m];
        let mut atoms = Vec::with_capacity(rows.len());
        let mut coords = Vec::with_capacity(rows.len());
        let mut ids = Vec::with_capacity(rows.len());
        let mut last_residue: Option<String> = None;
        let mut resindex = 0usize;
        for row in rows {
            let label_atom_v = required(row, label_atom.or(auth_atom), "_atom_site.label_atom_id")?;
            let label_comp_v = required(row, label_comp.or(auth_comp), "_atom_site.label_comp_id")?;
            let label_asym_v = value(row, label_asym.or(auth_asym)).unwrap_or(" ");
            let label_seq_v = value(row, label_seq.or(auth_seq)).unwrap_or("0");
            let atom_name = value(row, auth_atom).unwrap_or(label_atom_v);
            let resname = value(row, auth_comp).unwrap_or(label_comp_v);
            let chain_source = value(row, auth_asym).unwrap_or(label_asym_v);
            let resid_s = value(row, auth_seq).unwrap_or(label_seq_v);
            check_atom_str(atom_name, "atom name")?;
            check_atom_str(resname, "residue name")?;
            let chain = if chain_source == " " {
                ' '
            } else {
                chain_char(chain_source)?
            };
            if let Some(previous) = chain_sources.get(&chain) {
                if previous != chain_source {
                    return Err(CifHandlerError::ChainCollision {
                        first: previous.clone(),
                        second: chain_source.into(),
                        mapped: chain,
                    });
                }
            } else {
                chain_sources.insert(chain, chain_source.into());
            }
            let resid: i32 = parse_value(resid_s, "_atom_site.auth_seq_id")?;
            let residue_key = format!(
                "{}\x1f{}\x1f{}\x1f{}",
                chain_source,
                resid_s,
                value(row, ins).unwrap_or(""),
                resname
            );
            if last_residue
                .as_ref()
                .is_some_and(|last| last != &residue_key)
            {
                resindex += 1;
            }
            last_residue = Some(residue_key);

            let px = parse_value::<Float>(
                required(row, x, "_atom_site.Cartn_x")?,
                "_atom_site.Cartn_x",
            )? * ANGSTROM_TO_NM;
            let py = parse_value::<Float>(
                required(row, y, "_atom_site.Cartn_y")?,
                "_atom_site.Cartn_y",
            )? * ANGSTROM_TO_NM;
            let pz = parse_value::<Float>(
                required(row, z, "_atom_site.Cartn_z")?,
                "_atom_site.Cartn_z",
            )? * ANGSTROM_TO_NM;
            coords.push(Pos::new(px, py, pz));

            let occ = value(row, occupancy)
                .map(|s| parse_value(s, "_atom_site.occupancy"))
                .transpose()?
                .unwrap_or(1.0);
            let bf = value(row, bfactor)
                .map(|s| parse_value(s, "_atom_site.B_iso_or_equiv"))
                .transpose()?
                .unwrap_or(0.0);
            let mut atom = Atom::new()
                .with_name(atom_name)
                .with_resname(resname)
                .with_resid(resid)
                .with_resindex(resindex)
                .with_chain(chain)
                .with_occupancy(occ)
                .with_bfactor(bf);
            if let Some(s) = value(row, formal_charge) {
                atom = atom.with_formal_charge(parse_value(s, "_atom_site.pdbx_formal_charge")?);
            }
            let atomic_number = value(row, elem).map(atomic_number_from_symbol).unwrap_or(0);
            atom = if atomic_number != 0 {
                atom.with_atomic_number(atomic_number)
                    .with_mass(ELEMENT_MASS[atomic_number as usize] as Float)
            } else {
                atom.guess()
            };
            atoms.push(atom);
            ids.push(AtomIdentity {
                label_asym: label_asym_v.into(),
                label_seq: label_seq_v.into(),
                label_comp: label_comp_v.into(),
                label_atom: label_atom_v.into(),
                alt: value(row, alt).unwrap_or("").into(),
            });
        }
        if let Some(first) = &first_ids {
            if first != &ids {
                return Err(CifHandlerError::ModelMismatch { model: m });
            }
        } else {
            first_ids = Some(ids.clone());
        }
        let top = Topology {
            atoms: atoms.into_iter().collect(),
            ..Default::default()
        };
        frames.push_back((
            top,
            State {
                coords,
                time: 0.0,
                pbox: pbox.clone(),
                ..Default::default()
            },
        ));
    }

    // `_struct_conn` describes both bonds and non-bond contacts.  Only topology-like
    // connection types are imported.
    if let (Some(conn), Some(ids)) = (
        doc.loops
            .iter()
            .find(|l| l.tags.iter().any(|t| t.starts_with("_struct_conn."))),
        first_ids.as_ref(),
    ) {
        let atom_map: HashMap<String, usize> = ids
            .iter()
            .enumerate()
            .map(|(i, id)| (identity_key(id), i))
            .collect();
        let c = |name: &'static str| column(&conn.tags, name);
        let ty = c("_struct_conn.conn_type_id");
        let order = c("_struct_conn.pdbx_value_order");
        let sym1 = c("_struct_conn.ptnr1_symmetry");
        let sym2 = c("_struct_conn.ptnr2_symmetry");
        let mut resolved: HashMap<[usize; 2], BondOrder> = HashMap::new();
        for row in &conn.rows {
            let Some(kind) = value(row, ty) else { continue };
            if !matches!(
                kind.to_ascii_lowercase().as_str(),
                "covale" | "disulf" | "metalc"
            ) {
                continue;
            }
            if [value(row, sym1), value(row, sym2)]
                .into_iter()
                .flatten()
                .any(|sym| sym != "1_555")
            {
                continue;
            }
            let partner = |n: usize| -> Option<String> {
                let get = |suffix: &str| {
                    value(
                        row,
                        column(&conn.tags, &format!("_struct_conn.ptnr{n}_{suffix}")),
                    )
                };
                let id = AtomIdentity {
                    label_asym: get("label_asym_id")?.into(),
                    label_seq: get("label_seq_id")?.into(),
                    label_comp: get("label_comp_id")?.into(),
                    label_atom: get("label_atom_id")?.into(),
                    alt: get("label_alt_id").unwrap_or("").into(),
                };
                Some(identity_key(&id))
            };
            let (Some(k1), Some(k2)) = (partner(1), partner(2)) else {
                continue;
            };
            let (Some(&i), Some(&j)) = (atom_map.get(&k1), atom_map.get(&k2)) else {
                continue;
            };
            if i == j {
                continue;
            }
            let mut pair = [i, j];
            pair.sort();
            let new_order = bond_order(value(row, order))?;
            match resolved.get_mut(&pair) {
                None => {
                    resolved.insert(pair, new_order);
                }
                Some(old) if *old == new_order => {}
                Some(old) if *old == BondOrder::Unspecified => *old = new_order,
                Some(_) if new_order == BondOrder::Unspecified => {}
                Some(_) => return Err(CifHandlerError::ConflictingBondOrder(pair[0], pair[1])),
            }
        }
        let mut bonds: Vec<_> = resolved.into_iter().collect();
        bonds.sort_by_key(|(pair, _)| *pair);
        for (top, _) in &mut frames {
            top.bonds = bonds
                .iter()
                .map(|([i, j], order)| Bond::with_order(*i, *j, *order))
                .collect();
        }
    }
    Ok(frames)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const CHARGED_BONDED: &str = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
ATOM 1 N N LIG A 1 10.0 20.0 30.0 0.5 12.0 1
HETATM 2 O O LIG A 1 11.0 20.0 30.0 1.0 13.0 -1
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order
c1 covale LIG A 1 N LIG A 1 O doub
c2 hydrog LIG A 1 N LIG A 1 O sing
#
"#;

    fn memory_handler(text: &str) -> CifFileHandler {
        CifFileHandler::from_source(DynSource(Box::new(Cursor::new(text.as_bytes().to_vec()))))
            .unwrap()
    }

    #[test]
    fn reads_supplied_cif() {
        let mut h = CifFileHandler::open("tests/system_12A.cif").unwrap();
        let (top, state) = h.read().unwrap();
        assert_eq!(top.len(), 1628);
        assert_eq!(state.coords.len(), 1628);
        let first = top.get_atom(0).unwrap();
        assert_eq!(first.get_name(), "N");
        assert_eq!(first.get_resname(), "PRO");
        assert_eq!(first.get_resid(), 53);
        assert_eq!(first.get_chain(), 'A');
        assert_eq!(first.get_atomic_number(), 7);
        assert!((state.coords[0].x - 6.07).abs() < 1e-5);
        assert!((state.coords[0].y - 3.4326).abs() < 1e-5);
        assert!((state.coords[0].z + 16.6363).abs() < 1e-5);

        let residue_one: Vec<_> = top
            .iter_atoms()
            .filter(|a| matches!(a.get_chain(), 'E' | 'F' | 'G') && a.get_resid() == 1)
            .map(|a| (a.get_chain(), a.get_resindex()))
            .collect();
        let e = residue_one.iter().find(|(c, _)| *c == 'E').unwrap().1;
        let f = residue_one.iter().find(|(c, _)| *c == 'F').unwrap().1;
        let g = residue_one.iter().find(|(c, _)| *c == 'G').unwrap().1;
        assert_ne!(e, f);
        assert_ne!(f, g);
    }

    #[test]
    fn reads_formal_charges_and_bond_order() {
        let (top, state) = memory_handler(CHARGED_BONDED).read().unwrap();
        assert_eq!(top.len(), 2);
        assert_eq!(top.get_atom(0).unwrap().get_formal_charge(), Some(1));
        assert_eq!(top.get_atom(1).unwrap().get_formal_charge(), Some(-1));
        assert_eq!(
            top.bonds.len(),
            1,
            "hydrogen-bond contact is not a topology bond"
        );
        assert!(top.bonds.has_orders());
        assert_eq!(top.bonds.get(0).unwrap().pair(), [0, 1]);
        assert_eq!(top.bonds.get(0).unwrap().order(), BondOrder::Double);
        assert_eq!(state.coords[0], Pos::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn write_read_preserves_charge_and_order() {
        let (top, state) = memory_handler(CHARGED_BONDED).read().unwrap();
        let path = concat!(env!("OUT_DIR"), "/cif_charge_order_roundtrip.cif");
        let mut writer = CifFileHandler::create(path).unwrap();
        writer.write(&System::new(top, state).unwrap()).unwrap();
        drop(writer);

        let (got, _) = CifFileHandler::open(path).unwrap().read().unwrap();
        assert_eq!(got.get_atom(0).unwrap().get_formal_charge(), Some(1));
        assert_eq!(got.get_atom(1).unwrap().get_formal_charge(), Some(-1));
        assert_eq!(got.bonds.len(), 1);
        assert_eq!(got.bonds.get(0).unwrap().order(), BondOrder::Double);
    }

    #[test]
    fn missing_charge_and_unknown_order_stay_unset() {
        let text = CHARGED_BONDED
            .replace("_atom_site.pdbx_formal_charge\n", "")
            .replace(" 12.0 1\n", " 12.0\n")
            .replace(" 13.0 -1\n", " 13.0\n")
            .replace(
                "c1 covale LIG A 1 N LIG A 1 O doub",
                "c1 covale LIG A 1 N LIG A 1 O ?",
            );
        let (top, _) = memory_handler(&text).read().unwrap();
        assert_eq!(top.get_atom(0).unwrap().get_formal_charge(), None);
        assert!(!top.bonds.has_orders());
        assert_eq!(top.bonds.get(0).unwrap().order(), BondOrder::Unspecified);
    }

    #[test]
    fn quoted_and_multiline_values_are_tokenized() {
        let text = r#"
data_test
_struct.title
;
This ignored value has spaces and reserved words such as loop_.
;
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM C 'C 1' "LIG" A 1 1 2 3
"#;
        let (top, _) = memory_handler(text).read().unwrap();
        assert_eq!(top.get_atom(0).unwrap().get_name(), "C 1");
    }

    #[test]
    fn reads_models_and_cell_through_public_memory_dispatch() {
        let text = r#"
data_test
_cell.length_a 10
_cell.length_b 20
_cell.length_c 30
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM C C MOL A 1 1 2 3 1
ATOM C C MOL A 1 4 5 6 2
"#;
        let mut h = FileHandler::from_reader("cif", Cursor::new(text.as_bytes().to_vec())).unwrap();
        let top = h.read_topology().unwrap();
        let first = h.read_state().unwrap();
        let second = h.read_state().unwrap();
        assert_eq!(top.len(), 1);
        assert!((first.coords[0] - Pos::new(0.1, 0.2, 0.3)).norm() < 1e-6);
        assert!((second.coords[0] - Pos::new(0.4, 0.5, 0.6)).norm() < 1e-6);
        let (lengths, angles) = first.pbox.unwrap().to_vectors_angles();
        assert!((lengths[0] - 1.0).abs() < 1e-6);
        assert!((lengths[1] - 2.0).abs() < 1e-6);
        assert!((lengths[2] - 3.0).abs() < 1e-6);
        assert_eq!(angles, Vector3f::new(90.0, 90.0, 90.0));
        assert!(matches!(
            h.read_state().unwrap_err().kind(),
            FileFormatError::Eof
        ));
    }

    #[test]
    fn rejects_unrepresentable_quadruple_order() {
        let text = CHARGED_BONDED.replace("doub", "quad");
        let err = memory_handler(&text).read().unwrap_err();
        assert!(matches!(
            err,
            FileFormatError::Cif(CifHandlerError::UnsupportedBondOrder(_))
        ));
    }
}

fn cif_value(value: &str) -> Result<String, CifHandlerError> {
    let lower = value.to_ascii_lowercase();
    let bare = !value.is_empty()
        && !value.bytes().any(|b| b.is_ascii_whitespace() || b == b'#')
        && !value.starts_with('_')
        && !matches!(lower.as_str(), "loop_" | "stop_" | "global_")
        && !lower.starts_with("data_")
        && !lower.starts_with("save_");
    if bare {
        Ok(value.into())
    } else if !value.contains('\'') {
        Ok(format!("'{value}'"))
    } else if !value.contains('"') {
        Ok(format!("\"{value}\""))
    } else {
        Err(CifHandlerError::InvalidValue {
            field: "text output",
            value: value.into(),
        })
    }
}

impl CifFileHandler {
    pub(crate) fn from_source(src: DynSource) -> Result<Self, FileFormatError> {
        Ok(Self {
            reader: Some(BufReader::new(src)),
            writer: None,
            frames: None,
            stored_topology: None,
            stored_state: None,
            at_least_one_state_read: false,
        })
    }

    fn ensure_frames(&mut self) -> Result<(), CifHandlerError> {
        if self.frames.is_none() {
            let mut bytes = Vec::new();
            self.reader
                .as_mut()
                .ok_or(CifHandlerError::Empty)?
                .read_to_end(&mut bytes)?;
            self.frames = Some(parse_frames(&parse_document(String::from_utf8(bytes)?)?)?);
        }
        Ok(())
    }
}

impl FileFormatHandler for CifFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let file = File::open(fname).map_err(CifHandlerError::OpenRead)?;
        Self::from_source(DynSource(Box::new(file)))
    }

    fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        Ok(Self {
            reader: None,
            writer: Some(BufWriter::new(
                File::create(fname).map_err(CifHandlerError::OpenWrite)?,
            )),
            frames: None,
            stored_topology: None,
            stored_state: None,
            at_least_one_state_read: false,
        })
    }

    fn read(&mut self) -> Result<(Topology, State), FileFormatError> {
        self.ensure_frames()?;
        match self.frames.as_mut().and_then(VecDeque::pop_front) {
            Some(frame) => {
                self.at_least_one_state_read = true;
                Ok(frame)
            }
            None if self.at_least_one_state_read => Err(FileFormatError::Eof),
            None => Err(CifHandlerError::Empty.into()),
        }
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
        let w = self.writer.as_mut().ok_or(FileFormatError::NotWritable)?;
        writeln!(w, "data_molar\n#")?;
        if let Some(b) = data.get_box() {
            let (lengths, angles) = b.to_vectors_angles();
            writeln!(w, "_cell.length_a {:.6}", lengths[0] * NM_TO_ANGSTROM)?;
            writeln!(w, "_cell.length_b {:.6}", lengths[1] * NM_TO_ANGSTROM)?;
            writeln!(w, "_cell.length_c {:.6}", lengths[2] * NM_TO_ANGSTROM)?;
            writeln!(w, "_cell.angle_alpha {:.6}", angles[0])?;
            writeln!(w, "_cell.angle_beta {:.6}", angles[1])?;
            writeln!(w, "_cell.angle_gamma {:.6}\n#", angles[2])?;
        }
        writeln!(w, "loop_")?;
        for tag in [
            "group_PDB",
            "id",
            "type_symbol",
            "label_atom_id",
            "label_alt_id",
            "label_comp_id",
            "label_asym_id",
            "label_seq_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "occupancy",
            "B_iso_or_equiv",
            "pdbx_formal_charge",
            "auth_seq_id",
            "auth_comp_id",
            "auth_asym_id",
            "auth_atom_id",
            "pdbx_PDB_model_num",
        ] {
            writeln!(w, "_atom_site.{tag}")?;
        }
        for (i, (atom, pos)) in data.iter_atoms_dyn().zip(data.iter_pos_dyn()).enumerate() {
            let name = cif_value(atom.name())?;
            let resname = cif_value(atom.resname())?;
            let chain = cif_value(&atom.get_chain().to_string())?;
            let charge = atom
                .get_formal_charge()
                .map_or_else(|| "?".into(), |v| v.to_string());
            writeln!(
                w,
                "ATOM {} {} {} . {} {} {} {:.6} {:.6} {:.6} {:.6} {:.6} {} {} {} {} {} 1",
                i + 1,
                element_symbol(atom.get_atomic_number()),
                name,
                resname,
                chain,
                atom.get_resid(),
                pos.x * NM_TO_ANGSTROM,
                pos.y * NM_TO_ANGSTROM,
                pos.z * NM_TO_ANGSTROM,
                atom.get_occupancy(),
                atom.get_bfactor(),
                charge,
                atom.get_resid(),
                resname,
                chain,
                name
            )?;
        }
        if data.num_bonds() != 0 {
            writeln!(w, "#\nloop_")?;
            for tag in [
                "id",
                "conn_type_id",
                "ptnr1_label_comp_id",
                "ptnr1_label_asym_id",
                "ptnr1_label_seq_id",
                "ptnr1_label_atom_id",
                "ptnr2_label_comp_id",
                "ptnr2_label_asym_id",
                "ptnr2_label_seq_id",
                "ptnr2_label_atom_id",
                "pdbx_value_order",
            ] {
                writeln!(w, "_struct_conn.{tag}")?;
            }
            let atoms: Vec<_> = data.iter_atoms_dyn().collect();
            for (i, bond) in data.iter_bonds_dyn().enumerate() {
                let a = atoms[bond.i1()];
                let b = atoms[bond.i2()];
                let order = match bond.order() {
                    BondOrder::Unspecified => "?",
                    BondOrder::Single => "sing",
                    BondOrder::Double => "doub",
                    BondOrder::Triple => "trip",
                    BondOrder::Aromatic => "arom",
                };
                writeln!(
                    w,
                    "mol{} covale {} {} {} {} {} {} {} {} {}",
                    i + 1,
                    cif_value(a.resname())?,
                    cif_value(&a.get_chain().to_string())?,
                    a.get_resid(),
                    cif_value(a.name())?,
                    cif_value(b.resname())?,
                    cif_value(&b.get_chain().to_string())?,
                    b.get_resid(),
                    cif_value(b.name())?,
                    order
                )?;
            }
        }
        writeln!(w, "#")?;
        Ok(())
    }
}
