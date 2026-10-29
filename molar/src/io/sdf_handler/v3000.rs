//! V3000 connection-table reader. Coordinates, elements, masses, formal charges,
//! and ordinary bond orders are retained. Stereo, query, Sgroup, and collection
//! annotations have no representation in Molar's topology and are not retained.
use super::*;
use std::collections::{HashMap, HashSet};

fn bad(message: impl Into<String>) -> SdfHandlerError {
    SdfHandlerError::BadV3000(message.into())
}

/// Read a logical entry, joining continuation lines without inserting whitespace.
fn entry(buf: &mut impl BufRead) -> Result<String, SdfHandlerError> {
    let mut result = String::new();
    loop {
        let line = next_line(buf)?.ok_or_else(|| bad("unexpected end of file"))?;
        let body = line
            .trim_end()
            .strip_prefix("M  V30 ")
            .ok_or_else(|| bad(format!("expected M  V30 entry, got {:?}", line.trim_end())))?;
        if let Some(body) = body.strip_suffix('-') {
            result.push_str(body);
        } else {
            result.push_str(body);
            return Ok(result.trim().to_owned());
        }
    }
}

/// Keep quoted strings and parenthesized lists together, including doubled quotes.
fn tokens(line: &str) -> Result<Vec<String>, SdfHandlerError> {
    let mut out = Vec::new();
    let mut word = String::new();
    let mut quoted = false;
    let mut depth = 0usize;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' if quoted && chars.peek() == Some(&'"') => {
                word.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            '(' if !quoted => {
                depth += 1;
                word.push(c);
            }
            ')' if !quoted => {
                depth = depth
                    .checked_sub(1)
                    .ok_or_else(|| bad("unmatched closing parenthesis"))?;
                word.push(c);
            }
            c if c.is_whitespace() && !quoted && depth == 0 => {
                if !word.is_empty() {
                    out.push(std::mem::take(&mut word));
                }
            }
            _ => word.push(c),
        }
    }
    if quoted || depth != 0 {
        return Err(bad("unterminated string or list"));
    }
    if !word.is_empty() {
        out.push(word);
    }
    Ok(out)
}

fn integer(value: &str, field: &str) -> Result<usize, SdfHandlerError> {
    value
        .parse()
        .map_err(|source| SdfHandlerError::V3000Integer(field.to_owned(), source))
}

fn number(value: &str, field: &str) -> Result<Float, SdfHandlerError> {
    let value: Float = value
        .parse()
        .map_err(|source| SdfHandlerError::V3000Number(field.to_owned(), source))?;
    if !value.is_finite() {
        return Err(bad(format!("{field} must be finite")));
    }
    Ok(value)
}

fn expect(buf: &mut impl BufRead, expected: &str) -> Result<(), SdfHandlerError> {
    let line = entry(buf)?;
    if line != expected {
        return Err(bad(format!("expected {expected}, got {line:?}")));
    }
    Ok(())
}

fn skip_block(buf: &mut impl BufRead, name: &str) -> Result<(), SdfHandlerError> {
    let mut blocks = vec![name.to_owned()];
    while let Some(expected) = blocks.last() {
        let line = entry(buf)?;
        if let Some(name) = line.strip_prefix("BEGIN ") {
            blocks.push(name.to_owned());
        } else if let Some(name) = line.strip_prefix("END ") {
            if name != expected {
                return Err(bad(format!("expected END {expected}, got {line:?}")));
            }
            blocks.pop();
        } else {
            tokens(&line)?;
        }
    }
    Ok(())
}

pub(super) fn read(
    buf: &mut impl BufRead,
) -> Result<(Vec<Atom>, Vec<Pos>, Vec<Bond>), SdfHandlerError> {
    expect(buf, "BEGIN CTAB")?;
    let counts = tokens(&entry(buf)?)?;
    if counts.len() < 6 || counts[0] != "COUNTS" {
        return Err(bad("expected COUNTS na nb nsg n3d chiral"));
    }
    let natoms = integer(&counts[1], "atom count")?;
    let nbonds = integer(&counts[2], "bond count")?;
    integer(&counts[3], "Sgroup count")?;
    integer(&counts[4], "3D constraint count")?;
    integer(&counts[5], "chiral flag")?;
    if natoms == 0 {
        return Err(bad("atom count must be positive"));
    }
    // Grow from parsed records rather than allocating from untrusted counts.
    let mut atoms = Vec::new();
    let mut coords = Vec::new();
    let mut atom_ids = HashMap::new();
    expect(buf, "BEGIN ATOM")?;
    for index in 0..natoms {
        let fields = tokens(&entry(buf)?)?;
        if fields.len() < 6 {
            return Err(bad(format!("incomplete atom entry {}", index + 1)));
        }
        let id = integer(&fields[0], "atom ID")?;
        if id == 0 || atom_ids.insert(id, index).is_some() {
            return Err(bad(format!("invalid or duplicate atom ID {id}")));
        }
        let elem = fields[1].as_str();
        // Query atom lists and template atoms are not ordinary elements.
        let z = atomic_number_from_symbol(elem);
        if z == 0 {
            return Err(bad(format!("unsupported atom type {elem:?} at atom {id}")));
        }
        let mut atom = Atom::new()
            .with_name(elem)
            .with_resname("MOL")
            .with_resid(1)
            .with_chain('A')
            .with_atomic_number(z)
            .with_mass(ELEMENT_MASS[z as usize] as Float);
        coords.push(
            Pos::new(
                number(&fields[2], "x coordinate")?,
                number(&fields[3], "y coordinate")?,
                number(&fields[4], "z coordinate")?,
            ) * ANGSTROM_TO_NM,
        );
        integer(&fields[5], "atom mapping")?;
        for property in &fields[6..] {
            let (key, value) = property
                .split_once('=')
                .ok_or_else(|| bad(format!("invalid atom property {property:?}")))?;
            match key {
                "CHG" => {
                    atom.formal_charge = Some(value.parse().map_err(|source| {
                        SdfHandlerError::V3000Integer(format!("CHG at atom {id}"), source)
                    })?)
                }
                "MASS" => {
                    let mass = number(value, "MASS")?;
                    if mass <= 0.0 {
                        return Err(bad("MASS must be positive"));
                    }
                    atom.mass = mass;
                }
                "CFG" => {
                    integer(value, "atom CFG")?;
                }
                _ => (), // Annotations outside the topology data model.
            }
        }
        atoms.push(atom);
    }
    expect(buf, "END ATOM")?;
    let mut bonds = Vec::new();
    let mut bond_ids = HashSet::new();
    let mut line = entry(buf)?;
    if line == "BEGIN BOND" {
        for index in 0..nbonds {
            let fields = tokens(&entry(buf)?)?;
            if fields.len() < 4 {
                return Err(bad(format!("incomplete bond entry {}", index + 1)));
            }
            let id = integer(&fields[0], "bond ID")?;
            if id == 0 || !bond_ids.insert(id) {
                return Err(bad(format!("invalid or duplicate bond ID {id}")));
            }
            let order = match integer(&fields[1], "bond order")? {
                1 => BondOrder::Single,
                2 => BondOrder::Double,
                3 => BondOrder::Triple,
                4 => BondOrder::Aromatic,
                5..=8 => BondOrder::Unspecified,
                value => return Err(bad(format!("unsupported bond order {value} at bond {id}"))),
            };
            let endpoint = |s: &str| -> Result<usize, SdfHandlerError> {
                let atom = integer(s, "bond endpoint")?;
                atom_ids
                    .get(&atom)
                    .copied()
                    .ok_or_else(|| bad(format!("bond {id} references missing atom {atom}")))
            };
            let a = endpoint(&fields[2])?;
            let b = endpoint(&fields[3])?;
            if a == b {
                return Err(bad(format!(
                    "bond {id} connects atom {} to itself",
                    fields[2]
                )));
            }
            for property in &fields[4..] {
                let (key, value) = property
                    .split_once('=')
                    .ok_or_else(|| bad(format!("invalid bond property {property:?}")))?;
                match key {
                    "CFG" => {
                        integer(value, "bond CFG")?;
                    }
                    "ENDPTS" | "ATTACH" => {
                        return Err(bad(format!("multicenter bond {id} is not supported")));
                    }
                    _ => (),
                }
            }
            bonds.push(Bond::with_order(a, b, order));
        }
        expect(buf, "END BOND")?;
        line = entry(buf)?;
    } else if nbonds != 0 {
        return Err(bad(format!("expected BEGIN BOND, got {line:?}")));
    }
    while line != "END CTAB" {
        match line.as_str() {
            "BEGIN COLLECTION" => skip_block(buf, "COLLECTION")?,
            "BEGIN SGROUP" => skip_block(buf, "SGROUP")?,
            "BEGIN OBJ3D" => skip_block(buf, "OBJ3D")?,
            _ => return Err(bad(format!("unsupported entry {line:?}"))),
        }
        line = entry(buf)?;
    }
    // Require the molfile terminator before the common SDF data-field reader runs.
    loop {
        let line = next_line(buf)?.ok_or_else(|| bad("missing M  END"))?;
        match line.trim_end() {
            "M  END" => break,
            "M  V30 BEGIN COLLECTION" => skip_block(buf, "COLLECTION")?,
            line => return Err(bad(format!("expected M  END, got {line:?}"))),
        }
    }
    Ok((atoms, coords, bonds))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const MOL: &str = "v3000 test\n  molar\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 3 2 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 10 C 0 0 0 0 MASS=13\n\
M  V30 30 O 1.25 0 0 0 CHG=-1\n\
M  V30 20 \"Cl\" 0 1.5 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 15 2 10 30 CFG=3\n\
M  V30 35 1 10 20\n\
M  V30 END BOND\n\
M  V30 BEGIN COLLECTION\n\
M  V30 MDLV30/STERAC1 ATOMS=(2 10 -\n\
M  V30 30)\n\
M  V30 END COLLECTION\n\
M  V30 END CTAB\n\
M  END\n\
> <Catalog ID>\nTEST\n\n$$$$\n";

    fn reader(text: &str) -> SdfFileHandler {
        SdfFileHandler::from_source(DynSource(Box::new(Cursor::new(text.as_bytes().to_vec()))))
            .unwrap()
    }

    #[test]
    fn reads_coordinates_charges_masses_and_sparse_ids() {
        let mut h = reader(MOL);
        let (top, state) = h.read().unwrap();
        assert_eq!(top.atoms.len(), 3);
        assert_eq!(top.atoms.get(2).unwrap().get_atomic_number(), 17);
        assert_eq!(top.atoms.get(0).unwrap().get_mass(), 13.0);
        assert_eq!(top.atoms.get(1).unwrap().get_formal_charge(), Some(-1));
        assert!((state.coords[1].x - 0.125).abs() < 1e-6);
        assert_eq!(
            Bond::from(&top.bonds.get(0).unwrap()),
            Bond::with_order(0, 1, BondOrder::Double)
        );
        assert_eq!(
            Bond::from(&top.bonds.get(1).unwrap()),
            Bond::with_order(0, 2, BondOrder::Single)
        );
        assert!(matches!(h.read(), Err(FileFormatError::Eof)));
    }

    #[test]
    fn mixed_versions_and_separate_topology_state_reads() {
        let input = format!(
            "{}{MOL}{}",
            super::super::tests::ETHENE,
            super::super::tests::ETHENE
        );
        let mut h = reader(&input);
        assert_eq!(h.read().unwrap().0.atoms.len(), 6);
        assert_eq!(h.read_topology().unwrap().atoms.len(), 3);
        assert_eq!(h.read_state().unwrap().coords.len(), 3);
        assert_eq!(h.read().unwrap().0.atoms.len(), 6);
        assert!(matches!(h.read(), Err(FileFormatError::Eof)));
        assert_eq!(
            reader(MOL.split("M  END").next().unwrap())
                .read()
                .unwrap_err()
                .to_string(),
            "in sdf/mol format handler"
        );
        assert!(reader(MOL.split('>').next().unwrap()).read().is_ok());
    }

    #[test]
    fn reads_more_than_999_atoms_without_bonds() {
        let mut input = String::from(
            "large\n  molar\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1001 0 0 0 0\nM  V30 BEGIN ATOM\n",
        );
        for id in 1..=1001 {
            input.push_str(&format!("M  V30 {} C 0 0 0 0\n", id * 2));
        }
        input.push_str("M  V30 END ATOM\nM  V30 END CTAB\nM  END\n");
        let (top, state) = reader(&input).read().unwrap();
        assert_eq!(top.atoms.len(), 1001);
        assert_eq!(state.coords.len(), 1001);
        assert!(top.bonds.is_empty());
    }

    #[test]
    fn rejects_malformed_records_with_nested_parse_causes() {
        for (from, to) in [
            ("COUNTS 3 2", "COUNTS 4 2"),
            ("COUNTS 3 2", "COUNTS 3 1"),
            ("30 O", "10 O"),
            ("35 1", "15 1"),
            ("10 20\n", "10 999\n"),
            ("1.25", "NaN"),
            ("M  V30 END ATOM", "M  V30 END BOND"),
            ("M  V30 END COLLECTION", "M  V30 END CTAB"),
            ("CHG=-1", "CHG=bad"),
        ] {
            assert!(
                reader(&MOL.replace(from, to)).read().is_err(),
                "{from} -> {to}"
            );
        }
        let error = reader(&MOL.replace("1.25", "broken")).read().unwrap_err();
        let mut chain: &dyn std::error::Error = &error;
        while let Some(source) = chain.source() {
            chain = source;
        }
        assert_eq!(chain.to_string(), "invalid float literal");
    }

    #[test]
    fn joins_continuations_and_handles_quoted_properties() {
        let mut input = Cursor::new(
            b"M  V30 LABEL=\"abc-\nM  V30 def\" NAME=\"a \"\"quote\"\"\" LIST=(2 10 30)\n",
        );
        let line = entry(&mut input).unwrap();
        assert_eq!(
            tokens(&line).unwrap(),
            ["LABEL=abcdef", "NAME=a \"quote\"", "LIST=(2 10 30)"]
        );
        assert!(tokens("LABEL=\"unterminated").is_err());
        assert!(tokens("ATOMS=(2 1").is_err());
        assert!(tokens("ATOMS=2)").is_err());
    }
}
