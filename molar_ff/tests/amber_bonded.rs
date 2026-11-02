use molar::prelude::*;
use molar_ff::{FFType, amber::*};

#[test]
fn bundled_tables_and_ethanol() {
    for ff in [FFType::Gaff, FFType::Gaff2] {
        let _ = AmberParameters::builtin(ff);
        let sys = System::from_file("tests/data/gaff_ref/sdf/Ethanol.sdf").unwrap();
        let model = BondedModel::gaff(&sys, ff).unwrap();
        assert_eq!(model.bonds().len(), 8);
        assert_eq!(model.angles().len(), 13);
        assert_eq!(model.propers().len(), 12);
        assert!(model.impropers().is_empty());
        assert!(sys.iter_atoms().all(|a| a.get_type_name().is_none()));
        let p: Vec<_> = sys
            .state()
            .coords
            .iter()
            .map(|p| [p.x as f64, p.y as f64, p.z as f64])
            .collect();
        let evaluation = model.evaluate(&p).unwrap();
        assert!(evaluation.energy.total().is_finite());
        check_gradient(&model, &p);
    }
}
fn check_gradient(model: &BondedModel, positions: &[[f64; 3]]) {
    let actual = model.evaluate(positions).unwrap();
    for i in 0..positions.len() {
        for axis in 0..3 {
            let mut plus = positions.to_vec();
            let mut minus = positions.to_vec();
            let h = 1e-6;
            plus[i][axis] += h;
            minus[i][axis] -= h;
            let expected = -(model.evaluate(&plus).unwrap().energy.total()
                - model.evaluate(&minus).unwrap().energy.total())
                / (2.0 * h);
            let force = actual.forces[i][axis];
            assert!(
                (force - expected).abs() < 1e-3 + expected.abs() * 1e-6,
                "atom {i} axis {axis}: force={force} expected={expected}"
            );
        }
    }
}

// Also used by the independent Python/OpenMM reference generator. It deliberately
// exports resolved parameters, not the implementation's energy or forces.
#[test]
#[ignore = "writes reference-generator input; set MOLAR_AMBER_EXPORT"]
fn export_reference_input() {
    use serde_json::json;
    let mut records = vec![];
    for name in [
        "Ethanol",
        "Acetamide",
        "Caffeine",
        "Glycine",
        "Chloroacetic_acid",
        "Indomethacin_glucuronide",
        "Taurocholic_acid",
        "Phloretin",
        "Amprenavir",
    ] {
        for (ff, ff_name) in [(FFType::Gaff, "gaff"), (FFType::Gaff2, "gaff2")] {
            let sys = System::from_file(format!("tests/data/gaff_ref/sdf/{name}.sdf")).unwrap();
            let model = BondedModel::gaff(&sys, ff).unwrap();
            // Break exact symmetries to test signed derivatives and improper ordering.
            let positions: Vec<_> = sys
                .state()
                .coords
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    [
                        p.x as f64 + 0.003 * (i as f64).sin(),
                        p.y as f64 + 0.002 * (i as f64 * 0.7).cos(),
                        p.z as f64 + 0.004 * (i as f64 * 1.3).sin(),
                    ]
                })
                .collect();
            let torsions = |terms: &[BondedTerm<4, Vec<FourierTerm>>]| {
                terms.iter().map(|t|
                json!({"atoms":t.atoms,"terms":t.parameter.iter().map(|p|json!([p.k,p.phase,p.periodicity])).collect::<Vec<_>>(),
                    "method":format!("{:?}",t.resolution.method),"matched":t.resolution.matched_types,"penalty":t.resolution.penalty})).collect::<Vec<_>>()
            };
            records.push(json!({"name":name,"ff":ff_name,"positions":positions,"types":model.atom_types(),
                "atomic_numbers":sys.iter_atoms().map(|a|a.get_atomic_number()).collect::<Vec<_>>(),
                "formal_charges":sys.iter_atoms().map(|a|a.get_formal_charge().unwrap_or(0)).collect::<Vec<_>>(),
                "connectivity":sys.iter_bonds().map(|b|{let [a,c]=b.pair();[a,c,match b.order(){BondOrder::Single=>1,BondOrder::Double=>2,BondOrder::Triple=>3,BondOrder::Aromatic=>4,_=>0}]}).collect::<Vec<_>>(),
                "bonds":model.bonds().iter().map(|t|json!({"atoms":t.atoms,"k":t.parameter.k,"r":t.parameter.length,"method":format!("{:?}",t.resolution.method),"matched":t.resolution.matched_types,"penalty":t.resolution.penalty})).collect::<Vec<_>>(),
                "angles":model.angles().iter().map(|t|json!({"atoms":t.atoms,"k":t.parameter.k,"theta":t.parameter.angle,"method":format!("{:?}",t.resolution.method),"matched":t.resolution.matched_types,"penalty":t.resolution.penalty})).collect::<Vec<_>>(),
                "propers":torsions(model.propers()),"impropers":torsions(model.impropers())}));
        }
    }
    std::fs::write(
        std::env::var("MOLAR_AMBER_EXPORT").unwrap(),
        serde_json::to_string_pretty(&records).unwrap(),
    )
    .unwrap();
}

#[test]
#[ignore = "full automatic parameter coverage audit"]
fn automatic_corpus() {
    let mut failed = vec![];
    let mut methods = std::collections::BTreeMap::new();
    let mut count = 0;
    let mut files: Vec<_> = std::fs::read_dir("tests/data/gaff_ref/sdf")
        .unwrap()
        .map(|v| v.unwrap().path())
        .collect();
    files.sort();
    for path in files {
        for ff in [FFType::Gaff, FFType::Gaff2] {
            let sys = System::from_file(&path).unwrap();
            match BondedModel::gaff(&sys, ff) {
                Ok(model) => {
                    count += 1;
                    for (_, _, r) in model.estimated_terms() {
                        let entry = methods
                            .entry(format!("{:?}", r.method))
                            .or_insert((0, path.clone()));
                        entry.0 += 1;
                    }
                }
                Err(e) => failed.push((path.clone(), ff, format!("{e}"))),
            }
        }
    }
    eprintln!(
        "models={count}, methods={methods:?}, failures={}",
        failed.len()
    );
    for f in failed.iter().take(10) {
        eprintln!("{f:?}");
    }
    assert_eq!(count, 1190);
    assert_eq!(failed.len(), 4);
    for (path, _, error) in failed {
        assert!(matches!(
            path.file_stem().unwrap().to_str().unwrap(),
            "Hydrochloric_acid" | "Hydrofluoric_acid"
        ));
        assert!(error.contains("MissingParameter { kind: Bond"));
    }
}

#[test]
fn openmm_reference_energies_and_forces() {
    #[derive(serde::Deserialize)]
    struct Fixture {
        records: Vec<Record>,
    }
    #[derive(serde::Deserialize)]
    struct Record {
        name: String,
        ff: String,
        positions: Vec<[f64; 3]>,
        energies: [f64; 4],
        forces: Vec<[f64; 3]>,
        atomic_numbers: Vec<u8>,
        formal_charges: Vec<i32>,
        connectivity: Vec<[usize; 3]>,
    }
    let fixture: Fixture =
        serde_json::from_str(include_str!("data/amber_bonded_reference.json")).unwrap();
    for r in fixture.records {
        let mut top = Topology::default();
        for (&z, &charge) in r.atomic_numbers.iter().zip(&r.formal_charges) {
            let mut atom = Atom::new().with_name("atom");
            atom.atomic_number = z;
            atom.formal_charge = Some(charge);
            top.atoms.push(&atom);
        }
        for [a, b, order] in &r.connectivity {
            top.bonds.push(&Bond::with_order(
                *a,
                *b,
                match order {
                    1 => BondOrder::Single,
                    2 => BondOrder::Double,
                    3 => BondOrder::Triple,
                    4 => BondOrder::Aromatic,
                    _ => panic!("fixture bond order"),
                },
            ));
        }
        let sys = System::new(
            top,
            State {
                coords: r
                    .positions
                    .iter()
                    .map(|p| Pos::new(p[0] as _, p[1] as _, p[2] as _))
                    .collect(),
                ..Default::default()
            },
        )
        .unwrap();
        let ff = if r.ff == "gaff" {
            FFType::Gaff
        } else {
            FFType::Gaff2
        };
        let model = BondedModel::gaff(&sys, ff).unwrap();
        let result = model.evaluate(&r.positions).unwrap();
        let actual = [
            result.energy.bonds,
            result.energy.angles,
            result.energy.propers,
            result.energy.impropers,
        ];
        for (i, (a, b)) in actual.iter().zip(r.energies).enumerate() {
            assert!(
                (a - b).abs() < 1e-7,
                "{} {} term {i}: {a} != {b}",
                r.name,
                r.ff
            );
        }
        for (i, (a, b)) in result.forces.iter().zip(r.forces).enumerate() {
            for j in 0..3 {
                assert!(
                    (a[j] - b[j]).abs() < 1e-5,
                    "{} {} atom {i}: {a:?} != {b:?}",
                    r.name,
                    r.ff
                );
            }
        }
        // Force balance catches accidental double-counting and wrong center gradients.
        for axis in 0..3 {
            assert!(result.forces.iter().map(|f| f[axis]).sum::<f64>().abs() < 1e-7);
        }
    }
}

fn typed_system(types: &[&str], pairs: &[[usize; 2]]) -> System {
    let mut top = Topology::default();
    for t in types {
        let mut atom = Atom::new().with_name("C").guess();
        atom.type_name = Some(t.parse().unwrap());
        top.atoms.push(&atom);
    }
    for &[a, b] in pairs {
        top.bonds.push(&Bond::with_order(a, b, BondOrder::Single));
    }
    System::new(
        top,
        State {
            coords: vec![Pos::origin(); types.len()],
            ..Default::default()
        },
    )
    .unwrap()
}
#[test]
fn unresolved_terms_are_aggregated_with_source_indices_and_no_partial_model() {
    let sys = typed_system(&["zz", "zz", "zz"], &[[0, 1], [1, 2]]);
    let error = BondedModel::from_typed(
        &sys,
        AmberParameters::builtin(FFType::Gaff2),
        ParameterPolicy::Automatic,
    )
    .unwrap_err();
    let BondedError::Unresolved(terms) = error else {
        panic!("wrong error")
    };
    assert_eq!(terms.len(), 3);
    assert_eq!(terms.iter().filter(|t| t.kind == TermKind::Bond).count(), 2);
    let sys = typed_system(&["c3", "c3", "cl", "ha"], &[[0, 1], [2, 3]]);
    let selection = sys.select_bound("index 2:3").unwrap();
    let BondedError::Unresolved(terms) = BondedModel::from_typed(
        &selection,
        AmberParameters::builtin(FFType::Gaff2),
        ParameterPolicy::Automatic,
    )
    .unwrap_err() else {
        panic!("wrong error")
    };
    assert_eq!(terms.len(), 1);
    assert_eq!(terms[0].atoms, vec![2, 3]);
    let open = sys.select_bound("index 2").unwrap();
    assert!(matches!(
        BondedModel::from_typed(
            &open,
            AmberParameters::builtin(FFType::Gaff2),
            ParameterPolicy::Automatic
        ),
        Err(BondedError::OpenSelection(2, 3))
    ));
}

#[test]
fn bond_complete_selection_maps_indices_and_strict_rejects_estimates() {
    let sys = typed_system(
        &["cl", "ha", "c3", "c", "nt", "hn"],
        &[[0, 1], [2, 3], [3, 4], [4, 5]],
    );
    let selection = sys.select_bound("index 2:5").unwrap();
    let parameters = AmberParameters::builtin(FFType::Gaff2);
    let model =
        BondedModel::from_typed(&selection, parameters, ParameterPolicy::Automatic).unwrap();
    assert_eq!(model.atom_indices(), &[2, 3, 4, 5]);
    assert_eq!(model.bonds().len(), 3);
    let mut bond_atoms = model.bonds()[0].atoms;
    bond_atoms.sort_unstable();
    assert_eq!(bond_atoms, [0, 1]);
    assert_eq!(
        model.propers()[0].resolution.method,
        ResolutionMethod::Equivalent
    );
    let BondedError::Unresolved(missing) =
        BondedModel::from_typed(&selection, parameters, ParameterPolicy::Strict).unwrap_err()
    else {
        panic!("wrong error")
    };
    assert_eq!(missing.len(), 1);
    assert_eq!(missing[0].atoms, vec![2, 3, 4, 5]);
}
