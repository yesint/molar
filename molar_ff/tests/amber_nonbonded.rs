use molar::prelude::*;
use molar_ff::{FFType, amber::*};

fn typed(types: &[&str], bonds: &[[usize; 2]]) -> Topology {
    let mut top = Topology::default();
    for t in types {
        let mut a = Atom::new().with_name("C").guess();
        a.type_name = Some(t.parse().unwrap());
        top.atoms.push(&a);
    }
    for &[a, b] in bonds {
        top.bonds.push(&Bond::with_order(a, b, BondOrder::Single));
    }
    top
}
fn model(types: &[&str], bonds: &[[usize; 2]], charges: &[f64], ff: FFType) -> NonbondedModel {
    let parameters = AmberParameters::builtin(ff);
    let bonded =
        BondedModel::from_typed(&typed(types, bonds), parameters, ParameterPolicy::Automatic)
            .unwrap();
    NonbondedModel::from_bonded(&bonded, parameters, charges, NonbondedOptions::default()).unwrap()
}
#[test]
fn amber_lj_units_zero_sites_and_overrides() {
    let parameters = AmberParameters::builtin(FFType::Gaff);
    let (c, _) = parameters.lennard_jones("c3").unwrap();
    assert!((c.radius - 0.1908).abs() < 1e-12);
    assert!((c.epsilon - 0.1094 * 4.184).abs() < 1e-12);
    let (h, _) = parameters.lennard_jones("ho").unwrap();
    assert_eq!(h.epsilon, 0.0);
    assert_eq!(h.radius, 0.0);
    let modified = parameters
        .clone()
        .with_frcmod("override\nNONBON\nc3 2.0 0.2\n", "override")
        .unwrap();
    assert_eq!(modified.lennard_jones("c3").unwrap().0.radius, 0.2);
    assert!(
        parameters
            .clone()
            .with_frcmod("invalid\nNONBON\nc3 2 NaN\n", "bad")
            .is_err()
    );
    assert!(
        parameters
            .clone()
            .with_frcmod("invalid\nLJEDIT\n", "bad")
            .is_err()
    );
}
#[test]
fn cutoff_derivatives_continuity_and_moving_cell_neighbors() {
    let m = model(&["c3", "o"], &[], &[0.4, -0.4], FFType::Gaff2);
    for r in [0.25, 0.9, 1.0, 1.05, 1.19, 1.2, 1.21] {
        let mut p = [[0.0; 3], [r, 0.0, 0.0]];
        let e = m.evaluate(&p).unwrap();
        let h = 1e-6;
        p[1][0] += h;
        let plus = m.evaluate(&p).unwrap().energy.total();
        p[1][0] -= 2.0 * h;
        let minus = m.evaluate(&p).unwrap().energy.total();
        let numerical = -(plus - minus) / (2.0 * h);
        assert!(
            (e.forces[1][0] - numerical).abs() < 1e-3 + 1e-6 * numerical.abs(),
            "r={r}"
        );
        if r >= 1.2 {
            assert_eq!(e.energy.total(), 0.0);
            assert_eq!(e.forces, vec![[0.0; 3]; 2]);
        }
    }
    let near = m.evaluate(&[[0.0; 3], [1.2 - 1e-8, 0.0, 0.0]]).unwrap();
    assert!(near.energy.total().abs() < 1e-10);
    assert!(near.forces[0][0].abs() < 1e-5);
    // Crossing cells, including negative coordinates, must not lose a pair.
    let reference = m.evaluate(&[[0.0; 3], [0.9, 0.0, 0.0]]).unwrap();
    let translated = m.evaluate(&[[-1.3, 0.0, 0.0], [-0.4, 0.0, 0.0]]).unwrap();
    assert!((reference.energy.total() - translated.energy.total()).abs() < 1e-10);
}
#[test]
fn exceptions_use_shortest_paths_and_never_cut_one_four() {
    let chain = model(
        &["c3"; 4],
        &[[0, 1], [1, 2], [2, 3]],
        &[0.2, -0.2, 0.3, -0.3],
        FFType::Gaff2,
    );
    assert_eq!(chain.exceptions().len(), 6);
    let scaled: Vec<_> = chain
        .exceptions()
        .iter()
        .filter(|e| e.lj_scale != 0.0)
        .collect();
    assert_eq!(scaled.len(), 1);
    assert_eq!(scaled[0].atoms, [0, 3]);
    let p = [[0.0; 3], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
    let result = chain.evaluate(&p).unwrap();
    assert!((result.energy.coulomb - COULOMB_CONSTANT * 0.2 * (-0.3) / 1.2 / 1.5).abs() < 1e-12);
    let ring = model(
        &["c3"; 4],
        &[[0, 1], [1, 2], [2, 3], [3, 0]],
        &[0.2; 4],
        FFType::Gaff2,
    );
    assert!(ring.exceptions().iter().all(|e| e.lj_scale == 0.0));
    assert_eq!(ring.evaluate(&p).unwrap().energy.total(), 0.0);
    let six = model(
        &["c3"; 6],
        &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]],
        &[0.0; 6],
        FFType::Gaff2,
    );
    assert_eq!(
        six.exceptions()
            .iter()
            .filter(|e| e.lj_scale != 0.0)
            .count(),
        3
    );
}
#[test]
fn invalid_parameters_and_geometry_are_errors() {
    let parameters = AmberParameters::builtin(FFType::Gaff2);
    let b = BondedModel::from_typed(
        &typed(&["unknown"], &[]),
        parameters,
        ParameterPolicy::Automatic,
    )
    .unwrap();
    assert!(matches!(
        NonbondedModel::from_bonded(&b, parameters, &[0.0], NonbondedOptions::default()),
        Err(NonbondedError::MissingParameters(_))
    ));
    let b = BondedModel::from_typed(&typed(&["c3"], &[]), parameters, ParameterPolicy::Automatic)
        .unwrap();
    assert!(matches!(
        NonbondedModel::from_bonded(&b, parameters, &[f64::NAN], NonbondedOptions::default()),
        Err(NonbondedError::NonFiniteCharge(0))
    ));
    assert!(
        NonbondedModel::from_bonded(
            &b,
            parameters,
            &[0.0],
            NonbondedOptions {
                cutoff: 0.8,
                switch_distance: 1.0
            }
        )
        .is_err()
    );
    let m = model(&["c3"; 2], &[], &[0.0; 2], FFType::Gaff2);
    assert!(matches!(
        m.evaluate(&[[0.0; 3]; 2]),
        Err(EvaluationError::CoincidentNonbonded { .. })
    ));
}
#[test]
fn openmm_nonbonded_reference() {
    #[derive(serde::Deserialize)]
    struct Fixture {
        records: Vec<Record>,
    }
    #[derive(serde::Deserialize)]
    struct Record {
        name: String,
        ff: String,
        types: Vec<String>,
        bonds: Vec<[usize; 2]>,
        charges: Vec<f64>,
        positions: Vec<[f64; 3]>,
        energies: [f64; 2],
        forces: Vec<[f64; 3]>,
    }
    let fixture: Fixture =
        serde_json::from_str(include_str!("data/amber_nonbonded_reference.json")).unwrap();
    for r in fixture.records {
        let types: Vec<_> = r.types.iter().map(String::as_str).collect();
        let ff = if r.ff == "gaff" {
            FFType::Gaff
        } else {
            FFType::Gaff2
        };
        let m = model(&types, &r.bonds, &r.charges, ff);
        let actual = m.evaluate(&r.positions).unwrap();
        for (a, b) in [actual.energy.lennard_jones, actual.energy.coulomb]
            .iter()
            .zip(r.energies)
        {
            assert!(
                (a - b).abs() < 1e-7 + 1e-10 * b.abs(),
                "{} {}: {a} vs {b}",
                r.name,
                r.ff
            );
        }
        for (i, (a, b)) in actual.forces.iter().zip(r.forces).enumerate() {
            for j in 0..3 {
                assert!(
                    (a[j] - b[j]).abs() < 1e-6 + 1e-10 * b[j].abs(),
                    "{} {} atom {i}: {a:?} vs {b:?}",
                    r.name,
                    r.ff
                );
            }
        }
    }
}
