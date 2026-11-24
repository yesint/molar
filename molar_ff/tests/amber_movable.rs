use molar::prelude::*;
use molar_ff::{FFType, amber::*};

fn molecule() -> (System, AmberModel, Vec<[f64; 3]>) {
    let sys = System::from_file("tests/data/gaff_ref/sdf/Ethanol.sdf").unwrap();
    let model = AmberModel::gaff(&sys, FFType::Gaff2, NonbondedOptions::default()).unwrap();
    let p = sys
        .iter_pos()
        .map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .collect();
    (sys, model, p)
}
#[test]
fn masked_energy_has_constant_offset_and_exact_movable_forces() {
    let (_, model, p) = molecule();
    let mut mask = vec![false; p.len()];
    mask[2] = true;
    let mut offset: Option<f64> = None;
    for shift in [0., 0.01, -0.015] {
        let mut p = p.clone();
        p[2][0] += shift;
        let full = model.evaluate(&p).unwrap();
        let selected = model.evaluate_movable(&p, &mask).unwrap();
        let delta = full.energy.total() - selected.energy.total();
        if let Some(expected) = offset {
            assert!((delta - expected).abs() < 1e-8);
        } else {
            offset = Some(delta);
        }
        for k in 0..3 {
            assert!((full.forces[2][k] - selected.forces[2][k]).abs() < 1e-8);
        }
        for k in 0..3 {
            let mut plus = p.clone();
            let mut minus = p.clone();
            plus[2][k] += 1e-6;
            minus[2][k] -= 1e-6;
            let numerical = -(model.evaluate_movable(&plus, &mask).unwrap().energy.total()
                - model
                    .evaluate_movable(&minus, &mask)
                    .unwrap()
                    .energy
                    .total())
                / 2e-6;
            assert!((numerical - selected.forces[2][k]).abs() < 1e-3);
        }
    }
    let all = model.evaluate_movable(&p, &vec![true; p.len()]).unwrap();
    assert_eq!(all.energy, model.evaluate(&p).unwrap().energy);
    let none = model.evaluate_movable(&p, &vec![false; p.len()]).unwrap();
    assert_eq!(none.energy.total(), 0.);
    assert!(none.forces.iter().flatten().all(|&v| v == 0.));
    assert!(matches!(
        model.evaluate_movable(&p, &[]),
        Err(EvaluationError::MaskSize { .. })
    ));
}
#[test]
fn fixed_group_singularities_are_never_evaluated_during_minimization() {
    let (source, _, p) = molecule();
    let n = source.len();
    let mut top = Topology::default();
    let mut state = State::default();
    for copy in 0..2 {
        for a in source.iter_atoms() {
            top.atoms.push(&Atom::from(&a));
        }
        for point in source.iter_pos() {
            state.coords.push(*point);
        }
        for b in source.topology().bonds.iter() {
            let [a, c] = b.pair();
            top.bonds
                .push(&Bond::with_order(a + copy * n, c + copy * n, b.order()));
        }
    }
    let sys = System::new(top, state).unwrap();
    let model = AmberModel::gaff(&sys, FFType::Gaff2, NonbondedOptions::default()).unwrap();
    let mut positions = vec![[0.; 3]; n];
    positions.extend(p.iter().map(|p| [p[0] + 2., p[1], p[2]]));
    let mut mask = vec![false; n];
    mask.extend(vec![true; n]);
    assert!(model.evaluate(&positions).is_err());
    assert!(model.evaluate_movable(&positions, &mask).is_ok());
    positions[n + 2][1] += 0.02;
    let result = model
        .minimize(&mut positions, &mask, MinimizeOptions::default())
        .unwrap();
    assert!(result.final_energy.total() < result.initial_energy.total());
    assert_eq!(
        result.final_energy,
        model.evaluate_movable(&positions, &mask).unwrap().energy
    );
    assert!(positions[..n].iter().all(|p| *p == [0.; 3]));
    // The excluded group's coincident nonbonded pairs are also skipped when
    // testing the nonbonded evaluator independently of its degenerate bonds.
    assert!(model.nonbonded().evaluate(&positions).is_err());
    assert!(
        model
            .nonbonded()
            .evaluate_movable(&positions, &mask)
            .is_ok()
    );
}

#[test]
fn reusable_energy_scorer_matches_reference_after_large_and_rejected_moves() {
    let (mut sys, _, p) = molecule();
    for (i, mut atom) in sys.iter_atoms_mut().enumerate() {
        atom.set_charge(if i % 2 == 0 { 0.25 } else { -0.25 });
    }
    let model = AmberModel::gaff(&sys, FFType::Gaff2, NonbondedOptions::default()).unwrap();
    let mask: Vec<_> = (0..p.len()).map(|i| i == 0 || i == 3).collect();
    let mut scorer = model.movable_scorer(&p, &mask).unwrap();
    // Move one atom at a time as well as both together. This crosses cell and
    // cutoff boundaries and revisits earlier poses with cached fixed neighbors.
    for (dx, dz) in [
        (0., 0.),
        (0.012, 0.),
        (0.012, 0.02),
        (-1.35, 0.02),
        (1.8, -0.9),
        (0., 0.),
        (0.003, 0.),
    ] {
        let mut trial = p.clone();
        trial[0][0] += dx;
        trial[3][2] += dz;
        let reference = model.evaluate_movable(&trial, &mask).unwrap().energy;
        let cached = scorer.energy(&trial).unwrap();
        assert_eq!(cached.bonded, reference.bonded);
        for (a, b) in [
            (
                cached.nonbonded.lennard_jones,
                reference.nonbonded.lennard_jones,
            ),
            (cached.nonbonded.coulomb, reference.nonbonded.coulomb),
        ] {
            assert!((a - b).abs() <= 1e-12 * b.abs().max(1.));
        }
    }
    let mut invalid = p.clone();
    invalid[0][0] = f64::NAN;
    assert!(matches!(
        scorer.energy(&invalid),
        Err(EvaluationError::NonFinite(0))
    ));
    invalid = p.clone();
    invalid[1][0] += 2.;
    assert!(matches!(
        scorer.energy(&invalid),
        Err(EvaluationError::FixedAtomMoved(1))
    ));
    assert_eq!(
        scorer.energy(&p).unwrap().bonded,
        model.evaluate_movable(&p, &mask).unwrap().energy.bonded
    );
    assert!(matches!(
        scorer.energy(&p[1..]),
        Err(EvaluationError::Size { .. })
    ));
    let mut none = model.movable_scorer(&p, &vec![false; p.len()]).unwrap();
    assert_eq!(none.energy(&p).unwrap(), AmberEnergy::default());
}

#[test]
fn reused_minimizer_matches_fresh_runs_after_mask_changes_and_errors() {
    let (_, model, positions) = molecule();
    let mut worker = model.minimizer();
    for (step, pattern) in [0, 1, 2, 3, 1].into_iter().enumerate() {
        let mask: Vec<_> = (0..positions.len())
            .map(|i| match pattern {
                0 => true,
                1 => i % 2 == 0,
                2 => false,
                _ => i % 2 == 1,
            })
            .collect();
        let mut trial = positions.clone();
        trial[0][0] += 0.012 * step as f64;
        let mut reference = trial.clone();
        let opts = MinimizeOptions {
            max_steps: 8,
            history: step + 1,
            ..Default::default()
        };
        let expected = model.minimize(&mut reference, &mask, opts).unwrap();
        let actual = worker
            .minimize_observed(&mut trial, &mask, opts, |_, _, _| {})
            .unwrap();
        assert_eq!(actual.final_energy, expected.final_energy);
        assert_eq!(actual.evaluations, expected.evaluations);
        assert_eq!(actual.neighbor_rebuilds, expected.neighbor_rebuilds);
        assert_eq!(trial, reference);
        trial[0][0] = f64::NAN;
        assert!(
            worker
                .minimize_observed(&mut trial, &mask, opts, |_, _, _| {})
                .is_err()
        );
    }
}
