use molar::prelude::*;
use molar_ff::{ApplyCharges, ChargeModel, FFType, amber::*};

#[test]
fn cleanup_lowers_energy_keeps_fixed_atoms_and_reports_final_force() {
    let mut system = System::from_file("tests/data/gaff_ref/sdf/Ethanol.sdf").unwrap();
    system.apply_charges(ChargeModel::Espaloma).unwrap();
    let model = AmberModel::gaff(&system, FFType::Gaff2, NonbondedOptions::default()).unwrap();
    let mut p: Vec<[f64; 3]> = system
        .state()
        .coords
        .iter()
        .enumerate()
        .map(|(i, p)| {
            [
                p.x as f64 + 0.015 * (i as f64).sin(),
                p.y as f64 + 0.012 * (i as f64 * 0.7).cos(),
                p.z as f64 + 0.01 * (i as f64 * 1.3).sin(),
            ]
        })
        .collect();
    let fixed = p[0];
    let mut movable = vec![true; p.len()];
    movable[0] = false;
    let result = model
        .minimize(&mut p, &movable, MinimizeOptions::default())
        .unwrap();
    assert!(result.converged(), "{result:?}");
    assert!(result.final_energy.total() < result.initial_energy.total() - 1.0);
    assert_eq!(p[0], fixed);
    let final_eval = model.evaluate(&p).unwrap();
    assert!((result.final_energy.total() - final_eval.energy.total()).abs() < 1e-10);
    let actual = final_eval
        .forces
        .iter()
        .zip(movable)
        .filter(|(_, m)| *m)
        .map(|(v, _)| v[0].hypot(v[1]).hypot(v[2]))
        .fold(0.0, f64::max);
    assert!((actual - result.max_force).abs() < 1e-10);
}
#[test]
fn fixed_atoms_zero_budget_and_invalid_inputs() {
    let system = System::from_file("tests/data/gaff_ref/sdf/Ethanol.sdf").unwrap();
    let model = AmberModel::gaff(&system, FFType::Gaff2, NonbondedOptions::default()).unwrap();
    let mut p: Vec<[f64; 3]> = system
        .state()
        .coords
        .iter()
        .map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .collect();
    let before = p.clone();
    let n = p.len();
    let result = model
        .minimize(&mut p, &vec![false; n], MinimizeOptions::default())
        .unwrap();
    assert!(result.converged());
    assert_eq!(result.steps, 0);
    assert_eq!(p, before);
    let result = model
        .minimize(
            &mut p,
            &vec![true; n],
            MinimizeOptions {
                max_steps: 0,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(result.reason, StopReason::StepLimit);
    assert_eq!(p, before);
    assert!(
        model
            .minimize(&mut p, &[], MinimizeOptions::default())
            .is_err()
    );
    assert_eq!(p, before);
    p[0][0] = f64::NAN;
    assert!(matches!(
        model.minimize(&mut p, &vec![true; n], MinimizeOptions::default()),
        Err(MinimizeError::Evaluation(EvaluationError::NonFinite(0)))
    ));
}
