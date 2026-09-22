use crate::prelude::*;

use super::BondPerceptionError;

/// Options for distance-based connectivity perception.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConnectivityOptions {
    /// Added to the sum of the two covalent radii, in nanometers.
    pub tolerance: Float,
    /// Reject atom pairs closer than this distance, in nanometers.
    pub minimum_distance: Float,
    /// Periodic lattice directions used for the distance calculation.
    pub pbc: PbcDims,
    /// Remove the longest candidate edges from over-coordinated atoms.
    pub cleanup_overcoordination: bool,
}

impl Default for ConnectivityOptions {
    fn default() -> Self {
        Self {
            // The Open Babel ConnectTheDots tolerance is 0.45 angstrom.
            tolerance: 0.045,
            minimum_distance: 0.04,
            pbc: PBC_NONE,
            cleanup_overcoordination: true,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    i: usize,
    j: usize,
    distance: Float,
}

/// Perceive a connection table from elements and coordinates.
///
/// Returned atom indices are local to `data`, and all returned bond orders are
/// [`BondOrder::Unspecified`]. The function does not change `data`.
pub fn perceive_connectivity(
    data: &(impl AtomProvider + PosProvider + BoxProvider),
    options: &ConnectivityOptions,
) -> Result<BondStorage, BondPerceptionError> {
    validate_options(options)?;
    if data.len() < 2 {
        return Ok(BondStorage::default());
    }

    let mut radii = Vec::with_capacity(data.len());
    let mut coordination = Vec::with_capacity(data.len());
    for (i, atom) in data.iter_atoms().enumerate() {
        let z = atom.get_atomic_number();
        let Some(radius) = covalent_radius(z) else {
            return Err(BondPerceptionError::UnsupportedElement {
                atom: i,
                atomic_number: z,
            });
        };
        radii.push(radius);
        coordination.push(maximum_coordination(z));
    }

    let maximum_radius = radii.iter().copied().fold(0.0 as Float, Float::max);
    let cutoff = 2.0 * maximum_radius + options.tolerance;
    let nearby: Vec<(usize, usize, Float)> = if options.pbc.any() {
        let pbox = data
            .get_box()
            .ok_or(BondPerceptionError::MissingPeriodicBox)?;
        distance_search_single_pbc(cutoff, data.iter_pos(), 0..data.len(), pbox, options.pbc)?
    } else {
        distance_search_single(cutoff, data, 0..data.len())?
    };

    let mut candidates = nearby
        .into_iter()
        .filter_map(|(a, b, distance)| {
            let (i, j) = if a < b { (a, b) } else { (b, a) };
            (i != j
                && distance >= options.minimum_distance
                && distance <= radii[i] + radii[j] + options.tolerance)
                .then_some(Candidate { i, j, distance })
        })
        .collect::<Vec<_>>();

    candidates.sort_unstable_by_key(|edge| (edge.i, edge.j));
    candidates.dedup_by_key(|edge| (edge.i, edge.j));

    if options.cleanup_overcoordination {
        remove_overcoordination(&mut candidates, &coordination);
    }

    candidates.sort_unstable_by_key(|edge| (edge.i, edge.j));
    Ok(candidates
        .into_iter()
        .map(|edge| Bond::new(edge.i, edge.j))
        .collect())
}

fn validate_options(options: &ConnectivityOptions) -> Result<(), BondPerceptionError> {
    if !options.tolerance.is_finite() || options.tolerance < 0.0 {
        return Err(BondPerceptionError::InvalidTolerance(options.tolerance));
    }
    if !options.minimum_distance.is_finite() || options.minimum_distance < 0.0 {
        return Err(BondPerceptionError::InvalidMinimumDistance(
            options.minimum_distance,
        ));
    }
    Ok(())
}

/// Remove long candidate edges until each atom is within its coordination cap.
fn remove_overcoordination(candidates: &mut Vec<Candidate>, limits: &[usize]) {
    let mut degree = vec![0usize; limits.len()];
    for edge in candidates.iter() {
        degree[edge.i] += 1;
        degree[edge.j] += 1;
    }

    let mut longest_first = (0..candidates.len()).collect::<Vec<_>>();
    longest_first.sort_unstable_by(|&a, &b| {
        candidates[b]
            .distance
            .total_cmp(&candidates[a].distance)
            .then_with(|| candidates[b].i.cmp(&candidates[a].i))
            .then_with(|| candidates[b].j.cmp(&candidates[a].j))
    });

    let mut keep = vec![true; candidates.len()];
    for edge_index in longest_first {
        let edge = candidates[edge_index];
        if degree[edge.i] > limits[edge.i] || degree[edge.j] > limits[edge.j] {
            keep[edge_index] = false;
            degree[edge.i] -= 1;
            degree[edge.j] -= 1;
        }
    }

    let mut index = 0usize;
    candidates.retain(|_| {
        let retain = keep[index];
        index += 1;
        retain
    });
}

/// Single-bond covalent radii in nanometers, based on Cordero et al. values.
///
/// This first table covers elements H through Kr. Later parity work will extend
/// the table together with explicit coordination and metal rules.
pub(super) fn covalent_radius(z: u8) -> Option<Float> {
    const RADII: [Float; 37] = [
        0.000, 0.031, 0.028, 0.128, 0.096, 0.084, 0.076, 0.071, 0.066, 0.057, 0.058, 0.166, 0.141,
        0.121, 0.111, 0.107, 0.105, 0.102, 0.106, 0.203, 0.176, 0.170, 0.160, 0.153, 0.139, 0.139,
        0.132, 0.126, 0.124, 0.132, 0.122, 0.122, 0.120, 0.119, 0.120, 0.120, 0.116,
    ];
    (z != 0).then(|| RADII.get(z as usize).copied()).flatten()
}

fn maximum_coordination(z: u8) -> usize {
    match z {
        2 | 10 | 18 | 36 => 0,
        1 | 9 | 17 | 35 => 1,
        3 | 11 | 19 => 8,
        4 | 12 | 20 => 8,
        5 | 6 | 14 | 32 => 4,
        7 => 4,
        8 => 3,
        13 | 31 => 6,
        15 | 16 | 33 | 34 => 6,
        21..=30 => 8,
        _ => 8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn system(elements: &[u8], xyz_nm: &[[Float; 3]]) -> System {
        let mut top = Topology::default();
        for &z in elements {
            top.atoms.push(&Atom::new().with_atomic_number(z));
        }
        let mut state = State::default();
        state.coords = xyz_nm.iter().map(|p| Pos::new(p[0], p[1], p[2])).collect();
        System::new(top, state).unwrap()
    }

    #[test]
    fn methane_connectivity_has_unspecified_orders() {
        let molecule = system(
            &[6, 1, 1, 1, 1],
            &[
                [0.0, 0.0, 0.0],
                [0.063, 0.063, 0.063],
                [-0.063, -0.063, 0.063],
                [-0.063, 0.063, -0.063],
                [0.063, -0.063, -0.063],
            ],
        );
        let bonds = perceive_connectivity(&molecule, &ConnectivityOptions::default()).unwrap();

        assert_eq!(bonds.len(), 4);
        assert!(!bonds.has_orders());
        assert!(
            bonds
                .iter()
                .all(|bond| bond.order() == BondOrder::Unspecified)
        );
        assert!(bonds.iter().all(|bond| bond.contains(0)));
    }

    #[test]
    fn rejects_a_pair_that_is_outside_the_radius_tolerance() {
        let molecule = system(&[6, 6], &[[0.0, 0.0, 0.0], [0.20, 0.0, 0.0]]);
        let bonds = perceive_connectivity(&molecule, &ConnectivityOptions::default()).unwrap();
        assert!(bonds.is_empty());
    }

    #[test]
    fn cleanup_keeps_only_one_bond_to_hydrogen() {
        let molecule = system(
            &[1, 6, 6],
            &[[0.0, 0.0, 0.0], [-0.09, 0.0, 0.0], [0.10, 0.0, 0.0]],
        );
        let bonds = perceive_connectivity(&molecule, &ConnectivityOptions::default()).unwrap();
        let hydrogen_bonds = bonds
            .iter()
            .filter(|bond| bond.contains(0))
            .collect::<Vec<_>>();
        assert_eq!(hydrogen_bonds.len(), 1);
        assert_eq!(hydrogen_bonds[0].pair(), [0, 1]);
    }

    #[test]
    fn invalid_options_fail_before_search() {
        let molecule = system(&[6, 6], &[[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]]);
        let options = ConnectivityOptions {
            tolerance: -0.1,
            ..ConnectivityOptions::default()
        };
        assert!(matches!(
            perceive_connectivity(&molecule, &options),
            Err(BondPerceptionError::InvalidTolerance(_))
        ));
    }

    #[test]
    fn periodic_search_connects_atoms_across_a_box_face() {
        let mut molecule = system(&[6, 6], &[[0.03, 0.0, 0.0], [0.97, 0.0, 0.0]]);
        let mut state = molecule.state().clone();
        state.pbox = Some(PeriodicBox::from_matrix(Matrix3f::identity()).unwrap());
        molecule.set_state(state).unwrap();
        assert!(
            perceive_connectivity(&molecule, &ConnectivityOptions::default())
                .unwrap()
                .is_empty()
        );

        let options = ConnectivityOptions {
            pbc: PBC_FULL,
            ..ConnectivityOptions::default()
        };
        let bonds = perceive_connectivity(&molecule, &options).unwrap();
        assert_eq!(bonds.len(), 1);
        assert_eq!(bonds.get(0).unwrap().pair(), [0, 1]);
    }

    #[test]
    fn system_method_replaces_the_bond_table() {
        let mut molecule = system(&[6, 6], &[[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]]);
        let old_bonds = std::iter::once(Bond::with_order(0, 1, BondOrder::Double)).collect();
        molecule.set_bonds(old_bonds).unwrap();

        let count = molecule
            .perceive_connectivity(&ConnectivityOptions::default())
            .unwrap();
        assert_eq!(count, 1);
        assert_eq!(molecule.topology().bonds.len(), 1);
        assert_eq!(
            molecule.topology().bonds.get(0).unwrap().order(),
            BondOrder::Unspecified
        );
        assert!(!molecule.topology().bonds.has_orders());
    }
}
