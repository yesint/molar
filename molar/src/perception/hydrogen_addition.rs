//! Hydrogen addition (delivery step 7 of the bond-perception plan).
//!
//! [`plan_hydrogen_addition`] reads a system's implicit-hydrogen counts (the shortfall between
//! each atom's bond-order sum and its valence, from [`implicit_hydrogens`](super::implicit_hydrogens))
//! and places that many explicit hydrogens per atom using local geometry. It does not change
//! the system; the resulting [`HydrogenAddition`] is applied as one structural transaction with
//! [`System::add_hydrogens`], which appends the new atoms, coordinates, and bonds together.
//!
//! New hydrogens are appended at the end, so every existing atom index is unchanged — selections
//! and coordinates built before the addition stay valid. Hydrogen addition does not choose a
//! protonation state; it materializes the counts perception already selected.
//!
//! # Geometry
//! Each atom's local shape follows its number of electron domains — bonded neighbors, hydrogens
//! to add, and lone pairs. Two domains are linear, three trigonal, four tetrahedral. The
//! existing neighbors fix part of that shape and the new hydrogens fill the open directions at
//! the standard bond length. Where one bond can still rotate freely (a lone terminal atom) a
//! deterministic reference direction is used, so the result is reproducible.

use crate::prelude::*;

use super::implicit_hydrogens;

/// Options for [`plan_hydrogen_addition`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HydrogenOptions {
    /// When the system carries velocities/forces, give the new hydrogens zero entries. When
    /// false, adding hydrogens to a system that has velocities or forces is an error.
    pub zero_fill_dynamics: bool,
}

impl Default for HydrogenOptions {
    fn default() -> Self {
        Self {
            zero_fill_dynamics: true,
        }
    }
}

/// A planned, validated set of hydrogens to add: for each new hydrogen, its parent heavy atom
/// and its position. Applied with [`System::add_hydrogens`].
#[derive(Debug, Clone)]
pub struct HydrogenAddition {
    /// Atom count the plan was made for; the apply step rejects a changed system.
    source_atom_count: usize,
    /// Parent heavy-atom index for each new hydrogen (in the source indexing).
    parents: Vec<usize>,
    /// Position of each new hydrogen, parallel to `parents`.
    positions: Vec<Pos>,
    zero_fill_dynamics: bool,
}

impl HydrogenAddition {
    /// How many hydrogens this plan adds.
    pub fn len(&self) -> usize {
        self.parents.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parents.is_empty()
    }

    pub fn source_atom_count(&self) -> usize {
        self.source_atom_count
    }

    pub fn parents(&self) -> &[usize] {
        &self.parents
    }

    pub fn positions(&self) -> &[Pos] {
        &self.positions
    }

    pub fn zero_fill_dynamics(&self) -> bool {
        self.zero_fill_dynamics
    }

    /// Apply this plan to `system` (see [`System::add_hydrogens`]).
    pub fn apply_to(&self, system: &mut System) -> Result<(), BondPerceptionError> {
        system.add_hydrogens(self)
    }
}

/// Standard element-hydrogen bond length in nanometers.
fn hydrogen_bond_length(z: u8) -> Float {
    match z {
        7 => 0.101,          // N-H
        8 => 0.097,          // O-H
        16 => 0.134,         // S-H
        15 => 0.142,         // P-H
        5 => 0.119,          // B-H
        _ => 0.109,          // C-H and default
    }
}

/// Number of lone pairs on a heavy atom, from element and formal charge. Only the elements that
/// carry hydrogen in ordinary organic structures are modeled; others report 0.
fn lone_pairs(z: u8, fc: i32) -> i32 {
    let base = match z {
        6 => 0,           // C
        7 => 1,           // N
        8 | 16 => 2,      // O, S
        15 => 1,          // P
        9 | 17 | 35 | 53 => 3, // halogens
        _ => 0,
    };
    (base - fc).max(0)
}

/// A deterministic unit vector perpendicular to `v`.
fn any_perpendicular(v: &Vector3f) -> Vector3f {
    let seed = if v.x.abs() < 0.9 {
        Vector3f::new(1.0, 0.0, 0.0)
    } else {
        Vector3f::new(0.0, 1.0, 0.0)
    };
    (seed - v * v.dot(&seed)).normalize()
}

/// Rotate unit `v` by `angle` (radians) about unit axis `axis` (Rodrigues; `v` need not be
/// perpendicular to `axis`).
fn rotate_about(v: &Vector3f, axis: &Vector3f, angle: Float) -> Vector3f {
    let (s, c) = angle.sin_cos();
    v * c + axis.cross(v) * s + axis * (axis.dot(v) * (1.0 - c))
}

/// Unit directions for `n_h` new hydrogens on an atom whose existing bonds point along
/// `existing` (unit vectors) and whose total electron-domain count is `total` (neighbors +
/// hydrogens + lone pairs → 2 linear, 3 trigonal, 4 tetrahedral). Returns the open directions
/// that complete the shape; lone pairs, if any, take the ones beyond the `n_h` returned.
fn hydrogen_directions(existing: &[Vector3f], n_h: usize, total: usize) -> Vec<Vector3f> {
    // Angle a bond makes with the reference axis in each ideal shape.
    const TETRA_COMP: Float = -1.0 / 3.0; // cos(109.47°)
    const TETRA_PERP: Float = 0.942_809; // sin(109.47°)
    const TETRA_HALF_COS: Float = 0.577_35; // cos(54.735°)
    const TETRA_HALF_SIN: Float = 0.816_5; // sin(54.735°)
    const TRIG_COMP: Float = -0.5; // cos(120°)
    const TRIG_PERP: Float = 0.866_025; // sin(120°)
    let deg120 = 120.0_f64.to_radians() as Float;

    let mut open: Vec<Vector3f> = match (existing.len(), total) {
        // No existing bond: use canonical directions for the shape.
        (0, 1) => vec![Vector3f::new(1.0, 0.0, 0.0)],
        (0, 2) => vec![Vector3f::new(1.0, 0.0, 0.0), Vector3f::new(-1.0, 0.0, 0.0)],
        (0, 3) => vec![
            Vector3f::new(1.0, 0.0, 0.0),
            Vector3f::new(TRIG_COMP, TRIG_PERP, 0.0),
            Vector3f::new(TRIG_COMP, -TRIG_PERP, 0.0),
        ],
        (0, _) => {
            let s = 1.0 / (3.0 as Float).sqrt();
            vec![
                Vector3f::new(s, s, s),
                Vector3f::new(s, -s, -s),
                Vector3f::new(-s, s, -s),
                Vector3f::new(-s, -s, s),
            ]
        }
        // One existing bond d0: complete the shape around it.
        (1, 2) => vec![-existing[0]],
        (1, 3) => {
            let d0 = existing[0];
            let p = any_perpendicular(&d0);
            vec![d0 * TRIG_COMP + p * TRIG_PERP, d0 * TRIG_COMP - p * TRIG_PERP]
        }
        (1, _) => {
            let d0 = existing[0];
            let p0 = any_perpendicular(&d0);
            (0..3)
                .map(|k| {
                    let p = rotate_about(&p0, &d0, deg120 * k as Float);
                    d0 * TETRA_COMP + p * TETRA_PERP
                })
                .collect()
        }
        // Two existing bonds.
        (2, 3) => vec![-(existing[0] + existing[1]).normalize()],
        (2, _) => {
            let (d0, d1) = (existing[0], existing[1]);
            let bisector = -(d0 + d1);
            let bisector = safe_normalize(&bisector, &any_perpendicular(&d0));
            let normal = safe_normalize(&d0.cross(&d1), &any_perpendicular(&bisector));
            vec![
                bisector * TETRA_HALF_COS + normal * TETRA_HALF_SIN,
                bisector * TETRA_HALF_COS - normal * TETRA_HALF_SIN,
            ]
        }
        // Three existing bonds: one open direction opposite their sum.
        (3, _) => {
            let sum = existing[0] + existing[1] + existing[2];
            vec![safe_normalize(&-sum, &any_perpendicular(&existing[0]))]
        }
        // Fully coordinated already: nothing to add.
        _ => Vec::new(),
    };
    open.truncate(n_h);
    open
}

/// Normalize `v`, or return `fallback` when `v` is degenerate (near zero length).
fn safe_normalize(v: &Vector3f, fallback: &Vector3f) -> Vector3f {
    let norm = v.norm();
    if norm > 1e-6 { v / norm } else { *fallback }
}

/// Plan the explicit hydrogens for `system` from its implicit-hydrogen counts and geometry.
pub fn plan_hydrogen_addition(system: &System, options: &HydrogenOptions) -> HydrogenAddition {
    let n = system.len();
    let coords: Vec<Pos> = system.iter_pos().copied().collect();
    let z: Vec<u8> = system.iter_atoms().map(|a| a.get_atomic_number()).collect();
    let fc: Vec<i32> = system
        .iter_atoms()
        .map(|a| a.get_formal_charge().unwrap_or(0))
        .collect();

    let adj = BondAdjacency::build(n, system.topology().bonds.iter_pairs());
    let counts = implicit_hydrogens(system, &adj);

    let mut parents = Vec::new();
    let mut positions = Vec::new();
    for a in 0..n {
        let h = counts[a] as usize;
        if h == 0 {
            continue;
        }
        let neighbors = adj.neighbors(a);
        let existing: Vec<Vector3f> = neighbors
            .iter()
            .map(|nb| (coords[nb.atom()] - coords[a]).normalize())
            .collect();
        let total = neighbors.len() + h + lone_pairs(z[a], fc[a]).max(0) as usize;
        let bond = hydrogen_bond_length(z[a]);
        for dir in hydrogen_directions(&existing, h, total) {
            positions.push(coords[a] + dir * bond);
            parents.push(a);
        }
    }

    HydrogenAddition {
        source_atom_count: n,
        parents,
        positions,
        zero_fill_dynamics: options.zero_fill_dynamics,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a system from elements, coordinates (nm) and bonds with orders.
    fn system(z: &[u8], xyz: &[[Float; 3]], bonds: &[(usize, usize, BondOrder)]) -> System {
        let mut top = Topology::default();
        for &zi in z {
            top.atoms.push(&Atom::new().with_atomic_number(zi));
        }
        for &(i, j, o) in bonds {
            top.bonds.push(&Bond::with_order(i, j, o));
        }
        let mut state = State::default();
        state.coords = xyz.iter().map(|p| Pos::new(p[0], p[1], p[2])).collect();
        System::new(top, state).unwrap()
    }

    fn dist(a: Pos, b: Pos) -> Float {
        (a - b).norm()
    }

    fn angle(center: Pos, a: Pos, b: Pos) -> Float {
        let u = (a - center).normalize();
        let v = (b - center).normalize();
        u.dot(&v).clamp(-1.0, 1.0).acos().to_degrees()
    }

    #[test]
    fn methane_gets_four_tetrahedral_hydrogens() {
        let sys = system(&[6], &[[0.0, 0.0, 0.0]], &[]);
        let plan = plan_hydrogen_addition(&sys, &HydrogenOptions::default());
        assert_eq!(plan.len(), 4);
        for &p in plan.positions() {
            assert!((dist(p, Pos::new(0.0, 0.0, 0.0)) - 0.109).abs() < 1e-4, "C-H length");
        }
        // Tetrahedral: every pair of H subtends ~109.5°.
        for i in 0..4 {
            for j in (i + 1)..4 {
                let a = angle(Pos::new(0.0, 0.0, 0.0), plan.positions()[i], plan.positions()[j]);
                assert!((a - 109.47).abs() < 1.0, "H-C-H angle {a}");
            }
        }
    }

    #[test]
    fn carbonyl_carbon_completes_trigonal() {
        // Formaldehyde heavy atoms: C0=O1 (double). C needs two hydrogens at ~120°, in plane.
        let sys = system(
            &[6, 8],
            &[[0.0, 0.0, 0.0], [0.121, 0.0, 0.0]],
            &[(0, 1, BondOrder::Double)],
        );
        let plan = plan_hydrogen_addition(&sys, &HydrogenOptions::default());
        assert_eq!(plan.len(), 2, "H2C=O");
        let c = Pos::new(0.0, 0.0, 0.0);
        let o = Pos::new(0.121, 0.0, 0.0);
        for &h in plan.positions() {
            let a = angle(c, o, h);
            assert!((a - 120.0).abs() < 2.0, "O=C-H angle {a}");
        }
    }

    #[test]
    fn hydroxyl_oxygen_gets_one_hydrogen() {
        // Methanol heavy atoms: C0-O1 single. O needs one hydrogen at ~109° from the C.
        let sys = system(
            &[6, 8],
            &[[0.0, 0.0, 0.0], [0.143, 0.0, 0.0]],
            &[(0, 1, BondOrder::Single)],
        );
        let plan = plan_hydrogen_addition(&sys, &HydrogenOptions::default());
        // C also gets three hydrogens; check the oxygen's one.
        let o = Pos::new(0.143, 0.0, 0.0);
        let oh: Vec<Pos> = plan
            .parents()
            .iter()
            .zip(plan.positions())
            .filter(|&(&p, _)| p == 1)
            .map(|(_, &pos)| pos)
            .collect();
        assert_eq!(oh.len(), 1, "one O-H");
        assert!((dist(oh[0], o) - 0.097).abs() < 1e-4, "O-H length");
        let a = angle(o, Pos::new(0.0, 0.0, 0.0), oh[0]);
        assert!((a - 109.47).abs() < 2.0, "C-O-H angle {a}");
    }

    #[test]
    fn apply_appends_atoms_bonds_and_coordinates() {
        let mut sys = system(&[6], &[[0.0, 0.0, 0.0]], &[]);
        let plan = plan_hydrogen_addition(&sys, &HydrogenOptions::default());
        plan.apply_to(&mut sys).unwrap();
        assert_eq!(sys.len(), 5, "C + 4 H");
        assert_eq!(sys.state().coords.len(), 5);
        assert_eq!(BondProvider::num_bonds(&sys), 4, "four C-H bonds");
        // Every added atom is a hydrogen bonded to the carbon.
        for i in 1..5 {
            assert_eq!(sys.topology().atoms.get(i).unwrap().get_atomic_number(), 1);
        }
    }

    #[test]
    fn end_to_end_benzene_from_heavy_atoms() {
        // Six bare carbons on a hexagon: perceive orders from geometry, then add the hydrogens.
        let r = 0.139; // regular hexagon: side length equals circumradius
        let coords: Vec<[Float; 3]> = (0..6)
            .map(|k| {
                let t = std::f64::consts::PI / 3.0 * k as f64;
                [r * t.cos() as Float, r * t.sin() as Float, 0.0]
            })
            .collect();
        let bonds: Vec<(usize, usize, BondOrder)> =
            (0..6).map(|k| (k, (k + 1) % 6, BondOrder::Unspecified)).collect();
        let mut sys = system(&[6; 6], &coords, &bonds);

        let opts = BondOrderOptions {
            hydrogens: HydrogenPolicy::InferFromGeometry,
            total_charge: Some(0),
            ..BondOrderOptions::default()
        };
        let assignment = sys.assign_bond_orders(&opts).unwrap();
        sys.apply_bond_assignment(&assignment).unwrap();

        let plan = plan_hydrogen_addition(&sys, &HydrogenOptions::default());
        assert_eq!(plan.len(), 6, "benzene is C6H6");
        plan.apply_to(&mut sys).unwrap();

        assert_eq!(sys.len(), 12);
        assert_eq!(BondProvider::num_bonds(&sys), 12, "6 ring + 6 C-H bonds");
        // Each ring carbon ends up bonded to exactly one hydrogen.
        let adj = BondAdjacency::build(
            sys.len(),
            sys.topology().bonds.iter_pairs(),
        );
        for c in 0..6 {
            let h = adj
                .neighbors(c)
                .iter()
                .filter(|nb| sys.topology().atoms.get(nb.atom()).unwrap().get_atomic_number() == 1)
                .count();
            assert_eq!(h, 1, "carbon {c} has one hydrogen");
        }
    }

    #[test]
    fn apply_rejects_a_changed_system() {
        let sys = system(&[6], &[[0.0, 0.0, 0.0]], &[]);
        let plan = plan_hydrogen_addition(&sys, &HydrogenOptions::default());
        // A different, larger system.
        let mut other = system(&[6, 6], &[[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]], &[(0, 1, BondOrder::Single)]);
        assert!(matches!(
            plan.apply_to(&mut other),
            Err(BondPerceptionError::AtomCountChanged { .. })
        ));
    }
}
