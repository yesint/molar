use crate::prelude::*;

use super::BondPerceptionError;

/// A non-fatal condition found during chemical perception.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionWarning {
    /// Two or more assignments have the same best score.
    AmbiguousAssignment,
    /// Formal charges were selected without a requested total charge.
    ChargeWasNotConstrained,
    /// A residue template did not agree with the input graph.
    ResidueTemplateMismatch { residue: usize },
    /// The search hit its branch limit. The returned assignment is the best found so far. It
    /// can be non-optimal, and remaining ambiguity was not fully explored.
    SearchTruncated,
}

/// A complete, validated proposal for bond orders and atom valence data.
///
/// Bond entries use source bond order. Atom entries use source atom order.
/// `None` keeps the value that is already in the topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondAssignment {
    atom_count: usize,
    source_atomic_numbers: Vec<u8>,
    source_pairs: Vec<[usize; 2]>,
    bond_orders: Vec<Option<BondOrder>>,
    formal_charges: Vec<Option<i32>>,
    implicit_hydrogens: Vec<u8>,
    warnings: Vec<PerceptionWarning>,
}

impl BondAssignment {
    /// Make an assignment for `topology` from parallel result arrays.
    pub fn new(
        topology: &Topology,
        bond_orders: Vec<Option<BondOrder>>,
        formal_charges: Vec<Option<i32>>,
        implicit_hydrogens: Vec<u8>,
        warnings: Vec<PerceptionWarning>,
    ) -> Result<Self, BondPerceptionError> {
        check_length("bond-order", topology.bonds.len(), bond_orders.len())?;
        check_length("formal-charge", topology.atoms.len(), formal_charges.len())?;
        check_length(
            "implicit-hydrogen",
            topology.atoms.len(),
            implicit_hydrogens.len(),
        )?;
        Ok(Self {
            atom_count: topology.atoms.len(),
            source_atomic_numbers: topology
                .atoms
                .iter()
                .map(|atom| atom.get_atomic_number())
                .collect(),
            source_pairs: topology.bonds.iter_pairs().collect(),
            bond_orders,
            formal_charges,
            implicit_hydrogens,
            warnings,
        })
    }

    pub fn bond_orders(&self) -> &[Option<BondOrder>] {
        &self.bond_orders
    }

    pub fn formal_charges(&self) -> &[Option<i32>] {
        &self.formal_charges
    }

    pub fn implicit_hydrogens(&self) -> &[u8] {
        &self.implicit_hydrogens
    }

    pub fn warnings(&self) -> &[PerceptionWarning] {
        &self.warnings
    }

    /// Apply this result after verifying that the source graph did not change.
    ///
    /// All validation occurs before the first write. A validation error leaves
    /// `topology` unchanged.
    pub fn apply_to(&self, topology: &mut Topology) -> Result<(), BondPerceptionError> {
        if topology.atoms.len() != self.atom_count {
            return Err(BondPerceptionError::AtomCountChanged {
                expected: self.atom_count,
                actual: topology.atoms.len(),
            });
        }

        for (atom, (actual, expected)) in topology
            .atoms
            .iter()
            .map(|atom| atom.get_atomic_number())
            .zip(self.source_atomic_numbers.iter().copied())
            .enumerate()
        {
            if actual != expected {
                return Err(BondPerceptionError::AtomTableChanged {
                    atom,
                    expected,
                    actual,
                });
            }
        }

        let actual_pairs = topology.bonds.iter_pairs().collect::<Vec<_>>();
        if actual_pairs != self.source_pairs {
            let common = actual_pairs.len().min(self.source_pairs.len());
            let changed = (0..common)
                .find(|&i| actual_pairs[i] != self.source_pairs[i])
                .unwrap_or(common);
            return Err(BondPerceptionError::BondTableChanged {
                bond: changed,
                expected: self.source_pairs.get(changed).copied(),
                actual: actual_pairs.get(changed).copied(),
            });
        }

        for (i, order) in self.bond_orders.iter().copied().enumerate() {
            if let Some(order) = order {
                topology.bonds.set_order(i, order);
            }
        }
        for (i, charge) in self.formal_charges.iter().copied().enumerate() {
            if let Some(charge) = charge {
                topology.atoms.get_mut(i).unwrap().set_formal_charge(charge);
            }
        }
        Ok(())
    }
}

fn check_length(
    field: &'static str,
    expected: usize,
    actual: usize,
) -> Result<(), BondPerceptionError> {
    if actual == expected {
        Ok(())
    } else {
        Err(BondPerceptionError::AssignmentLength {
            field,
            expected,
            actual,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn topology() -> Topology {
        let mut top = Topology::default();
        top.atoms.push(&Atom::new().with_atomic_number(6));
        top.atoms.push(&Atom::new().with_atomic_number(8));
        top.bonds.push(&Bond::new(0, 1));
        top
    }

    #[test]
    fn applies_orders_and_formal_charges() {
        let mut top = topology();
        let assignment = BondAssignment::new(
            &top,
            vec![Some(BondOrder::Double)],
            vec![Some(1), Some(-1)],
            vec![1, 0],
            vec![],
        )
        .unwrap();

        assignment.apply_to(&mut top).unwrap();
        assert_eq!(top.bonds.get(0).unwrap().order(), BondOrder::Double);
        assert_eq!(top.atoms.get(0).unwrap().get_formal_charge(), Some(1));
        assert_eq!(top.atoms.get(1).unwrap().get_formal_charge(), Some(-1));
    }

    #[test]
    fn stale_assignment_does_not_partly_update_topology() {
        let mut top = topology();
        let assignment = BondAssignment::new(
            &top,
            vec![Some(BondOrder::Double)],
            vec![Some(1), Some(-1)],
            vec![1, 0],
            vec![],
        )
        .unwrap();
        top.bonds = std::iter::once(Bond::new(1, 0)).collect();

        assert!(matches!(
            assignment.apply_to(&mut top),
            Err(BondPerceptionError::BondTableChanged { .. })
        ));
        assert_eq!(top.bonds.get(0).unwrap().order(), BondOrder::Unspecified);
        assert_eq!(top.atoms.get(0).unwrap().get_formal_charge(), None);
    }

    #[test]
    fn constructor_rejects_mismatched_arrays() {
        let top = topology();
        assert!(matches!(
            BondAssignment::new(&top, vec![], vec![None, None], vec![0, 0], vec![]),
            Err(BondPerceptionError::AssignmentLength {
                field: "bond-order",
                ..
            })
        ));
    }

    #[test]
    fn atom_change_rejects_assignment_before_writes() {
        let mut top = topology();
        let assignment = BondAssignment::new(
            &top,
            vec![Some(BondOrder::Double)],
            vec![Some(1), Some(-1)],
            vec![1, 0],
            vec![],
        )
        .unwrap();
        top.atoms.get_mut(0).unwrap().set_atomic_number(7);

        assert!(matches!(
            assignment.apply_to(&mut top),
            Err(BondPerceptionError::AtomTableChanged { atom: 0, .. })
        ));
        assert_eq!(top.bonds.get(0).unwrap().order(), BondOrder::Unspecified);
        assert_eq!(top.atoms.get(0).unwrap().get_formal_charge(), None);
    }
}
