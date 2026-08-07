//! The OpenSMILES **default-valence** model, which decides whether an atom can be written as a
//! bare organic-subset symbol or needs the bracket form.
//!
//! ⚠ This is deliberately **not** [`perception::base_valence`](crate::perception) /
//! `target_valence`. Those are molar's charge-adjusted *chemistry perception* — one valence per
//! element (N is 3, S is 2) adjusted by formal charge, which is what ring/aromaticity/implicit-H
//! perception needs. The table here is the *normative spec* table from OpenSMILES §3.1.5, where
//! an element may have several valences (N 3 or 5, S 2, 4 or 6) and formal charge plays no part
//! — a charged atom is bracketed instead, and a bracket atom carries no implicit hydrogens at
//! all. Both are correct in their own domain; collapsing them into one table would break either
//! aromaticity perception or SMILES round-tripping.

/// Default valences of the SMILES organic subset, in ascending order. An element absent from
/// this table can never be written bare — it always brackets.
fn default_valences(z: u8) -> Option<&'static [u8]> {
    Some(match z {
        5 => &[3],           // B
        6 => &[4],           // C
        7 => &[3, 5],        // N
        8 => &[2],           // O
        15 => &[3, 5],       // P
        16 => &[2, 4, 6],    // S
        9 | 17 | 35 | 53 => &[1], // F, Cl, Br, I
        _ => return None,
    })
}

/// How many hydrogens a **bare** organic-subset atom is understood to carry, given the summed
/// order of its explicitly written bonds: the smallest default valence that is at least
/// `bond_order_sum`, minus that sum.
///
/// `None` for elements outside the organic subset — they must be bracketed, where the hydrogen
/// count is stated outright rather than implied. A bond-order sum exceeding every default
/// valence yields `Some(0)`: the atom is hypervalent for the notation, so no hydrogens are
/// implied and (since the count will not match a real hydrogen) it ends up bracketed anyway.
pub(super) fn implied_h(z: u8, bond_order_sum: u32) -> Option<u8> {
    let valences = default_valences(z)?;
    let v = valences
        .iter()
        .copied()
        .find(|&v| u32::from(v) >= bond_order_sum)
        .map_or(0, |v| u32::from(v) - bond_order_sum);
    Some(v as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn organic_subset_implied_counts() {
        assert_eq!(implied_h(6, 0), Some(4), "methane C");
        assert_eq!(implied_h(6, 3), Some(1), "aromatic-ring CH written Kekulé (1 + 2)");
        assert_eq!(implied_h(6, 4), Some(0), "fully substituted C");
        assert_eq!(implied_h(7, 3), Some(0), "amine N with three bonds");
        assert_eq!(implied_h(7, 0), Some(3), "ammonia N");
        assert_eq!(implied_h(8, 1), Some(1), "hydroxyl O");
        assert_eq!(implied_h(8, 2), Some(0), "carbonyl / ether O");
    }

    /// The multi-valence elements are the whole reason this table is separate from the
    /// perception one: N jumps 3 → 5 and S 2 → 4 → 6 rather than being a single number.
    #[test]
    fn multiple_valences_step_up() {
        assert_eq!(implied_h(7, 4), Some(1), "N past 3 steps to the 5 valence");
        assert_eq!(implied_h(7, 5), Some(0), "nitro / azo N at valence 5");
        assert_eq!(implied_h(16, 2), Some(0), "thioether S");
        assert_eq!(implied_h(16, 3), Some(1), "S past 2 steps to the 4 valence");
        assert_eq!(implied_h(16, 6), Some(0), "sulfone S at valence 6");
    }

    #[test]
    fn non_subset_elements_have_no_implied_hydrogens() {
        for z in [1, 11, 14, 26, 34] {
            // H, Na, Si, Fe, Se
            assert_eq!(implied_h(z, 1), None, "Z={z} must be bracketed");
        }
    }

    /// Hypervalent for the notation: no hydrogens are implied, which is what pushes the atom
    /// into the bracket form where the true count can be stated.
    #[test]
    fn beyond_the_largest_valence_implies_no_hydrogens() {
        assert_eq!(implied_h(6, 5), Some(0));
        assert_eq!(implied_h(16, 8), Some(0));
    }
}
