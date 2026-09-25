//! Bond-order and formal-charge perception (delivery step 3 of the bond-perception plan).
//!
//! [`assign_bond_orders`] takes a connection table with known connectivity — orders may be
//! absent, partly present, or aromatic — and returns a validated [`BondAssignment`] with a
//! bond order and a formal charge for every solved bond and atom. It does not change the
//! topology; apply the result with [`BondAssignment::apply_to`].
//!
//! # Model
//! Each element has an ordered set of candidate *valence states*. A state is a total
//! bond-order sum (the valence), the formal charge that goes with it, and a penalty; a lower
//! penalty is a more common state. The solver assigns an integer order to every bond so that
//! each atom's incident-order sum matches one of its states, and it minimizes the summed
//! penalty over the whole fragment. An optional total charge constrains the result.
//!
//! # Hydrogen
//! Two policies (see [`HydrogenPolicy`]):
//! - [`HydrogenPolicy::AllExplicit`] (default): every hydrogen is an explicit atom, so an
//!   atom's bond-order sum must reach a valence state exactly, with no implicit hydrogen. This
//!   is the right choice for structures read with hydrogen (SDF, most ligands).
//! - [`HydrogenPolicy::InferFromGeometry`]: hydrogens may be missing. Bond lengths bound each
//!   heavy-heavy bond's order (a long bond is single, a short one may be double or triple), and
//!   the shortfall between an atom's bond-order sum and its valence becomes implicit hydrogen.
//!   This lets a hydrogen-free structure be perceived; geometry is what separates, e.g., benzene
//!   from cyclohexane, which share a heavy-atom graph.
//!
//! Elements the table does not model impose no valence constraint and take single bonds only.
//! Residue templates for polymers are a later step in the plan.

use crate::prelude::*;

use super::{
    BondAssignment, BondPerceptionError, PerceptionWarning, connected_components, kekulize,
};

/// How the solver treats bond orders already present on the input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputOrders {
    /// Keep every concrete (non-[`BondOrder::Unspecified`]) input order and solve only the
    /// unspecified bonds. Aromatic input is kekulized first and the result is kept.
    PreserveKnown,
    /// Treat every non-aromatic input order as a free variable. Aromatic input is still
    /// kekulized first and the resulting single/double bonds are kept.
    ReassignAll,
}

/// Explicit work limits for the per-fragment search. Reaching a limit produces an error (when
/// no solution was found) or a [`PerceptionWarning::SearchTruncated`] (when one was), never a
/// silent partial assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SearchLimits {
    /// Maximum number of search-tree nodes visited per bonded fragment.
    pub max_branches: u64,
}

impl Default for SearchLimits {
    fn default() -> Self {
        Self {
            max_branches: 5_000_000,
        }
    }
}

/// How the solver accounts for hydrogen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HydrogenPolicy {
    /// Every hydrogen is an explicit atom. Each atom's bond-order sum must reach a full valence
    /// state exactly; no implicit hydrogen is inferred. This is the default and is the right
    /// choice for structures read with hydrogen (SDF, most ligands).
    AllExplicit,
    /// Hydrogens may be missing. The shortfall between an atom's bond-order sum and its valence
    /// becomes implicit hydrogen, and bond lengths bound each heavy-heavy bond's order (a long
    /// bond stays single, a short one may be double or triple). Requires coordinates.
    InferFromGeometry,
}

/// Options for [`assign_bond_orders`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondOrderOptions {
    /// How to treat orders already on the input bonds.
    pub input_orders: InputOrders,
    /// How to account for hydrogen (see [`HydrogenPolicy`]).
    pub hydrogens: HydrogenPolicy,
    /// Pin recognized functional groups (nitro, carboxyl, azide, …) to their canonical bond
    /// orders before the search. On by default.
    pub use_functional_groups: bool,
    /// Pin standard biopolymer residues (amino-acid backbones and side chains) to their
    /// canonical bond orders before the search. On by default.
    pub use_residue_templates: bool,
    /// Optional net formal charge of the molecule. Only supported when the molecule is a
    /// single bonded fragment; with more than one fragment it is an error.
    pub total_charge: Option<i32>,
    /// Search work limits.
    pub limits: SearchLimits,
}

impl Default for BondOrderOptions {
    fn default() -> Self {
        Self {
            input_orders: InputOrders::ReassignAll,
            hydrogens: HydrogenPolicy::AllExplicit,
            use_functional_groups: true,
            use_residue_templates: true,
            total_charge: None,
            limits: SearchLimits::default(),
        }
    }
}

/// One candidate valence state of an element: its total bond-order sum, its formal charge, and
/// a preference penalty (lower is more common).
#[derive(Debug, Clone, Copy)]
struct ValenceState {
    valence: u8,
    formal_charge: i32,
    penalty: u32,
}

const fn vs(valence: u8, formal_charge: i32, penalty: u32) -> ValenceState {
    ValenceState {
        valence,
        formal_charge,
        penalty,
    }
}

// Ordered candidate valence states per element. Penalties are chemistry-derived, not copied
// from Antechamber; step 4 of the plan (the Antechamber reference harness) tunes them.
const STATES_H: &[ValenceState] = &[vs(1, 0, 0)];
const STATES_B: &[ValenceState] = &[vs(3, 0, 0), vs(4, -1, 32)];
const STATES_C: &[ValenceState] = &[vs(4, 0, 0), vs(3, 1, 32), vs(3, -1, 32)];
const STATES_N: &[ValenceState] = &[vs(3, 0, 0), vs(4, 1, 1), vs(2, -1, 3)];
const STATES_O: &[ValenceState] = &[vs(2, 0, 0), vs(1, -1, 1), vs(3, 1, 4)];
const STATES_P: &[ValenceState] = &[vs(3, 0, 0), vs(5, 0, 0), vs(4, 1, 2), vs(2, -1, 4)];
// Divalent sulfur is by far the common case. Hypervalent sulfur (a sulfoxide at valence 4, a
// sulfone/sulfate at valence 6) is real but must be dispreferred, or a ring sulfur that should
// stay divalent and aromatic (thiophene, thiazole) is instead given two double bonds.
const STATES_S: &[ValenceState] =
    &[vs(2, 0, 0), vs(4, 0, 2), vs(6, 0, 3), vs(1, -1, 1), vs(3, 1, 3)];
const STATES_HALOGEN: &[ValenceState] = &[vs(1, 0, 0)];

/// Ordered candidate valence states for the common organic elements. `None` marks an element
/// the solver does not model: it imposes no valence constraint and takes single bonds only.
fn valence_states(z: u8) -> Option<&'static [ValenceState]> {
    Some(match z {
        1 => STATES_H,
        5 => STATES_B,
        6 => STATES_C,
        7 => STATES_N,
        8 => STATES_O,
        15 => STATES_P,
        16 => STATES_S,
        9 | 17 | 35 | 53 => STATES_HALOGEN,
        _ => return None,
    })
}

/// The largest bond order an element accepts. Hydrogen, halogens, and unmodeled elements take
/// single bonds only; oxygen up to double; everything else up to triple.
fn elem_max_order(z: u8) -> u8 {
    if valence_states(z).is_none() {
        return 1;
    }
    match z {
        1 | 9 | 17 | 35 | 53 => 1,
        8 => 2,
        _ => 3,
    }
}

/// Geometry-implied bond order from length: the most plausible order for a heavy-heavy bond of
/// these elements at `length_nm`, capped by each element's chemistry. The single-bond
/// reference is the covalent-radius sum; a double bond sits near 0.87 of it and a triple near
/// 0.78, so the midpoints 0.935 and 0.825 separate the orders. Returns 1 when a covalent radius
/// is unknown. Used only in [`HydrogenPolicy::InferFromGeometry`].
fn geometric_order(zi: u8, zj: u8, length_nm: Float) -> u8 {
    let cap = elem_max_order(zi).min(elem_max_order(zj));
    if cap <= 1 {
        return cap;
    }
    let (Some(ri), Some(rj)) = (
        super::connectivity_perception::covalent_radius(zi),
        super::connectivity_perception::covalent_radius(zj),
    ) else {
        return 1;
    };
    let single = ri + rj;
    let order = if length_nm < 0.825 * single {
        3
    } else if length_nm < 0.935 * single {
        2
    } else {
        1
    };
    order.min(cap)
}

fn order_to_int(o: BondOrder) -> Option<u8> {
    match o {
        BondOrder::Single => Some(1),
        BondOrder::Double => Some(2),
        BondOrder::Triple => Some(3),
        BondOrder::Unspecified | BondOrder::Aromatic => None,
    }
}

fn int_to_order(o: u8) -> BondOrder {
    match o {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        _ => unreachable!("bond order out of range: {o}"),
    }
}

/// Perceive bond orders and formal charges for `top`, returning a validated [`BondAssignment`].
///
/// `coords` (one position per atom) is required by [`HydrogenPolicy::InferFromGeometry`], where
/// bond lengths bound each order, and ignored otherwise. This does not change `top`.
pub fn assign_bond_orders(
    top: &Topology,
    coords: Option<&[Pos]>,
    options: &BondOrderOptions,
) -> Result<BondAssignment, BondPerceptionError> {
    let n = top.atoms.len();
    let m = top.bonds.len();

    let infer = options.hydrogens == HydrogenPolicy::InferFromGeometry;
    let coords = if infer {
        match coords {
            Some(c) if c.len() == n => Some(c),
            Some(_) => return Err(BondPerceptionError::CoordinateCountMismatch),
            None => return Err(BondPerceptionError::MissingCoordinates),
        }
    } else {
        None
    };

    let z: Vec<u8> = top.atoms.iter().map(|a| a.get_atomic_number()).collect();
    let input_fc: Vec<i32> = top
        .atoms
        .iter()
        .map(|a| a.get_formal_charge().unwrap_or(0))
        .collect();
    let pairs: Vec<[usize; 2]> = top.bonds.iter_pairs().collect();
    let input_orders: Vec<BondOrder> = top.bonds.iter().map(|b| b.order()).collect();

    // A throwaway adjacency, built without mutating `top` (this is a non-mutating calculation).
    let adj = BondAdjacency::build(n, top.bonds.iter_pairs());

    // Aromatic input becomes a concrete Kekulé structure before the valence search.
    let kekulized = if input_orders.contains(&BondOrder::Aromatic) {
        kekulize(&z, &input_fc, &input_orders, &adj)
            .map_err(|e| BondPerceptionError::Kekulization(e.to_string()))?
    } else {
        input_orders.clone()
    };

    // Canonical bond orders pinned before the search: standard residues first (most specific),
    // then functional groups. These take priority over the input orders and the search.
    let mut template = if options.use_functional_groups {
        super::functional_groups::functional_group_orders(&z, &adj, coords)
    } else {
        vec![None; m]
    };
    if options.use_residue_templates {
        let residue = super::residue_templates::residue_template_orders(top, &pairs);
        for b in 0..m {
            if residue[b].is_some() {
                template[b] = residue[b];
            }
        }
    }

    // Per bond: a fixed integer order, or `None` for a free variable. The domain cap is the
    // chemical limit, tightened to the geometry-implied order when inferring from coordinates.
    let mut fixed: Vec<Option<u8>> = vec![None; m];
    let dom_max: Vec<u8> = (0..m)
        .map(|b| {
            let [i, j] = pairs[b];
            let chem = elem_max_order(z[i]).min(elem_max_order(z[j]));
            match coords {
                Some(c) => geometric_order(z[i], z[j], (c[i] - c[j]).norm()).min(chem),
                None => chem,
            }
        })
        .collect();
    for b in 0..m {
        let [i, j] = pairs[b];
        // A functional-group template pins this bond outright.
        if let Some(o) = template[b] {
            fixed[b] = order_to_int(o);
            continue;
        }
        // A bond to hydrogen, a halogen, or an unmodeled element can only be single.
        if elem_max_order(z[i]) == 1 || elem_max_order(z[j]) == 1 {
            fixed[b] = Some(1);
            continue;
        }
        let was_aromatic = input_orders[b] == BondOrder::Aromatic;
        match options.input_orders {
            InputOrders::PreserveKnown => {
                if let Some(o) = order_to_int(kekulized[b]) {
                    fixed[b] = Some(o);
                }
            }
            InputOrders::ReassignAll => {
                // A kekulized aromatic bond is already resolved; keep it. Every other order is
                // a free variable.
                if was_aromatic {
                    fixed[b] = order_to_int(kekulized[b]);
                }
            }
        }
    }

    // Carbons and oxygens the geometry shows as sp2/sp: a heavy-atom bond short enough for a
    // multiple order. A saturated neutral state contradicts that geometry (see
    // `geometry_contradiction`), which the search counts as a tie-break (only when inferring
    // hydrogens).
    let sp2_geometry: Vec<bool> = (0..n)
        .map(|a| {
            infer
                && matches!(z[a], 6 | 8)
                && adj.neighbors(a).iter().any(|nb| dom_max[nb.bond()] >= 2)
        })
        .collect();
    let degree: Vec<i32> = (0..n).map(|a| adj.neighbors(a).len() as i32).collect();
    // Two-coordinate C and N the geometry shows as bent. Two π bonds (C=C=C, HN=C=O, R-C≡N)
    // make an atom sp, which is linear; on a bent atom that is a contradiction (see
    // `Search::linear_violations`). Only C and N can carry two π bonds among the modeled
    // second-row elements.
    let bent: Vec<bool> = match coords {
        Some(c) => (0..n)
            .map(|a| {
                let nbs = adj.neighbors(a);
                matches!(z[a], 6 | 7) && nbs.len() == 2 && {
                    let u = c[nbs[0].atom()] - c[a];
                    let v = c[nbs[1].atom()] - c[a];
                    u.dot(&v) / (u.norm() * v.norm()) > SP_MIN_ANGLE_DEG.to_radians().cos()
                }
            })
            .collect(),
        None => vec![false; n],
    };

    let comp = connected_components(&adj);
    let n_comp = comp.iter().copied().max().map_or(0, |c| c + 1);
    let mut comp_bonds: Vec<Vec<usize>> = vec![Vec::new(); n_comp];
    let mut comp_atoms: Vec<Vec<usize>> = vec![Vec::new(); n_comp];
    for b in 0..m {
        comp_bonds[comp[pairs[b][0]]].push(b);
    }
    for a in 0..n {
        comp_atoms[comp[a]].push(a);
    }
    let bonded_comp_count = comp_bonds.iter().filter(|b| !b.is_empty()).count();

    if options.total_charge.is_some() && bonded_comp_count > 1 {
        return Err(BondPerceptionError::TotalChargeWithMultipleComponents);
    }

    // Rings and ring membership, for the cumulene guard and the tautomer tie-break (only needed
    // when inferring hydrogens).
    let rings = if infer { super::sssr(&adj) } else { Vec::new() };
    let mut bond_in_ring = vec![false; m];
    let mut atom_in_ring = vec![false; n];
    for ring in &rings {
        for &b in &ring.bonds {
            bond_in_ring[b] = true;
        }
        for &a in &ring.atoms {
            atom_in_ring[a] = true;
        }
    }
    // Order of every bond that is not searched (template/input pins and propagation), 0 while
    // unknown; lets a leaf read whole rings, not only its own free bonds.
    let mut known_order = vec![0u8; m];

    let mut bond_orders: Vec<Option<BondOrder>> = vec![None; m];
    let mut formal_charges: Vec<Option<i32>> = vec![None; n];
    // Implicit-hydrogen count per atom: 0 everywhere under `AllExplicit`, the per-atom shortfall
    // between valence and bond-order sum under `InferFromGeometry`.
    let mut implicit_hydrogens = vec![0u8; n];
    let mut warnings: Vec<PerceptionWarning> = Vec::new();

    // Global-sized scratch, reused across fragments and reset after each.
    let mut base = vec![0i32; n];
    let mut rem_min = vec![0i32; n];
    let mut rem_max = vec![0i32; n];
    let mut assigned_sum = vec![0i32; n];
    // Per-bond order cap the search uses, narrowed per fragment by propagation. Fragments are
    // bond-disjoint, so one shared array is safe.
    let mut dom_max_eff = dom_max.clone();

    for c in 0..n_comp {
        if comp_bonds[c].is_empty() {
            continue; // an isolated atom: leave its order/charge unchanged
        }
        let comp_at = &comp_atoms[c];

        // Fixed input bonds contribute to the base sum; the rest start as free.
        let mut free_bonds: Vec<usize> = Vec::new();
        for &b in &comp_bonds[c] {
            let [i, j] = pairs[b];
            match fixed[b] {
                Some(o) => {
                    base[i] += o as i32;
                    base[j] += o as i32;
                }
                None => free_bonds.push(b),
            }
        }

        // Propagate valence bounds to fix the forced bonds and narrow the rest before searching.
        let Some((dmin, dmax_p)) =
            propagate_domains(&free_bonds, &pairs, &z, &base, &dom_max, infer)
        else {
            return Err(BondPerceptionError::NoValidAssignment { atom: comp_at[0] });
        };
        let mut all_vars: Vec<usize> = Vec::new();
        let mut prop_fixed: Vec<(usize, u8)> = Vec::new();
        for (p, &b) in free_bonds.iter().enumerate() {
            let [i, j] = pairs[b];
            if dmax_p[p] == dmin[p] {
                // Only one order is possible: treat it as fixed.
                base[i] += dmax_p[p] as i32;
                base[j] += dmax_p[p] as i32;
                prop_fixed.push((b, dmax_p[p]));
            } else {
                dom_max_eff[b] = dmax_p[p];
                all_vars.push(b);
            }
        }
        for &b in &comp_bonds[c] {
            if let Some(o) = fixed[b] {
                known_order[b] = o;
            }
        }
        for &(b, o) in &prop_fixed {
            known_order[b] = o;
        }

        // A whole-molecule charge target couples the free bonds, so when one is set they are
        // solved together — still small after propagation, so still fast on the drug-sized
        // molecules where a target is used. Without a target the free bonds are split into
        // independent clusters (typically one per aromatic ring) that share no atoms, and each
        // is solved on its own atoms; this is what keeps a whole protein tractable.
        let active: std::collections::HashSet<usize> =
            all_vars.iter().flat_map(|&b| pairs[b]).collect();
        let constrained = bonded_comp_count == 1 && options.total_charge.is_some();
        // Inferring hydrogen, an atom whose bonds are all fixed still has a choice: an amine N
        // with one single bond is R-NH2 or R-NH3+. With a charge target that choice cannot be
        // made up front. Instead the clusters stay independent, each is solved for its best
        // assignment *per net charge*, and a small dynamic program over those tables and the
        // fixed atoms' own states then meets the target (see `combine_per_net`). (With explicit
        // hydrogen a fixed atom's valence, and so its charge, is determined, and resolving it up
        // front loses nothing.)
        let joint_states = constrained && infer;
        let clusters = if constrained && !joint_states {
            if all_vars.is_empty() {
                Vec::new()
            } else {
                vec![all_vars.clone()]
            }
        } else {
            cluster_free_bonds(&all_vars, &pairs)
        };
        let single_cluster = constrained && !joint_states;

        // Atoms with no free bond have a determined valence: resolve them directly.
        let mut comp_charge = 0i32;
        for &a in comp_at {
            if active.contains(&a) || joint_states {
                continue;
            }
            match fixed_atom_state(z[a], base[a], infer) {
                FixedState::Wildcard => {}
                FixedState::Infeasible => {
                    return Err(BondPerceptionError::NoValidAssignment { atom: a });
                }
                FixedState::Feasible { fc, implicit_h } => {
                    comp_charge += fc;
                    if fc != input_fc[a] {
                        formal_charges[a] = Some(fc);
                    }
                    implicit_hydrogens[a] = implicit_h;
                }
            }
        }
        let fixed_charge = comp_charge;

        let mut solved_orders: Vec<(usize, u8)> = Vec::new();
        let mut ambiguous_any = false;
        let mut truncated_any = false;
        // Joint mode only: one per-net-charge table per cluster, with that cluster's variables.
        let mut cluster_tables: Vec<(Vec<usize>, NetTable)> = Vec::new();
        for cluster in &clusters {
            let mut cl_atoms: Vec<usize> = Vec::new();
            let mut seen = std::collections::HashSet::new();
            for &b in cluster {
                for a in pairs[b] {
                    if seen.insert(a) {
                        cl_atoms.push(a);
                    }
                }
            }
            cl_atoms.sort_unstable();
            for &b in cluster {
                let [i, j] = pairs[b];
                let dm = dom_max_eff[b] as i32;
                rem_min[i] += 1;
                rem_max[i] += dm;
                rem_min[j] += 1;
                rem_max[j] += dm;
            }
            let cl_vars = connected_order(&cl_atoms, cluster, &pairs);
            let tautomer = if infer {
                TautomerScope::new(&cl_vars, &cl_atoms, &rings, &adj, &z, &known_order)
            } else {
                TautomerScope::default()
            };
            // A whole-molecule charge target can be honored only when a single cluster carries
            // the remaining freedom; otherwise each cluster's charge follows from its own atoms.
            let target = if single_cluster {
                options.total_charge.map(|t| t - fixed_charge)
            } else {
                None
            };
            let mut search = Search {
                pairs: &pairs,
                z: &z,
                dom_max: &dom_max_eff,
                bond_in_ring: &bond_in_ring,
                atom_in_ring: &atom_in_ring,
                rings: &rings,
                adj: &adj,
                known_order: &known_order,
                tautomer,
                sp2_geometry: &sp2_geometry,
                bent: &bent,
                degree: &degree,
                infer,
                vars: &cl_vars,
                atoms: &cl_atoms,
                base: &base,
                rem_min: &mut rem_min,
                rem_max: &mut rem_max,
                assigned_sum: &mut assigned_sum,
                order_out: vec![0u8; cl_vars.len()],
                total_charge: target,
                max_branches: options.limits.max_branches,
                branches: 0,
                best_cost: None,
                best_orders: Vec::new(),
                best_fc: Vec::new(),
                best_ih: Vec::new(),
                per_net: joint_states.then(NetTable::new),
                ambiguous: false,
                truncated: false,
            };
            search.run(0);
            let per_net = search.per_net.take();
            let best_cost = search.best_cost;
            let best_orders = search.best_orders.clone();
            let best_fc = search.best_fc.clone();
            let best_ih = search.best_ih.clone();
            ambiguous_any |= search.ambiguous;
            truncated_any |= search.truncated;
            let truncated = search.truncated;
            drop(search);

            for &a in &cl_atoms {
                rem_min[a] = 0;
                rem_max[a] = 0;
                assigned_sum[a] = 0;
            }

            if let Some(table) = per_net {
                if table.is_empty() {
                    let atom = cl_atoms[0];
                    return Err(if truncated {
                        BondPerceptionError::SearchLimitExceeded { atom }
                    } else {
                        BondPerceptionError::NoValidAssignment { atom }
                    });
                }
                cluster_tables.push((cl_vars, table));
                continue;
            }
            if best_cost.is_none() {
                let atom = cl_atoms[0];
                return Err(if truncated {
                    BondPerceptionError::SearchLimitExceeded { atom }
                } else {
                    BondPerceptionError::NoValidAssignment { atom }
                });
            }
            for (pos, &b) in cl_vars.iter().enumerate() {
                solved_orders.push((b, best_orders[pos]));
            }
            for &(a, fc) in &best_fc {
                comp_charge += fc;
                if fc != input_fc[a] {
                    formal_charges[a] = Some(fc);
                }
            }
            for &(a, h) in &best_ih {
                implicit_hydrogens[a] = h;
            }
        }

        if joint_states {
            // Every atom off the clusters is a one-atom table of its own states.
            let mut tables: Vec<(Vec<usize>, NetTable)> = cluster_tables;
            for &a in comp_at {
                if active.contains(&a) || valence_states(z[a]).is_none() {
                    continue;
                }
                let mut table = NetTable::new();
                for (fc, cost, ih) in atom_state_options(z[a], base[a], sp2_geometry[a], degree[a])
                {
                    let entry = NetEntry {
                        cost,
                        orders: Vec::new(),
                        fcs: vec![(a, fc)],
                        ihs: vec![(a, ih)],
                        ambiguous: false,
                    };
                    offer_net(&mut table, fc, entry);
                }
                if table.is_empty() {
                    return Err(BondPerceptionError::NoValidAssignment { atom: a });
                }
                tables.push((Vec::new(), table));
            }
            let target = options.total_charge.expect("joint mode has a target");
            let Some((picks, ambiguous)) = combine_per_net(&tables, target) else {
                return Err(BondPerceptionError::NoValidAssignment { atom: comp_at[0] });
            };
            ambiguous_any |= ambiguous;
            for ((vars, _), entry) in tables.iter().zip(picks) {
                for (pos, &b) in vars.iter().enumerate() {
                    solved_orders.push((b, entry.orders[pos]));
                }
                for &(a, fc) in &entry.fcs {
                    comp_charge += fc;
                    if fc != input_fc[a] {
                        formal_charges[a] = Some(fc);
                    }
                }
                for &(a, h) in &entry.ihs {
                    implicit_hydrogens[a] = h;
                }
            }
        }

        // Write out every bond whose solved order differs from the input (`None` = unchanged):
        // the searched bonds, the bonds propagation fixed, and the originally-fixed bonds.
        let mut emit = |b: usize, order: u8| {
            let solved = int_to_order(order);
            if solved != input_orders[b] {
                bond_orders[b] = Some(solved);
            }
        };
        for &(b, order) in &solved_orders {
            emit(b, order);
        }
        for &(b, order) in &prop_fixed {
            emit(b, order);
        }
        for &b in &comp_bonds[c] {
            if let Some(o) = fixed[b] {
                emit(b, o);
            }
        }

        if ambiguous_any {
            warnings.push(PerceptionWarning::AmbiguousAssignment);
        }
        if options.total_charge.is_none() && comp_charge != 0 {
            warnings.push(PerceptionWarning::ChargeWasNotConstrained);
        }
        if truncated_any {
            warnings.push(PerceptionWarning::SearchTruncated);
        }

        // Restore the scratch this fragment touched.
        for &a in comp_at {
            base[a] = 0;
            rem_min[a] = 0;
            rem_max[a] = 0;
            assigned_sum[a] = 0;
        }
    }

    warnings.dedup();
    BondAssignment::new(top, bond_orders, formal_charges, implicit_hydrogens, warnings)
}

/// Reorder `vars` (free bond indices) into a connected traversal: a BFS over the fragment's
/// atoms that emits each atom's incident free bonds when the atom is first reached. Adjacent
/// decisions then share atoms, so an atom closes (all its bonds assigned) early and the
/// feasibility bound prunes as soon as its valence is fixed. Deterministic.
fn connected_order(atoms: &[usize], vars: &[usize], pairs: &[[usize; 2]]) -> Vec<usize> {
    use std::collections::{HashMap, HashSet, VecDeque};
    let mut inc: HashMap<usize, Vec<usize>> = HashMap::new();
    for &b in vars {
        let [i, j] = pairs[b];
        inc.entry(i).or_default().push(b);
        inc.entry(j).or_default().push(b);
    }
    for v in inc.values_mut() {
        v.sort_unstable();
    }
    let mut ordered = Vec::with_capacity(vars.len());
    let mut emitted: HashSet<usize> = HashSet::new();
    let mut visited: HashSet<usize> = HashSet::new();
    let mut queue = VecDeque::new();
    for &start in atoms {
        // `atoms` is ascending, so this is deterministic.
        if !visited.insert(start) {
            continue;
        }
        queue.push_back(start);
        while let Some(a) = queue.pop_front() {
            let Some(incident) = inc.get(&a) else {
                continue;
            };
            for &b in incident {
                if emitted.insert(b) {
                    ordered.push(b);
                }
                let [i, j] = pairs[b];
                let other = if i == a { j } else { i };
                if visited.insert(other) {
                    queue.push_back(other);
                }
            }
        }
    }
    ordered
}

/// Smallest bond angle, in degrees, still read as linear (sp). Real sp centers sit at 170-180°;
/// sp2 at ~120°, so 150° separates them with margin for crystal-structure noise.
const SP_MIN_ANGLE_DEG: Float = 150.0;

/// Penalty for an atom carrying more than one ring double bond (a cumulated aromatic ring, as
/// opposed to a Kekulé structure with one double per atom). Large so it dominates the state
/// penalties, but applied only when inferring hydrogens (see [`Search::record_leaf`]).
const RING_CUMULENE_PENALTY: u32 = 100;

/// Bounds-consistency propagation over a fragment's free bonds. Given the fixed-bond
/// contribution per atom (`base`) and each free bond's order cap (`dom_cap`, indexed by global
/// bond), it repeatedly tightens each free bond's `[min, max]` order range from the valence its
/// two endpoints can still reach, to a fixpoint.
///
/// This is what makes a large fragment (a whole protein, one connected piece once connectivity
/// is perceived) tractable: a valence-saturated atom — a backbone or aliphatic carbon carrying
/// its hydrogens — forces every incident bond to single, and that propagates, so only genuinely
/// free clusters (aromatic rings) are left for the search.
///
/// Returns the per-free-bond `(min, max)` order ranges (parallel to `free_bonds`), or `None` if
/// some atom has no reachable valence state (the fragment cannot be solved).
fn propagate_domains(
    free_bonds: &[usize],
    pairs: &[[usize; 2]],
    z: &[u8],
    base: &[i32],
    dom_cap: &[u8],
    infer: bool,
) -> Option<(Vec<u8>, Vec<u8>)> {
    use std::collections::HashMap;
    let mut dmin = vec![1u8; free_bonds.len()];
    let mut dmax: Vec<u8> = free_bonds.iter().map(|&b| dom_cap[b]).collect();
    let mut incident: HashMap<usize, Vec<usize>> = HashMap::new();
    for (p, &b) in free_bonds.iter().enumerate() {
        let [i, j] = pairs[b];
        incident.entry(i).or_default().push(p);
        incident.entry(j).or_default().push(p);
    }

    loop {
        let mut changed = false;
        for (&a, positions) in &incident {
            let Some(states) = valence_states(z[a]) else {
                continue; // unmodeled element: no valence constraint
            };
            let lo = base[a] + positions.iter().map(|&p| dmin[p] as i32).sum::<i32>();
            let hi = base[a] + positions.iter().map(|&p| dmax[p] as i32).sum::<i32>();
            let mut vmin = i32::MAX;
            let mut vmax = i32::MIN;
            for st in states {
                let v = st.valence as i32;
                let fits = if infer { v >= lo } else { v >= lo && v <= hi };
                if fits {
                    vmin = vmin.min(v);
                    vmax = vmax.max(v);
                }
            }
            if vmax == i32::MIN {
                return None; // no reachable valence state
            }
            // Inferring hydrogens only bounds the sum from above (hydrogen fills any shortfall),
            // and the sum can never exceed `hi`.
            let vmax = if infer { vmax.min(hi) } else { vmax };
            for &p in positions {
                // Raising this bond alone (others at their minimum) must keep the sum <= vmax.
                let new_max = (dmin[p] as i32 + (vmax - lo)).clamp(dmin[p] as i32, dmax[p] as i32);
                if (new_max as u8) < dmax[p] {
                    dmax[p] = new_max as u8;
                    changed = true;
                }
                // With explicit hydrogen the sum must also reach vmin, so lowering this bond
                // alone (others at their maximum) is bounded from below.
                if !infer {
                    let new_min =
                        (dmax[p] as i32 - (hi - vmin)).clamp(dmin[p] as i32, dmax[p] as i32);
                    if (new_min as u8) > dmin[p] {
                        dmin[p] = new_min as u8;
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    Some((dmin, dmax))
}

/// Group the free bonds into connected clusters — maximal sets that share atoms. Distinct
/// clusters have no atom in common, so each is an independent bond-order sub-problem. Bonds are
/// returned in a deterministic order (clusters seeded by ascending position in `vars`).
fn cluster_free_bonds(vars: &[usize], pairs: &[[usize; 2]]) -> Vec<Vec<usize>> {
    use std::collections::{HashMap, VecDeque};
    let mut incident: HashMap<usize, Vec<usize>> = HashMap::new();
    for (p, &b) in vars.iter().enumerate() {
        let [i, j] = pairs[b];
        incident.entry(i).or_default().push(p);
        incident.entry(j).or_default().push(p);
    }
    let mut cluster_of = vec![usize::MAX; vars.len()];
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for start in 0..vars.len() {
        if cluster_of[start] != usize::MAX {
            continue;
        }
        let id = clusters.len();
        cluster_of[start] = id;
        let mut members = Vec::new();
        let mut queue = VecDeque::from([start]);
        while let Some(p) = queue.pop_front() {
            members.push(vars[p]);
            for a in pairs[vars[p]] {
                for &q in &incident[&a] {
                    if cluster_of[q] == usize::MAX {
                        cluster_of[q] = id;
                        queue.push_back(q);
                    }
                }
            }
        }
        clusters.push(members);
    }
    clusters
}

/// The resolution of an atom whose incident bonds are all determined (`bosum` known).
enum FixedState {
    /// An element the tables do not model: leave its charge unchanged, add no hydrogen.
    Wildcard,
    /// No valence state matches the determined bond-order sum.
    Infeasible,
    /// A formal charge and implicit-hydrogen count.
    Feasible { fc: i32, implicit_h: u8 },
}

/// Resolve an atom whose bond-order sum is fixed at `bosum`: pick its minimum-penalty valence
/// state (smallest valence, charge nearest neutral on a tie) and read off the charge and the
/// implicit-hydrogen shortfall. Mirrors [`Search::record_leaf`] for the single-atom case.
fn fixed_atom_state(z: u8, bosum: i32, infer: bool) -> FixedState {
    let Some(states) = valence_states(z) else {
        return FixedState::Wildcard;
    };
    let fits = |v: i32| if infer { v >= bosum } else { v == bosum };
    let Some(min_pen) = states
        .iter()
        .filter(|s| fits(s.valence as i32))
        .map(|s| s.penalty)
        .min()
    else {
        return FixedState::Infeasible;
    };
    let valence = states
        .iter()
        .filter(|s| fits(s.valence as i32) && s.penalty == min_pen)
        .map(|s| s.valence as i32)
        .min()
        .unwrap();
    let fc = states
        .iter()
        .filter(|s| s.valence as i32 == valence && s.penalty == min_pen)
        .map(|s| s.formal_charge)
        .min_by_key(|c| c.unsigned_abs())
        .unwrap();
    FixedState::Feasible {
        fc,
        implicit_h: (valence - bosum) as u8,
    }
}

/// The result of [`Search::bound`]: a penalty lower bound and two net-charge ranges (see there).
struct Bound {
    penalty: u32,
    full_lo: i32,
    full_hi: i32,
    tight_lo: i32,
    tight_hi: i32,
}

/// One connected fragment's branch-and-bound over the free bond orders.
struct Search<'a> {
    pairs: &'a [[usize; 2]],
    z: &'a [u8],
    dom_max: &'a [u8],
    /// Whether each bond lies on a ring (used only when inferring hydrogens).
    bond_in_ring: &'a [bool],
    /// Per global atom: whether it lies on a ring.
    atom_in_ring: &'a [bool],
    /// The fragment's SSSR rings (only when inferring).
    rings: &'a [super::RingData],
    adj: &'a BondAdjacency,
    /// Order of each bond outside the search, 0 if unknown (see `known_order`).
    known_order: &'a [u8],
    /// What this cluster's leaves score for the tautomer tie-break.
    tautomer: TautomerScope,
    /// Per global atom: a C or O whose geometry rules out saturation (used only when inferring).
    sp2_geometry: &'a [bool],
    /// Per global atom: a two-coordinate C or N whose bond angle is not linear.
    bent: &'a [bool],
    /// Number of bonds per global atom.
    degree: &'a [i32],
    /// When true, an atom's bond-order sum may fall short of its valence, the shortfall being
    /// implicit hydrogen; otherwise the sum must reach a valence state exactly.
    infer: bool,
    /// Global bond indices of the free (unfixed) bonds, in connected-walk order.
    vars: &'a [usize],
    /// Global atom indices in this fragment.
    atoms: &'a [usize],
    /// Bond-order sum contributed by the fixed bonds, per global atom.
    base: &'a [i32],
    /// Sum of the minimum order (1) of each atom's still-unassigned incident free bonds.
    rem_min: &'a mut [i32],
    /// Sum of the maximum order of each atom's still-unassigned incident free bonds.
    rem_max: &'a mut [i32],
    /// Order already committed on assigned free bonds, per global atom.
    assigned_sum: &'a mut [i32],
    /// The order chosen for each free bond, parallel to `vars`.
    order_out: Vec<u8>,
    total_charge: Option<i32>,
    max_branches: u64,
    branches: u64,
    /// Best objective found: `(summed penalty, net-charge magnitude, geometry contradictions,
    /// total implicit H, non-aromatic rings, terminal imines, non-amide N)`. The penalty (state
    /// penalties plus any ring-cumulene penalty) is primary; the charge magnitude breaks ties
    /// toward the least charge-separated form. The contradiction count then rejects what the
    /// geometry rules out: a C or O left saturated although its geometry is sp2/sp (an sp3 CG in
    /// an indole, an O-H on a lactam C=O), and a bent atom given two π bonds (HN=C=O read into
    /// formamide). The implicit-hydrogen count breaks the remaining ties toward the most-saturated
    /// structure the geometry allows. The contradiction count is a tie-break, not a penalty, so
    /// noisy geometry can never force a charged form. The last three terms only separate
    /// tautomers, which tie on everything before them (see [`Search::tautomer_terms`]).
    best_cost: Option<(u32, u32, u32, u32, u32, u32, u32)>,
    best_orders: Vec<u8>,
    best_fc: Vec<(usize, i32)>,
    best_ih: Vec<(usize, u8)>,
    /// Joint mode (inferring hydrogen under a charge target): instead of one best leaf, the best
    /// leaf for every net charge the cluster can take. `best_cost` is then unused.
    per_net: Option<NetTable>,
    ambiguous: bool,
    truncated: bool,
}

impl Search<'_> {
    fn run(&mut self, k: usize) {
        self.branches += 1;
        if self.branches > self.max_branches {
            self.truncated = true;
            return;
        }

        // Feasibility + admissible lower bound on the primary (penalty) objective, plus the
        // achievable net-charge ranges over the fragment.
        let Some(bound) = self.bound() else {
            return;
        };
        let lb = bound.penalty;
        // With a fixed total charge, prune as soon as the atoms still open cannot bring the net
        // charge to the target. This is the key cut on large fragments: without it the search
        // descends to a full leaf before the leaf-level charge check rejects it.
        if let Some(target) = self.total_charge
            && (target < bound.full_lo || target > bound.full_hi)
        {
            return;
        }

        // Per-net mode: a branch can only land on a net charge in `[full_lo, full_hi]`. Once every
        // such charge already has a leaf strictly cheaper than this branch's penalty bound, no
        // completion can improve any entry.
        if let Some(table) = &self.per_net
            && !table.is_empty()
            && (bound.full_lo..=bound.full_hi)
                .all(|net| table.get(&net).is_some_and(|e| e.cost.0 < lb))
        {
            return;
        }

        if let Some((best_penalty, best_dev, ..)) = self.best_cost {
            // The penalty lower bound alone rules this branch out.
            if lb > best_penalty {
                return;
            }
            // For an equal-penalty branch, bound the secondary cost (net-charge magnitude) from
            // the charge the minimum-penalty states can still reach: with a target every tying
            // leaf has magnitude 0, otherwise the smallest reachable magnitude. Once a solution
            // is known (ambiguity seen) and this branch cannot beat it on penalty or charge,
            // prune it — this stops the tie explosion on large poly-aromatic, multiply-charged
            // fragments (whole proteins). Skipped when inferring hydrogens, where the
            // implicit-hydrogen tie-break could still improve such a branch, and those fragments
            // are small.
            let dev_lb = match self.total_charge {
                Some(_) => 0,
                None if bound.tight_lo <= 0 && bound.tight_hi >= 0 => 0,
                None => bound.tight_lo.unsigned_abs().min(bound.tight_hi.unsigned_abs()),
            };
            if !self.infer && lb == best_penalty && dev_lb >= best_dev && self.ambiguous {
                return;
            }
        }

        if k == self.vars.len() {
            self.record_leaf();
            return;
        }

        let b = self.vars[k];
        let [i, j] = self.pairs[b];
        let dm = self.dom_max[b];
        self.rem_min[i] -= 1;
        self.rem_max[i] -= dm as i32;
        self.rem_min[j] -= 1;
        self.rem_max[j] -= dm as i32;
        for o in 1..=dm {
            self.assigned_sum[i] += o as i32;
            self.assigned_sum[j] += o as i32;
            self.order_out[k] = o;
            self.run(k + 1);
            self.assigned_sum[i] -= o as i32;
            self.assigned_sum[j] -= o as i32;
        }
        self.rem_min[i] += 1;
        self.rem_max[i] += dm as i32;
        self.rem_min[j] += 1;
        self.rem_max[j] += dm as i32;
    }

    /// Whether a candidate valence state with `valence` fits an atom whose bond-order sum is
    /// bounded to `[lo, hi]`: exactly in range when hydrogens are explicit, or at least `lo`
    /// when inferring (implicit hydrogen makes up any shortfall).
    fn state_fits(&self, valence: i32, lo: i32, hi: i32) -> bool {
        if self.infer {
            valence >= lo
        } else {
            valence >= lo && valence <= hi
        }
    }

    /// For the current partial assignment: the least summed penalty any completion could reach
    /// (an admissible lower bound), the range of net formal charge any completion allows
    /// (`full_*`, for target feasibility), and the tighter range that only the **minimum-penalty**
    /// states allow (`tight_*`). A completion that ties the penalty lower bound must use each
    /// atom's minimum-penalty state, so its net charge lies in the tight range — which is what
    /// bounds the charge tie-break. `None` if some atom already has no reachable valence state.
    fn bound(&self) -> Option<Bound> {
        let mut lb = 0u32;
        let (mut full_lo, mut full_hi) = (0i32, 0i32);
        let (mut tight_lo, mut tight_hi) = (0i32, 0i32);
        for &a in self.atoms {
            let Some(states) = valence_states(self.z[a]) else {
                continue; // unmodeled element: no constraint, no penalty, no charge
            };
            let lo = self.base[a] + self.assigned_sum[a] + self.rem_min[a];
            let hi = self.base[a] + self.assigned_sum[a] + self.rem_max[a];
            let mut best_pen: Option<u32> = None;
            let (mut min_fc, mut max_fc) = (i32::MAX, i32::MIN);
            for st in states {
                if self.state_fits(st.valence as i32, lo, hi) {
                    best_pen = Some(best_pen.map_or(st.penalty, |p| p.min(st.penalty)));
                    min_fc = min_fc.min(st.formal_charge);
                    max_fc = max_fc.max(st.formal_charge);
                }
            }
            let best_pen = best_pen?;
            lb += best_pen;
            full_lo += min_fc;
            full_hi += max_fc;
            // Charge range over only the minimum-penalty states in the window.
            let (mut tmin, mut tmax) = (i32::MAX, i32::MIN);
            for st in states {
                if st.penalty == best_pen && self.state_fits(st.valence as i32, lo, hi) {
                    tmin = tmin.min(st.formal_charge);
                    tmax = tmax.max(st.formal_charge);
                }
            }
            tight_lo += tmin;
            tight_hi += tmax;
        }
        Some(Bound {
            penalty: lb,
            full_lo,
            full_hi,
            tight_lo,
            tight_hi,
        })
    }

    /// The number of ring double (or higher) bonds incident to each atom, for the cumulene
    /// guard. Fixed bonds are single when inferring, so only the free bonds contribute.
    fn ring_double_counts(&self) -> std::collections::HashMap<usize, u32> {
        let mut counts = std::collections::HashMap::new();
        for (pos, &b) in self.vars.iter().enumerate() {
            if self.order_out[pos] >= 2 && self.bond_in_ring[b] {
                let [i, j] = self.pairs[b];
                *counts.entry(i).or_insert(0) += 1;
                *counts.entry(j).or_insert(0) += 1;
            }
        }
        counts
    }

    fn record_leaf(&mut self) {
        let ring_doubles = if self.infer {
            self.ring_double_counts()
        } else {
            std::collections::HashMap::new()
        };
        if self.per_net.is_some() {
            self.record_leaf_per_net(&ring_doubles);
            return;
        }
        let target = self.total_charge.unwrap_or(0);

        let mut penalty = 0u32;
        let mut net = 0i32;
        let mut fcs: Vec<(usize, i32)> = Vec::with_capacity(self.atoms.len());
        let mut ihs: Vec<(usize, u8)> = Vec::with_capacity(self.atoms.len());
        let mut charge_choice_tie = false;
        let mut contradictions = 0u32;

        for &a in self.atoms {
            let Some(states) = valence_states(self.z[a]) else {
                continue;
            };
            let bosum = self.base[a] + self.assigned_sum[a];
            // Candidate states: valence == bosum (explicit) or valence >= bosum (inferring).
            let min_pen = states
                .iter()
                .filter(|s| self.state_fits(s.valence as i32, bosum, bosum))
                .map(|s| s.penalty)
                .min();
            // `bound()` proved a candidate exists before this leaf was reached.
            let min_pen =
                min_pen.expect("leaf reached with an atom that has no candidate valence state");
            // Among min-penalty candidates prefer the smallest valence (fewest implicit
            // hydrogens); then the formal charge is chosen to steer the net toward the target.
            let valence = states
                .iter()
                .filter(|s| self.state_fits(s.valence as i32, bosum, bosum) && s.penalty == min_pen)
                .map(|s| s.valence as i32)
                .min()
                .unwrap();
            let mut fc_options: Vec<i32> = states
                .iter()
                .filter(|s| s.valence as i32 == valence && s.penalty == min_pen)
                .map(|s| s.formal_charge)
                .collect();
            fc_options.dedup();
            if fc_options.len() > 1 {
                charge_choice_tie = true;
            }
            let fc = *fc_options
                .iter()
                .min_by_key(|&&c| (target - (net + c)).abs())
                .expect("every atom has at least one candidate charge");
            net += fc;

            penalty += min_pen;
            if let Some(&count) = ring_doubles.get(&a)
                && count >= 2
            {
                penalty += (count - 1) * RING_CUMULENE_PENALTY;
            }
            if geometry_contradiction(
                self.z[a],
                self.sp2_geometry[a],
                valence,
                fc,
                bosum,
                self.degree[a],
            ) {
                contradictions += 1;
            }
            fcs.push((a, fc));
            ihs.push((a, (valence - bosum) as u8));
        }

        let tautomer = self.tautomer_terms();
        self.offer(
            penalty,
            charge_dev_of(self.total_charge, net),
            contradictions + self.linear_violations(),
            tautomer,
            fcs,
            ihs,
            charge_choice_tie,
        );
    }

    /// Bent atoms given two π bonds at this leaf: a geometry contradiction that depends on the
    /// bond orders, not on the atom's state, so it is counted per leaf.
    fn linear_violations(&self) -> u32 {
        self.atoms
            .iter()
            .filter(|&&a| {
                self.bent[a]
                    && self
                        .adj
                        .neighbors(a)
                        .iter()
                        .map(|nb| self.leaf_order(nb.bond()) as u32 - 1)
                        .sum::<u32>()
                        >= 2
            })
            .count() as u32
    }

    /// Order of bond `b` at this leaf: searched, or known from outside the search.
    fn leaf_order(&self, b: usize) -> u8 {
        match self.tautomer.var_pos.get(&b) {
            Some(&pos) => self.order_out[pos],
            None => self.known_order[b].max(1),
        }
    }

    /// The tautomer tie-break, `(non-aromatic rings, terminal imines, non-amide N)`, all to be
    /// minimized. Two tautomers carry the same hydrogens, charges and penalty, so everything
    /// before these in the cost ties; chemistry then prefers, in this order,
    /// 1. the form that keeps more rings aromatic (molar's own Hückel test) — purine N7/N9-H
    ///    over N1/N3-H, which breaks a ring;
    /// 2. the amino over the imino form: no C=N to an N without another heavy neighbour —
    ///    cytosine and guanine keep their NH2, not a ring N-H plus an exocyclic C=NH;
    /// 3. the lactam: an N-H next to a C=O rather than elsewhere — guanine N1-H over N3-H,
    ///    3H-quinazolin-4-one over the 1H form. This comes after (2), which it would otherwise
    ///    overrule: imino-oxo cytosine has one more N-H beside its C=O than the real amino form.
    ///
    /// Measured over the rings and nitrogens this cluster can change ([`TautomerScope`]).
    fn tautomer_terms(&self) -> (u32, u32, u32) {
        if !self.infer {
            return (0, 0, 0);
        }
        let as_order = |o: u8| match o {
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            _ => BondOrder::Single,
        };
        let non_aromatic = self
            .tautomer
            .ring_ids
            .iter()
            .filter(|&&r| {
                !super::ring_is_huckel_aromatic(
                    &self.rings[r],
                    |b| as_order(self.leaf_order(b)),
                    self.adj,
                    self.z,
                    self.atom_in_ring,
                )
            })
            .count() as u32;
        let non_amide = self
            .tautomer
            .nitrogens
            .iter()
            .filter(|&&nn| {
                let nbs = self.adj.neighbors(nn);
                nbs.iter().all(|nb| self.leaf_order(nb.bond()) == 1)
                    && !nbs.iter().any(|nb| {
                        self.z[nb.atom()] == 6
                            && self.adj.neighbors(nb.atom()).iter().any(|cb| {
                                self.z[cb.atom()] == 8
                                    && self.degree[cb.atom()] == 1
                                    && self.leaf_order(cb.bond()) == 2
                            })
                    })
            })
            .count() as u32;
        let terminal_imines = self
            .tautomer
            .nitrogens
            .iter()
            .filter(|&&nn| {
                let nbs = self.adj.neighbors(nn);
                nbs.iter().filter(|nb| self.z[nb.atom()] != 1).count() == 1
                    && nbs
                        .iter()
                        .any(|nb| self.z[nb.atom()] == 6 && self.leaf_order(nb.bond()) == 2)
            })
            .count() as u32;
        (non_aromatic, terminal_imines, non_amide)
    }

    /// Per-net leaf: for this bond assignment, the cheapest choice of atom states for every
    /// reachable net charge, by a dynamic program over the running charge minimizing
    /// `(penalty, contradictions, implicit H)`, with the leaf's tautomer terms added. Each result
    /// is offered to `per_net`.
    fn record_leaf_per_net(&mut self, ring_doubles: &std::collections::HashMap<usize, u32>) {
        let mut cumulene = 0u32;
        let mut units: Vec<Vec<(i32, Cost, usize, u8)>> = Vec::new();
        for &a in self.atoms {
            if valence_states(self.z[a]).is_none() {
                continue;
            }
            if let Some(&count) = ring_doubles.get(&a)
                && count >= 2
            {
                cumulene += (count - 1) * RING_CUMULENE_PENALTY;
            }
            let bosum = self.base[a] + self.assigned_sum[a];
            let opts = atom_state_options(self.z[a], bosum, self.sp2_geometry[a], self.degree[a]);
            units.push(
                opts.into_iter()
                    .map(|(fc, cost, ih)| (fc, cost, a, ih))
                    .collect(),
            );
        }
        let (non_aromatic, terminal_imines, non_amide) = self.tautomer_terms();
        let linear = self.linear_violations();
        for (net, cost, ways, picks) in charge_dp(&units) {
            let fcs = picks.iter().map(|&(fc, a, _)| (a, fc)).collect();
            let ihs = picks.iter().map(|&(_, a, ih)| (a, ih)).collect();
            let entry = NetEntry {
                cost: add_cost(
                    cost,
                    (
                        cumulene,
                        linear,
                        0,
                        non_aromatic,
                        terminal_imines,
                        non_amide,
                    ),
                ),
                orders: self.order_out.clone(),
                fcs,
                ihs,
                ambiguous: ways > 1,
            };
            offer_net(self.per_net.as_mut().expect("per-net mode"), net, entry);
        }
    }

    /// Compare a completed leaf with the best so far and keep it if it is better.
    fn offer(
        &mut self,
        penalty: u32,
        charge_dev: Option<u32>,
        contradictions: u32,
        (non_aromatic, terminal_imines, non_amide): (u32, u32, u32),
        fcs: Vec<(usize, i32)>,
        ihs: Vec<(usize, u8)>,
        charge_choice_tie: bool,
    ) {
        let Some(charge_dev) = charge_dev else {
            return; // missed the charge target
        };
        let implicit_total: u32 = ihs.iter().map(|&(_, h)| h as u32).sum();
        let cost = (
            penalty,
            charge_dev,
            contradictions,
            implicit_total,
            non_aromatic,
            terminal_imines,
            non_amide,
        );

        match self.best_cost {
            Some(best) if cost > best => {}
            Some(best) if cost == best => {
                if self.order_out != self.best_orders || charge_choice_tie {
                    self.ambiguous = true;
                }
            }
            _ => {
                self.best_cost = Some(cost);
                self.best_orders = self.order_out.clone();
                self.best_fc = fcs;
                self.best_ih = ihs;
                self.ambiguous = charge_choice_tie;
            }
        }
    }
}

/// Per-net cost: `(penalty, geometry contradictions, implicit H, non-aromatic rings, terminal
/// imines, non-amide N)`, compared lexicographically — the search's cost without the charge term, which the table
/// key replaces.
type Cost = (u32, u32, u32, u32, u32, u32);

fn add_cost(a: Cost, b: Cost) -> Cost {
    (
        a.0 + b.0,
        a.1 + b.1,
        a.2 + b.2,
        a.3 + b.3,
        a.4 + b.4,
        a.5 + b.5,
    )
}

/// What one cluster's leaves score for the tautomer tie-break: the rings and nitrogens whose
/// every bond is either searched in this cluster or known from outside the search, and which the
/// cluster touches. A ring or N that also depends on another cluster's bonds is left out, so each
/// cluster's score is exact over what it measures.
#[derive(Default)]
struct TautomerScope {
    /// Position of each of the cluster's variables in its walk order, by bond index.
    var_pos: std::collections::HashMap<usize, usize>,
    /// Indices into the fragment's ring list.
    ring_ids: Vec<usize>,
    nitrogens: Vec<usize>,
}

impl TautomerScope {
    fn new(
        vars: &[usize],
        atoms: &[usize],
        rings: &[super::RingData],
        adj: &BondAdjacency,
        z: &[u8],
        known_order: &[u8],
    ) -> Self {
        let var_pos: std::collections::HashMap<usize, usize> =
            vars.iter().enumerate().map(|(p, &b)| (b, p)).collect();
        let known = |b: usize| known_order[b] > 0 || var_pos.contains_key(&b);
        let settled = |a: usize| adj.neighbors(a).iter().all(|nb| known(nb.bond()));
        let touches = |a: usize| atoms.binary_search(&a).is_ok();
        let ring_ids = (0..rings.len())
            .filter(|&r| {
                rings[r].bonds.iter().any(|b| var_pos.contains_key(b))
                    && rings[r].atoms.iter().all(|&a| settled(a))
            })
            .collect();
        // An amide N needs its own bonds and its carbon neighbours' bonds settled.
        let mut nitrogens: Vec<usize> = atoms
            .iter()
            .flat_map(|&a| std::iter::once(a).chain(adj.neighbors(a).iter().map(|nb| nb.atom())))
            .filter(|&a| z[a] == 7)
            .filter(|&a| {
                settled(a)
                    && adj
                        .neighbors(a)
                        .iter()
                        .all(|nb| z[nb.atom()] != 6 || settled(nb.atom()))
                    && (touches(a) || adj.neighbors(a).iter().any(|nb| touches(nb.atom())))
            })
            .collect();
        nitrogens.sort_unstable();
        nitrogens.dedup();
        Self {
            var_pos,
            ring_ids,
            nitrogens,
        }
    }
}

/// A cluster's (or fixed atom's) best assignment at one net charge.
#[derive(Clone)]
struct NetEntry {
    cost: Cost,
    /// Orders for the cluster's variables, in its walk order (empty for a fixed atom).
    orders: Vec<u8>,
    fcs: Vec<(usize, i32)>,
    ihs: Vec<(usize, u8)>,
    /// Another assignment ties this one.
    ambiguous: bool,
}

type NetTable = std::collections::BTreeMap<i32, NetEntry>;

/// Keep `entry` at `net` if it beats what is there; a tie marks the entry ambiguous.
fn offer_net(table: &mut NetTable, net: i32, entry: NetEntry) {
    match table.get_mut(&net) {
        Some(e) if entry.cost > e.cost => {}
        Some(e) if entry.cost == e.cost => e.ambiguous = true,
        _ => {
            table.insert(net, entry);
        }
    }
}

/// Whether a state contradicts the geometry: a C or O with a bond short enough for a multiple
/// order (`sp2_geometry`), left saturated — all bonds single, neutral, the rest implicit H. That
/// is an sp3 CG in an indole or a C-OH on a 0.123 nm lactam C=O. N is exempt: a pyrrole or
/// aniline N-H is saturated with short bonds. A charged state is exempt too (a phenolate's short
/// C-O⁻), and since this is only a tie-break it never forces a charge either.
fn geometry_contradiction(
    z: u8,
    sp2_geometry: bool,
    valence: i32,
    fc: i32,
    bosum: i32,
    degree: i32,
) -> bool {
    let neutral_valence = match z {
        6 => 4,
        8 => 2,
        _ => return false,
    };
    sp2_geometry && fc == 0 && bosum == degree && valence == neutral_valence
}

/// Every valence state an atom with bond-order sum `bosum` can take when inferring hydrogen, as
/// `(formal charge, cost, implicit H)`.
fn atom_state_options(z: u8, bosum: i32, sp2_geometry: bool, degree: i32) -> Vec<(i32, Cost, u8)> {
    let Some(states) = valence_states(z) else {
        return Vec::new();
    };
    states
        .iter()
        .filter(|st| st.valence as i32 >= bosum)
        .map(|st| {
            let valence = st.valence as i32;
            let contradiction =
                geometry_contradiction(z, sp2_geometry, valence, st.formal_charge, bosum, degree)
                    as u32;
            let ih = (valence - bosum) as u8;
            (
                st.formal_charge,
                (st.penalty, contradiction, ih as u32, 0, 0, 0),
                ih,
            )
        })
        .collect()
}

/// Pick one option per unit to reach each reachable net charge at least cost. Returns, per net
/// charge, `(net, cost, optimal ways capped at 2, picks)` with picks as `(fc, atom, ih)`.
fn charge_dp(units: &[Vec<(i32, Cost, usize, u8)>]) -> Vec<(i32, Cost, u8, Vec<(i32, usize, u8)>)> {
    use std::collections::BTreeMap;
    // net -> (cost, ways, back-pointer into `steps`)
    let mut layer: BTreeMap<i32, (Cost, u8, usize)> =
        BTreeMap::from([(0, ((0, 0, 0, 0, 0, 0), 1, usize::MAX))]);
    let mut steps: Vec<((i32, usize, u8), usize)> = Vec::new();
    for opts in units {
        let mut next: BTreeMap<i32, (Cost, u8, usize)> = BTreeMap::new();
        for (&net, &(cost, ways, back)) in &layer {
            for &(fc, c, a, ih) in opts {
                let total = add_cost(cost, c);
                match next.get_mut(&(net + fc)) {
                    Some(e) if total > e.0 => {}
                    Some(e) if total == e.0 => e.1 = (e.1 + ways).min(2),
                    _ => {
                        steps.push(((fc, a, ih), back));
                        next.insert(net + fc, (total, ways, steps.len() - 1));
                    }
                }
            }
        }
        layer = next;
    }
    layer
        .into_iter()
        .map(|(net, (cost, ways, mut at))| {
            let mut picks = Vec::new();
            while at != usize::MAX {
                picks.push(steps[at].0);
                at = steps[at].1;
            }
            picks.reverse();
            (net, cost, ways, picks)
        })
        .collect()
}

/// Meet the charge `target` by picking one entry from each table, minimizing the summed cost.
/// Returns the picks (parallel to `tables`) and whether the choice is ambiguous: another
/// combination ties, or a picked entry was itself tied.
fn combine_per_net(
    tables: &[(Vec<usize>, NetTable)],
    target: i32,
) -> Option<(Vec<NetEntry>, bool)> {
    use std::collections::BTreeMap;
    let mut layer: BTreeMap<i32, (Cost, u8, usize)> =
        BTreeMap::from([(0, ((0, 0, 0, 0, 0, 0), 1, usize::MAX))]);
    // (unit, net charge chosen in that unit, back-pointer)
    let mut steps: Vec<(usize, i32, usize)> = Vec::new();
    for (u, (_, table)) in tables.iter().enumerate() {
        let mut next: BTreeMap<i32, (Cost, u8, usize)> = BTreeMap::new();
        for (&net, &(cost, ways, back)) in &layer {
            for (&q, e) in table {
                let total = add_cost(cost, e.cost);
                match next.get_mut(&(net + q)) {
                    Some(x) if total > x.0 => {}
                    Some(x) if total == x.0 => x.1 = (x.1 + ways).min(2),
                    _ => {
                        steps.push((u, q, back));
                        next.insert(net + q, (total, ways, steps.len() - 1));
                    }
                }
            }
        }
        layer = next;
    }
    let &(_, ways, mut at) = layer.get(&target)?;
    let mut picks: Vec<Option<NetEntry>> = vec![None; tables.len()];
    let mut ambiguous = ways > 1;
    while at != usize::MAX {
        let (u, q, back) = steps[at];
        let e = &tables[u].1[&q];
        ambiguous |= e.ambiguous;
        picks[u] = Some(e.clone());
        at = back;
    }
    Some((
        picks
            .into_iter()
            .map(|p| p.expect("every table picked once"))
            .collect(),
        ambiguous,
    ))
}

/// The charge tie-break term: 0 on target (or `None` when the target is missed), else the net
/// magnitude when no target is set.
fn charge_dev_of(target: Option<i32>, net: i32) -> Option<u32> {
    match target {
        Some(t) => (net == t).then_some(0),
        None => Some(net.unsigned_abs()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn topo(z: &[u8], bonds: &[(usize, usize, BondOrder)]) -> Topology {
        let mut top = Topology::default();
        for &zi in z {
            top.atoms.push(&Atom::new().with_atomic_number(zi));
        }
        for &(i, j, o) in bonds {
            top.bonds.push(&Bond::with_order(i, j, o));
        }
        top
    }

    fn solve(top: &Topology, options: &BondOrderOptions) -> BondAssignment {
        assign_bond_orders(top, None, options).unwrap()
    }

    /// The final order per bond after applying the assignment onto the input.
    fn final_orders(top: &Topology, a: &BondAssignment) -> Vec<BondOrder> {
        let mut t = top.clone();
        a.apply_to(&mut t).unwrap();
        t.bonds.iter().map(|b| b.order()).collect()
    }

    fn final_charges(top: &Topology, a: &BondAssignment) -> Vec<i32> {
        let mut t = top.clone();
        a.apply_to(&mut t).unwrap();
        t.atoms.iter().map(|x| x.get_formal_charge().unwrap_or(0)).collect()
    }

    use BondOrder::{Double as D, Single as S, Triple as T, Unspecified as U};

    #[test]
    fn methane_all_single() {
        let top = topo(&[6, 1, 1, 1, 1], &[(0, 1, U), (0, 2, U), (0, 3, U), (0, 4, U)]);
        let a = solve(&top, &BondOrderOptions::default());
        assert_eq!(final_orders(&top, &a), vec![S, S, S, S]);
        assert_eq!(final_charges(&top, &a), vec![0, 0, 0, 0, 0]);
        assert!(a.warnings().is_empty());
    }

    #[test]
    fn carbon_dioxide_two_double_bonds() {
        // O=C=O
        let top = topo(&[8, 6, 8], &[(1, 0, U), (1, 2, U)]);
        let a = solve(&top, &BondOrderOptions::default());
        assert_eq!(final_orders(&top, &a), vec![D, D]);
        assert_eq!(final_charges(&top, &a), vec![0, 0, 0]);
    }

    #[test]
    fn ethene_and_acetylene() {
        // H2C=CH2
        let ethene = topo(
            &[6, 6, 1, 1, 1, 1],
            &[(0, 1, U), (0, 2, U), (0, 3, U), (1, 4, U), (1, 5, U)],
        );
        let a = solve(&ethene, &BondOrderOptions::default());
        assert_eq!(final_orders(&ethene, &a)[0], D);

        // HC#CH
        let acetylene = topo(&[6, 6, 1, 1], &[(0, 1, U), (0, 2, U), (1, 3, U)]);
        let a = solve(&acetylene, &BondOrderOptions::default());
        assert_eq!(final_orders(&acetylene, &a)[0], T);
    }

    #[test]
    fn formate_anion_is_a_carboxylate() {
        // H-C(=O)-O(-), net charge -1, two symmetric oxygens -> ambiguous
        let top = topo(&[6, 1, 8, 8], &[(0, 1, U), (0, 2, U), (0, 3, U)]);
        let opts = BondOrderOptions {
            total_charge: Some(-1),
            ..BondOrderOptions::default()
        };
        let a = solve(&top, &opts);
        let orders = final_orders(&top, &a);
        // one C=O, one C-O(-)
        assert_eq!(orders[0], S); // C-H
        let co = [orders[1], orders[2]];
        assert!(co.contains(&S) && co.contains(&D), "one single, one double: {co:?}");
        assert_eq!(final_charges(&top, &a).iter().sum::<i32>(), -1);
        // The carboxyl template pins a canonical form, so the symmetric choice is not flagged.
        assert!(a.warnings().is_empty());
    }

    #[test]
    fn formate_without_charge_constraint_warns() {
        let top = topo(&[6, 1, 8, 8], &[(0, 1, U), (0, 2, U), (0, 3, U)]);
        let a = solve(&top, &BondOrderOptions::default());
        assert_eq!(final_charges(&top, &a).iter().sum::<i32>(), -1);
        assert!(a.warnings().contains(&PerceptionWarning::ChargeWasNotConstrained));
    }

    #[test]
    fn formamide_is_neutral_amide() {
        // H-C(=O)-NH2
        let top = topo(
            &[6, 8, 7, 1, 1, 1],
            &[(0, 1, U), (0, 2, U), (0, 3, U), (2, 4, U), (2, 5, U)],
        );
        let a = solve(&top, &BondOrderOptions::default());
        let orders = final_orders(&top, &a);
        assert_eq!(orders[0], D, "C=O"); // C-O
        assert_eq!(orders[1], S, "C-H");
        assert_eq!(orders[2], S, "C-N");
        assert_eq!(final_charges(&top, &a), vec![0, 0, 0, 0, 0, 0]);
        assert!(a.warnings().is_empty());
    }

    #[test]
    fn nitromethane_is_net_neutral() {
        // H3C-N(=O)-O(-), N is +1, one O is -1 -> net 0
        let top = topo(
            &[6, 1, 1, 1, 7, 8, 8],
            &[(0, 1, U), (0, 2, U), (0, 3, U), (0, 4, U), (4, 5, U), (4, 6, U)],
        );
        let a = solve(&top, &BondOrderOptions::default());
        let charges = final_charges(&top, &a);
        assert_eq!(charges[4], 1, "N is +1");
        assert_eq!(charges.iter().sum::<i32>(), 0);
        // one N=O, one N-O
        let no = [final_orders(&top, &a)[4], final_orders(&top, &a)[5]];
        assert!(no.contains(&S) && no.contains(&D));
        // The nitro template pins a canonical form: deterministic, net-neutral, no warnings.
        assert!(a.warnings().is_empty());
    }

    #[test]
    fn nitrobenzene_aromatic_ring_plus_nitro() {
        // Ring C0..C5 (C0 bears the nitro, C1..C5 bear H), N=11, O=12/13. All orders unknown.
        let z = [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 7, 8, 8];
        let bonds = [
            (0, 1, U), (1, 2, U), (2, 3, U), (3, 4, U), (4, 5, U), (5, 0, U),
            (1, 6, U), (2, 7, U), (3, 8, U), (4, 9, U), (5, 10, U),
            (0, 11, U), (11, 12, U), (11, 13, U),
        ];
        let top = topo(&z, &bonds);
        let opts = BondOrderOptions {
            total_charge: Some(0),
            ..BondOrderOptions::default()
        };
        let a = solve(&top, &opts);
        let charges = final_charges(&top, &a);
        assert_eq!(charges[11], 1, "nitro N is +1");
        assert_eq!(charges[12] + charges[13], -1, "exactly one nitro O is -1");
        assert_eq!(charges.iter().sum::<i32>(), 0, "molecule is net neutral");
        // The ring alternates to a valid Kekulé structure: three doubles among the six ring bonds.
        let orders = final_orders(&top, &a);
        let ring_doubles = (0..6).filter(|&b| orders[b] == D).count();
        assert_eq!(ring_doubles, 3, "benzene ring Kekulizes to three double bonds");
    }

    #[test]
    fn preserve_known_keeps_a_given_double_bond() {
        // Give the C=O as concrete; only the C-N (unspecified) is solved.
        let top = topo(
            &[6, 8, 7, 1, 1, 1],
            &[(0, 1, D), (0, 2, U), (0, 3, U), (2, 4, U), (2, 5, U)],
        );
        let opts = BondOrderOptions {
            input_orders: InputOrders::PreserveKnown,
            ..BondOrderOptions::default()
        };
        let a = solve(&top, &opts);
        // The preserved C=O bond is unchanged -> None.
        assert_eq!(a.bond_orders()[0], None);
        assert_eq!(final_orders(&top, &a)[0], D);
        assert_eq!(final_orders(&top, &a)[2], S); // C-N
    }

    #[test]
    fn aromatic_input_is_kekulized_then_solved() {
        // Benzene ring given as aromatic bonds, with explicit ring H.
        let mut bonds = Vec::new();
        for k in 0..6 {
            bonds.push((k, (k + 1) % 6, BondOrder::Aromatic));
        }
        for k in 0..6 {
            bonds.push((k, 6 + k, U)); // ring C - H
        }
        let z = [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1];
        let top = topo(&z, &bonds);
        let a = solve(&top, &BondOrderOptions::default());
        let ring: Vec<BondOrder> = final_orders(&top, &a)[..6].to_vec();
        assert_eq!(ring.iter().filter(|&&o| o == D).count(), 3, "three double bonds");
        assert_eq!(ring.iter().filter(|&&o| o == S).count(), 3, "three single bonds");
        assert!(ring.iter().all(|&o| o != BondOrder::Aromatic), "no aromatic bonds left");
        assert_eq!(final_charges(&top, &a), vec![0; 12]);
    }

    #[test]
    fn ammonium_is_a_cation() {
        // NH4+ : four N-H bonds, N is +1
        let top = topo(&[7, 1, 1, 1, 1], &[(0, 1, U), (0, 2, U), (0, 3, U), (0, 4, U)]);
        let opts = BondOrderOptions {
            total_charge: Some(1),
            ..BondOrderOptions::default()
        };
        let a = solve(&top, &opts);
        assert_eq!(final_charges(&top, &a)[0], 1);
        assert_eq!(final_orders(&top, &a), vec![S, S, S, S]);
    }

    #[test]
    fn total_charge_with_two_fragments_is_rejected() {
        // Two separate H2 molecules.
        let top = topo(&[1, 1, 1, 1], &[(0, 1, U), (2, 3, U)]);
        let opts = BondOrderOptions {
            total_charge: Some(0),
            ..BondOrderOptions::default()
        };
        assert!(matches!(
            assign_bond_orders(&top, None, &opts),
            Err(BondPerceptionError::TotalChargeWithMultipleComponents)
        ));
    }

    #[test]
    fn impossible_valence_fails_without_a_solution() {
        // A single naked carbon with one hydrogen cannot reach any valence state
        // (valence-complete assumption): C would need bosum 4 but has only one bond.
        let top = topo(&[6, 1], &[(0, 1, U)]);
        assert!(matches!(
            assign_bond_orders(&top, None, &BondOrderOptions::default()),
            Err(BondPerceptionError::NoValidAssignment { atom: 0 })
        ));
    }

    // --- Geometry-informed implicit-hydrogen mode --------------------------------------------

    /// Solve a hydrogen-free structure from connectivity + coordinates (nm) in
    /// `InferFromGeometry` mode. Every bond is unspecified; positions decide the orders.
    fn solve_geom(z: &[u8], bonds: &[(usize, usize)], coords: &[[Float; 3]]) -> BondAssignment {
        let bond_rows: Vec<(usize, usize, BondOrder)> =
            bonds.iter().map(|&(i, j)| (i, j, U)).collect();
        let top = topo(z, &bond_rows);
        let pos: Vec<Pos> = coords.iter().map(|c| Pos::new(c[0], c[1], c[2])).collect();
        let opts = BondOrderOptions {
            hydrogens: HydrogenPolicy::InferFromGeometry,
            ..BondOrderOptions::default()
        };
        assign_bond_orders(&top, Some(&pos), &opts).unwrap()
    }

    /// Like [`solve_geom`], with a total-charge target.
    fn solve_geom_charged(
        z: &[u8],
        bonds: &[(usize, usize)],
        coords: &[[Float; 3]],
        total_charge: i32,
    ) -> Result<BondAssignment, BondPerceptionError> {
        let bond_rows: Vec<(usize, usize, BondOrder)> =
            bonds.iter().map(|&(i, j)| (i, j, U)).collect();
        let top = topo(z, &bond_rows);
        let pos: Vec<Pos> = coords.iter().map(|c| Pos::new(c[0], c[1], c[2])).collect();
        let opts = BondOrderOptions {
            hydrogens: HydrogenPolicy::InferFromGeometry,
            total_charge: Some(total_charge),
            ..BondOrderOptions::default()
        };
        assign_bond_orders(&top, Some(&pos), &opts)
    }

    /// Heavy-atom ethylamine with a +1 target is ethylammonium. Every bond is single by length,
    /// so no bond is free: the charge must come from choosing the nitrogen's state, which the
    /// solver used to fix as neutral before the search and then report no valid assignment.
    #[test]
    fn geometry_with_charge_target_protonates_a_fixed_amine() {
        let coords = [[0.0, 0.0, 0.0], [0.153, 0.0, 0.0], [0.200, 0.140, 0.0]];
        let a = solve_geom_charged(&[6, 6, 7], &[(0, 1), (1, 2)], &coords, 1).unwrap();
        assert_eq!(a.implicit_hydrogens(), &[3, 2, 3], "CH3-CH2-NH3+");
        assert_eq!(a.formal_charges()[2], Some(1));
    }

    /// The same with a free cluster present: benzylamine + 1. The ring is solved on its own and
    /// the +1 still lands on the amine, not in the ring.
    #[test]
    fn geometry_with_charge_target_combines_ring_and_fixed_atoms() {
        let r = 0.139;
        let mut coords: Vec<[Float; 3]> = (0..6)
            .map(|k| {
                let t = std::f64::consts::PI / 3.0 * k as f64;
                [r * t.cos() as Float, r * t.sin() as Float, 0.0]
            })
            .collect();
        coords.push([0.290, 0.0, 0.0]); // CH2, 0.151 nm from ring C0
        coords.push([0.340, 0.138, 0.0]); // N, 0.147 nm from CH2
        let mut bonds: Vec<(usize, usize)> = (0..6).map(|k| (k, (k + 1) % 6)).collect();
        bonds.extend([(0, 6), (6, 7)]);
        let a = solve_geom_charged(&[6, 6, 6, 6, 6, 6, 6, 7], &bonds, &coords, 1).unwrap();
        assert_eq!(
            a.implicit_hydrogens(),
            &[0, 1, 1, 1, 1, 1, 2, 3],
            "C6H5-CH2-NH3+"
        );
        assert_eq!(a.formal_charges()[7], Some(1));
    }

    /// Heavy-atom formamide. The 0.135 nm amide C-N is short enough to pass the double-bond
    /// cap, and HN=C=O needs fewer implicit hydrogens, but its carbon would be sp while the
    /// N-C-O angle is 124°: the linear-geometry rule keeps H2N-CH=O.
    #[test]
    fn geometry_rejects_two_pi_bonds_on_a_bent_atom() {
        let t = (124.0 as Float).to_radians();
        let coords = [
            [0.135, 0.0, 0.0],                       // N
            [0.0, 0.0, 0.0],                         // C
            [0.123 * t.cos(), 0.123 * t.sin(), 0.0], // O
        ];
        let a = solve_geom(&[7, 6, 8], &[(0, 1), (1, 2)], &coords);
        assert_eq!(a.implicit_hydrogens(), &[2, 1, 0], "H2N-CH=O, not HN=C=O");
        assert_eq!(a.bond_orders()[0], Some(S));
        assert_eq!(a.bond_orders()[1], Some(D));
    }

    /// ...while a linear atom keeps its two π bonds: acetonitrile's C≡N stays triple.
    #[test]
    fn geometry_keeps_a_linear_nitrile() {
        let coords = [[-0.146, 0.0, 0.0], [0.0, 0.0, 0.0], [0.116, 0.0, 0.0]];
        let a = solve_geom(&[6, 6, 7], &[(0, 1), (1, 2)], &coords);
        assert_eq!(a.bond_orders()[1], Some(T));
        assert_eq!(a.implicit_hydrogens(), &[3, 0, 0], "CH3-C≡N");
    }

    #[test]
    fn geometry_infers_order_and_implicit_h_for_two_carbons() {
        // Bond length alone distinguishes ethane / ethene / ethyne (heavy atoms only).
        for (len, order, h) in [(0.154, S, 3u8), (0.134, D, 2), (0.120, T, 1)] {
            let a = solve_geom(&[6, 6], &[(0, 1)], &[[0.0, 0.0, 0.0], [len, 0.0, 0.0]]);
            assert_eq!(a.bond_orders()[0], Some(order), "len {len}");
            assert_eq!(a.implicit_hydrogens(), &[h, h], "len {len}");
        }
    }

    #[test]
    fn geometry_infers_carbonyl_vs_hydroxyl_oxygen() {
        // Formaldehyde: short C=O -> O takes no hydrogen, C takes two.
        let carbonyl = solve_geom(&[6, 8], &[(0, 1)], &[[0.0, 0.0, 0.0], [0.121, 0.0, 0.0]]);
        assert_eq!(carbonyl.bond_orders()[0], Some(D));
        assert_eq!(carbonyl.implicit_hydrogens(), &[2, 0]);

        // Methanol: long C-O -> single bond, O takes one hydrogen, C takes three.
        let hydroxyl = solve_geom(&[6, 8], &[(0, 1)], &[[0.0, 0.0, 0.0], [0.143, 0.0, 0.0]]);
        assert_eq!(hydroxyl.bond_orders()[0], Some(S));
        assert_eq!(hydroxyl.implicit_hydrogens(), &[3, 1]);
    }

    #[test]
    fn geometry_infers_benzene_not_cumulene_or_cyclohexane() {
        // Six carbons on a regular hexagon with aromatic ~0.139 nm sides, no hydrogens.
        let r = 0.139; // regular hexagon: side length equals circumradius
        let coords: Vec<[Float; 3]> = (0..6)
            .map(|k| {
                let t = std::f64::consts::PI / 3.0 * k as f64;
                [r * t.cos() as Float, r * t.sin() as Float, 0.0]
            })
            .collect();
        let bonds: Vec<(usize, usize)> = (0..6).map(|k| (k, (k + 1) % 6)).collect();
        let a = solve_geom(&[6; 6], &bonds, &coords);

        let orders: Vec<BondOrder> = a.bond_orders().iter().map(|o| o.unwrap()).collect();
        assert_eq!(orders.iter().filter(|&&o| o == D).count(), 3, "three double bonds");
        assert_eq!(orders.iter().filter(|&&o| o == S).count(), 3, "three single bonds");
        // Every ring carbon takes exactly one implicit hydrogen (benzene C6H6), which rules out
        // both cumulene (0 H) and cyclohexane (2 H each).
        assert_eq!(a.implicit_hydrogens(), &[1, 1, 1, 1, 1, 1]);
        // The two Kekulé structures are equally good, so the assignment is reported ambiguous.
        assert!(a.warnings().contains(&PerceptionWarning::AmbiguousAssignment));
    }

    #[test]
    fn infer_mode_requires_coordinates() {
        let top = topo(&[6, 6], &[(0, 1, U)]);
        let opts = BondOrderOptions {
            hydrogens: HydrogenPolicy::InferFromGeometry,
            ..BondOrderOptions::default()
        };
        assert!(matches!(
            assign_bond_orders(&top, None, &opts),
            Err(BondPerceptionError::MissingCoordinates)
        ));
    }
}
