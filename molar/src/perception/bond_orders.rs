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

    // Canonical bond orders for recognized functional groups; these take priority over the
    // input orders and the search.
    let template = if options.use_functional_groups {
        super::functional_groups::functional_group_orders(&z, &adj)
    } else {
        vec![None; m]
    };

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

    // Ring membership per bond, for the cumulene guard (only needed when inferring hydrogens).
    let bond_in_ring: Vec<bool> = if infer {
        let mut r = vec![false; m];
        for ring in super::sssr(&adj) {
            for b in ring.bonds {
                r[b] = true;
            }
        }
        r
    } else {
        vec![false; m]
    };

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

    for c in 0..n_comp {
        if comp_bonds[c].is_empty() {
            continue; // an isolated atom: leave its order/charge unchanged
        }
        let atoms = &comp_atoms[c];

        let mut vars: Vec<usize> = Vec::new();
        for &b in &comp_bonds[c] {
            let [i, j] = pairs[b];
            match fixed[b] {
                Some(o) => {
                    base[i] += o as i32;
                    base[j] += o as i32;
                }
                None => {
                    let dm = dom_max[b] as i32;
                    vars.push(b);
                    rem_min[i] += 1;
                    rem_max[i] += dm;
                    rem_min[j] += 1;
                    rem_max[j] += dm;
                }
            }
        }
        // Order the free bonds by a connected walk so that each atom's incident bonds are
        // decided together. An atom whose bonds are all assigned has a fixed valence, which the
        // feasibility bound then checks immediately — this keeps the search near-linear on
        // sparse molecules instead of scattering decisions and pruning only at the leaves.
        vars = connected_order(atoms, &vars, &pairs);

        let target = if bonded_comp_count == 1 {
            options.total_charge
        } else {
            None
        };

        let mut search = Search {
            pairs: &pairs,
            z: &z,
            dom_max: &dom_max,
            bond_in_ring: &bond_in_ring,
            infer,
            vars: &vars,
            atoms,
            base: &base,
            rem_min: &mut rem_min,
            rem_max: &mut rem_max,
            assigned_sum: &mut assigned_sum,
            order_out: vec![0u8; vars.len()],
            total_charge: target,
            max_branches: options.limits.max_branches,
            branches: 0,
            best_cost: None,
            best_orders: Vec::new(),
            best_fc: Vec::new(),
            best_ih: Vec::new(),
            ambiguous: false,
            truncated: false,
        };
        search.run(0);

        let best_cost = search.best_cost;
        let best_orders = search.best_orders.clone();
        let best_fc = search.best_fc.clone();
        let best_ih = search.best_ih.clone();
        let ambiguous = search.ambiguous;
        let truncated = search.truncated;
        drop(search);

        // Restore the scratch this fragment touched (the search left the counters at their
        // pre-search values; zero them for the next fragment).
        for &a in atoms {
            base[a] = 0;
            rem_min[a] = 0;
            rem_max[a] = 0;
            assigned_sum[a] = 0;
        }

        let Some(_) = best_cost else {
            let atom = atoms[0];
            return Err(if truncated {
                BondPerceptionError::SearchLimitExceeded { atom }
            } else {
                BondPerceptionError::NoValidAssignment { atom }
            });
        };

        // Write out every bond whose solved order differs from the input (`None` = unchanged).
        for (pos, &b) in vars.iter().enumerate() {
            let solved = int_to_order(best_orders[pos]);
            if solved != input_orders[b] {
                bond_orders[b] = Some(solved);
            }
        }
        for &b in &comp_bonds[c] {
            if let Some(o) = fixed[b] {
                let solved = int_to_order(o);
                if solved != input_orders[b] {
                    bond_orders[b] = Some(solved);
                }
            }
        }

        let mut comp_charge = 0i32;
        for &(a, fc) in &best_fc {
            comp_charge += fc;
            if fc != input_fc[a] {
                formal_charges[a] = Some(fc);
            }
        }
        for &(a, h) in &best_ih {
            implicit_hydrogens[a] = h;
        }

        if ambiguous {
            warnings.push(PerceptionWarning::AmbiguousAssignment);
        }
        if target.is_none() && comp_charge != 0 {
            warnings.push(PerceptionWarning::ChargeWasNotConstrained);
        }
        if truncated {
            warnings.push(PerceptionWarning::SearchTruncated);
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

/// Penalty for an atom carrying more than one ring double bond (a cumulated aromatic ring, as
/// opposed to a Kekulé structure with one double per atom). Large so it dominates the state
/// penalties, but applied only when inferring hydrogens (see [`Search::record_leaf`]).
const RING_CUMULENE_PENALTY: u32 = 100;

/// One connected fragment's branch-and-bound over the free bond orders.
struct Search<'a> {
    pairs: &'a [[usize; 2]],
    z: &'a [u8],
    dom_max: &'a [u8],
    /// Whether each bond lies on a ring (used only when inferring hydrogens).
    bond_in_ring: &'a [bool],
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
    /// Best objective found: `(summed penalty, net-charge magnitude, total implicit H)`. The
    /// penalty (state penalties plus any ring-cumulene penalty) is primary; the charge
    /// magnitude breaks ties toward the least charge-separated form; the implicit-hydrogen
    /// count breaks the remaining ties toward the most-saturated structure the geometry allows.
    best_cost: Option<(u32, u32, u32)>,
    best_orders: Vec<u8>,
    best_fc: Vec<(usize, i32)>,
    best_ih: Vec<(usize, u8)>,
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
        // achievable net-charge range over the fragment.
        let Some((lb, charge_lo, charge_hi)) = self.bound() else {
            return;
        };
        // With a fixed total charge, prune as soon as the atoms still open cannot bring the net
        // charge to the target. This is the key cut on large fragments: without it the search
        // descends to a full leaf before the leaf-level charge check rejects it.
        if let Some(target) = self.total_charge
            && (target < charge_lo || target > charge_hi)
        {
            return;
        }

        if let Some((best_penalty, best_dev, _)) = self.best_cost {
            // The penalty lower bound alone rules this branch out.
            if lb > best_penalty {
                return;
            }
            // Equal-penalty branches are kept only while they might still improve the result: to
            // find a strictly better secondary cost, or the first distinct equal-cost solution
            // that proves the assignment ambiguous. Once the charge cost is at its floor and
            // ambiguity is known, pruning them stops the tie explosion on big symmetric and
            // poly-aromatic fragments. This tie-prune is skipped when inferring hydrogens, where
            // the implicit-hydrogen tie-break can still improve an equal-penalty, floor-charge
            // branch; those fragments are small enough not to need it.
            if !self.infer && lb == best_penalty && best_dev == 0 && self.ambiguous {
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
    /// (an admissible lower bound), together with the range of net formal charge the still-open
    /// atoms allow. `None` if some atom already has no reachable valence state.
    fn bound(&self) -> Option<(u32, i32, i32)> {
        let mut lb = 0u32;
        let mut charge_lo = 0i32;
        let mut charge_hi = 0i32;
        for &a in self.atoms {
            let Some(states) = valence_states(self.z[a]) else {
                continue; // unmodeled element: no constraint, no penalty, no charge
            };
            let lo = self.base[a] + self.assigned_sum[a] + self.rem_min[a];
            let hi = self.base[a] + self.assigned_sum[a] + self.rem_max[a];
            let mut best_pen: Option<u32> = None;
            let mut min_fc = i32::MAX;
            let mut max_fc = i32::MIN;
            for st in states {
                if self.state_fits(st.valence as i32, lo, hi) {
                    best_pen = Some(best_pen.map_or(st.penalty, |p| p.min(st.penalty)));
                    min_fc = min_fc.min(st.formal_charge);
                    max_fc = max_fc.max(st.formal_charge);
                }
            }
            lb += best_pen?;
            charge_lo += min_fc;
            charge_hi += max_fc;
        }
        Some((lb, charge_lo, charge_hi))
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
        let target = self.total_charge.unwrap_or(0);

        let mut penalty = 0u32;
        let mut net = 0i32;
        let mut fcs: Vec<(usize, i32)> = Vec::with_capacity(self.atoms.len());
        let mut ihs: Vec<(usize, u8)> = Vec::with_capacity(self.atoms.len());
        let mut charge_choice_tie = false;

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
            fcs.push((a, fc));
            ihs.push((a, (valence - bosum) as u8));
        }

        let charge_dev = match self.total_charge {
            Some(t) => {
                if net != t {
                    return;
                }
                0
            }
            None => net.unsigned_abs(),
        };
        let implicit_total: u32 = ihs.iter().map(|&(_, h)| h as u32).sum();
        let cost = (penalty, charge_dev, implicit_total);

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
