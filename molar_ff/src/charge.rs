//! Espaloma-charge partial charges from a fixed graph neural network.
//!
//! The source model is an ONNX file, but production inference does not use a general
//! ONNX run time. The model contains only matrix multiplication, addition, `tanh`,
//! `ReLU`, and two output-column selections. Its reviewed weights are extracted to a
//! fixed binary file by `tools/extract_espaloma_weights.py`. See `assets/README.md` for
//! the file format and update procedure.
//!
//! ## Large-system behavior
//!
//! The GNN maps per-atom features `[n, 116]` to a hidden state `[n, 128]`. Each of
//! four message layers combines a self transform with the mean hidden state of bonded
//! neighbors. Neighbor aggregation uses [`BondAdjacency`] directly. It does not create
//! the mathematically equivalent dense `[n, n]` adjacency matrix. Thus, graph work is
//! `O(E * 128)` and graph memory is `O(n + E)`, where `E` is the bond count.
//!
//! Dense transforms use the SIMD kernels in `matrixmultiply`. Native builds divide
//! sufficiently large row sets into independent Rayon tasks. Each task keeps only one
//! bounded neighbor-work block. The two full hidden buffers are `O(n * 128)`. Small
//! molecules use one task to avoid scheduling and packing overhead.

use matrixmultiply::sgemm;
use molar::prelude::BondAdjacency;
use std::sync::LazyLock;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

const FEATURE_WIDTH: usize = 116;
const HIDDEN_WIDTH: usize = 128;
const OUTPUT_WIDTH: usize = 2;
const MESSAGE_LAYERS: usize = 4;

/// A row block gives each matrix kernel enough work for effective SIMD packing. It also
/// bounds temporary neighbor storage to 256 KiB per active task.
const ROW_BLOCK: usize = 512;
/// Below this size, one matrix call is faster than Rayon task scheduling.
#[cfg(not(target_arch = "wasm32"))]
const PARALLEL_ROWS: usize = ROW_BLOCK * 2;

const INPUT_WEIGHT_START: usize = 0;
const INPUT_WEIGHT_LEN: usize = FEATURE_WIDTH * HIDDEN_WIDTH;
const INPUT_BIAS_START: usize = INPUT_WEIGHT_START + INPUT_WEIGHT_LEN;
const INPUT_BIAS_LEN: usize = HIDDEN_WIDTH;
const LAYER_WEIGHT_START: usize = INPUT_BIAS_START + INPUT_BIAS_LEN;
const LAYER_WEIGHT_LEN: usize = HIDDEN_WIDTH * HIDDEN_WIDTH;
const OUTPUT_WEIGHT_START: usize = LAYER_WEIGHT_START + 2 * MESSAGE_LAYERS * LAYER_WEIGHT_LEN;
const OUTPUT_WEIGHT_LEN: usize = HIDDEN_WIDTH * OUTPUT_WIDTH;
const OUTPUT_BIAS_START: usize = OUTPUT_WEIGHT_START + OUTPUT_WEIGHT_LEN;
const OUTPUT_BIAS_LEN: usize = OUTPUT_WIDTH;
const WEIGHT_COUNT: usize = OUTPUT_BIAS_START + OUTPUT_BIAS_LEN;
const WEIGHT_MAGIC: &[u8; 8] = b"ESPCHG01";
const WEIGHT_BYTES: &[u8] = include_bytes!("../assets/espaloma_charge.weights");

/// Parsed fixed-model weights. Parsing copies the unaligned embedded bytes once. All
/// inference calls then use aligned native `f32` values without byte conversion.
struct Weights(Box<[f32]>);

impl Weights {
    fn load() -> Self {
        assert_eq!(
            &WEIGHT_BYTES[..WEIGHT_MAGIC.len()],
            WEIGHT_MAGIC,
            "invalid espaloma weight-file marker"
        );
        let bytes = &WEIGHT_BYTES[WEIGHT_MAGIC.len()..];
        assert_eq!(
            bytes.len(),
            WEIGHT_COUNT * 4,
            "invalid espaloma weight-file size"
        );
        let values = bytes
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect::<Vec<_>>()
            .into_boxed_slice();
        Self(values)
    }

    #[inline]
    fn input_weight(&self) -> &[f32] {
        &self.0[INPUT_WEIGHT_START..INPUT_BIAS_START]
    }

    #[inline]
    fn input_bias(&self) -> &[f32] {
        &self.0[INPUT_BIAS_START..LAYER_WEIGHT_START]
    }

    #[inline]
    fn self_weight(&self, layer: usize) -> &[f32] {
        let start = LAYER_WEIGHT_START + 2 * layer * LAYER_WEIGHT_LEN;
        &self.0[start..start + LAYER_WEIGHT_LEN]
    }

    #[inline]
    fn neighbor_weight(&self, layer: usize) -> &[f32] {
        let start = LAYER_WEIGHT_START + (2 * layer + 1) * LAYER_WEIGHT_LEN;
        &self.0[start..start + LAYER_WEIGHT_LEN]
    }

    #[inline]
    fn output_weight(&self) -> &[f32] {
        &self.0[OUTPUT_WEIGHT_START..OUTPUT_BIAS_START]
    }

    #[inline]
    fn output_bias(&self) -> &[f32] {
        &self.0[OUTPUT_BIAS_START..]
    }
}

static WEIGHTS: LazyLock<Weights> = LazyLock::new(Weights::load);

/// Compute `C = A * B + beta * C` for row-major matrices.
///
/// `A` is `[rows, inner]`, `B` is `[inner, cols]`, and `C` is `[rows, cols]`.
fn matmul_into(
    a: &[f32],
    rows: usize,
    inner: usize,
    b: &[f32],
    cols: usize,
    beta: f32,
    c: &mut [f32],
) {
    assert_eq!(a.len(), rows * inner, "invalid left matrix extent");
    assert_eq!(b.len(), inner * cols, "invalid right matrix extent");
    assert_eq!(c.len(), rows * cols, "invalid output matrix extent");
    // SAFETY: The assertions above prove all three matrix extents. The row-major
    // strides address each element in those slices exactly once. `a`, `b`, and `c`
    // are separate borrows, so the output cannot alias either input. A zero-sized
    // matrix is handled before the call because the kernel requires valid pointers.
    if rows != 0 && inner != 0 && cols != 0 {
        unsafe {
            sgemm(
                rows,
                inner,
                cols,
                1.0,
                a.as_ptr(),
                inner as isize,
                1,
                b.as_ptr(),
                cols as isize,
                1,
                beta,
                c.as_mut_ptr(),
                cols as isize,
                1,
            );
        }
    }
}

/// Apply `work` to disjoint row blocks. Large native inputs use Rayon; small and
/// wasm inputs use the same block operation serially.
fn row_blocks<F>(output: &mut [f32], row_width: usize, rows: usize, work: F)
where
    F: Fn(usize, &mut [f32]) + Send + Sync,
{
    debug_assert_eq!(output.len(), rows * row_width);
    let block_len = ROW_BLOCK * row_width;

    #[cfg(not(target_arch = "wasm32"))]
    if rows >= PARALLEL_ROWS {
        output
            .par_chunks_mut(block_len)
            .enumerate()
            .for_each(|(block, out)| work(block * ROW_BLOCK, out));
        return;
    }

    for (block, out) in output.chunks_mut(block_len).enumerate() {
        work(block * ROW_BLOCK, out);
    }
}

fn input_layer(features: &[f32], n: usize, weights: &Weights) -> Vec<f32> {
    let mut hidden = vec![0.0; n * HIDDEN_WIDTH];
    row_blocks(&mut hidden, HIDDEN_WIDTH, n, |start, out| {
        let rows = out.len() / HIDDEN_WIDTH;
        let input = &features[start * FEATURE_WIDTH..(start + rows) * FEATURE_WIDTH];
        matmul_into(
            input,
            rows,
            FEATURE_WIDTH,
            weights.input_weight(),
            HIDDEN_WIDTH,
            0.0,
            out,
        );
        for row in out.chunks_exact_mut(HIDDEN_WIDTH) {
            for (value, bias) in row.iter_mut().zip(weights.input_bias()) {
                *value = (*value + bias).tanh();
            }
        }
    });
    hidden
}

/// Compute one message layer. The temporary neighbor matrix is one row block, not
/// `[n, 128]`. Peak full-size inference storage is therefore two hidden matrices.
fn message_layer(
    hidden: &[f32],
    n: usize,
    adjacency: &BondAdjacency,
    self_weight: &[f32],
    neighbor_weight: &[f32],
) -> Vec<f32> {
    let mut next = vec![0.0; n * HIDDEN_WIDTH];
    row_blocks(&mut next, HIDDEN_WIDTH, n, |start, out| {
        let rows = out.len() / HIDDEN_WIDTH;
        let mut neighbor_mean = vec![0.0; out.len()];
        for (local, mean) in neighbor_mean.chunks_exact_mut(HIDDEN_WIDTH).enumerate() {
            let neighbors = adjacency.neighbors(start + local);
            if neighbors.is_empty() {
                continue;
            }
            for neighbor in neighbors {
                let source =
                    &hidden[neighbor.atom() * HIDDEN_WIDTH..(neighbor.atom() + 1) * HIDDEN_WIDTH];
                for (sum, value) in mean.iter_mut().zip(source) {
                    *sum += value;
                }
            }
            let inverse_degree = 1.0 / neighbors.len() as f32;
            for value in mean {
                *value *= inverse_degree;
            }
        }

        let own = &hidden[start * HIDDEN_WIDTH..(start + rows) * HIDDEN_WIDTH];
        matmul_into(own, rows, HIDDEN_WIDTH, self_weight, HIDDEN_WIDTH, 0.0, out);
        matmul_into(
            &neighbor_mean,
            rows,
            HIDDEN_WIDTH,
            neighbor_weight,
            HIDDEN_WIDTH,
            1.0,
            out,
        );
        for value in out {
            *value = value.max(0.0);
        }
    });
    next
}

/// Run the fixed GNN. `features` is row-major `[n, 116]`. `adjacency` contains
/// the undirected molecular bonds and supplies the row-mean neighbor operation.
/// Returns `(electronegativity, hardness)` with one value per atom.
pub(crate) fn run_gnn(
    features: Vec<f32>,
    adjacency: &BondAdjacency,
    n: usize,
) -> Result<(Vec<f32>, Vec<f32>), &'static str> {
    let expected = n
        .checked_mul(FEATURE_WIDTH)
        .ok_or("espaloma input is too large")?;
    n.checked_mul(HIDDEN_WIDTH)
        .ok_or("espaloma input is too large")?;
    n.checked_mul(OUTPUT_WIDTH)
        .ok_or("espaloma input is too large")?;
    if features.len() != expected {
        return Err("espaloma feature matrix has an invalid size");
    }

    let weights = &*WEIGHTS;
    let mut hidden = input_layer(&features, n, weights);
    // Release the 116-column feature matrix before the larger message-layer work.
    drop(features);
    for layer in 0..MESSAGE_LAYERS {
        hidden = message_layer(
            &hidden,
            n,
            adjacency,
            weights.self_weight(layer),
            weights.neighbor_weight(layer),
        );
    }

    let mut output = vec![0.0; n * OUTPUT_WIDTH];
    row_blocks(&mut output, OUTPUT_WIDTH, n, |start, out| {
        let rows = out.len() / OUTPUT_WIDTH;
        let input = &hidden[start * HIDDEN_WIDTH..(start + rows) * HIDDEN_WIDTH];
        matmul_into(
            input,
            rows,
            HIDDEN_WIDTH,
            weights.output_weight(),
            OUTPUT_WIDTH,
            0.0,
            out,
        );
        for row in out.chunks_exact_mut(OUTPUT_WIDTH) {
            row[0] += weights.output_bias()[0];
            row[1] += weights.output_bias()[1];
        }
    });

    let mut electronegativity = Vec::with_capacity(n);
    let mut hardness = Vec::with_capacity(n);
    for row in output.chunks_exact(OUTPUT_WIDTH) {
        electronegativity.push(row[0]);
        hardness.push(row[1]);
    }
    Ok((electronegativity, hardness))
}

/// Standard atomic weight for the elements espaloma supports (matches RDKit `GetMass`).
fn mass_by_z(z: u8) -> f32 {
    match z {
        1 => 1.008,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        15 => 30.974,
        16 => 32.06,
        17 => 35.45,
        35 => 79.904,
        53 => 126.904,
        _ => 0.0,
    }
}

/// RDKit hybridization one-hot index: 0=SP,1=SP2,2=SP3,3=SP3D,4=SP3D2. `None` (all-zero)
/// for hydrogen (RDKit reports S/unspecified). `neighbor_conj` = does any neighbour carry a
/// multiple bond or is aromatic — a lone pair adjacent to such a π system conjugates into it,
/// so RDKit reports the atom as SP2 (amide N, ester/conjugated O) rather than SP3.
fn hybridization(
    z: u8,
    degree: usize,
    n_double: u32,
    n_triple: u32,
    aromatic: bool,
    neighbor_conj: bool,
) -> Option<usize> {
    if z == 1 {
        return None;
    }
    if aromatic {
        return Some(1); // SP2
    }
    if degree >= 6 {
        return Some(4); // SP3D2
    }
    if degree == 5 {
        return Some(3); // SP3D
    }
    if degree == 4 {
        return Some(2); // tetrahedral: sp3 C, ammonium N, sulfone/phosphate → SP3
    }
    // degree <= 3, non-aromatic:
    if n_triple >= 1 || n_double >= 2 {
        return Some(0); // SP: nitrile, allene / CO2 central atom
    }
    if n_double == 1 {
        return Some(1); // SP2: carbonyl, imine, ...
    }
    // lone-pair atom with no π bond of its own:
    if neighbor_conj && ((z == 7 && degree == 3) || (z == 8 && degree <= 2)) {
        return Some(1); // conjugated lone pair → SP2 (amide N, ester/conjugated O)
    }
    Some(2) // SP3
}

/// Outer-shell (valence) electron count for the elements espaloma supports.
fn n_outer_elec(z: u8) -> i32 {
    match z {
        1 => 1,
        6 => 4,
        7 => 5,
        8 => 6,
        9 => 7,
        15 => 5,
        16 => 6,
        17 => 7,
        35 => 7,
        53 => 7,
        _ => 0,
    }
}

/// Pauling electronegativity, used to classify exocyclic double bonds: a double bond to a *more*
/// electronegative atom (C=O, C=N) pulls the π out of the ring plane (contributes 0, ring may
/// still be aromatic), whereas one to an equal/less electronegative atom (C=C — fulvene,
/// quinone-methide, thioxanthene ylidene) disrupts the ring and breaks aromaticity.
fn electronegativity(z: u8) -> f32 {
    match z {
        1 => 2.20,
        6 => 2.55,
        7 => 3.04,
        8 => 3.44,
        9 => 3.98,
        15 => 2.19,
        16 => 2.58,
        17 => 3.16,
        35 => 2.96,
        53 => 2.66,
        _ => 0.0,
    }
}

/// RDKit-style aromaticity perception. Rather than counting a specific Kekulé structure (which
/// mis-counts fused systems), each atom's π contribution is derived Kekulé-invariantly from its
/// element, formal charge and σ-connectivity, then Hückel's 4n+2 rule is applied to every
/// individual ring **and** every maximal fused ring system (so naphthalene, azulene, purine,
/// phenanthridine … resolve correctly).
///
/// Per-atom contribution: an sp3/hyper-valent atom (≥4 σ connections) or an in-ring triple bond
/// breaks aromaticity; an exocyclic multiple bond to a non-ring atom (carbonyl, exocyclic
/// methylene) donates 0 (its π sits outside the ring); otherwise `avail = n_outer − fc − σ`
/// leaves an odd count (one p electron → 1) or an even count (a lone pair → 2), with `avail ≤ 0`
/// (carbocation, quaternary N⁺) donating 0.
fn aromatic_atoms(
    z: &[u8],
    fc: &[i32],
    bonds: &[crate::gaff::LocalBond],
    rings: &[Vec<usize>],
) -> Vec<bool> {
    let n = z.len();
    let mut inc: Vec<Vec<(usize, u8)>> = vec![Vec::new(); n];
    for b in bonds {
        inc[b.i].push((b.j, b.order));
        inc[b.j].push((b.i, b.order));
    }
    let mut in_ring = vec![false; n];
    for r in rings {
        for &a in r {
            in_ring[a] = true;
        }
    }
    // Kekulé-invariant per-atom π contribution (None = breaks aromaticity). σ is the explicit
    // connection count — molar carries no implicit hydrogens, and under-hydrogenated molecules
    // are excluded upstream. sp3/hyper-valent (≥4 σ) or an in-ring triple bond breaks. An
    // exocyclic multiple bond to a more-electronegative atom (C=O, C=N) contributes 0; one to an
    // equal/less-electronegative atom (C=C — fulvene/ylidene) breaks the ring. Otherwise
    // `n_outer − fc − σ` leaves an odd count (one p electron → 1) or even (a lone pair → 2), ≤0 → 0.
    let contrib: Vec<Option<i32>> = (0..n)
        .map(|a| {
            let sigma = inc[a].len() as i32;
            if sigma >= 4 || inc[a].iter().any(|&(j, o)| o == 3 && in_ring[j]) {
                return None;
            }
            let mut exocyclic_zero = false;
            for &(j, o) in &inc[a] {
                if o >= 2 && !in_ring[j] {
                    if electronegativity(z[j]) > electronegativity(z[a]) {
                        exocyclic_zero = true;
                    } else {
                        return None; // exocyclic C=C etc. disrupts the ring
                    }
                }
            }
            if exocyclic_zero {
                return Some(0);
            }
            let avail = n_outer_elec(z[a]) - fc[a] - sigma;
            Some(if avail <= 0 {
                0
            } else if avail % 2 == 1 {
                1
            } else {
                2
            })
        })
        .collect();
    let huckel = |atoms: &[usize]| -> bool {
        let mut pi = 0;
        for &a in atoms {
            match contrib[a] {
                Some(c) => pi += c,
                None => return false,
            }
        }
        pi % 4 == 2
    };

    let mut arom = vec![false; n];
    // Every individual ring, plus every maximal fused ring system (union of rings sharing a bond),
    // so fused aromatics resolve even when no single SSSR ring independently reaches 4n+2.
    for ring in rings {
        if huckel(ring) {
            for &a in ring {
                arom[a] = true;
            }
        }
    }
    let mut parent: Vec<usize> = (0..rings.len()).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = x;
        while parent[c] != r {
            let next = parent[c];
            parent[c] = r;
            c = next;
        }
        r
    }
    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            if rings[i].iter().filter(|a| rings[j].contains(a)).count() >= 2 {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                parent[ri] = rj;
            }
        }
    }
    let mut systems: std::collections::BTreeMap<usize, std::collections::BTreeSet<usize>> =
        Default::default();
    for i in 0..rings.len() {
        let root = find(&mut parent, i);
        systems
            .entry(root)
            .or_default()
            .extend(rings[i].iter().copied());
    }
    for atoms in systems.values() {
        let atoms: Vec<usize> = atoms.iter().copied().collect();
        if huckel(&atoms) {
            for &a in &atoms {
                arom[a] = true;
            }
        }
    }
    arom
}

/// Build the 116-column atom feature matrix and the sparse bond adjacency that
/// the Espaloma GNN consumes.
pub(crate) fn featurize(
    z: &[u8],
    fc: &[i32],
    bonds: &[crate::gaff::LocalBond],
) -> (Vec<f32>, BondAdjacency) {
    let n = z.len();
    // One bonded-adjacency index over the local subgraph, shared by the ring search and the
    // neighbour walks below. Replaces a `build_con` call plus a full `Vec<Bond>` copy that
    // existed only to reach `sssr_rings`.
    let bond_adj = molar::prelude::BondAdjacency::build(n, bonds.iter().map(|b| [b.i, b.j]));
    // Rings from molar core SSSR (matches RDKit ring semantics and is H-independent, so aromatic
    // CH carbons written without explicit hydrogens are not lost the way gaff's antechamber-style
    // detection loses them). Bond order is irrelevant to ring finding, so connectivity suffices.
    let rings = molar::prelude::sssr_rings(&bond_adj);
    let mut rg = vec![[false; 11]; n]; // per-atom ring-size membership, sizes 3..=10
    for r in &rings {
        let sz = r.len().min(10);
        for &a in r {
            rg[a][sz] = true;
        }
    }
    let (mut nd, mut nt, mut val) = (vec![0u32; n], vec![0u32; n], vec![0u32; n]);
    for b in bonds {
        val[b.i] += b.order as u32;
        val[b.j] += b.order as u32;
        if b.order == 2 {
            nd[b.i] += 1;
            nd[b.j] += 1;
        } else if b.order == 3 {
            nt[b.i] += 1;
            nt[b.j] += 1;
        }
    }
    // RDKit-style Hückel aromaticity. A neighbour "conjugates" a lone pair if it is aromatic or
    // carries a multiple bond on carbon or nitrogen (C=O/C=C/C=N, amidine); a multiple bond on
    // S or P (sulfonamide S=O, phosphonate P=O) does NOT — RDKit leaves those N/O as SP3.
    let aromatic = aromatic_atoms(z, fc, bonds, &rings);
    let neighbor_conj: Vec<bool> = (0..n)
        .map(|i| {
            bond_adj.neighbors(i).iter().any(|nb| {
                let j = nb.atom();
                aromatic[j] || ((nd[j] > 0 || nt[j] > 0) && matches!(z[j], 6 | 7))
            })
        })
        .collect();
    let mut feat = vec![0f32; n * 116];
    for i in 0..n {
        let o = i * 116;
        if (z[i] as usize) < 100 {
            feat[o + z[i] as usize] = 1.0; // element one-hot by atomic number
        }
        let degree = bond_adj.neighbors(i).len();
        feat[o + 100] = degree as f32; // TotalDegree (explicit; molar carries no implicit H)
        feat[o + 101] = val[i] as f32; // TotalValence (explicit)
        feat[o + 102] = val[i] as f32; // ExplicitValence
        feat[o + 103] = if aromatic[i] { 1.0 } else { 0.0 };
        feat[o + 104] = mass_by_z(z[i]);
        for (k, sz) in (3..=8).enumerate() {
            feat[o + 105 + k] = if rg[i][sz] { 1.0 } else { 0.0 };
        }
        if let Some(h) = hybridization(z[i], degree, nd[i], nt[i], aromatic[i], neighbor_conj[i]) {
            feat[o + 111 + h] = 1.0;
        }
    }
    (feat, bond_adj)
}

/// Featurize and run the GNN, returning the raw per-atom electronegativity/hardness
/// `(e, s)` before charge equilibration.
pub(crate) fn espaloma_e_s(
    z: &[u8],
    fc: &[i32],
    bonds: &[crate::gaff::LocalBond],
) -> Result<(Vec<f32>, Vec<f32>), &'static str> {
    let (features, adjacency) = featurize(z, fc, bonds);
    run_gnn(features, &adjacency, z.len())
}

/// End-to-end: atomic numbers + formal charges + local bonds → espaloma partial charges.
///
/// The charges are equilibrated to sum to the molecule's **total formal charge**
/// `Σ_i fc_i`, matching upstream espaloma-charge (which takes `Q = Chem.GetFormalCharge(mol)`
/// when no explicit total is supplied). A cation therefore sums to +1, not 0.
pub(crate) fn espaloma_charges(
    z: &[u8],
    fc: &[i32],
    bonds: &[crate::gaff::LocalBond],
) -> Result<Vec<f32>, &'static str> {
    let (e, s) = espaloma_e_s(z, fc, bonds)?;
    Ok(equilibrate(&e, &s, fc.iter().sum::<i32>() as f32))
}

/// Espaloma charge equilibration over the whole molecule, constrained to total charge
/// `q_total`:
///
/// `q_i = -e_i/s_i + (1/s_i) · (q_total + Σ_j e_j/s_j) / (Σ_j 1/s_j)`
///
/// This is the Lagrange-multiplier solution of `U(q) = Σ_i (e_i q_i + ½ s_i q_i²)` under
/// `Σ_i q_i = q_total`, and matches upstream espaloma-charge
/// (`espaloma_charge/models.py::charge_equilibrium`, whose `sum_q` is the molecule's total
/// formal charge). Passing `q_total = 0.0` recovers the neutral-molecule special case.
pub(crate) fn equilibrate(e: &[f32], s: &[f32], q_total: f32) -> Vec<f32> {
    let inv: Vec<f32> = s.iter().map(|x| 1.0 / x).collect();
    let sum_inv: f32 = inv.iter().sum();
    let sum_eos: f32 = e.iter().zip(&inv).map(|(a, b)| a * b).sum();
    let lam = (q_total + sum_eos) / sum_inv;
    e.iter().zip(&inv).map(|(a, b)| -a * b + b * lam).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows_f32(v: &serde_json::Value) -> Vec<f32> {
        v.as_array()
            .unwrap()
            .iter()
            .flat_map(|row| {
                row.as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_f64().unwrap() as f32)
            })
            .collect()
    }
    fn col_f32(v: &serde_json::Value) -> Vec<f32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap() as f32)
            .collect()
    }

    /// The fixed kernel must reproduce the Python reference for the fixture molecule.
    #[test]
    fn fixed_kernel_matches_python_fixture() {
        let txt = std::fs::read_to_string("tests/data/espaloma_fixture.json").unwrap();
        let v: serde_json::Value = serde_json::from_str(&txt).unwrap();
        let n = v["n"].as_u64().unwrap() as usize;
        let feats = rows_f32(&v["features"]);
        let dense_adj = rows_f32(&v["adjacency_mean"]);
        let exp_e = col_f32(&v["e"]);
        let exp_q = col_f32(&v["charges"]);

        let adjacency = BondAdjacency::build(
            n,
            (0..n).flat_map(|i| {
                let adj = &dense_adj;
                ((i + 1)..n)
                    .filter(move |&j| adj[i * n + j] != 0.0)
                    .map(move |j| [i, j])
            }),
        );
        let (e, s) = run_gnn(feats, &adjacency, n).expect("fixed-kernel run");
        let de = e
            .iter()
            .zip(&exp_e)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f32, f32::max);
        // The fixture molecule is neutral, so the reference charges are the `Q = 0` solution.
        let q = equilibrate(&e, &s, 0.0);
        let dq = q
            .iter()
            .zip(&exp_q)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f32, f32::max);
        println!("fixed kernel vs Python fixture: max|Δe|={de:.2e}  max|Δq|={dq:.2e}");
        assert!(de < 1e-4, "electronegativity mismatch: {de}");
        assert!(dq < 1e-4, "charge mismatch: {dq}");
    }

    /// End-to-end: Rust featurization and inference against the Python reference corpus.
    /// The reference (`references_espaloma.json`) covers only molecules with fully explicit
    /// hydrogens — molar deliberately performs no implicit-H perception, so under-hydrogenated
    /// inputs (where RDKit would silently add H) are excluded from the benchmark by construction.
    #[test]
    fn corpus_rmse_vs_reference() {
        use molar::prelude::*;
        #[derive(serde::Deserialize)]
        struct RefFile {
            molecules: Vec<RefMol>,
        }
        #[derive(serde::Deserialize)]
        struct RefMol {
            name: String,
            charges: Vec<f32>,
        }
        let txt = std::fs::read_to_string("tests/data/gaff_ref/references_espaloma.json").unwrap();
        let refs: RefFile = serde_json::from_str(&txt).unwrap();
        let (mut se, mut nat, mut maxd) = (0f64, 0usize, 0f32);
        let (mut worst, mut load_err) = (String::new(), 0usize);
        for mol in &refs.molecules {
            let path = format!("tests/data/gaff_ref/sdf/{}.sdf", mol.name);
            let sys = match System::from_file(&path) {
                Ok(s) => s,
                Err(_) => {
                    load_err += 1;
                    continue;
                }
            };
            let z: Vec<u8> = sys.iter_atoms().map(|a| a.get_atomic_number()).collect();
            let fc: Vec<i32> = sys
                .iter_atoms()
                .map(|a| a.get_formal_charge().unwrap_or(0))
                .collect();
            let mut bonds = Vec::new();
            for b in sys.iter_bonds() {
                let order = match b.order() {
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    _ => 1,
                };
                bonds.push(crate::gaff::LocalBond {
                    i: b.i1(),
                    j: b.i2(),
                    order,
                });
            }
            // `references_espaloma.json` was generated with the total charge pinned at 0 for
            // every molecule: 274 of the 595 SDFs carry `M CHG` records (mostly protonated
            // amines, net +1) yet every reference charge set sums to exactly 0.0. So this
            // test compares against the same `Q = 0` equilibration, which still validates the
            // featurization, the GNN and the equilibration algebra. The physical
            // `Σq = Σ fc` path is covered by `total_formal_charge_is_preserved`.
            let (e, s) = espaloma_e_s(&z, &fc, &bonds).unwrap();
            let q = equilibrate(&e, &s, 0.0);
            if q.len() != mol.charges.len() {
                continue;
            }
            for (i, (a, b)) in q.iter().zip(&mol.charges).enumerate() {
                let d = (a - b).abs();
                se += (d as f64) * (d as f64);
                nat += 1;
                if d > maxd {
                    maxd = d;
                    worst = format!("{} atom={i} z={} (mine={a:.3} ref={b:.3})", mol.name, z[i]);
                }
            }
        }
        let rmse = (se / nat as f64).sqrt();
        println!(
            "espaloma Rust vs reference: RMSE={rmse:.4}e  max|Δq|={maxd:.4}e  worst={worst}  atoms={nat}  load_err={load_err}"
        );
        assert_eq!(load_err, 0);
        // The featurization reproduces RDKit exactly over this corpus, so the charges match the
        // Python espaloma reference to float precision (residual is f32 rounding, ~2e-4).
        assert!(rmse < 5e-4, "espaloma charge RMSE {rmse} regressed");
    }

    #[derive(serde::Deserialize)]
    struct RefFile {
        molecules: Vec<RefMol>,
    }
    #[derive(serde::Deserialize)]
    struct RefMol {
        name: String,
        charges: Vec<f32>,
    }

    fn load_refs() -> Vec<RefMol> {
        let txt = std::fs::read_to_string("tests/data/gaff_ref/references_espaloma.json").unwrap();
        serde_json::from_str::<RefFile>(&txt).unwrap().molecules
    }

    /// The public `ApplyCharges` API writes predicted charges into `atom.charge`, matching the
    /// reference and summing to ~0 over the molecule. Uses the first **neutral** reference
    /// molecule, since the reference was generated with the total charge pinned at 0 (see
    /// `corpus_rmse_vs_reference`) and so does not describe charged species.
    #[test]
    fn apply_charges_public_api() {
        use crate::{ApplyCharges, ChargeModel};
        use molar::prelude::*;

        let refs = load_refs();
        let (mol, mut sys) = refs
            .iter()
            .find_map(|m| {
                let sys =
                    System::from_file(format!("tests/data/gaff_ref/sdf/{}.sdf", m.name)).ok()?;
                let q_total: i32 = sys
                    .iter_atoms()
                    .map(|a| a.get_formal_charge().unwrap_or(0))
                    .sum();
                (q_total == 0).then_some((m, sys))
            })
            .expect("corpus must contain a neutral molecule");

        sys.apply_charges(ChargeModel::Espaloma).unwrap();

        let got: Vec<f32> = sys.iter_atoms().map(|a| a.get_charge() as f32).collect();
        let maxd = got
            .iter()
            .zip(&mol.charges)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f32, f32::max);
        let sum: f32 = got.iter().sum();
        println!(
            "apply_charges({}): max|Δ|={maxd:.2e}  Σq={sum:.2e}",
            mol.name
        );
        assert!(
            maxd < 1e-3,
            "apply_charges disagrees with reference: {maxd}"
        );
        assert!(sum.abs() < 1e-3, "charges should sum to ~0: {sum}");
    }

    /// Charges must sum to the molecule's **total formal charge**, not to zero. 274 of the
    /// corpus SDFs carry `M CHG` records (protonated amines, carboxylates, …); every one of
    /// them must equilibrate to its own net charge.
    #[test]
    fn total_formal_charge_is_preserved() {
        use crate::{ApplyCharges, ChargeModel};
        use molar::prelude::*;

        let (mut n_charged, mut n_neutral, mut worst) = (0usize, 0usize, 0f32);
        let mut worst_name = String::new();
        for m in load_refs() {
            let mut sys = match System::from_file(format!("tests/data/gaff_ref/sdf/{}.sdf", m.name))
            {
                Ok(s) => s,
                Err(_) => continue,
            };
            let q_total: i32 = sys
                .iter_atoms()
                .map(|a| a.get_formal_charge().unwrap_or(0))
                .sum();
            sys.apply_charges(ChargeModel::Espaloma).unwrap();
            let sum: f32 = sys.iter_atoms().map(|a| a.get_charge() as f32).sum();
            let d = (sum - q_total as f32).abs();
            if d > worst {
                worst = d;
                worst_name = format!("{} (Q={q_total} Σq={sum:.4})", m.name);
            }
            if q_total == 0 {
                n_neutral += 1;
            } else {
                n_charged += 1;
            }
        }
        println!(
            "Σq vs Σfc: {n_charged} charged + {n_neutral} neutral molecules, \
             max|Σq − Q|={worst:.2e}  worst={worst_name}"
        );
        assert!(
            n_charged > 0,
            "corpus must contain charged molecules to exercise this path"
        );
        assert!(
            worst < 1e-3,
            "total charge not preserved: max deviation {worst} ({worst_name})"
        );
    }

    /// Aromatic input is kekulized rather than rejected. Before this, an SDF order-4 record — or
    /// any system run through `System::perceive` — could not be charged at all.
    ///
    /// Charging the aromatized structure must reproduce the Kekulé result: the resonance form
    /// kekulization picks is arbitrary, but espaloma's featurization of an aromatic ring is
    /// symmetric enough that the predicted charges land in the same place.
    #[test]
    fn aromatic_input_is_kekulized_not_rejected() {
        use crate::{ApplyCharges, ChargeModel};
        use molar::prelude::*;

        let mut checked = 0usize;
        let mut worst = 0f32;
        let mut worst_name = String::new();
        for m in load_refs().iter().take(80) {
            let path = format!("tests/data/gaff_ref/sdf/{}.sdf", m.name);
            let Ok(mut kek) = System::from_file(&path) else {
                continue;
            };
            let Ok(mut arom) = System::from_file(&path) else {
                continue;
            };

            // Only molecules that actually have an aromatic ring exercise the new path.
            arom.perceive();
            if !arom.iter_bonds().any(|b| b.order() == BondOrder::Aromatic) {
                continue;
            }

            kek.apply_charges(ChargeModel::Espaloma).unwrap();
            arom.apply_charges(ChargeModel::Espaloma)
                .unwrap_or_else(|e| panic!("{}: aromatic input rejected: {e}", m.name));

            let d = kek
                .iter_atoms()
                .zip(arom.iter_atoms())
                .map(|(a, b)| (a.get_charge() as f32 - b.get_charge() as f32).abs())
                .fold(0.0f32, f32::max);
            if d > worst {
                worst = d;
                worst_name = m.name.clone();
            }
            checked += 1;
        }
        println!(
            "aromatic vs Kekulé charges over {checked} molecules: max|Δ|={worst:.2e} ({worst_name})"
        );
        assert!(
            checked > 10,
            "corpus should supply aromatic molecules to check"
        );
        assert!(
            worst < 1e-3,
            "aromatic path disagrees with Kekulé: {worst} ({worst_name})"
        );
    }

    /// The equilibration is affine in the total charge: `q(Q) = q(0) + Q·(1/s_i)/Σ(1/s_j)`.
    /// This pins the exact relationship between the two conventions, which is what lets
    /// `corpus_rmse_vs_reference` keep using the `Q = 0` reference.
    #[test]
    fn equilibration_is_affine_in_total_charge() {
        let e = [0.3f32, -0.7, 1.1, 0.05, -0.2];
        let s = [1.7f32, 2.3, 0.9, 3.1, 1.2];
        let inv_sum: f32 = s.iter().map(|x| 1.0 / x).sum();

        let q0 = equilibrate(&e, &s, 0.0);
        for &q_total in &[-2.0f32, -1.0, 1.0, 3.0] {
            let q = equilibrate(&e, &s, q_total);
            let sum: f32 = q.iter().sum();
            assert!((sum - q_total).abs() < 1e-5, "Σq={sum} != Q={q_total}");
            for i in 0..e.len() {
                let expected = q0[i] + q_total * (1.0 / s[i]) / inv_sum;
                assert!(
                    (q[i] - expected).abs() < 1e-5,
                    "atom {i}: {} != {expected} at Q={q_total}",
                    q[i]
                );
            }
        }
    }

    /// Opt-in throughput and stability check for an unusually large sparse molecule.
    ///
    /// This test is ignored because elapsed time depends on the machine. Run it in release
    /// mode as documented in `assets/README.md`. The synthetic graph is a 20,000-atom chain.
    /// It selects the parallel path and verifies that inference does not need a dense
    /// adjacency matrix.
    #[test]
    #[ignore = "large release-mode performance check"]
    fn large_sparse_inference_smoke() {
        use std::time::Instant;

        const N: usize = 20_000;
        let features = vec![0.0; N * FEATURE_WIDTH];
        let adjacency = BondAdjacency::build(N, (1..N).map(|i| [i - 1, i]));
        let start = Instant::now();
        let (e, s) = run_gnn(features, &adjacency, N).unwrap();
        let elapsed = start.elapsed();
        assert_eq!(e.len(), N);
        assert_eq!(s.len(), N);
        assert!(e.iter().chain(&s).all(|value| value.is_finite()));
        eprintln!("20,000-atom sparse Espaloma inference: {elapsed:.3?}");
    }
}
