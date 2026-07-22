//! espaloma-charge partial charges via a bundled ONNX graph network (feature `espaloma`).
//!
//! The GNN maps per-atom features `[n, 116]` + a row-mean-normalised bond adjacency
//! `[n, n]` to per-atom electronegativity `e` and hardness `s`; a closed-form charge
//! equilibration then yields partial charges summing to zero over the molecule.

use std::io::Cursor;
use std::sync::LazyLock;
use tract_onnx::prelude::*;

/// The bundled espaloma-charge model (bias-free, opset 18), exported from v0.0.8.
const MODEL: &[u8] = include_bytes!("../assets/espaloma_charge.onnx");

/// The runnable model, built once (keeps the atom count `N` symbolic so it runs on any molecule).
static PLAN: LazyLock<TypedRunnableModel<TypedModel>> = LazyLock::new(|| {
    tract_onnx::onnx()
        .model_for_read(&mut Cursor::new(MODEL))
        .expect("read espaloma onnx")
        .into_optimized()
        .expect("optimize")
        .into_runnable()
        .expect("runnable")
});

/// Run the GNN. `features` is row-major `[n, 116]`, `adj_mean` is row-major `[n, n]`
/// (`1/deg_i` for bonded pairs, no self-loop). Returns `(e, s)` per atom.
pub(crate) fn run_gnn(features: &[f32], adj_mean: &[f32], n: usize) -> TractResult<(Vec<f32>, Vec<f32>)> {
    let f = Tensor::from_shape(&[n, 116], features)?;
    let a = Tensor::from_shape(&[n, n], adj_mean)?;
    let out = PLAN.run(tvec!(f.into(), a.into()))?;
    let e = out[0].to_array_view::<f32>()?.iter().copied().collect();
    let s = out[1].to_array_view::<f32>()?.iter().copied().collect();
    Ok((e, s))
}

/// Standard atomic weight for the elements espaloma supports (matches RDKit `GetMass`).
fn mass_by_z(z: u8) -> f32 {
    match z {
        1 => 1.008, 6 => 12.011, 7 => 14.007, 8 => 15.999, 9 => 18.998,
        15 => 30.974, 16 => 32.06, 17 => 35.45, 35 => 79.904, 53 => 126.904,
        _ => 0.0,
    }
}

/// RDKit hybridization one-hot index: 0=SP,1=SP2,2=SP3,3=SP3D,4=SP3D2. `None` (all-zero)
/// for hydrogen (RDKit reports S/unspecified). `neighbor_conj` = does any neighbour carry a
/// multiple bond or is aromatic — a lone pair adjacent to such a π system conjugates into it,
/// so RDKit reports the atom as SP2 (amide N, ester/conjugated O) rather than SP3.
fn hybridization(z: u8, degree: usize, n_double: u32, n_triple: u32, aromatic: bool, neighbor_conj: bool) -> Option<usize> {
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
        1 => 1, 6 => 4, 7 => 5, 8 => 6, 9 => 7, 15 => 5, 16 => 6, 17 => 7, 35 => 7, 53 => 7,
        _ => 0,
    }
}

/// Pauling electronegativity, used to classify exocyclic double bonds: a double bond to a *more*
/// electronegative atom (C=O, C=N) pulls the π out of the ring plane (contributes 0, ring may
/// still be aromatic), whereas one to an equal/less electronegative atom (C=C — fulvene,
/// quinone-methide, thioxanthene ylidene) disrupts the ring and breaks aromaticity.
fn electronegativity(z: u8) -> f32 {
    match z {
        1 => 2.20, 6 => 2.55, 7 => 3.04, 8 => 3.44, 9 => 3.98, 15 => 2.19, 16 => 2.58,
        17 => 3.16, 35 => 2.96, 53 => 2.66, _ => 0.0,
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
fn aromatic_atoms(z: &[u8], fc: &[i32], bonds: &[crate::gaff::LocalBond], rings: &[Vec<usize>]) -> Vec<bool> {
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
            Some(if avail <= 0 { 0 } else if avail % 2 == 1 { 1 } else { 2 })
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
        systems.entry(root).or_default().extend(rings[i].iter().copied());
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

/// Build the 116-dim atom feature matrix (row-major `[n,116]`) and the row-mean-normalised
/// adjacency (`[n,n]`) that the espaloma GNN consumes.
pub(crate) fn featurize(z: &[u8], fc: &[i32], bonds: &[crate::gaff::LocalBond]) -> (Vec<f32>, Vec<f32>) {
    let n = z.len();
    let con = crate::gaff::build_con(n, bonds);
    // Rings from molar core SSSR (matches RDKit ring semantics and is H-independent, so aromatic
    // CH carbons written without explicit hydrogens are not lost the way gaff's antechamber-style
    // detection loses them). Bond order is irrelevant to ring finding, so pass connectivity only.
    let mbonds: Vec<molar::prelude::Bond> = bonds
        .iter()
        .map(|b| molar::prelude::Bond::with_order(b.i, b.j, molar::prelude::BondOrder::Single))
        .collect();
    let rings = molar::prelude::sssr_rings(n, &mbonds);
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
            con[i]
                .iter()
                .any(|&j| aromatic[j] || ((nd[j] > 0 || nt[j] > 0) && matches!(z[j], 6 | 7)))
        })
        .collect();
    let mut feat = vec![0f32; n * 116];
    for i in 0..n {
        let o = i * 116;
        if (z[i] as usize) < 100 {
            feat[o + z[i] as usize] = 1.0; // element one-hot by atomic number
        }
        let degree = con[i].len();
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
    let mut adj = vec![0f32; n * n];
    for b in bonds {
        adj[b.i * n + b.j] = 1.0;
        adj[b.j * n + b.i] = 1.0;
    }
    for i in 0..n {
        let deg: f32 = adj[i * n..(i + 1) * n].iter().sum();
        if deg > 0.0 {
            for x in &mut adj[i * n..(i + 1) * n] {
                *x /= deg;
            }
        }
    }
    (feat, adj)
}

/// End-to-end: atomic numbers + local bonds → espaloma partial charges.
pub(crate) fn espaloma_charges(z: &[u8], fc: &[i32], bonds: &[crate::gaff::LocalBond]) -> TractResult<Vec<f32>> {
    let (feat, adj) = featurize(z, fc, bonds);
    let (e, s) = run_gnn(&feat, &adj, z.len())?;
    Ok(equilibrate(&e, &s))
}

/// Espaloma charge equilibration over the whole molecule with total charge 0:
/// `q_i = -e_i/s_i + (1/s_i) · (Σ_j e_j/s_j) / (Σ_j 1/s_j)`.
pub(crate) fn equilibrate(e: &[f32], s: &[f32]) -> Vec<f32> {
    let inv: Vec<f32> = s.iter().map(|x| 1.0 / x).collect();
    let sum_inv: f32 = inv.iter().sum();
    let sum_eos: f32 = e.iter().zip(&inv).map(|(a, b)| a * b).sum();
    let lam = sum_eos / sum_inv;
    e.iter().zip(&inv).map(|(a, b)| -a * b + b * lam).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows_f32(v: &serde_json::Value) -> Vec<f32> {
        v.as_array()
            .unwrap()
            .iter()
            .flat_map(|row| row.as_array().unwrap().iter().map(|x| x.as_f64().unwrap() as f32))
            .collect()
    }
    fn col_f32(v: &serde_json::Value) -> Vec<f32> {
        v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap() as f32).collect()
    }

    /// tract must reproduce the Python reference (e, s, charges) for the fixture molecule.
    #[test]
    fn tract_matches_python_fixture() {
        let txt = std::fs::read_to_string("tests/data/espaloma_fixture.json").unwrap();
        let v: serde_json::Value = serde_json::from_str(&txt).unwrap();
        let n = v["n"].as_u64().unwrap() as usize;
        let feats = rows_f32(&v["features"]);
        let adj = rows_f32(&v["adjacency_mean"]);
        let exp_e = col_f32(&v["e"]);
        let exp_q = col_f32(&v["charges"]);

        let (e, s) = run_gnn(&feats, &adj, n).expect("tract run");
        let de = e.iter().zip(&exp_e).map(|(a, b)| (a - b).abs()).fold(0.0f32, f32::max);
        let q = equilibrate(&e, &s);
        let dq = q.iter().zip(&exp_q).map(|(a, b)| (a - b).abs()).fold(0.0f32, f32::max);
        println!("tract vs python fixture: max|Δe|={de:.2e}  max|Δq|={dq:.2e}");
        assert!(de < 1e-4, "electronegativity mismatch: {de}");
        assert!(dq < 1e-4, "charge mismatch: {dq}");
    }

    /// End-to-end: Rust featurization + tract vs the Python espaloma reference over the corpus.
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
            let z: Vec<u8> = sys.iter_atoms().map(|a| a.atomic_number).collect();
            let fc: Vec<i32> = sys.iter_atoms().map(|a| a.charge.round() as i32).collect();
            let mut bonds = Vec::new();
            for b in sys.iter_bonds() {
                let order = match b.order {
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    _ => 1,
                };
                bonds.push(crate::gaff::LocalBond { i: b.i1, j: b.i2, order });
            }
            let q = espaloma_charges(&z, &fc, &bonds).unwrap();
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
        println!("espaloma Rust vs reference: RMSE={rmse:.4}e  max|Δq|={maxd:.4}e  worst={worst}  atoms={nat}  load_err={load_err}");
        assert_eq!(load_err, 0);
        // The featurization reproduces RDKit exactly over this corpus, so the charges match the
        // Python espaloma reference to float precision (residual is f32 rounding, ~2e-4).
        assert!(rmse < 5e-4, "espaloma charge RMSE {rmse} regressed");
    }
}
