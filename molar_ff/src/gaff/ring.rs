//! Ring perception: detect all chordless simple rings (size 3..=10) and record per-atom
//! ring-size counts.
//!
//! Detects **all chordless simple rings of size 3..=10** (this ring set differs from a
//! minimum-cycle SSSR on fused/bridged systems), and records per-atom ring-size counts
//! `rg[k]`.

/// An atom is ring-eligible (same test used both for ring start atoms and for atoms
/// traversed during cycle detection): C with >2 neighbours, N, P, O/S with !=1
/// neighbour. Everything else is skipped.
fn eligible(z: u8, connum: usize) -> bool {
    match z {
        6 => connum > 2,
        7 => true,
        15 => true,
        8 => connum != 1,
        16 => connum != 1,
        _ => false,
    }
}

/// Detect all chordless simple rings (size 3..=10). `con[a]` lists the neighbours of
/// atom `a` in input-bond order. Returns each ring as a sorted list of member atoms.
pub fn detect_rings(z: &[u8], con: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let n = z.len();
    let mut raw: Vec<Vec<usize>> = Vec::new();
    for i in 0..n {
        if !eligible(z[i], con[i].len()) {
            continue;
        }
        let mut path: Vec<usize> = Vec::new();
        walk(i, con, z, &mut path, &mut raw);
    }
    purify(con, raw)
}

/// Recursive path-DFS for cycle detection: traverses at most the first 4 neighbours
/// (`con[0..3]`), records a ring when the current path (length 2..=9) can be closed back
/// to its origin through one of the origin's first 4 neighbours.
fn walk(cur: usize, con: &[Vec<usize>], z: &[u8], path: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    path.push(cur);
    let sn = path.len();
    if sn <= 10 {
        let a0 = path[0];
        for i in 0..con[cur].len().min(4) {
            let start = con[cur][i];
            if !eligible(z[start], con[start].len()) {
                continue;
            }
            if path.contains(&start) {
                continue;
            }
            // ring closure: path has 2..=9 atoms and `start` neighbours the origin
            if (2..=9).contains(&sn) && con[a0].iter().take(4).any(|&x| x == start) {
                let mut r = path.clone();
                r.push(start);
                out.push(r);
            }
            walk(start, con, z, path, out);
        }
    }
    path.pop();
}

/// Dedup identical rings, then drop any ring containing a chord (a member with exactly
/// 3 of its neighbours inside the ring — e.g. a fused-ring envelope).
fn purify(con: &[Vec<usize>], raw: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut unique: Vec<Vec<usize>> = Vec::new();
    for mut r in raw {
        r.sort_unstable();
        if !unique.iter().any(|u| *u == r) {
            unique.push(r);
        }
    }
    unique
        .into_iter()
        .filter(|r| {
            for &m in r {
                let cnt = con[m].iter().filter(|&&nb| r.contains(&nb)).count();
                if cnt == 3 {
                    return false;
                }
            }
            true
        })
        .collect()
}

/// Per-atom ring-size counts. `rg[a][0]` = total rings atom `a` is in;
/// `rg[a][k]` = number of k-membered rings (k = 3..=10).
pub fn ring_property(n: usize, rings: &[Vec<usize>]) -> Vec<[u16; 11]> {
    let mut rg = vec![[0u16; 11]; n];
    for r in rings {
        let sz = r.len();
        for &m in r {
            rg[m][0] += 1;
            if sz <= 10 {
                rg[m][sz] += 1;
            }
        }
    }
    rg
}
