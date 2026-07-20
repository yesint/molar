//! Aromaticity / ring-class perception.
//!
//! Assigns the five GAFF ring classes AR1..AR5 (as per-atom counts, since a fused atom
//! can belong to rings of several classes), plus the per-atom electron-withdrawing flag
//! `ewd` and the "not in any classified ring" flag `nr`. This is NOT a Hückel 4n+2 test
//! — it is an element+connectivity heuristic.

use super::LocalBond;

/// Result of aromaticity perception. Indexed by atom.
pub struct Arom {
    /// `ar[a][k]` = number of AR`k` rings atom `a` belongs to (k = 1..=5; index 0 unused).
    pub ar: Vec<[u16; 6]>,
    /// Electron-withdrawing flag: 1 = withdrawing, 0 = other.
    pub ewd: Vec<i8>,
    /// True if the atom is in no classified ring.
    pub nr: Vec<bool>,
}

/// Per-atom π/hybridisation score from element + coordination only.
fn init_arom(z: u8, connum: usize) -> i32 {
    match z {
        6 => match connum {
            3 => 2,
            4 => -2,
            _ => 0,
        },
        7 => {
            if connum <= 3 {
                2
            } else {
                0
            }
        }
        8 => {
            if connum == 2 {
                1
            } else {
                0
            }
        }
        15 => match connum {
            2 => 2,
            3 => 1,
            c if c >= 4 => -1,
            _ => 0,
        },
        16 => {
            if connum == 2 {
                1
            } else if connum >= 3 {
                -1
            } else {
                0
            }
        }
        _ => 0,
    }
}

fn ewd_flag(z: u8) -> i8 {
    match z {
        7 | 8 | 16 | 9 | 17 | 35 | 53 => 1,
        _ => 0,
    }
}

pub fn aromatic(
    z: &[u8],
    con: &[Vec<usize>],
    bonds: &[LocalBond],
    rings: &[Vec<usize>],
    rg: &[[u16; 11]],
) -> Arom {
    let n = z.len();
    let initarom: Vec<i32> = (0..n).map(|i| init_arom(z[i], con[i].len())).collect();
    let ewd: Vec<i8> = (0..n).map(|i| ewd_flag(z[i])).collect();
    let mut ar = vec![[0u16; 6]; n];
    let mut nr = vec![true; n];

    for r in rings {
        let num = r.len() as i32;
        let tmpint: i32 = r.iter().map(|&m| initarom[m]).sum();

        // AR5: pure aliphatic ring (all sp3 C).
        if tmpint == -2 * num {
            for &m in r {
                ar[m][5] += 1;
            }
            continue;
        }

        // AR4: ring containing an sp3/hypervalent atom (negative π/hybridisation score).
        if r.iter().any(|&m| initarom[m] < 0) {
            for &m in r {
                ar[m][4] += 1;
            }
            continue;
        }

        // AR3: planar ring with an exocyclic double bond to a non-ring atom.
        if (num..=2 * num).contains(&tmpint) {
            let mut found = false;
            for b in bonds {
                let mut index = 0;
                if r.contains(&b.i) && rg[b.j][0] == 0 {
                    index += 1;
                }
                if r.contains(&b.j) && rg[b.i][0] == 0 {
                    index += 1;
                }
                if index == 1 && (b.order == 2 || b.order == 8) {
                    found = true;
                    break;
                }
            }
            if found {
                for &m in r {
                    ar[m][3] += 1;
                }
                continue;
            }
        }

        // AR1: pure aromatic 6-ring (benzene/pyridine); every ring N/P must carry a π bond.
        if tmpint == 12 && num == 6 {
            let mut bad = false;
            for &m in r {
                if z[m] == 7 || z[m] == 15 {
                    let has_pi = bonds.iter().any(|b| {
                        (b.i == m || b.j == m) && (b.order == 8 || b.order == 2 || b.order == 10)
                    });
                    if !has_pi {
                        bad = true;
                    }
                }
            }
            if !bad {
                for &m in r {
                    ar[m][1] += 1;
                }
                continue;
            }
        }

        // AR2: other planar/aromatic ring (5-membered aromatics, non-benzene 6-rings, …).
        if tmpint >= num + 3 {
            for &m in r {
                ar[m][2] += 1;
            }
            continue;
        }

        // Otherwise: bucket into AR4 ("other rings").
        for &m in r {
            ar[m][4] += 1;
        }
    }

    for i in 0..n {
        if (1..=5).any(|k| ar[i][k] > 0) {
            nr[i] = false;
        }
    }

    Arom { ar, ewd, nr }
}
