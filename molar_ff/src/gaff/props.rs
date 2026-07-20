//! Per-atom property precompute — the bond/neighbour counts the rule matcher needs.
//!
//! With Kekulé input (bond orders 1/2/3 only) the lowercase/uppercase bond counts
//! coincide: `sb == SB` (order-1 bonds), `db == DB` (order-2), `tb == TB` (order-3);
//! `AB` and `DL` are always 0 and are not stored.

use super::LocalBond;

/// Neighbour traversal is capped at six atoms: several per-atom counting loops only ever
/// look at the first six neighbours. This constant is that cap.
const MAX_CON: usize = 6;

/// Precomputed per-atom properties, indexed by local atom id.
pub(crate) struct Props {
    /// Coordination number (number of bonds).
    pub connum: Vec<usize>,
    /// Number of directly attached hydrogens (capped at 6 neighbours).
    pub nh: Vec<usize>,
    /// Number of electron-withdrawing neighbours (`ewd == 1`), capped at 6.
    pub ewd_neigh: Vec<usize>,
    /// Single-bond count per atom (`sb == SB`).
    pub sb: Vec<u16>,
    /// Double-bond count per atom (`db == DB`).
    pub db: Vec<u16>,
    /// Triple-bond count per atom (`tb == TB`).
    pub tb: Vec<u16>,
}

/// Compute [`Props`] from atomic numbers, neighbour lists, bond list and the EW flags.
pub(crate) fn compute(
    z: &[u8],
    con: &[Vec<usize>],
    bonds: &[LocalBond],
    ewd: &[i8],
) -> Props {
    let n = z.len();

    let connum: Vec<usize> = con.iter().map(|c| c.len()).collect();

    let nh: Vec<usize> = (0..n)
        .map(|i| {
            con[i]
                .iter()
                .take(MAX_CON)
                .filter(|&&j| z[j] == 1)
                .count()
        })
        .collect();

    let ewd_neigh: Vec<usize> = (0..n)
        .map(|i| {
            con[i]
                .iter()
                .take(MAX_CON)
                .filter(|&&j| ewd[j] == 1)
                .count()
        })
        .collect();

    let mut sb = vec![0u16; n];
    let mut db = vec![0u16; n];
    let mut tb = vec![0u16; n];
    for b in bonds {
        match b.order {
            1 => {
                sb[b.i] += 1;
                sb[b.j] += 1;
            }
            2 => {
                db[b.i] += 1;
                db[b.j] += 1;
            }
            3 => {
                tb[b.i] += 1;
                tb[b.j] += 1;
            }
            _ => {} // aromatic/other codes do not occur with Kekulé input
        }
    }

    Props { connum, nh, ewd_neigh, sb, db, tb }
}
