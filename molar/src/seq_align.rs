//Needleman Wunsch global alignment with affine gaps that mirrors 
// rust bio’s pairwise aligner logic.
// Made by ChatGPT 5.2

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentOperation {
    Match,
    Subst,
    Ins, // gap in y (consume x)
    Del, // gap in x (consume y)
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    pub operations: Vec<AlignmentOperation>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layer {
    S,
    I,
    D,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TbS {
    Start,
    DiagMatch,
    DiagSubst,
    FromI,
    FromD,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TbGap {
    FromS, // opened a gap from S
    FromGap, // extended an existing gap (I->I or D->D)
}

/// Global alignment (Needleman–Wunsch) with affine gap penalties:
/// - gap length k costs: gap_open + gap_extend * k
/// - scoring is generic via `score_fn(&T,&T) -> i32`
///
/// Complexity: O(m*n) time, O(m*n) memory (traceback stored).
pub fn global_align_affine<T, F>(
    x: &[T],
    y: &[T],
    gap_open: i32,
    gap_extend: i32,
    score_fn: F,
) -> Alignment
where
    T: PartialEq,
    F: Fn(&T, &T) -> i32,
{
    // A "negative infinity" that is safe to add to without overflowing.
    // Using /4 gives headroom even if you add many penalties.
    const MIN_SCORE: i32 = i32::MIN / 4;

    #[inline]
    fn idx(i: usize, j: usize, n: usize) -> usize {
        // row-major, (m+1) x (n+1)
        i * (n + 1) + j
    }

    #[inline]
    fn add(a: i32, b: i32) -> i32 {
        // If a is effectively -inf, keep it -inf.
        // Otherwise saturating add avoids underflow/overflow surprises.
        if a <= MIN_SCORE / 2 {
            MIN_SCORE
        } else {
            a.saturating_add(b)
        }
    }

    let m = x.len();
    let n = y.len();
    let size = (m + 1) * (n + 1);

    let mut s = vec![MIN_SCORE; size];
    let mut i_mat = vec![MIN_SCORE; size];
    let mut d_mat = vec![MIN_SCORE; size];

    let mut tb_s = vec![TbS::Start; size];
    let mut tb_i = vec![TbGap::FromS; size];
    let mut tb_d = vec![TbGap::FromS; size];

    // (0,0)
    s[idx(0, 0, n)] = 0;
    tb_s[idx(0, 0, n)] = TbS::Start;

    // Initialize first column (j = 0): must consume x via Ins operations (gap in y).
    for ii in 1..=m {
        let pos = idx(ii, 0, n);
        i_mat[pos] = gap_open + gap_extend * (ii as i32); // gap length ii
        tb_i[pos] = if ii == 1 { TbGap::FromS } else { TbGap::FromGap };

        s[pos] = i_mat[pos];
        tb_s[pos] = TbS::FromI;
        // d_mat[pos] stays -inf
    }

    // Initialize first row (i = 0): must consume y via Del operations (gap in x).
    for jj in 1..=n {
        let pos = idx(0, jj, n);
        d_mat[pos] = gap_open + gap_extend * (jj as i32); // gap length jj
        tb_d[pos] = if jj == 1 { TbGap::FromS } else { TbGap::FromGap };

        s[pos] = d_mat[pos];
        tb_s[pos] = TbS::FromD;
        // i_mat[pos] stays -inf
    }

    // Fill DP.
    for ii in 1..=m {
        for jj in 1..=n {
            let pos = idx(ii, jj, n);

            // I(ii,jj): x[ii-1] aligned to a gap (Ins) => come from (ii-1,jj)
            let from_i = add(i_mat[idx(ii - 1, jj, n)], gap_extend);
            let from_s_open = add(s[idx(ii - 1, jj, n)], gap_open + gap_extend);

            // strict '>' like rust-bio: ties prefer opening from S
            if from_i > from_s_open {
                i_mat[pos] = from_i;
                tb_i[pos] = TbGap::FromGap;
            } else {
                i_mat[pos] = from_s_open;
                tb_i[pos] = TbGap::FromS;
            }

            // D(ii,jj): y[jj-1] aligned to a gap (Del) => come from (ii,jj-1)
            let from_d = add(d_mat[idx(ii, jj - 1, n)], gap_extend);
            let from_s_open = add(s[idx(ii, jj - 1, n)], gap_open + gap_extend);

            if from_d > from_s_open {
                d_mat[pos] = from_d;
                tb_d[pos] = TbGap::FromGap;
            } else {
                d_mat[pos] = from_s_open;
                tb_d[pos] = TbGap::FromS;
            }

            // Diagonal: match/subst
            let diag = add(
                s[idx(ii - 1, jj - 1, n)],
                score_fn(&x[ii - 1], &y[jj - 1]),
            );

            // S(ii,jj) = max(diag, I(ii,jj), D(ii,jj)) with strict '>' priority:
            // diag wins ties over I, and I wins ties over D (same comparison order as rust-bio).
            let mut best = diag;
            tb_s[pos] = if x[ii - 1] == y[jj - 1] {
                TbS::DiagMatch
            } else {
                TbS::DiagSubst
            };

            if i_mat[pos] > best {
                best = i_mat[pos];
                tb_s[pos] = TbS::FromI;
            }
            if d_mat[pos] > best {
                best = d_mat[pos];
                tb_s[pos] = TbS::FromD;
            }

            s[pos] = best;
        }
    }

    // Traceback from (m,n), starting in S layer (global).
    let mut ops = Vec::with_capacity(m + n);
    let mut ii = m;
    let mut jj = n;
    let mut layer = Layer::S;

    while !(ii == 0 && jj == 0 && layer == Layer::S) {
        match layer {
            Layer::S => match tb_s[idx(ii, jj, n)] {
                TbS::Start => break,
                TbS::DiagMatch => {
                    ops.push(AlignmentOperation::Match);
                    ii -= 1;
                    jj -= 1;
                    layer = Layer::S;
                }
                TbS::DiagSubst => {
                    ops.push(AlignmentOperation::Subst);
                    ii -= 1;
                    jj -= 1;
                    layer = Layer::S;
                }
                TbS::FromI => {
                    layer = Layer::I;
                }
                TbS::FromD => {
                    layer = Layer::D;
                }
            },
            Layer::I => {
                // An insertion consumes one element of x (ii-1) while jj stays.
                ops.push(AlignmentOperation::Ins);
                let prev = tb_i[idx(ii, jj, n)];
                ii -= 1;
                layer = match prev {
                    TbGap::FromGap => Layer::I,
                    TbGap::FromS => Layer::S,
                };
            }
            Layer::D => {
                // A deletion consumes one element of y (jj-1) while ii stays.
                ops.push(AlignmentOperation::Del);
                let prev = tb_d[idx(ii, jj, n)];
                jj -= 1;
                layer = match prev {
                    TbGap::FromGap => Layer::D,
                    TbGap::FromS => Layer::S,
                };
            }
        }
    }

    ops.reverse();

    Alignment {
        score: s[idx(m, n, n)],
        operations: ops,
    }
}
