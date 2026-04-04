//! Ungapped and gapped extension for BLAST.

use crate::matrix::ScoringMatrix;

/// Result of ungapped extension.
#[derive(Debug, Clone)]
pub struct UngappedHit {
    pub score: i32,
    pub q_start: usize,
    pub q_end: usize,   // exclusive
    pub s_start: usize,
    pub s_end: usize,   // exclusive
}

/// Perform ungapped extension from a seed position.
/// Extends left and right from (q_pos, s_pos).
/// Stops when score drops > x_drop below best seen score.
///
/// Returns the best-scoring ungapped alignment.
pub fn ungapped_extend(
    query: &[u8],
    subject: &[u8],
    q_pos: usize,
    s_pos: usize,
    matrix: &ScoringMatrix,
    x_drop: i32,
) -> UngappedHit {
    let mut best = 0i32;
    let mut best_q_start = q_pos;
    let mut best_s_start = s_pos;

    // Extend left
    let mut score = 0i32;
    let mut qi = q_pos as isize - 1;
    let mut si = s_pos as isize - 1;
    while qi >= 0 && si >= 0 {
        score += matrix.score(query[qi as usize], subject[si as usize]);
        if score > best {
            best = score;
            best_q_start = qi as usize;
            best_s_start = si as usize;
        } else if best - score > x_drop {
            break;
        }
        qi -= 1;
        si -= 1;
    }
    let left_score = best;
    let q_start = best_q_start;
    let s_start = best_s_start;

    // Extend right from q_pos, s_pos
    let mut qi = q_pos;
    let mut si = s_pos;
    score = 0;
    best = 0;
    let mut best_q_end = q_pos;
    let mut best_s_end = s_pos;
    while qi < query.len() && si < subject.len() {
        score += matrix.score(query[qi], subject[si]);
        if score > best {
            best = score;
            best_q_end = qi + 1;
            best_s_end = si + 1;
        } else if best - score > x_drop {
            break;
        }
        qi += 1;
        si += 1;
    }
    let right_score = best;

    // Center residue score
    let center_score = if q_pos < query.len() && s_pos < subject.len() {
        matrix.score(query[q_pos], subject[s_pos])
    } else {
        0
    };

    UngappedHit {
        score: left_score + center_score + right_score,
        q_start,
        q_end: best_q_end,
        s_start,
        s_end: best_s_end,
    }
}

/// Result of gapped extension (full alignment).
#[derive(Debug, Clone)]
pub struct GappedHit {
    pub score: i32,
    pub q_start: usize,
    pub q_end: usize,   // exclusive
    pub s_start: usize,
    pub s_end: usize,   // exclusive
    pub num_identities: usize,
    pub num_gaps: usize,
    /// Alignment strings: query, midline, subject
    pub query_aln: Vec<u8>,
    pub midline: Vec<u8>,
    pub subject_aln: Vec<u8>,
}

/// Affine gap DP with X-dropoff.
/// Uses two DP vectors: best[] (current best score ending with a match/mismatch at column j)
/// and best_gap[] (ending with a gap in the query at column j).
///
/// This implements the standard BLAST gapped extension starting from a "seed" hit.
/// We extend from a center point (q_center, s_center) bidirectionally.
#[allow(clippy::too_many_arguments)]
pub fn gapped_extend(
    query: &[u8],
    subject: &[u8],
    q_center: usize,
    s_center: usize,
    matrix: &ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> GappedHit {
    // Extend right
    let right = extend_one_direction(
        &query[q_center..],
        &subject[s_center..],
        matrix,
        gap_open,
        gap_extend,
        x_drop,
    );

    // Extend left (reversed)
    let query_rev: Vec<u8> = query[..q_center].iter().rev().cloned().collect();
    let subject_rev: Vec<u8> = subject[..s_center].iter().rev().cloned().collect();
    let left = extend_one_direction(
        &query_rev,
        &subject_rev,
        matrix,
        gap_open,
        gap_extend,
        x_drop,
    );

    // Combine
    let total_score = left.score + right.score;
    let q_start = q_center - left.query_len;
    let q_end = q_center + right.query_len;
    let s_start = s_center - left.subject_len;
    let s_end = s_center + right.subject_len;

    // Build alignment strings
    let mut query_aln = Vec::new();
    let mut midline = Vec::new();
    let mut subject_aln = Vec::new();
    // Left part (reversed)
    let mut left_q_aln: Vec<u8> = left.query_aln.iter().rev().cloned().collect();
    let mut left_m_aln: Vec<u8> = left.midline.iter().rev().cloned().collect();
    let mut left_s_aln: Vec<u8> = left.subject_aln.iter().rev().cloned().collect();
    query_aln.append(&mut left_q_aln);
    midline.append(&mut left_m_aln);
    subject_aln.append(&mut left_s_aln);

    // Right part
    query_aln.extend_from_slice(&right.query_aln);
    midline.extend_from_slice(&right.midline);
    subject_aln.extend_from_slice(&right.subject_aln);

    let num_identities = midline.iter().filter(|&&c| c == b'|').count();
    let num_gaps = query_aln.iter().filter(|&&c| c == b'-').count()
        + subject_aln.iter().filter(|&&c| c == b'-').count();

    GappedHit {
        score: total_score,
        q_start,
        q_end,
        s_start,
        s_end,
        num_identities,
        num_gaps,
        query_aln,
        midline,
        subject_aln,
    }
}

struct DirectionResult {
    score: i32,
    query_len: usize,
    subject_len: usize,
    query_aln: Vec<u8>,
    midline: Vec<u8>,
    subject_aln: Vec<u8>,
}

/// Extend in one direction using affine gap DP.
fn extend_one_direction(
    query: &[u8],
    subject: &[u8],
    matrix: &ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> DirectionResult {
    let qlen = query.len();
    let slen = subject.len();

    if qlen == 0 || slen == 0 {
        return DirectionResult {
            score: 0,
            query_len: 0,
            subject_len: 0,
            query_aln: vec![],
            midline: vec![],
            subject_aln: vec![],
        };
    }

    // Banded X-drop DP with traceback.
    // Only columns within the active band [j_lo, j_hi] are computed per row.
    // The band shrinks from the left when cells drop below X-drop threshold
    // and extends rightward only as far as cells remain viable.
    // This reduces effective complexity from O(m*n) to O(m*W) where W is band width.
    const NEG_INF: i32 = i32::MIN / 2;
    let rows = qlen + 1;
    let cols = slen + 1;

    // Precompute query score profile
    let profile: Vec<[i32; 28]> = query.iter().map(|&q| {
        let mut row = [matrix.min_score; 28];
        for r in 0u8..28 {
            row[r as usize] = matrix.score(q, r);
        }
        row
    }).collect();

    // Use two-row rolling storage for h, e, f (saves memory for long sequences).
    // Traceback still needs full matrix.
    let mut h_prev = vec![NEG_INF; cols];
    let mut h_curr = vec![NEG_INF; cols];
    let mut e_curr = vec![NEG_INF; cols];
    let mut f_prev = vec![NEG_INF; cols];
    let mut f_curr = vec![NEG_INF; cols];
    let mut tb = vec![0u8; rows * cols];

    h_prev[0] = 0;

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    // Active band: [j_lo, j_hi) — only iterate these columns per row
    let mut j_lo: usize = 1;
    let mut j_hi: usize = cols; // start unbounded, will shrink

    for i in 1..rows {
        let q_profile = &profile[i - 1];

        // Reset current row
        h_curr.fill(NEG_INF);
        e_curr.fill(NEG_INF);
        f_curr.fill(NEG_INF);

        let mut row_has_valid = false;
        let mut new_j_lo = j_hi; // will be set to first valid j
        let mut new_j_hi = j_lo; // will be set to last valid j + 1

        // Extend j_lo leftward by 1 to allow gaps from diagonal, but not below 1
        let j_start = if j_lo > 1 { j_lo - 1 } else { 1 };
        // Extend j_hi rightward by 1 to allow gap opening
        let j_end = (j_hi + 1).min(cols);

        for j in j_start..j_end {
            let diag = h_prev[j - 1];
            let s = q_profile[subject[j - 1] as usize % 28];
            let match_score = if diag == NEG_INF { NEG_INF } else { diag + s };

            // E[i][j] = max(H[i][j-1] - gap_open - gap_extend, E[i][j-1] - gap_extend)
            let h_left = if j > 0 { h_curr[j - 1] } else { NEG_INF };
            let e_left = if j > 0 { e_curr[j - 1] } else { NEG_INF };
            let from_h_e = if h_left == NEG_INF { NEG_INF } else { h_left - gap_open - gap_extend };
            let from_e_e = if e_left == NEG_INF { NEG_INF } else { e_left - gap_extend };
            e_curr[j] = from_h_e.max(from_e_e);

            // F[i][j] = max(H[i-1][j] - gap_open - gap_extend, F[i-1][j] - gap_extend)
            let h_up = h_prev[j];
            let f_up = f_prev[j];
            let from_h_f = if h_up == NEG_INF { NEG_INF } else { h_up - gap_open - gap_extend };
            let from_f_f = if f_up == NEG_INF { NEG_INF } else { f_up - gap_extend };
            f_curr[j] = from_h_f.max(from_f_f);

            let cell_score = match_score.max(e_curr[j]).max(f_curr[j]);

            let tb_idx = i * cols + j;
            if cell_score < best_score - x_drop {
                // Pruned
                tb[tb_idx] = 0;
            } else {
                h_curr[j] = cell_score;
                if cell_score == match_score {
                    tb[tb_idx] = 1;
                } else if cell_score == e_curr[j] {
                    tb[tb_idx] = 3;
                } else {
                    tb[tb_idx] = 2;
                }

                if cell_score > best_score {
                    best_score = cell_score;
                    best_i = i;
                    best_j = j;
                }

                row_has_valid = true;
                if j < new_j_lo { new_j_lo = j; }
                if j + 1 > new_j_hi { new_j_hi = j + 1; }
            }
        }

        if !row_has_valid {
            break; // Entire row pruned — no point continuing
        }

        j_lo = new_j_lo;
        j_hi = new_j_hi;

        // Swap rows
        std::mem::swap(&mut h_prev, &mut h_curr);
        std::mem::swap(&mut f_prev, &mut f_curr);
    }

    // Traceback from (best_i, best_j)
    let mut query_aln = Vec::new();
    let mut midline = Vec::new();
    let mut subject_aln = Vec::new();

    let mut i = best_i;
    let mut j = best_j;
    while i > 0 || j > 0 {
        let idx = i * cols + j;
        match tb[idx] {
            1 => {
                let qi = query[i - 1];
                let si = subject[j - 1];
                query_aln.push(qi);
                subject_aln.push(si);
                midline.push(if qi == si { b'|' } else { b' ' });
                i -= 1;
                j -= 1;
            }
            2 => {
                query_aln.push(query[i - 1]);
                midline.push(b' ');
                subject_aln.push(b'-');
                i -= 1;
            }
            3 => {
                query_aln.push(b'-');
                midline.push(b' ');
                subject_aln.push(subject[j - 1]);
                j -= 1;
            }
            _ => break,
        }
    }

    query_aln.reverse();
    midline.reverse();
    subject_aln.reverse();

    DirectionResult {
        score: best_score,
        query_len: best_i,
        subject_len: best_j,
        query_aln,
        midline,
        subject_aln,
    }
}

/// Nucleotide ungapped extension with match/mismatch scoring.
pub fn ungapped_extend_nucleotide(
    query: &[u8],
    subject: &[u8],
    q_pos: usize,
    s_pos: usize,
    match_score: i32,
    mismatch: i32,
    x_drop: i32,
) -> UngappedHit {
    // Use a simple scoring: treat query/subject as ASCII nucleotides
    let score_fn = |a: u8, b: u8| -> i32 {
        if a.eq_ignore_ascii_case(&b) {
            match_score
        } else {
            mismatch
        }
    };

    let mut best = 0i32;
    let mut best_q_start = q_pos;
    let mut best_s_start = s_pos;

    // Extend left
    let mut score = 0i32;
    let mut qi = q_pos as isize - 1;
    let mut si = s_pos as isize - 1;
    while qi >= 0 && si >= 0 {
        score += score_fn(query[qi as usize], subject[si as usize]);
        if score > best {
            best = score;
            best_q_start = qi as usize;
            best_s_start = si as usize;
        } else if best - score > x_drop {
            break;
        }
        qi -= 1;
        si -= 1;
    }
    let left_score = best;
    let q_start = best_q_start;
    let s_start = best_s_start;

    // Extend right
    let mut score = 0i32;
    best = 0;
    let mut best_q_end = q_pos;
    let mut best_s_end = s_pos;
    let mut qi = q_pos;
    let mut si = s_pos;
    while qi < query.len() && si < subject.len() {
        score += score_fn(query[qi], subject[si]);
        if score > best {
            best = score;
            best_q_end = qi + 1;
            best_s_end = si + 1;
        } else if best - score > x_drop {
            break;
        }
        qi += 1;
        si += 1;
    }
    let right_score = best;

    let center = if q_pos < query.len() && s_pos < subject.len() {
        score_fn(query[q_pos], subject[s_pos])
    } else { 0 };

    UngappedHit {
        score: left_score + center + right_score,
        q_start,
        q_end: best_q_end,
        s_start,
        s_end: best_s_end,
    }
}
