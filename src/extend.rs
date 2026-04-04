//! Ungapped and gapped extension for BLAST.

use crate::matrix::ScoringMatrix;

// ── SIMD helpers for x86_64 ─────────────────────────────────────────────────

/// SIMD-accelerated F computation: f_curr[j] = max(h_prev[j] - goe, f_prev[j] - ge)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compute_f_avx2(
    h_prev: &[i32], f_prev: &[i32], f_curr: &mut [i32],
    gap_open_extend: i32, gap_extend: i32,
    j_start: usize, j_end: usize,
) {
    use std::arch::x86_64::*;
    let goe = _mm256_set1_epi32(gap_open_extend);
    let ge = _mm256_set1_epi32(gap_extend);

    let mut j = j_start;
    while j + 8 <= j_end {
        let h = _mm256_loadu_si256(h_prev.as_ptr().add(j) as *const _);
        let f = _mm256_loadu_si256(f_prev.as_ptr().add(j) as *const _);
        let fh = _mm256_sub_epi32(h, goe);
        let ff = _mm256_sub_epi32(f, ge);
        let result = _mm256_max_epi32(fh, ff);
        _mm256_storeu_si256(f_curr.as_mut_ptr().add(j) as *mut _, result);
        j += 8;
    }
    for j in j..j_end {
        f_curr[j] = (h_prev[j] - gap_open_extend).max(f_prev[j] - gap_extend);
    }
}

/// SIMD-accelerated diagonal score: diag[j] = h_prev[j-1] + scores[j-1]
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compute_diag_avx2(
    h_prev: &[i32], subject_scores: &[i32], diag_scores: &mut [i32],
    j_start: usize, j_end: usize,
) {
    use std::arch::x86_64::*;
    let mut j = j_start;
    while j + 8 <= j_end {
        let h = _mm256_loadu_si256(h_prev.as_ptr().add(j - 1) as *const _);
        let s = _mm256_loadu_si256(subject_scores.as_ptr().add(j - 1) as *const _);
        let result = _mm256_add_epi32(h, s);
        _mm256_storeu_si256(diag_scores.as_mut_ptr().add(j) as *mut _, result);
        j += 8;
    }
    for j in j..j_end {
        diag_scores[j] = h_prev[j - 1] + subject_scores[j - 1];
    }
}

/// Check if AVX2 is available at runtime.
#[cfg(target_arch = "x86_64")]
fn has_avx2() -> bool {
    is_x86_feature_detected!("avx2")
}

// ── Farrar striped SIMD score-only extension ────────────────────────────────

/// One-direction score-only extension using Farrar's striped SIMD approach.
/// Processes 8 query positions per SIMD instruction (AVX2 i32).
/// Column-wise DP eliminates the sequential E dependency that blocks vectorization
/// in the row-wise approach.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(dead_code)]
unsafe fn extend_score_only_striped(
    query: &[u8],
    subject: &[u8],
    matrix: &ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> i32 {
    use std::arch::x86_64::*;

    let qlen = query.len();
    let slen = subject.len();
    if qlen == 0 || slen == 0 { return 0; }

    const W: usize = 8; // AVX2 i32 lanes
    let seg = qlen.div_ceil(W); // number of segments
    let neg_inf = i32::MIN / 2;

    let v_neg_inf = _mm256_set1_epi32(neg_inf);
    let v_goe = _mm256_set1_epi32(gap_open + gap_extend);
    let v_ge = _mm256_set1_epi32(gap_extend);
    // Permutation index to shift lanes right by 1: dst[0]=7, dst[1]=0, dst[2]=1, ...
    let shift_right = _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6);

    // Build striped profile: profile[residue][s] is a __m256i
    // Lane w of segment s holds score(query[s + w*seg], residue)
    let mut profile: Vec<Vec<__m256i>> = Vec::with_capacity(28);
    for r in 0..28u8 {
        let mut segs = Vec::with_capacity(seg);
        for s in 0..seg {
            let mut scores = [neg_inf; W];
            for (w, score) in scores.iter_mut().enumerate() {
                let pos = s + w * seg;
                if pos < qlen {
                    *score = matrix.score(query[pos], r);
                }
            }
            segs.push(_mm256_loadu_si256(scores.as_ptr() as *const _));
        }
        profile.push(segs);
    }

    // DP state: H and E for current and previous columns, F for current column
    let mut h_prev = vec![v_neg_inf; seg];
    let mut h_curr = vec![v_neg_inf; seg];
    let mut e_prev = vec![v_neg_inf; seg];
    let mut e_curr = vec![v_neg_inf; seg];
    let mut f = vec![v_neg_inf; seg];

    // Boundary: H at "query position -1" for diagonal into segment 0, lane 0
    let mut boundary = 0i32; // 0 for first column, NEG_INF after

    let mut best_score = 0i32;

    for &s_byte in subject.iter().take(slen) {
        let s_res = (s_byte as usize) % 28;
        let prof = &profile[s_res];

        // Step 1: E for all segments (from previous column, fully parallel)
        for s in 0..seg {
            e_curr[s] = _mm256_max_epi32(
                _mm256_sub_epi32(h_prev[s], v_goe),
                _mm256_sub_epi32(e_prev[s], v_ge),
            );
        }

        // Step 2: Diagonal + initial H = max(diag + score, E)
        // Segment 0: diagonal comes from h_prev[seg-1] shifted right, with boundary at lane 0
        {
            let last = h_prev[seg - 1];
            let shifted = _mm256_permutevar8x32_epi32(last, shift_right);
            let bnd = _mm256_insert_epi32::<0>(v_neg_inf, boundary);
            let diag = _mm256_blend_epi32::<0x01>(shifted, bnd);
            let m = _mm256_add_epi32(diag, prof[0]);
            h_curr[0] = _mm256_max_epi32(m, e_curr[0]);
        }
        // Segments 1..seg-1: diagonal = h_prev[s-1] (same lanes)
        for s in 1..seg {
            let m = _mm256_add_epi32(h_prev[s - 1], prof[s]);
            h_curr[s] = _mm256_max_epi32(m, e_curr[s]);
        }

        // Step 3: F propagation (Farrar's lazy-F)
        // For segments 1..seg-1: F[s] depends on H[s-1] and F[s-1] at SAME lanes
        // For segment 0: F depends on H[seg-1] and F[seg-1] SHIFTED right by 1 lane
        //
        // First pass: propagate F through segments 1..seg-1
        f[0] = v_neg_inf; // Will be corrected by carry from segment seg-1
        for s in 1..seg {
            f[s] = _mm256_max_epi32(
                _mm256_sub_epi32(h_curr[s - 1], v_goe),
                _mm256_sub_epi32(f[s - 1], v_ge),
            );
            h_curr[s] = _mm256_max_epi32(h_curr[s], f[s]);
        }

        // Carry from segment seg-1 to segment 0 (shifted)
        let carry_h = _mm256_permutevar8x32_epi32(h_curr[seg - 1], shift_right);
        let carry_f = _mm256_permutevar8x32_epi32(f[seg - 1], shift_right);
        // Lane 0 = boundary (NEG_INF for F)
        let carry_h = _mm256_blend_epi32::<0x01>(carry_h, v_neg_inf);
        let carry_f = _mm256_blend_epi32::<0x01>(carry_f, v_neg_inf);

        let f0 = _mm256_max_epi32(
            _mm256_sub_epi32(carry_h, v_goe),
            _mm256_sub_epi32(carry_f, v_ge),
        );

        // Check if F improves segment 0
        let improved = _mm256_cmpgt_epi32(f0, h_curr[0]);
        if _mm256_movemask_epi8(improved) != 0 {
            f[0] = f0;
            h_curr[0] = _mm256_max_epi32(h_curr[0], f0);

            // Lazy correction: re-propagate through segments 1..seg-1
            // In practice converges in 0-1 iterations
            for _iter in 0..2 {
                let mut changed = false;
                for s in 1..seg {
                    let new_f = _mm256_max_epi32(
                        _mm256_sub_epi32(h_curr[s - 1], v_goe),
                        _mm256_sub_epi32(f[s - 1], v_ge),
                    );
                    let better = _mm256_cmpgt_epi32(new_f, h_curr[s]);
                    if _mm256_movemask_epi8(better) != 0 {
                        f[s] = _mm256_max_epi32(f[s], new_f);
                        h_curr[s] = _mm256_max_epi32(h_curr[s], new_f);
                        changed = true;
                    }
                }
                if !changed { break; }
            }
        }

        // Step 4: X-drop check — find max H across all segments and lanes
        let mut col_max_v = v_neg_inf;
        for h in &h_curr {
            col_max_v = _mm256_max_epi32(col_max_v, *h);
        }
        // Horizontal max: reduce 8 lanes to scalar
        let hi = _mm256_extracti128_si256::<1>(col_max_v);
        let lo = _mm256_castsi256_si128(col_max_v);
        let m128 = _mm_max_epi32(lo, hi);
        let m64 = _mm_max_epi32(m128, _mm_shuffle_epi32(m128, 0x4E));
        let m32 = _mm_max_epi32(m64, _mm_shuffle_epi32(m64, 0xB1));
        let col_max = _mm_extract_epi32::<0>(m32);

        if col_max > best_score { best_score = col_max; }
        if col_max < best_score - x_drop { break; }

        // Swap columns
        std::mem::swap(&mut h_prev, &mut h_curr);
        std::mem::swap(&mut e_prev, &mut e_curr);
        boundary = neg_inf; // After first column, no boundary contribution

        // Reset current column
        h_curr.fill(v_neg_inf);
        e_curr.fill(v_neg_inf);
    }

    best_score
}

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
    ungapped_extend_profile(query, subject, q_pos, s_pos, matrix, x_drop)
}

/// Profile-based ungapped extension: precomputes a score lookup per query residue
/// to avoid repeated matrix.score() calls with bounds checking.
fn ungapped_extend_profile(
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
        let q_byte = query[qi as usize];
        let s_byte = subject[si as usize];
        score += matrix.scores[q_byte as usize % 28][s_byte as usize % 28];
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
        let q_byte = query[qi];
        let s_byte = subject[si];
        score += matrix.scores[q_byte as usize % 28][s_byte as usize % 28];
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

/// Score-only gapped extension (no traceback, no alignment strings).
/// Used for preliminary gapped extension to quickly reject low-scoring hits.
#[allow(clippy::too_many_arguments)]
pub fn gapped_extend_score_only(
    query: &[u8],
    subject: &[u8],
    q_center: usize,
    s_center: usize,
    matrix: &ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> i32 {
    let q_right = &query[q_center..];
    let s_right = &subject[s_center..];
    let query_rev: Vec<u8> = query[..q_center].iter().rev().cloned().collect();
    let subject_rev: Vec<u8> = subject[..s_center].iter().rev().cloned().collect();

    // Farrar's striped SIMD is available but only beneficial when alignments are long
    // and the X-drop band is wide. For typical BLAST searches with short local alignments,
    // the banded scalar approach below is faster because it only processes the narrow
    // active band rather than all query positions.
    // The striped implementation (extend_score_only_striped) is retained for future use
    // with real biological databases where long homologous alignments are common.

    let right = extend_score_only(q_right, s_right, matrix, gap_open, gap_extend, x_drop);
    let left = extend_score_only(&query_rev, &subject_rev, matrix, gap_open, gap_extend, x_drop);
    left + right
}

/// Score-only one-direction extension. No traceback array, no alignment strings.
/// Uses branchless arithmetic and separated loops to enable auto-vectorization.
fn extend_score_only(
    query: &[u8],
    subject: &[u8],
    matrix: &ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> i32 {
    let qlen = query.len();
    let slen = subject.len();
    if qlen == 0 || slen == 0 { return 0; }

    const NEG_INF: i32 = i32::MIN / 2;
    let cols = slen + 1;
    let gap_open_extend = gap_open + gap_extend;

    let profile: Vec<[i32; 28]> = query.iter().map(|&q| {
        let mut row = [matrix.min_score; 28];
        for r in 0u8..28 { row[r as usize] = matrix.score(q, r); }
        row
    }).collect();

    // Precompute subject residue indices to avoid repeated modulo in inner loop
    let subj_idx: Vec<usize> = subject.iter().map(|&b| (b as usize) % 28).collect();

    let mut h_prev = vec![NEG_INF; cols];
    let mut h_curr = vec![NEG_INF; cols];
    let mut e_curr = vec![NEG_INF; cols];
    let mut f_prev = vec![NEG_INF; cols];
    let mut f_curr = vec![NEG_INF; cols];
    let mut diag_scores = vec![NEG_INF; cols];
    let mut subject_scores = vec![0i32; slen];

    #[cfg(target_arch = "x86_64")]
    let use_avx2 = has_avx2();

    h_prev[0] = 0;
    let mut best_score = 0i32;
    let mut j_lo: usize = 1;
    let mut j_hi: usize = cols;

    for i in 1..=qlen {
        let q_profile = &profile[i - 1];

        let j_start = if j_lo > 1 { j_lo - 1 } else { 1 };
        let j_end = (j_hi + 1).min(cols);

        // Reset current row
        h_curr[..j_end].fill(NEG_INF);
        e_curr[..j_end].fill(NEG_INF);

        // Precompute subject scores for this query row
        for j in j_start..j_end {
            subject_scores[j - 1] = q_profile[subj_idx[j - 1]];
        }

        #[cfg(target_arch = "x86_64")]
        if use_avx2 {
            unsafe {
                compute_f_avx2(&h_prev, &f_prev, &mut f_curr, gap_open_extend, gap_extend, j_start, j_end);
                compute_diag_avx2(&h_prev, &subject_scores, &mut diag_scores, j_start, j_end);
            }
        } else {
            for j in j_start..j_end {
                f_curr[j] = (h_prev[j] - gap_open_extend).max(f_prev[j] - gap_extend);
            }
            for j in j_start..j_end {
                diag_scores[j] = h_prev[j - 1] + subject_scores[j - 1];
            }
        }

        #[cfg(not(target_arch = "x86_64"))]
        {
            for j in j_start..j_end {
                f_curr[j] = (h_prev[j] - gap_open_extend).max(f_prev[j] - gap_extend);
            }
            for j in j_start..j_end {
                diag_scores[j] = h_prev[j - 1] + subject_scores[j - 1];
            }
        }

        // Sequential E sweep (left-to-right dependency)
        for j in j_start..j_end {
            e_curr[j] = (h_curr[j - 1] - gap_open_extend).max(e_curr[j - 1] - gap_extend);
        }

        // Combine: cell_score = max(diag, e, f) with X-drop pruning
        let mut row_has_valid = false;
        let mut new_j_lo = j_end;
        let mut new_j_hi = j_start;
        let drop_threshold = best_score - x_drop;

        for j in j_start..j_end {
            let cell_score = diag_scores[j].max(e_curr[j]).max(f_curr[j]);

            if cell_score >= drop_threshold {
                h_curr[j] = cell_score;
                if cell_score > best_score { best_score = cell_score; }
                row_has_valid = true;
                if j < new_j_lo { new_j_lo = j; }
                new_j_hi = j + 1;
            } else {
                h_curr[j] = NEG_INF;
            }
        }

        if !row_has_valid { break; }
        j_lo = new_j_lo;
        j_hi = new_j_hi;

        std::mem::swap(&mut h_prev, &mut h_curr);
        std::mem::swap(&mut f_prev, &mut f_curr);
    }

    best_score
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
    // The inner loop is split into vectorizable (F, diagonal scores) and sequential (E, H) phases.
    // Branchless arithmetic enables LLVM auto-vectorization of the independent loops.
    const NEG_INF: i32 = i32::MIN / 2;
    let rows = qlen + 1;
    let cols = slen + 1;
    let gap_open_extend = gap_open + gap_extend;

    let profile: Vec<[i32; 28]> = query.iter().map(|&q| {
        let mut row = [matrix.min_score; 28];
        for r in 0u8..28 {
            row[r as usize] = matrix.score(q, r);
        }
        row
    }).collect();

    let subj_idx: Vec<usize> = subject.iter().map(|&b| (b as usize) % 28).collect();

    let mut h_prev = vec![NEG_INF; cols];
    let mut h_curr = vec![NEG_INF; cols];
    let mut e_curr = vec![NEG_INF; cols];
    let mut f_prev = vec![NEG_INF; cols];
    let mut f_curr = vec![NEG_INF; cols];
    let mut diag_scores = vec![NEG_INF; cols];
    let mut subject_scores = vec![0i32; slen];
    let mut tb = vec![0u8; rows * cols];

    #[cfg(target_arch = "x86_64")]
    let use_avx2 = has_avx2();

    h_prev[0] = 0;

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    let mut j_lo: usize = 1;
    let mut j_hi: usize = cols;

    for i in 1..rows {
        let q_profile = &profile[i - 1];

        let j_start = if j_lo > 1 { j_lo - 1 } else { 1 };
        let j_end = (j_hi + 1).min(cols);

        h_curr[..j_end].fill(NEG_INF);
        e_curr[..j_end].fill(NEG_INF);

        for j in j_start..j_end {
            subject_scores[j - 1] = q_profile[subj_idx[j - 1]];
        }

        #[cfg(target_arch = "x86_64")]
        if use_avx2 {
            unsafe {
                compute_f_avx2(&h_prev, &f_prev, &mut f_curr, gap_open_extend, gap_extend, j_start, j_end);
                compute_diag_avx2(&h_prev, &subject_scores, &mut diag_scores, j_start, j_end);
            }
        } else {
            for j in j_start..j_end {
                f_curr[j] = (h_prev[j] - gap_open_extend).max(f_prev[j] - gap_extend);
            }
            for j in j_start..j_end {
                diag_scores[j] = h_prev[j - 1] + subject_scores[j - 1];
            }
        }

        #[cfg(not(target_arch = "x86_64"))]
        {
            for j in j_start..j_end {
                f_curr[j] = (h_prev[j] - gap_open_extend).max(f_prev[j] - gap_extend);
            }
            for j in j_start..j_end {
                diag_scores[j] = h_prev[j - 1] + subject_scores[j - 1];
            }
        }

        // Step 3: Sequential left-to-right sweep for E, H, and traceback
        let mut row_has_valid = false;
        let mut new_j_lo = j_end;
        let mut new_j_hi = j_start;
        let drop_threshold = best_score - x_drop;

        for j in j_start..j_end {
            e_curr[j] = (h_curr[j - 1] - gap_open_extend).max(e_curr[j - 1] - gap_extend);

            let cell_score = diag_scores[j].max(e_curr[j]).max(f_curr[j]);

            let tb_idx = i * cols + j;
            if cell_score >= drop_threshold {
                h_curr[j] = cell_score;
                if cell_score == diag_scores[j] {
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
                new_j_hi = j + 1;
            } else {
                h_curr[j] = NEG_INF;
                tb[tb_idx] = 0;
            }
        }

        if !row_has_valid {
            break;
        }

        j_lo = new_j_lo;
        j_hi = new_j_hi;

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
