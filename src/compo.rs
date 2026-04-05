//! Composition-based statistics adjustment.
//!
//! Adjusts E-values to account for amino acid composition bias between
//! query and subject sequences, using the method of Schaffer et al. (2001).
//!
//! The core idea: find λ' such that
//!     `Σᵢⱼ qᵢ · rⱼ · exp(λ' · s(i,j)) = 1`
//! where qᵢ and rⱼ are the per-sequence amino acid frequencies and s is the
//! scoring matrix. Then use λ' instead of λ for the E-value computation.

use crate::matrix::ScoringMatrix;

/// BLOSUM62 background amino acid frequencies (Henikoff & Henikoff 1992).
/// Indexed by Ncbistdaa code (0–27). Codes without a standard AA get 0.0.
pub static BACKGROUND_FREQ: [f64; 28] = [
    0.000, // 0  -  gap
    0.074, // 1  A
    0.050, // 2  B  (average N/D)
    0.025, // 3  C
    0.054, // 4  D
    0.054, // 5  E
    0.047, // 6  F
    0.074, // 7  G
    0.026, // 8  H
    0.068, // 9  I
    0.058, // 10 K
    0.099, // 11 L
    0.025, // 12 M
    0.045, // 13 N
    0.039, // 14 P
    0.034, // 15 Q
    0.052, // 16 R
    0.057, // 17 S
    0.051, // 18 T
    0.073, // 19 V
    0.013, // 20 W
    0.050, // 21 X  (unknown – use mean)
    0.032, // 22 Y
    0.050, // 23 Z  (average Q/E)
    0.000, // 24 U
    0.000, // 25 *
    0.000, // 26 O
    0.000, // 27 J
];

/// Compute the amino acid composition (as frequencies) of a Ncbistdaa-encoded sequence.
/// Only codes 1–22 (standard amino acids) contribute.
pub fn composition_ncbistdaa(seq: &[u8]) -> [f64; 28] {
    let mut counts = [0u64; 28];
    let mut total = 0u64;
    for &c in seq {
        let idx = c as usize;
        if (1..=22).contains(&idx) {
            counts[idx] += 1;
            total += 1;
        }
    }
    let mut freq = [0.0f64; 28];
    if total > 0 {
        let t = total as f64;
        for i in 0..28 { freq[i] = counts[i] as f64 / t; }
    } else {
        freq.copy_from_slice(&BACKGROUND_FREQ);
    }
    freq
}

/// Compute the amino acid composition of an ASCII protein sequence.
pub fn composition_ascii(seq: &[u8]) -> [f64; 28] {
    // Map ASCII → Ncbistdaa, then delegate
    let ncbi: Vec<u8> = seq.iter().map(|&c| ascii_to_ncbistdaa(c)).collect();
    composition_ncbistdaa(&ncbi)
}

fn ascii_to_ncbistdaa(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 1,  b'B' => 2,  b'C' => 3,  b'D' => 4,  b'E' => 5,
        b'F' => 6,  b'G' => 7,  b'H' => 8,  b'I' => 9,  b'K' => 10,
        b'L' => 11, b'M' => 12, b'N' => 13, b'P' => 14, b'Q' => 15,
        b'R' => 16, b'S' => 17, b'T' => 18, b'V' => 19, b'W' => 20,
        b'X' => 21, b'Y' => 22, b'Z' => 23, _ => 21,
    }
}

/// Compute the expected score Σᵢⱼ pᵢ·qⱼ·s[i][j] under background frequencies.
#[allow(clippy::needless_range_loop)]
fn expected_score(pq: &[f64; 28], rq: &[f64; 28], matrix: &ScoringMatrix) -> f64 {
    let mut mu = 0.0f64;
    for i in 1..23usize {
        for j in 1..23usize {
            mu += pq[i] * rq[j] * matrix.score(i as u8, j as u8) as f64;
        }
    }
    mu
}

/// Find λ' via bisection such that `Σᵢⱼ pᵢ·rⱼ·exp(λ'·s(i,j)) = 1`.
///
/// Returns `None` if no positive λ' exists (e.g., the expected score is non-negative,
/// meaning composition correction is not applicable).
#[allow(clippy::needless_range_loop)]
pub fn find_adjusted_lambda(
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
) -> Option<f64> {
    // Guard: if expected score >= 0, correction is not meaningful.
    if expected_score(q, r, matrix) >= 0.0 { return None; }

    // Σᵢⱼ pᵢ·rⱼ·exp(λ·s[i][j]) should equal 1 at λ=0 for a valid matrix.
    // We search in [0, 2*lambda_standard].
    let eval_sum = |lam: f64| -> f64 {
        let mut sum = 0.0f64;
        for i in 1..23usize {
            for j in 1..23usize {
                let s = matrix.score(i as u8, j as u8) as f64;
                sum += q[i] * r[j] * (lam * s).exp();
            }
        }
        sum
    };

    let lo_val = eval_sum(0.0);
    if lo_val < 1.0 {
        // No positive solution; return None to skip adjustment.
        return None;
    }

    let mut lo = 0.0f64;
    let mut hi = lambda_standard * 4.0;

    // Verify hi bracket gives sum < 1
    if eval_sum(hi) >= 1.0 {
        return None; // Can't bracket — very unusual; skip adjustment
    }

    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if eval_sum(mid) > 1.0 { lo = mid; } else { hi = mid; }
        if hi - lo < 1e-10 { break; }
    }

    let lambda_prime = (lo + hi) / 2.0;
    Some(lambda_prime)
}

/// Adjust an E-value using per-sequence composition correction (mode 1).
///
/// Returns the adjusted E-value, or the original if correction is inapplicable.
#[allow(clippy::too_many_arguments)]
pub fn adjust_evalue(
    raw_evalue: f64,
    score: i32,
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
    k: f64,
    eff_query_len: usize,
    eff_db_len: u64,
) -> f64 {
    match find_adjusted_lambda(q, r, matrix, lambda_standard) {
        None => raw_evalue,
        Some(lambda_prime) => {
            (eff_query_len as f64) * (eff_db_len as f64) * k
                * (-(lambda_prime * score as f64)).exp()
        }
    }
}

/// Adjust E-value with mode selection:
///   0 = no adjustment (returns raw_evalue)
///   1 = unconditional composition-based λ adjustment (original method)
///   2 = conditional: only apply if expected score diverges significantly from standard
///   3 = unconditional adjustment (same as 1 but always applied, even if expected_score ≥ 0)
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub fn adjust_evalue_with_mode(
    raw_evalue: f64,
    score: i32,
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
    k: f64,
    eff_query_len: usize,
    eff_db_len: u64,
    mode: u8,
) -> f64 {
    match mode {
        0 => raw_evalue,
        1 => adjust_evalue(raw_evalue, score, q, r, matrix, lambda_standard, k, eff_query_len, eff_db_len),
        2 => {
            // Conditional: only apply when composition differs significantly
            // from background. Threshold: expected score < -0.2 * lambda_standard
            let mu_bg = expected_score_with_bg(matrix, lambda_standard);
            let mu_actual = expected_score(q, r, matrix);
            let threshold = -0.2 * lambda_standard;
            if (mu_actual - mu_bg).abs() > threshold.abs() {
                adjust_evalue(raw_evalue, score, q, r, matrix, lambda_standard, k, eff_query_len, eff_db_len)
            } else {
                raw_evalue
            }
        }
        3 => {
            // Unconditional: force find adjusted lambda even if expected score ≥ 0
            let eval_sum = |lam: f64| -> f64 {
                let mut sum = 0.0f64;
                for i in 1..23usize {
                    for j in 1..23usize {
                        let s = matrix.score(i as u8, j as u8) as f64;
                        sum += q[i] * r[j] * (lam * s).exp();
                    }
                }
                sum
            };

            let mut lo = 0.0f64;
            let mut hi = lambda_standard * 4.0;
            if eval_sum(hi) >= 1.0 { return raw_evalue; }
            if eval_sum(lo) < 1.0 { return raw_evalue; }

            for _ in 0..60 {
                let mid = (lo + hi) / 2.0;
                if eval_sum(mid) > 1.0 { lo = mid; } else { hi = mid; }
                if hi - lo < 1e-10 { break; }
            }
            let lambda_prime = (lo + hi) / 2.0;
            (eff_query_len as f64) * (eff_db_len as f64) * k
                * (-(lambda_prime * score as f64)).exp()
        }
        _ => raw_evalue,
    }
}

/// Expected score under background frequencies (for mode 2 comparison).
#[allow(clippy::needless_range_loop)]
fn expected_score_with_bg(matrix: &ScoringMatrix, _lambda: f64) -> f64 {
    let mut mu = 0.0f64;
    for i in 1..23usize {
        for j in 1..23usize {
            mu += BACKGROUND_FREQ[i] * BACKGROUND_FREQ[j] * matrix.score(i as u8, j as u8) as f64;
        }
    }
    mu
}
