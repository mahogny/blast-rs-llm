//! Position-Specific Scoring Matrix (PSSM) and PSI-BLAST iteration.
//!
//! ## PSSM construction
//!
//! Given a set of aligned sequences from a previous BLAST round:
//!   1. For each query position, tally residue frequencies from the alignment.
//!   2. Apply pseudocounts using BLOSUM62 background frequencies.
//!   3. Convert to log-odds scores scaled by 1/λ.
//!
//! ## PSI-BLAST
//!
//! Run `blast_search` for the first iteration, then iteratively rebuild the
//! PSSM from hits below the inclusion threshold and re-search with PSSM scoring.

use std::collections::HashSet;
use crate::matrix::{ScoringMatrix, MatrixType};
use crate::stats::{KarlinAltschul, GapPenalty, lookup_ka_params};
use crate::lookup::ProteinLookup;
use crate::extend::{ungapped_extend, gapped_extend};
use crate::hsp::{Hsp, SearchResult};
use crate::compo::BACKGROUND_FREQ;
use crate::search::{SearchParams, neighbor_threshold};
use blast_db::BlastDb;
use rayon::prelude::*;

// ─── PSSM data structure ─────────────────────────────────────────────────────

/// A position-specific scoring matrix for a protein query.
///
/// `scores[i][j]` is the score for Ncbistdaa code `j` at query position `i`.
#[derive(Debug, Clone)]
pub struct Pssm {
    pub query_len: usize,
    /// scores[position][ncbistdaa_code 0..28]
    pub scores: Vec<[i32; 28]>,
    /// Lambda used for scaling (BLOSUM62 ungapped λ ≈ 0.3176).
    pub lambda: f64,
}

impl Pssm {
    /// Build an identity PSSM from a query sequence (equivalent to the scoring matrix).
    pub fn from_query(query: &[u8], matrix: &ScoringMatrix, lambda: f64) -> Self {
        let scores: Vec<[i32; 28]> = query.iter().map(|&qa| {
            let mut row = [0i32; 28];
            for j in 0..28usize {
                row[j] = matrix.score(qa, j as u8);
            }
            row
        }).collect();
        Pssm { query_len: query.len(), scores, lambda }
    }

    /// Look up the score for residue `aa` (Ncbistdaa) at query position `pos`.
    #[inline]
    pub fn score(&self, pos: usize, aa: u8) -> i32 {
        if pos < self.scores.len() && (aa as usize) < 28 {
            self.scores[pos][aa as usize]
        } else {
            0
        }
    }

    /// Minimum score across all positions (used for ungapped extension X-drop).
    pub fn min_score(&self) -> i32 {
        self.scores.iter()
            .flat_map(|row| row.iter().copied())
            .min()
            .unwrap_or(-4)
    }
}

// ─── PSSM construction from alignments ───────────────────────────────────────

/// Pseudocount weight (β). NCBI BLAST uses ~10.
const PSEUDOCOUNT: f64 = 10.0;

/// Build a PSSM from a set of search results.
///
/// Only HSPs with E-value ≤ `inclusion_evalue` contribute.
pub fn build_pssm(
    query: &[u8],
    results: &[SearchResult],
    inclusion_evalue: f64,
    matrix: MatrixType,
    lambda: f64,
) -> Pssm {
    let query_len = query.len();
    // freq_counts[pos][ncbistdaa_code] = weighted count of observations
    let mut freq_counts: Vec<[f64; 28]> = vec![[0.0f64; 28]; query_len];
    let mut obs_total:   Vec<f64>       = vec![0.0f64;       query_len];

    for result in results {
        for hsp in &result.hsps {
            if hsp.evalue > inclusion_evalue { continue; }

            // Walk the gapped alignment and tally residues at each query position.
            let mut q_pos = hsp.query_start;
            let q_aln = &hsp.query_aln;
            let s_aln = &hsp.subject_aln;

            for (qa, sa) in q_aln.iter().zip(s_aln.iter()) {
                if *qa == b'-' { continue; } // gap in query → skip
                if q_pos >= query_len { break; }
                if *sa != b'-' {
                    // Subject residue at this query position
                    let ncbi = ascii_to_ncbistdaa(*sa);
                    if ncbi >= 1 && ncbi <= 22 {
                        freq_counts[q_pos][ncbi as usize] += 1.0;
                        obs_total[q_pos] += 1.0;
                    }
                }
                q_pos += 1;
            }
        }
    }

    // Always include the query itself
    for (i, &qa) in query.iter().enumerate() {
        let ncbi = qa as usize; // query is already Ncbistdaa
        if ncbi < 28 {
            freq_counts[i][ncbi] += 1.0;
            obs_total[i] += 1.0;
        }
    }

    // Build PSSM scores with pseudocounts
    let scale = 1.0 / lambda; // bits per nat
    let mat = ScoringMatrix::from_type(matrix);

    let scores: Vec<[i32; 28]> = (0..query_len).map(|i| {
        let total = obs_total[i].max(1.0);
        let mut row = [0i32; 28];
        for j in 0..28usize {
            let bg = BACKGROUND_FREQ[j];
            if bg <= 0.0 { continue; }
            // Pseudocount-adjusted frequency
            let f_adj = (freq_counts[i][j] + PSEUDOCOUNT * bg) / (total + PSEUDOCOUNT);
            // Log-odds scaled by 1/λ
            let log_odds = (f_adj / bg).ln() * scale;
            row[j] = log_odds.round() as i32;
        }
        // For codes without background (gap, *, …), use the matrix value against query aa
        let qa = query[i];
        for j in [0, 25, 26, 27usize] {
            row[j] = mat.score(qa, j as u8);
        }
        row
    }).collect();

    Pssm { query_len, scores, lambda }
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

// ─── PSSM-based search ───────────────────────────────────────────────────────

/// Run a single BLAST iteration using a PSSM instead of a substitution matrix.
pub fn search_with_pssm(
    db: &BlastDb,
    query: &[u8], // Ncbistdaa
    pssm: &Pssm,
    params: &SearchParams,
) -> Vec<SearchResult> {
    let mat = ScoringMatrix::from_type(params.matrix);
    let gap = GapPenalty::new(params.gap_open, params.gap_extend);

    let ka = match lookup_ka_params(params.matrix, gap) {
        Some(k) => k,
        None => lookup_ka_params(MatrixType::Blosum62, GapPenalty::blosum62_default()).unwrap(),
    };

    let db_len = db.volume_length();
    let num_seqs = db.num_sequences() as u64;
    let (eff_q, eff_db) = ka.effective_lengths(query.len(), db_len, num_seqs);

    // Build the neighbor lookup using the standard matrix (not the PSSM)
    // so that we seed on common residue types; PSSM scoring is used during extension.
    let threshold = neighbor_threshold(params.matrix, params.word_size);
    let lookup = ProteinLookup::build(query, params.word_size, &mat, threshold);

    let oids: Vec<u32> = (0..db.num_sequences()).collect();

    let mut results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
        let subject = match db.get_sequence_protein_raw(oid) {
            Ok(s) => s.to_vec(),
            Err(_) => return None,
        };
        if subject.is_empty() { return None; }

        let hsps = search_one_pssm(query, &subject, pssm, &lookup, &ka, params, eff_q, eff_db);
        if hsps.is_empty() { return None; }

        let header = db.get_header(oid).unwrap_or_default();
        Some(SearchResult {
            subject_oid: oid,
            subject_title: header.title,
            subject_accession: header.accession,
            subject_len: subject.len(),
            hsps,
        })
    }).collect();

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    results.truncate(params.max_target_seqs);
    results
}

fn search_one_pssm(
    query: &[u8],
    subject: &[u8],
    pssm: &Pssm,
    lookup: &ProteinLookup,
    ka: &KarlinAltschul,
    params: &SearchParams,
    eff_q: usize,
    eff_db: u64,
) -> Vec<Hsp> {
    let slen = subject.len();
    let ws = lookup.word_size;
    if slen < ws { return vec![]; }

    let mat_for_extend = pssm_to_scoring_matrix(pssm, query);

    let diag_offset = query.len();
    let mut diag_hit = vec![false; query.len() + slen + 1];
    let mut ungapped_hits = Vec::new();

    for s_pos in 0..=(slen - ws) {
        let word = &subject[s_pos..s_pos + ws];
        if let Some(q_positions) = lookup.lookup(word) {
            for &q_pos in q_positions {
                let q_pos = q_pos as usize;
                let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
                if diag < diag_hit.len() && diag_hit[diag] { continue; }

                let hit = ungapped_extend(query, subject, q_pos, s_pos, &mat_for_extend, params.x_drop_ungapped);
                let cutoff = params.ungapped_cutoff.max(1);
                if hit.score >= cutoff {
                    if diag < diag_hit.len() { diag_hit[diag] = true; }
                    ungapped_hits.push(hit);
                }
            }
        }
    }

    if ungapped_hits.is_empty() { return vec![]; }
    ungapped_hits.sort_by(|a, b| b.score.cmp(&a.score));

    let mut hsps = Vec::new();
    let mut covered = vec![false; query.len()];

    for uh in ungapped_hits {
        let center_q = (uh.q_start + uh.q_end) / 2;
        if center_q < covered.len() && covered[center_q] { continue; }
        let center_s = (uh.s_start + uh.s_end) / 2;

        let gh = gapped_extend(query, subject, center_q, center_s, &mat_for_extend,
                               params.gap_open, params.gap_extend, params.x_drop_gapped);
        if gh.score <= 0 { continue; }

        let evalue = ka.evalue(gh.score, eff_q, eff_db);
        if evalue > params.evalue_threshold { continue; }

        for i in gh.q_start.min(covered.len())..gh.q_end.min(covered.len()) { covered[i] = true; }

        let query_aln = blast_db::sequence::decode_protein(&gh.query_aln);
        let subject_aln = blast_db::sequence::decode_protein(&gh.subject_aln);

        hsps.push(Hsp {
            score: gh.score, bit_score: ka.bit_score(gh.score), evalue,
            query_start: gh.q_start, query_end: gh.q_end,
            subject_start: gh.s_start, subject_end: gh.s_end,
            num_identities: gh.num_identities, num_gaps: gh.num_gaps,
            alignment_length: gh.query_aln.len(),
            query_aln, midline: gh.midline, subject_aln,
            query_frame: 0, subject_frame: 0,
        });
    }

    hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    hsps
}

/// Build a ScoringMatrix from PSSM scores for use in the existing extend functions.
/// The matrix encodes query position in the row dimension by using the query residue
/// as a proxy; this is an approximation — the full PSSM search ideally needs a
/// position-aware extension kernel. For now we build an "average" PSSM row matrix.
///
/// For a more correct implementation, extend.rs would need to accept a PSSM directly.
/// As a practical approximation: use column averages.
fn pssm_to_scoring_matrix(pssm: &Pssm, query: &[u8]) -> ScoringMatrix {
    // Build a per-query-residue average: for each query aa qa, average PSSM scores across
    // all positions where query[i] == qa.
    let mut sums = [[0i64; 28]; 28];
    let mut cnts = [0u64; 28];
    for (i, &qa) in query.iter().enumerate() {
        let qa = qa as usize;
        if qa < 28 {
            cnts[qa] += 1;
            for j in 0..28usize {
                sums[qa][j] += pssm.scores[i][j] as i64;
            }
        }
    }
    let mut scores = [[0i32; 28]; 28];
    for i in 0..28 {
        for j in 0..28 {
            scores[i][j] = if cnts[i] > 0 { (sums[i][j] / cnts[i] as i64) as i32 } else { 0 };
        }
    }
    let min_score = scores.iter().flat_map(|r| r.iter().copied()).min().unwrap_or(-4);
    ScoringMatrix { name: MatrixType::Blosum62, scores, min_score }
}

// ─── PSI-BLAST iteration ─────────────────────────────────────────────────────

/// Run PSI-BLAST for up to `num_iterations` rounds.
///
/// Returns the final hit list and the PSSM from the last iteration.
/// `inclusion_evalue` controls which hits are incorporated into the PSSM (default: 0.001).
pub fn psiblast_search(
    db: &BlastDb,
    query: &[u8], // Ncbistdaa
    params: &SearchParams,
    num_iterations: u32,
    inclusion_evalue: f64,
) -> (Vec<SearchResult>, Pssm) {
    use crate::search::blast_search;

    let mat = ScoringMatrix::from_type(params.matrix);
    let gap = GapPenalty::new(params.gap_open, params.gap_extend);
    let lambda = lookup_ka_params(params.matrix, gap)
        .map(|ka| ka.lambda)
        .unwrap_or(0.267);

    // Iteration 0: standard blastp
    let mut results = blast_search(db, query, params);
    let mut pssm = build_pssm(query, &results, inclusion_evalue, params.matrix, lambda);

    let mut prev_hit_set: HashSet<u32> = results.iter().map(|r| r.subject_oid).collect();

    for _iter in 1..num_iterations {
        let new_results = search_with_pssm(db, query, &pssm, params);
        let new_hit_set: HashSet<u32> = new_results.iter().map(|r| r.subject_oid).collect();

        // Converged?
        if new_hit_set == prev_hit_set { break; }
        prev_hit_set = new_hit_set;

        pssm = build_pssm(query, &new_results, inclusion_evalue, params.matrix, lambda);
        results = new_results;
    }

    (results, pssm)
}

// ─── PSSM checkpoint save/load ──────────────────────────────────────────────

/// Ncbistdaa code to ASCII
fn ncbistdaa_to_ascii(code: u8) -> char {
    match code {
        1 => 'A', 2 => 'B', 3 => 'C', 4 => 'D', 5 => 'E', 6 => 'F', 7 => 'G',
        8 => 'H', 9 => 'I', 10 => 'K', 11 => 'L', 12 => 'M', 13 => 'N', 14 => 'P',
        15 => 'Q', 16 => 'R', 17 => 'S', 18 => 'T', 19 => 'V', 20 => 'W',
        21 => 'X', 22 => 'Y', 23 => 'Z', 25 => '*', _ => 'X',
    }
}

impl Pssm {
    /// Write the PSSM in ASCII format compatible with NCBI's `-out_ascii_pssm` output.
    pub fn write_ascii<W: std::io::Write>(&self, out: &mut W, query: &[u8]) -> std::io::Result<()> {
        // Header line with residue codes
        let aa_order: &[u8] = &[1,16,13,4,3,15,5,7,8,9,11,10,12,6,14,17,18,20,22,19,2,23,21,25];
        write!(out, "       ")?;
        for &code in aa_order {
            write!(out, " {:>4}", ncbistdaa_to_ascii(code))?;
        }
        writeln!(out)?;

        for i in 0..self.query_len {
            let aa = if i < query.len() {
                ncbistdaa_to_ascii(query[i])
            } else { 'X' };
            write!(out, "{:>5} {}", i + 1, aa)?;
            for &code in aa_order {
                write!(out, " {:>4}", self.scores[i][code as usize])?;
            }
            writeln!(out)?;
        }
        writeln!(out)?;
        writeln!(out, "Lambda: {:.4}", self.lambda)?;
        Ok(())
    }

    /// Write binary PSSM checkpoint (simple format: query_len, lambda, then scores).
    pub fn write_checkpoint<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<()> {
        use std::io::Write;
        let len_bytes = (self.query_len as u32).to_le_bytes();
        out.write_all(&len_bytes)?;
        let lambda_bytes = self.lambda.to_le_bytes();
        out.write_all(&lambda_bytes)?;
        for row in &self.scores {
            for &val in row.iter() {
                out.write_all(&val.to_le_bytes())?;
            }
        }
        Ok(())
    }

    /// Read binary PSSM checkpoint.
    pub fn read_checkpoint<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        use std::io::Read;
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        reader.read_exact(&mut buf4)?;
        let query_len = u32::from_le_bytes(buf4) as usize;

        reader.read_exact(&mut buf8)?;
        let lambda = f64::from_le_bytes(buf8);

        let mut scores = Vec::with_capacity(query_len);
        for _ in 0..query_len {
            let mut row = [0i32; 28];
            for j in 0..28 {
                reader.read_exact(&mut buf4)?;
                row[j] = i32::from_le_bytes(buf4);
            }
            scores.push(row);
        }

        Ok(Pssm { query_len, scores, lambda })
    }
}
