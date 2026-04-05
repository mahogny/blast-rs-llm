//! Main BLAST search pipeline.

use std::io::BufRead;
use rayon::prelude::*;
use crate::db::BlastDb;

use crate::matrix::{ScoringMatrix, MatrixType, nt_to_2bit};
use crate::stats::{KarlinAltschul, GapPenalty, lookup_ka_params, blastn_ka_params};
use crate::lookup::{ProteinLookup, NucleotideLookup, DiscontiguousLookup};
use crate::extend::{ungapped_extend, gapped_extend, gapped_extend_score_only, ungapped_extend_nucleotide};
use crate::hsp::{Hsp, SearchResult};
use crate::translate::{six_frame_translate_with_code, strip_stops, TranslatedFrame, reverse_complement};
use crate::compo::{composition_ncbistdaa, adjust_evalue};

/// Parameters for a BLAST search.
#[derive(Debug, Clone)]
pub struct SearchParams {
    pub word_size: usize,
    pub matrix: MatrixType,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue_threshold: f64,
    pub max_target_seqs: usize,
    pub x_drop_ungapped: i32,
    pub x_drop_gapped: i32,
    pub ungapped_cutoff: i32,
    /// Minimum score to keep after gapped extension (before E-value filter)
    pub min_score: i32,
    /// Nucleotide match/mismatch (only used for blastn/blastx/tblastn/tblastx)
    pub match_score: i32,
    pub mismatch: i32,
    pub num_threads: usize,
    /// Apply SEG masking to protein queries/subjects.
    pub filter_low_complexity: bool,
    /// Apply composition-based statistics adjustment.
    pub comp_adjust: bool,
    /// Query strand: "both", "plus", or "minus" (for nucleotide searches).
    pub strand: String,
    /// Genetic code for query translation (1=standard, 2=vert mito, etc.)
    pub query_gencode: u8,
    /// Genetic code for database translation.
    pub db_gencode: u8,
    /// Maximum HSPs per subject sequence (None = no limit).
    pub max_hsps: Option<usize>,
    /// Culling limit (limits hits aligning to same region; None = disabled).
    pub culling_limit: Option<usize>,
    /// Use 2-hit algorithm for seeding (default true for protein)
    pub two_hit: bool,
    /// Window size for 2-hit algorithm (default 40)
    pub two_hit_window: usize,
    /// Final X-dropoff for gapped extension (higher than initial, default 25 for protein)
    pub x_drop_final: i32,
    /// Soft masking: mask for seeding/lookup only, use unmasked for extension
    pub soft_masking: bool,
    /// Lowercase masking: treat lowercase letters in query as masked
    pub lcase_masking: bool,
}

impl Default for SearchParams {
    fn default() -> Self { SearchParams::blastp_defaults() }
}

impl SearchParams {
    pub fn blastp_defaults() -> Self {
        SearchParams {
            word_size: 3,
            matrix: MatrixType::Blosum62,
            gap_open: 11, gap_extend: 1,
            evalue_threshold: 10.0, max_target_seqs: 500,
            x_drop_ungapped: 7, x_drop_gapped: 15,
            ungapped_cutoff: 0, min_score: 0,
            match_score: 2, mismatch: -3,
            num_threads: 0,
            filter_low_complexity: true,
            comp_adjust: true,
            strand: "both".to_string(),
            query_gencode: 1,
            db_gencode: 1,
            max_hsps: None,
            culling_limit: None,
            two_hit: true,
            two_hit_window: 40,
            x_drop_final: 25,
            soft_masking: true,
            lcase_masking: false,
        }
    }

    pub fn blastn_defaults() -> Self {
        SearchParams {
            word_size: 11,
            matrix: MatrixType::Blosum62, // not used for blastn
            gap_open: 5, gap_extend: 2,
            evalue_threshold: 10.0, max_target_seqs: 500,
            x_drop_ungapped: 20, x_drop_gapped: 30,
            ungapped_cutoff: 0, min_score: 0,
            match_score: 2, mismatch: -3,
            num_threads: 0,
            filter_low_complexity: true,
            comp_adjust: false, // composition adjustment not standard for blastn
            strand: "both".to_string(),
            query_gencode: 1,
            db_gencode: 1,
            max_hsps: None,
            culling_limit: None,
            two_hit: false,
            two_hit_window: 0,
            x_drop_final: 100,
            soft_masking: false,
            lcase_masking: false,
        }
    }

    /// Defaults for BLASTX (translate nucleotide query, search protein DB).
    pub fn blastx_defaults() -> Self {
        let mut p = Self::blastp_defaults();
        p.filter_low_complexity = true;
        p
    }

    /// Defaults for TBLASTN (protein query, translate nucleotide DB).
    pub fn tblastn_defaults() -> Self {
        Self::blastp_defaults()
    }

    /// Defaults for TBLASTX (translate both query and DB).
    pub fn tblastx_defaults() -> Self {
        let mut p = Self::blastp_defaults();
        p.gap_open = 0; p.gap_extend = 0; // TBLASTX is ungapped
        p
    }

    /// Short alias for [`blastp_defaults`](Self::blastp_defaults).
    pub fn blastp() -> Self { Self::blastp_defaults() }
    /// Short alias for [`blastn_defaults`](Self::blastn_defaults).
    pub fn blastn() -> Self { Self::blastn_defaults() }
    /// Short alias for [`blastx_defaults`](Self::blastx_defaults).
    pub fn blastx() -> Self { Self::blastx_defaults() }
    /// Short alias for [`tblastn_defaults`](Self::tblastn_defaults).
    pub fn tblastn() -> Self { Self::tblastn_defaults() }
    /// Short alias for [`tblastx_defaults`](Self::tblastx_defaults).
    pub fn tblastx() -> Self { Self::tblastx_defaults() }

    // ── Builder setters ──────────────────────────────────────────────────────

    pub fn evalue(mut self, v: f64) -> Self { self.evalue_threshold = v; self }
    pub fn max_target_seqs(mut self, v: usize) -> Self { self.max_target_seqs = v; self }
    pub fn matrix(mut self, v: MatrixType) -> Self { self.matrix = v; self }
    pub fn num_threads(mut self, v: usize) -> Self { self.num_threads = v; self }
    pub fn word_size(mut self, v: usize) -> Self { self.word_size = v; self }
    pub fn gap_open(mut self, v: i32) -> Self { self.gap_open = v; self }
    pub fn gap_extend(mut self, v: i32) -> Self { self.gap_extend = v; self }
    pub fn filter_low_complexity(mut self, v: bool) -> Self { self.filter_low_complexity = v; self }
    pub fn comp_adjust(mut self, v: bool) -> Self { self.comp_adjust = v; self }
    pub fn match_score(mut self, v: i32) -> Self { self.match_score = v; self }
    pub fn mismatch(mut self, v: i32) -> Self { self.mismatch = v; self }
    pub fn strand(mut self, v: &str) -> Self { self.strand = v.to_string(); self }
    pub fn query_gencode(mut self, v: u8) -> Self { self.query_gencode = v; self }
    pub fn db_gencode(mut self, v: u8) -> Self { self.db_gencode = v; self }
    pub fn max_hsps(mut self, v: Option<usize>) -> Self { self.max_hsps = v; self }
    pub fn culling_limit(mut self, v: Option<usize>) -> Self { self.culling_limit = v; self }
    pub fn two_hit(mut self, v: bool) -> Self { self.two_hit = v; self }
    pub fn two_hit_window(mut self, v: usize) -> Self { self.two_hit_window = v; self }
    pub fn x_drop_final(mut self, v: i32) -> Self { self.x_drop_final = v; self }
    pub fn soft_masking(mut self, v: bool) -> Self { self.soft_masking = v; self }
    pub fn lcase_masking(mut self, v: bool) -> Self { self.lcase_masking = v; self }
}

/// Neighbor word score threshold for protein lookup.
/// BLAST uses T=11 for BLOSUM62, word_size=3.
/// Ungapped Karlin-Altschul parameters (lambda, K) for gap_trigger computation.
/// NCBI uses kbp_std (ungapped params) not kbp_gap for the gap_trigger formula.
fn ungapped_ka_params(matrix: MatrixType) -> (f64, f64) {
    match matrix {
        MatrixType::Blosum62 => (0.3176, 0.134),
        MatrixType::Blosum45 => (0.2291, 0.0924),
        MatrixType::Blosum50 => (0.2326, 0.0450),
        MatrixType::Blosum80 => (0.3430, 0.177),
        MatrixType::Blosum90 => (0.3346, 0.190),
        MatrixType::Pam30    => (0.3400, 0.283),
        MatrixType::Pam70    => (0.3345, 0.310),
        MatrixType::Pam250   => (0.2252, 0.0868),
    }
}

pub(crate) fn neighbor_threshold(matrix: MatrixType, word_size: usize) -> i32 {
    match (matrix, word_size) {
        (MatrixType::Blosum62, 3) => 11,
        (MatrixType::Blosum62, 2) => 8,
        (MatrixType::Blosum45, 3) => 14,
        (MatrixType::Blosum50, 3) => 13,
        (MatrixType::Blosum80, 3) => 25,
        (MatrixType::Blosum90, 3) => 27,
        (MatrixType::Pam30, 3)    => 10,
        (MatrixType::Pam70, 3)    => 11,
        (MatrixType::Pam250, 3)   => 11,
        _ => 11,
    }
}


/// Run a protein BLAST search (blastp).
pub fn blast_search(
    db: &BlastDb,
    query: &[u8],   // Ncbistdaa encoded query
    params: &SearchParams,
) -> Vec<SearchResult> {
    // For soft masking: mask query for lookup, but use unmasked for extension
    let query_masked: Vec<u8> = if params.filter_low_complexity {
        let mut q = query.to_vec();
        crate::mask::apply_seg_ncbistdaa(&mut q);
        q
    } else {
        query.to_vec()
    };

    // For extension, use unmasked query if soft_masking is enabled
    let query_for_lookup = &query_masked;
    let query_for_extend: &[u8] = if params.soft_masking && params.filter_low_complexity {
        query // original unmasked
    } else {
        &query_masked
    };

    let matrix = ScoringMatrix::from_type(params.matrix);
    let gap = GapPenalty::new(params.gap_open, params.gap_extend);

    let ka = match lookup_ka_params(params.matrix, gap) {
        Some(k) => k,
        None => {
            eprintln!("Warning: no KA params for this matrix/gap combination, using defaults");
            lookup_ka_params(MatrixType::Blosum62, GapPenalty::blosum62_default()).unwrap()
        }
    };

    let db_len = db.volume_length();
    let num_seqs = db.num_sequences() as u64;
    let (eff_query_len, eff_db_len) = ka.effective_lengths(query_for_lookup.len(), db_len, num_seqs);

    let threshold = neighbor_threshold(params.matrix, params.word_size);
    let lookup = ProteinLookup::build(query_for_lookup, params.word_size, &matrix, threshold);

    // NCBI-style two-threshold approach:
    // 1. ungapped_cutoff: minimum ungapped score to keep a hit (low — keeps weak hits)
    // 2. gap_trigger: minimum ungapped score to attempt gapped extension (higher — saves time)
    // Hits below gap_trigger are kept but not gapped-extended; they may still appear in output
    // if they pass the evalue threshold based on ungapped score alone.
    // NCBI uses UNGAPPED KA parameters (kbp_std) for gap_trigger, not gapped params.
    // Ungapped params have higher lambda and K, giving a lower gap_trigger threshold.
    let (ug_lambda, ug_k) = ungapped_ka_params(params.matrix);
    let gap_trigger_bits = 22.0f64;
    let gap_trigger_score = ((gap_trigger_bits * std::f64::consts::LN_2 + ug_k.ln()) / ug_lambda) as i32;
    // Match NCBI: ungapped hits must meet gap_trigger to be saved for gapped extension.
    // Cap at the maximum achievable score for the query (NCBI's cutoff_score_max),
    // so short queries aren't filtered too aggressively.
    let max_self_score: i32 = query_for_extend.iter()
        .map(|&q| matrix.scores[q as usize % 28][q as usize % 28])
        .sum();
    let capped_trigger = gap_trigger_score.max(1).min(max_self_score);
    let effective_cutoff = if params.ungapped_cutoff > 0 {
        params.ungapped_cutoff
    } else {
        capped_trigger
    };

    let query_comp = if params.comp_adjust { Some(composition_ncbistdaa(query_for_extend)) } else { None };

    // Precompute query score profile for fast ungapped extension:
    // profile[q_residue][s_residue] = score, accessed as profile[query[i]][subject[j]]
    let score_profile: Vec<[i32; 28]> = query_for_extend.iter().map(|&q| {
        let mut row = [0i32; 28];
        for r in 0..28u8 {
            row[r as usize] = matrix.scores[q as usize % 28][r as usize];
        }
        row
    }).collect();

    // Process each OID (parallelized with rayon)
    // Thread-local scratch buffers avoid per-subject allocation
    use std::cell::RefCell;
    thread_local! {
        static SCRATCH: RefCell<SearchScratch> = RefCell::new(SearchScratch::new());
    }

    let oids: Vec<u32> = (0..db.num_sequences()).collect();

    let mut results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
        let subject = match db.get_sequence_protein_raw(oid) {
            Ok(s) => s,
            Err(_) => return None,
        };
        if subject.len() < lookup.word_size { return None; }

        let mut hsps = SCRATCH.with(|scratch| {
            let mut scratch = scratch.borrow_mut();
            if scratch.diag_mask == 0 {
                scratch.init_diags(query_for_extend.len());
            }
            search_one_protein_scratch(
                query_for_extend, subject, &lookup, &matrix, &ka, params,
                eff_query_len, eff_db_len, &score_profile, &mut scratch,
                effective_cutoff, capped_trigger,
            )
        });

        // Composition adjustment
        if let Some(ref qc) = query_comp {
            let sc = composition_ncbistdaa(subject);
            for hsp in &mut hsps {
                hsp.evalue = adjust_evalue(
                    hsp.evalue, hsp.score, qc, &sc, &matrix,
                    ka.lambda, ka.k, eff_query_len, eff_db_len,
                );
            }
            hsps.retain(|h| h.evalue <= params.evalue_threshold);
        }

        if hsps.is_empty() { return None; }

        if let Some(max) = params.max_hsps {
            hsps.truncate(max);
        }

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

/// Diagonal state: packed last_hit + flag in a single u32.
/// Matches NCBI's DiagStruct (last_hit:31, flag:1).
/// Uses a global offset counter to avoid zeroing the array between subjects.
#[derive(Clone, Copy)]
struct DiagEntry(u32);

impl DiagEntry {
    const FLAG_BIT: u32 = 1 << 31;

    #[inline]
    fn new(last_hit: i32, flag: bool) -> Self {
        let v = (last_hit as u32) & !Self::FLAG_BIT;
        DiagEntry(if flag { v | Self::FLAG_BIT } else { v })
    }

    #[inline]
    fn last_hit(self) -> i32 { (self.0 & !Self::FLAG_BIT) as i32 }

    #[inline]
    fn flag(self) -> bool { self.0 & Self::FLAG_BIT != 0 }
}

/// Thread-local scratch buffers reused across subjects.
/// Uses NCBI's offset trick: instead of zeroing the diagonal array between subjects,
/// we increment a global offset so stale entries are recognized by their old offset values.
struct SearchScratch {
    diag_array: Vec<DiagEntry>,
    diag_mask: usize,
    /// Offset added to subject positions to distinguish entries from different subjects.
    /// Incremented by (max_subject_len + window + 1) after each subject.
    diag_offset: i32,
    ungapped_hits: Vec<(i32, usize, usize, usize, usize)>,
    covered_query: Vec<bool>,
}

impl SearchScratch {
    fn new() -> Self {
        SearchScratch {
            diag_array: Vec::new(),
            diag_mask: 0,
            diag_offset: 0,
            ungapped_hits: Vec::new(),
            covered_query: Vec::new(),
        }
    }

    /// Initialize diagonal array. Size is rounded to power of 2 for mask-based indexing.
    /// The offset trick avoids zeroing between subjects.
    fn init_diags(&mut self, query_len: usize) {
        // Smallest power of 2 >= query_len
        let mut size = 1;
        while size < query_len { size <<= 1; }
        if self.diag_array.len() < size {
            self.diag_array.resize(size, DiagEntry(0));
        }
        self.diag_mask = size - 1;
        self.diag_offset = 0;
    }

    /// Advance offset after processing one subject to invalidate stale diagonal entries.
    #[inline]
    fn advance_offset(&mut self, subject_len: usize, window: usize) {
        self.diag_offset += (subject_len + window + 1) as i32;
        // Wrap to prevent overflow
        if self.diag_offset > i32::MAX / 2 {
            self.diag_offset = 0;
            self.diag_array.fill(DiagEntry(0));
        }
    }

    fn reset_covered(&mut self, len: usize) {
        self.covered_query.resize(len, false);
        self.covered_query[..len].fill(false);
    }
}

/// Optimized protein search using NCBI-style architecture:
/// 1. Batch scan: collect all (q_off, s_off) hit pairs (cache-friendly)
/// 2. Process hits: diagonal tracking with offset trick (no per-subject zeroing)
/// 3. Ungapped extension only for two-hit or single-hit candidates
#[allow(clippy::too_many_arguments)]
fn search_one_protein_scratch(
    query: &[u8],
    subject: &[u8],
    lookup: &ProteinLookup,
    matrix: &ScoringMatrix,
    ka: &KarlinAltschul,
    params: &SearchParams,
    eff_query_len: usize,
    eff_db_len: u64,
    score_profile: &[[i32; 28]],
    scratch: &mut SearchScratch,
    ungapped_cutoff: i32,
    gap_trigger: i32,
) -> Vec<Hsp> {
    let slen = subject.len();
    let ws = lookup.word_size;
    if slen < ws { return vec![]; }

    let cutoff = ungapped_cutoff;
    let x_drop = params.x_drop_ungapped;
    let diag_mask = scratch.diag_mask;
    let diag_offset = scratch.diag_offset;
    let window = params.two_hit_window as i32;
    // NCBI-style shift-based rolling hash: (code << charsize | residue) & mask
    let charsize = lookup.charsize;
    let mask = lookup.mask;

    scratch.ungapped_hits.clear();

    // Dispatch to specialized scanning functions (NCBI splits these too).
    // Separate functions keep instruction footprint small for icache.
    if params.two_hit {
        scan_two_hit(
            subject, lookup, score_profile, x_drop, cutoff,
            charsize, mask, ws, diag_mask, diag_offset, window,
            &mut scratch.diag_array, &mut scratch.ungapped_hits,
            query,
        );
    } else {
        scan_one_hit(
            subject, lookup, score_profile, x_drop, cutoff,
            charsize, mask, ws, diag_mask, diag_offset,
            &mut scratch.diag_array, &mut scratch.ungapped_hits,
            query,
        );
    }

    scratch.advance_offset(slen, params.two_hit_window);

    if scratch.ungapped_hits.is_empty() { return vec![]; }

    scratch.ungapped_hits.sort_unstable_by(|a, b| b.0.cmp(&a.0));

    // Gapped extension with two-stage approach (reuse covered_query buffer)
    scratch.reset_covered(query.len());
    let mut hsps = Vec::new();

    for &(_, qs, qe, ss, se) in &scratch.ungapped_hits {
        let center_q = (qs + qe) / 2;
        if center_q < query.len() && scratch.covered_query[center_q] { continue; }
        let center_s = (ss + se) / 2;

        // Match NCBI: all saved hits (score >= gap_trigger) get gapped extension.
        // No additional threshold check needed here — effective_cutoff already
        // filtered at save time.

        // Stage 1: preliminary score-only extension (fast rejection)
        let prelim_score = gapped_extend_score_only(
            query, subject, center_q, center_s, matrix,
            params.gap_open, params.gap_extend, params.x_drop_gapped,
        );
        // Match NCBI: check raw score against gap_trigger (same as ungapped cutoff)
        if prelim_score < gap_trigger { continue; }

        // Stage 2: full extension with traceback
        let final_x_drop = params.x_drop_final.max(params.x_drop_gapped);
        let gh = gapped_extend(
            query, subject, center_q, center_s, matrix,
            params.gap_open, params.gap_extend, final_x_drop,
        );
        // Match NCBI: gapped result must also meet cutoff_score
        if gh.score < gap_trigger { continue; }
        let evalue = ka.evalue(gh.score, eff_query_len, eff_db_len);
        if evalue > params.evalue_threshold { continue; }

        let dominated = hsps.iter().any(|existing: &Hsp| {
            let ov_start = gh.q_start.max(existing.query_start);
            let ov_end = gh.q_end.min(existing.query_end);
            if ov_end <= ov_start { return false; }
            let overlap = ov_end - ov_start;
            let shorter = (gh.q_end - gh.q_start).min(existing.query_end - existing.query_start);
            overlap * 2 > shorter
        });
        if dominated { continue; }

        let bit_score = ka.bit_score(gh.score);
        for i in gh.q_start.min(query.len())..gh.q_end.min(query.len()) {
            scratch.covered_query[i] = true;
        }

        let query_aln = crate::db::sequence::decode_protein(&gh.query_aln);
        let subject_aln = crate::db::sequence::decode_protein(&gh.subject_aln);

        hsps.push(Hsp {
            score: gh.score, bit_score, evalue,
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

/// Two-hit scanning: faithful port of NCBI's s_BlastAaWordFinder_TwoHit.
/// Kept as a separate function for icache efficiency.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn scan_two_hit(
    subject: &[u8],
    lookup: &crate::lookup::ProteinLookup,
    score_profile: &[[i32; 28]],
    x_drop: i32,
    cutoff: i32,
    charsize: u32,
    mask: u32,
    ws: usize,
    diag_mask: usize,
    diag_offset: i32,
    window: i32,
    diag_array: &mut [DiagEntry],
    ungapped_hits: &mut Vec<(i32, usize, usize, usize, usize)>,
    query: &[u8],
) {
    let slen = subject.len();
    let num_words = slen - ws + 1;
    let ws_i32 = ws as i32;

    unsafe {
        let subj_ptr = subject.as_ptr();
        let diag_arr = diag_array.as_mut_ptr();
        let pv_ptr = lookup.presence.as_ptr();
        let bb_ptr = lookup.backbone.as_ptr();
        let ovfl_ptr = lookup.overflow.as_ptr();
        let prof_ptr = score_profile.as_ptr();

        let mut code = 0u32;
        for i in 0..ws {
            code = (code << charsize) | (*subj_ptr.add(i) as u32 & 0x1F);
        }

        for s_pos in 0..num_words {
            if s_pos > 0 {
                code = ((code << charsize) | (*subj_ptr.add(s_pos + ws - 1) as u32 & 0x1F)) & mask;
            }

            let pv_word = *pv_ptr.add((code as usize) >> 6);
            if pv_word & (1u64 << (code & 63)) == 0 { continue; }

            let cell = &*bb_ptr.add(code as usize);
            let n = cell.num_used as usize;
            if n == 0 { continue; }

            let (hits_ptr, hits_len) = if n <= 3 {
                (cell.data.as_ptr(), n)
            } else {
                (ovfl_ptr.add(cell.data[0] as usize), n)
            };

            let s_off_i32 = s_pos as i32;
            let s_with_offset = s_off_i32 + diag_offset;

            // NCBI-style hit processing: minimize branches using combined conditions.
            // The flag+last_hit are packed in a single u32, enabling branchless comparisons.
            for h in 0..hits_len {
                let q_pos = *hits_ptr.add(h) as usize;
                let diag_coord = q_pos.wrapping_sub(s_pos) & diag_mask;
                let dc = diag_arr.add(diag_coord);
                let entry = *dc;

                // Combined check: if flag is set and we haven't passed the last_hit, skip.
                // If flag is set but we HAVE passed, start new hit.
                if entry.flag() {
                    if s_with_offset < entry.last_hit() { continue; }
                    *dc = DiagEntry::new(s_with_offset, false);
                    continue;
                }

                // Not flagged: check two-hit window
                let diff = s_off_i32 - (entry.last_hit() - diag_offset);

                // Record new hit if outside window or first hit on diagonal
                if (diff as u32) >= (window as u32) {
                    // Covers diff >= window AND diff < 0 (negative wraps to large u32)
                    *dc = DiagEntry::new(s_with_offset, false);
                    continue;
                }

                if diff < ws_i32 { continue; }

                // Two hits within window — extend
                let hit = ungapped_extend_unsafe(
                    query, subject, q_pos, s_pos, prof_ptr, x_drop,
                );
                if hit.0 >= cutoff {
                    ungapped_hits.push(hit);
                    *dc = DiagEntry::new(hit.4 as i32 + diag_offset, true);
                } else {
                    *dc = DiagEntry::new(s_with_offset, false);
                }
            }
        }
    }
}

/// One-hit scanning: faithful port of NCBI's s_BlastAaWordFinder_OneHit.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn scan_one_hit(
    subject: &[u8],
    lookup: &crate::lookup::ProteinLookup,
    score_profile: &[[i32; 28]],
    x_drop: i32,
    cutoff: i32,
    charsize: u32,
    mask: u32,
    ws: usize,
    diag_mask: usize,
    diag_offset: i32,
    diag_array: &mut [DiagEntry],
    ungapped_hits: &mut Vec<(i32, usize, usize, usize, usize)>,
    query: &[u8],
) {
    let slen = subject.len();
    let num_words = slen - ws + 1;

    unsafe {
        let subj_ptr = subject.as_ptr();
        let diag_arr = diag_array.as_mut_ptr();
        let pv_ptr = lookup.presence.as_ptr();
        let bb_ptr = lookup.backbone.as_ptr();
        let ovfl_ptr = lookup.overflow.as_ptr();
        let prof_ptr = score_profile.as_ptr();

        let mut code = 0u32;
        for i in 0..ws {
            code = (code << charsize) | (*subj_ptr.add(i) as u32 & 0x1F);
        }

        for s_pos in 0..num_words {
            if s_pos > 0 {
                code = ((code << charsize) | (*subj_ptr.add(s_pos + ws - 1) as u32 & 0x1F)) & mask;
            }

            let pv_word = *pv_ptr.add((code as usize) >> 6);
            if pv_word & (1u64 << (code & 63)) == 0 { continue; }

            let cell = &*bb_ptr.add(code as usize);
            let n = cell.num_used as usize;
            if n == 0 { continue; }

            let (hits_ptr, hits_len) = if n <= 3 {
                (cell.data.as_ptr(), n)
            } else {
                (ovfl_ptr.add(cell.data[0] as usize), n)
            };

            let s_with_offset = s_pos as i32 + diag_offset;

            for h in 0..hits_len {
                let q_pos = *hits_ptr.add(h) as usize;
                let diag_coord = q_pos.wrapping_sub(s_pos) & diag_mask;
                let entry = *diag_arr.add(diag_coord);

                if entry.flag() && s_with_offset < entry.last_hit() { continue; }

                let hit = ungapped_extend_unsafe(
                    query, subject, q_pos, s_pos, prof_ptr, x_drop,
                );
                if hit.0 >= cutoff {
                    *diag_arr.add(diag_coord) = DiagEntry::new(hit.4 as i32 + diag_offset, true);
                    ungapped_hits.push(hit);
                }
            }
        }
    }
}

/// Raw pointer ungapped extension — faithful port of NCBI's s_BlastAaExtendRight/Left.
/// All bounds are pre-validated by the caller; inner loops use pure pointer arithmetic.
#[inline(always)]
unsafe fn ungapped_extend_unsafe(
    query: &[u8], subject: &[u8], q_pos: usize, s_pos: usize,
    prof_ptr: *const [i32; 28], x_drop: i32,
) -> (i32, usize, usize, usize, usize) {
    let qlen = query.len();
    let slen = subject.len();
    let s_ptr = subject.as_ptr();

    // Word scan: find best start within word_size=3 positions
    let ws = 3usize;
    let scan_end = ws.min(qlen - q_pos).min(slen - s_pos);
    let mut scan_score = 0i32;
    let mut best_scan = 0i32;
    let mut right_d = 0usize;
    for i in 0..scan_end {
        let row = &*prof_ptr.add(q_pos + i);
        scan_score += row[*s_ptr.add(s_pos + i) as usize];
        if scan_score > best_scan {
            best_scan = scan_score;
            right_d = i + 1;
        }
    }

    let eq = q_pos + right_d;
    let es = s_pos + right_d;

    // Extend left (NCBI s_BlastAaExtendLeft)
    let n_left = eq.min(es);
    let mut score = 0i32;
    let mut maxscore = 0i32;
    let mut best_left = 0usize;
    for i in 1..=n_left {
        let row = &*prof_ptr.add(eq - i);
        score += row[*s_ptr.add(es - i) as usize];
        if score > maxscore {
            maxscore = score;
            best_left = i;
        }
        if maxscore - score >= x_drop { break; }
    }
    let left_score = maxscore;

    // Extend right (NCBI s_BlastAaExtendRight)
    let n_right = (qlen - eq).min(slen - es);
    score = 0;
    maxscore = 0;
    let mut best_right = 0usize;
    for i in 0..n_right {
        let row = &*prof_ptr.add(eq + i);
        score += row[*s_ptr.add(es + i) as usize];
        if score > maxscore {
            maxscore = score;
            best_right = i + 1;
        }
        if score <= 0 || maxscore - score >= x_drop { break; }
    }

    (left_score + maxscore, eq - best_left, eq + best_right, es - best_left, es + best_right)
}

/// Inline ungapped extension matching NCBI's s_BlastAaExtendTwoHit/Right/Left.
/// Uses unsafe pointer access to match NCBI's performance (no bounds checks in hot loop).
/// Returns (score, q_start, q_end, s_start, s_end).
#[inline]
#[allow(dead_code)]
fn ungapped_extend_inline(
    query: &[u8], subject: &[u8], q_pos: usize, s_pos: usize,
    score_profile: &[[i32; 28]], x_drop: i32,
) -> (i32, usize, usize, usize, usize) {
    let slen = subject.len();
    let qlen = query.len();

    // NCBI word-scan: find best starting position within word (lines 1104-1121)
    let ws = 3usize;
    let mut scan_score = 0i32;
    let mut best_scan = 0i32;
    let mut right_d = 0usize;
    let scan_end = ws.min(qlen - q_pos).min(slen - s_pos);
    for i in 0..scan_end {
        scan_score += unsafe {
            *score_profile.get_unchecked(q_pos + i).get_unchecked(
                *subject.get_unchecked(s_pos + i) as usize
            )
        };
        if scan_score > best_scan {
            best_scan = scan_score;
            right_d = i + 1;
        }
    }

    let ext_q = q_pos + right_d;
    let ext_s = s_pos + right_d;

    // Extend left (NCBI s_BlastAaExtendLeft pattern: pointer arithmetic, no bounds check)
    let n_left = ext_q.min(ext_s);
    let mut best = 0i32;
    let mut score = 0i32;
    let mut best_i = 0usize;
    for i in 1..=n_left {
        score += unsafe {
            *score_profile.get_unchecked(ext_q - i).get_unchecked(
                *subject.get_unchecked(ext_s - i) as usize
            )
        };
        if score > best {
            best = score;
            best_i = i;
        } else if best - score > x_drop {
            break;
        }
    }
    let left_score = best;
    let best_q_start = ext_q - best_i;
    let best_s_start = ext_s - best_i;

    // Extend right (NCBI s_BlastAaExtendRight pattern)
    let n_right = (qlen - ext_q).min(slen - ext_s);
    score = 0;
    best = 0;
    best_i = 0;
    for i in 0..n_right {
        score += unsafe {
            *score_profile.get_unchecked(ext_q + i).get_unchecked(
                *subject.get_unchecked(ext_s + i) as usize
            )
        };
        if score > best {
            best = score;
            best_i = i + 1;
        }
        if score <= 0 || best - score >= x_drop {
            break;
        }
    }
    let best_q_end = ext_q + best_i;
    let best_s_end = ext_s + best_i;

    (left_score + best, best_q_start, best_q_end, best_s_start, best_s_end)
}

/// Optimized protein search for one subject (legacy, used by search_one_protein_fast).
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
fn search_one_protein_fast(
    query: &[u8],
    subject: &[u8],
    lookup: &ProteinLookup,
    matrix: &ScoringMatrix,
    ka: &KarlinAltschul,
    params: &SearchParams,
    eff_query_len: usize,
    eff_db_len: u64,
    score_profile: &[[i32; 28]],
) -> Vec<Hsp> {
    let slen = subject.len();
    let ws = lookup.word_size;
    if slen < ws { return vec![]; }

    // Rolling word codes for subject
    let mask = lookup.mask;
    let num_words = slen - ws + 1;

    // Quick scan: count total hits before doing extensions.
    // If no hits at all, skip this subject immediately.
    let mut total_hits = 0u32;
    let mut code = 0u32;
    for &r in subject.iter().take(ws) {
        code = (code << 5) | (r as u32 & 0x1F);
    }
    if !lookup.get_hits(code).is_empty() { total_hits += 1; }
    for &r in subject.iter().skip(ws) {
        code = ((code << 5) | (r as u32 & 0x1F)) & mask;
        if !lookup.get_hits(code).is_empty() { total_hits += 1; }
    }
    if total_hits == 0 { return vec![]; }

    let diag_offset = query.len();
    let diag_len = query.len() + slen + 1;
    let cutoff = params.ungapped_cutoff.max(1);

    // Ungapped extension using score profile (avoids matrix.score() bounds checks)
    let ungapped_extend_fast = |q_pos: usize, s_pos: usize| -> (i32, usize, usize, usize, usize) {
        let mut best = 0i32;
        let mut best_q_start = q_pos;
        let mut best_s_start = s_pos;

        // Extend left
        let mut score = 0i32;
        let mut qi = q_pos;
        let mut si = s_pos;
        while qi > 0 && si > 0 {
            qi -= 1;
            si -= 1;
            score += score_profile[qi][subject[si] as usize];
            if score > best {
                best = score;
                best_q_start = qi;
                best_s_start = si;
            } else if best - score > params.x_drop_ungapped {
                break;
            }
        }
        let left_score = best;

        // Extend right
        let mut qi = q_pos;
        let mut si = s_pos;
        score = 0;
        best = 0;
        let mut best_q_end = q_pos;
        let mut best_s_end = s_pos;
        while qi < query.len() && si < slen {
            score += score_profile[qi][subject[si] as usize];
            if score > best {
                best = score;
                best_q_end = qi + 1;
                best_s_end = si + 1;
            } else if best - score > params.x_drop_ungapped {
                break;
            }
            qi += 1;
            si += 1;
        }

        let center = if q_pos < query.len() && s_pos < slen {
            score_profile[q_pos][subject[s_pos] as usize]
        } else { 0 };

        (left_score + center + best, best_q_start, best_q_end, best_s_start, best_s_end)
    };

    let mut ungapped_hits: Vec<(i32, usize, usize, usize, usize)> = Vec::new();

    if params.two_hit {
        let mut diag_last: Vec<i32> = vec![i32::MIN; diag_len];
        let mut diag_extended: Vec<bool> = vec![false; diag_len];

        let mut code = 0u32;
        for &r in subject.iter().take(ws) {
            code = (code << 5) | (r as u32 & 0x1F);
        }
        // Process first word
        for &q_pos in lookup.get_hits(code) {
            let q_pos = q_pos as usize;
            let diag = (0isize - q_pos as isize + diag_offset as isize) as usize;
            if diag < diag_len { diag_last[diag] = 0; }
        }
        // Process remaining words with rolling code
        for s_pos in 1..num_words {
            code = ((code << 5) | (subject[s_pos + ws - 1] as u32 & 0x1F)) & mask;
            let positions = lookup.get_hits(code);
            if positions.is_empty() { continue; }
            for &q_pos in positions {
                let q_pos = q_pos as usize;
                let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
                if diag >= diag_len { continue; }
                if diag_extended[diag] { continue; }

                let prev = diag_last[diag];
                if prev == i32::MIN {
                    diag_last[diag] = s_pos as i32;
                    continue;
                }
                let dist = (s_pos as i32) - prev;
                if dist > params.two_hit_window as i32 {
                    diag_last[diag] = s_pos as i32;
                    continue;
                }

                let (score, qs, qe, ss, se) = ungapped_extend_fast(q_pos, s_pos);
                if score >= cutoff {
                    diag_extended[diag] = true;
                    ungapped_hits.push((score, qs, qe, ss, se));
                }
            }
        }
    } else {
        let mut diag_hit: Vec<bool> = vec![false; diag_len];

        let mut code = 0u32;
        for &r in subject.iter().take(ws) {
            code = (code << 5) | (r as u32 & 0x1F);
        }
        // Inline the scan + extend loop with rolling codes
        {
            let positions = lookup.get_hits(code);
            for &q_pos in positions {
                let q_pos = q_pos as usize;
                let diag = (0isize - q_pos as isize + diag_offset as isize) as usize;
                if diag < diag_len && diag_hit[diag] { continue; }
                let (score, qs, qe, ss, se) = ungapped_extend_fast(q_pos, 0);
                if score >= cutoff {
                    if diag < diag_len { diag_hit[diag] = true; }
                    ungapped_hits.push((score, qs, qe, ss, se));
                }
            }
        }
        for s_pos in 1..num_words {
            code = ((code << 5) | (subject[s_pos + ws - 1] as u32 & 0x1F)) & mask;
            let positions = lookup.get_hits(code);
            if positions.is_empty() { continue; }
            for &q_pos in positions {
                let q_pos = q_pos as usize;
                let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
                if diag < diag_len && diag_hit[diag] { continue; }
                let (score, qs, qe, ss, se) = ungapped_extend_fast(q_pos, s_pos);
                if score >= cutoff {
                    if diag < diag_len { diag_hit[diag] = true; }
                    ungapped_hits.push((score, qs, qe, ss, se));
                }
            }
        }
    }

    if ungapped_hits.is_empty() { return vec![]; }

    ungapped_hits.sort_by(|a, b| b.0.cmp(&a.0));

    // Gapped extension with two-stage approach
    let mut hsps = Vec::new();
    let mut covered_query: Vec<bool> = vec![false; query.len()];

    for &(_, qs, qe, ss, se) in &ungapped_hits {
        let center_q = (qs + qe) / 2;
        if center_q < covered_query.len() && covered_query[center_q] {
            continue;
        }
        let center_s = (ss + se) / 2;

        // Stage 1: score-only preliminary check
        let prelim_score = gapped_extend_score_only(
            query, subject, center_q, center_s, matrix,
            params.gap_open, params.gap_extend, params.x_drop_gapped,
        );
        if prelim_score <= 0 { continue; }
        let prelim_evalue = ka.evalue(prelim_score, eff_query_len, eff_db_len);
        if prelim_evalue > params.evalue_threshold { continue; }

        // Stage 2: full extension with traceback
        let final_x_drop = params.x_drop_final.max(params.x_drop_gapped);
        let gh = gapped_extend(
            query, subject, center_q, center_s, matrix,
            params.gap_open, params.gap_extend, final_x_drop,
        );
        if gh.score <= 0 { continue; }
        let evalue = ka.evalue(gh.score, eff_query_len, eff_db_len);
        if evalue > params.evalue_threshold { continue; }

        let dominated = hsps.iter().any(|existing: &Hsp| {
            let ov_start = gh.q_start.max(existing.query_start);
            let ov_end = gh.q_end.min(existing.query_end);
            if ov_end <= ov_start { return false; }
            let overlap = ov_end - ov_start;
            let shorter = (gh.q_end - gh.q_start).min(existing.query_end - existing.query_start);
            overlap * 2 > shorter
        });
        if dominated { continue; }

        let bit_score = ka.bit_score(gh.score);
        let cover_start = gh.q_start.min(covered_query.len());
        let cover_end = gh.q_end.min(covered_query.len());
        for item in covered_query.iter_mut().take(cover_end).skip(cover_start) {
            *item = true;
        }

        let query_aln = crate::db::sequence::decode_protein(&gh.query_aln);
        let subject_aln = crate::db::sequence::decode_protein(&gh.subject_aln);

        hsps.push(Hsp {
            score: gh.score, bit_score, evalue,
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

/// Search one protein subject sequence (original version, used by tblastn/tblastx/psiblast).
#[allow(clippy::too_many_arguments)]
fn search_one_protein(
    query: &[u8],
    subject: &[u8],
    lookup: &ProteinLookup,
    matrix: &ScoringMatrix,
    ka: &KarlinAltschul,
    params: &SearchParams,
    eff_query_len: usize,
    eff_db_len: u64,
) -> Vec<Hsp> {
    let slen = subject.len();
    let ws = lookup.word_size;
    if slen < ws { return vec![]; }

    // Precompute rolling word codes for the subject to avoid per-position re-encoding
    let mask = lookup.mask;
    let mut subject_codes: Vec<u32> = Vec::with_capacity(slen - ws + 1);
    {
        let mut code = 0u32;
        for &r in subject.iter().take(ws) {
            code = (code << 5) | (r as u32 & 0x1F);
        }
        subject_codes.push(code);
        for &r in subject.iter().skip(ws) {
            code = ((code << 5) | (r as u32 & 0x1F)) & mask;
            subject_codes.push(code);
        }
    }

    let diag_offset = query.len();
    let diag_len = query.len() + slen + 1;
    let mut ungapped_hits = Vec::new();

    // Match NCBI: compute gap_trigger using UNGAPPED KA params, capped at max self-score
    let (ug_lam, ug_k) = ungapped_ka_params(params.matrix);
    let gap_trigger_bits = 22.0f64;
    let gap_trigger = ((gap_trigger_bits * std::f64::consts::LN_2 + ug_k.ln()) / ug_lam) as i32;
    let max_self_score: i32 = query.iter()
        .map(|&q| matrix.scores[q as usize % 28][q as usize % 28])
        .sum();
    let cutoff = if params.ungapped_cutoff > 0 {
        params.ungapped_cutoff
    } else {
        gap_trigger.max(1).min(max_self_score)
    };

    if params.two_hit {
        let mut diag_last: Vec<i32> = vec![i32::MIN; diag_len];
        let mut diag_extended: Vec<bool> = vec![false; diag_len];

        for (s_pos, &word_code) in subject_codes.iter().enumerate() {
            let positions = lookup.get_hits(word_code);
            if !positions.is_empty() {
                let q_positions = positions;
                for &q_pos in q_positions {
                    let q_pos = q_pos as usize;
                    let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
                    if diag >= diag_extended.len() { continue; }
                    if diag_extended[diag] { continue; }

                    let prev = diag_last[diag];
                    if prev == i32::MIN {
                        diag_last[diag] = s_pos as i32;
                        continue;
                    }

                    let dist = (s_pos as i32) - prev;
                    if dist > params.two_hit_window as i32 {
                        diag_last[diag] = s_pos as i32;
                        continue;
                    }

                    let hit = ungapped_extend(query, subject, q_pos, s_pos, matrix, params.x_drop_ungapped);
                    if hit.score >= cutoff {
                        diag_extended[diag] = true;
                        ungapped_hits.push(hit);
                    }
                }
            }
        }
    } else {
        // Single-hit mode (original behavior)
        let mut diag_hit: Vec<bool> = vec![false; diag_len];

        for (s_pos, &word_code) in subject_codes.iter().enumerate() {
            let positions = lookup.get_hits(word_code);
            if !positions.is_empty() {
                for &q_pos in positions {
                    let q_pos = q_pos as usize;
                    let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
                    if diag < diag_hit.len() && diag_hit[diag] {
                        continue;
                    }

                    let hit = ungapped_extend(query, subject, q_pos, s_pos, matrix, params.x_drop_ungapped);
                    if hit.score >= cutoff {
                        if diag < diag_hit.len() {
                            diag_hit[diag] = true;
                        }
                        ungapped_hits.push(hit);
                    }
                }
            }
        }
    }

    if ungapped_hits.is_empty() { return vec![]; }

    // Sort ungapped hits by score descending, deduplicate overlapping ones
    ungapped_hits.sort_by(|a, b| b.score.cmp(&a.score));

    // Two-stage gapped extension:
    // Stage 1: Score-only with low X-drop (fast rejection of most candidates)
    // Stage 2: Full extension with traceback only for hits that pass stage 1
    let mut hsps = Vec::new();
    let mut covered_query: Vec<bool> = vec![false; query.len()];

    for uh in ungapped_hits {
        let center_q = (uh.q_start + uh.q_end) / 2;
        if center_q < covered_query.len() && covered_query[center_q] {
            continue;
        }

        let center_s = (uh.s_start + uh.s_end) / 2;

        // Stage 1: quick score-only check (match NCBI: raw score threshold)
        let prelim_score = gapped_extend_score_only(
            query, subject, center_q, center_s, matrix,
            params.gap_open, params.gap_extend, params.x_drop_gapped,
        );
        if prelim_score < cutoff { continue; }

        // Stage 2: full extension with traceback
        let final_x_drop = params.x_drop_final.max(params.x_drop_gapped);
        let gh = gapped_extend(
            query, subject, center_q, center_s, matrix,
            params.gap_open, params.gap_extend, final_x_drop,
        );

        if gh.score < cutoff { continue; }

        let evalue = ka.evalue(gh.score, eff_query_len, eff_db_len);
        if evalue > params.evalue_threshold { continue; }

        // Check if the gapped result substantially overlaps an existing HSP
        // (query region overlap > 50% of the shorter alignment)
        let dominated = hsps.iter().any(|existing: &Hsp| {
            let ov_start = gh.q_start.max(existing.query_start);
            let ov_end = gh.q_end.min(existing.query_end);
            if ov_end <= ov_start { return false; }
            let overlap = ov_end - ov_start;
            let shorter = (gh.q_end - gh.q_start).min(existing.query_end - existing.query_start);
            overlap * 2 > shorter
        });
        if dominated { continue; }

        let bit_score = ka.bit_score(gh.score);

        // Mark covered region
        let cover_start = gh.q_start.min(covered_query.len());
        let cover_end = gh.q_end.min(covered_query.len());
        for item in covered_query.iter_mut().take(cover_end).skip(cover_start) {
            *item = true;
        }

        // Convert alignment to ASCII for display
        let query_aln = crate::db::sequence::decode_protein(&gh.query_aln);
        let subject_aln = crate::db::sequence::decode_protein(&gh.subject_aln);

        hsps.push(Hsp {
            score: gh.score, bit_score, evalue,
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

/// Run a nucleotide BLAST search (blastn).
///
/// Respects `params.strand`: "plus" searches forward only, "minus" searches
/// reverse complement only, "both" (default) searches both strands.
pub fn blastn_search(
    db: &BlastDb,
    query: &[u8],   // ASCII nucleotide query
    params: &SearchParams,
) -> Vec<SearchResult> {
    // Determine which query orientations to search based on strand selection.
    // Each entry: (query bytes, query_frame: 0=forward, 0=revcomp but we use 1/-1 convention)
    let queries_to_search: Vec<(Vec<u8>, i32)> = match params.strand.as_str() {
        "plus"  => vec![(query.to_vec(), 1)],
        "minus" => vec![(reverse_complement(query), -1)],
        _       => vec![(query.to_vec(), 1), (reverse_complement(query), -1)],
    };

    let ka = blastn_ka_params(
        params.match_score, params.mismatch,
        params.gap_open, params.gap_extend,
    );

    let db_len = db.volume_length();
    let num_seqs = db.num_sequences() as u64;
    let (eff_query_len, eff_db_len) = ka.effective_lengths(query.len(), db_len, num_seqs);

    use std::collections::HashMap;
    let mut merged: HashMap<u32, SearchResult> = HashMap::new();

    for (q_seq, q_frame) in &queries_to_search {
        let lookup = NucleotideLookup::build(q_seq, params.word_size);

        let oids: Vec<u32> = (0..db.num_sequences()).collect();

        let strand_results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
            let subject = match db.get_sequence_nucleotide(oid) {
                Ok(s) => s,
                Err(_) => return None,
            };
            if subject.is_empty() { return None; }

            let mut hsps = search_one_nucleotide(
                q_seq,
                &subject,
                &lookup,
                &ka,
                params,
                eff_query_len,
                eff_db_len,
            );

            // Set query_frame on HSPs so output can distinguish strands
            for hsp in &mut hsps {
                hsp.query_frame = *q_frame;
            }

            // Apply max_hsps limit per subject
            if let Some(max) = params.max_hsps {
                hsps.truncate(max);
            }

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

        for result in strand_results {
            merged.entry(result.subject_oid)
                .and_modify(|e| e.hsps.extend(result.hsps.iter().cloned()))
                .or_insert(result);
        }
    }

    let mut results: Vec<SearchResult> = merged.into_values().map(|mut r| {
        r.hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
        if let Some(max) = params.max_hsps {
            r.hsps.truncate(max);
        }
        r
    }).collect();

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    results.truncate(params.max_target_seqs);
    results
}

fn search_one_nucleotide(
    query: &[u8],
    subject: &[u8],
    lookup: &NucleotideLookup,
    ka: &KarlinAltschul,
    params: &SearchParams,
    eff_query_len: usize,
    eff_db_len: u64,
) -> Vec<Hsp> {
    let diag_offset = query.len();
    let diag_len = query.len() + subject.len() + 1;
    let mut ungapped_hits = Vec::new();

    if params.two_hit {
        // 2-hit mode for nucleotide
        let mut diag_last: Vec<i32> = vec![i32::MIN; diag_len];
        let mut diag_extended: Vec<bool> = vec![false; diag_len];

        for (q_pos, s_pos) in lookup.scan_subject(subject) {
            let q_pos = q_pos as usize;
            let s_pos = s_pos as usize;
            let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
            if diag >= diag_extended.len() { continue; }
            if diag_extended[diag] { continue; }

            let prev = diag_last[diag];
            if prev == i32::MIN {
                diag_last[diag] = s_pos as i32;
                continue;
            }

            let dist = (s_pos as i32) - prev;
            if dist > params.two_hit_window as i32 {
                diag_last[diag] = s_pos as i32;
                continue;
            }

            // 2nd hit within window -> extend
            let hit = ungapped_extend_nucleotide(
                query, subject,
                q_pos, s_pos,
                params.match_score, params.mismatch,
                params.x_drop_ungapped,
            );
            if hit.score > 0 {
                diag_extended[diag] = true;
                ungapped_hits.push(hit);
            }
        }
    } else {
        // Single-hit mode (original behavior)
        let mut diag_hit: Vec<bool> = vec![false; diag_len];

        for (q_pos, s_pos) in lookup.scan_subject(subject) {
            let q_pos = q_pos as usize;
            let s_pos = s_pos as usize;
            let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
            if diag < diag_hit.len() && diag_hit[diag] { continue; }

            let hit = ungapped_extend_nucleotide(
                query, subject,
                q_pos, s_pos,
                params.match_score, params.mismatch,
                params.x_drop_ungapped,
            );
            if hit.score > 0 {
                if diag < diag_hit.len() { diag_hit[diag] = true; }
                ungapped_hits.push(hit);
            }
        }
    }

    if ungapped_hits.is_empty() { return vec![]; }
    ungapped_hits.sort_by(|a, b| b.score.cmp(&a.score));

    // Gapped extension with nucleotide scoring
    // Convert query/subject to 2-bit encoding for the scoring matrix (indices 0-3)
    let nt_matrix = build_nt_matrix(params.match_score, params.mismatch);
    let query_2bit: Vec<u8> = query.iter().map(|&b| nt_to_2bit(b)).collect();
    let subject_2bit: Vec<u8> = subject.iter().map(|&b| nt_to_2bit(b)).collect();
    let mut hsps = Vec::new();
    let mut covered: Vec<bool> = vec![false; query.len()];

    for uh in ungapped_hits {
        let center_q = (uh.q_start + uh.q_end) / 2;
        if center_q < covered.len() && covered[center_q] { continue; }
        let center_s = (uh.s_start + uh.s_end) / 2;

        let gh = gapped_extend(
            &query_2bit, &subject_2bit,
            center_q, center_s,
            &nt_matrix,
            params.gap_open, params.gap_extend,
            params.x_drop_gapped,
        );
        if gh.score <= 0 { continue; }

        let evalue = ka.evalue(gh.score, eff_query_len, eff_db_len);
        if evalue > params.evalue_threshold { continue; }

        let bit_score = ka.bit_score(gh.score);

        for i in gh.q_start.min(covered.len())..gh.q_end.min(covered.len()) {
            covered[i] = true;
        }

        // Convert 2-bit alignment bytes back to ASCII nucleotides
        let query_aln = gh.query_aln.iter().map(|&b| bit2_to_ascii(b)).collect();
        let subject_aln = gh.subject_aln.iter().map(|&b| bit2_to_ascii(b)).collect();

        hsps.push(Hsp {
            score: gh.score, bit_score, evalue,
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

/// Convert a 2-bit encoded nucleotide (0-3) back to ASCII, preserving gap characters.
fn bit2_to_ascii(b: u8) -> u8 {
    match b {
        0 => b'A', 1 => b'C', 2 => b'G', 3 => b'T',
        b'-' => b'-',
        _ => b'N',
    }
}

/// Build a simple nucleotide scoring matrix (28×28, Ncbistdaa-compatible indices unused;
/// here we use ASCII-indexed scoring for nucleotides).
fn build_nt_matrix(match_score: i32, mismatch: i32) -> ScoringMatrix {
    let mut scores = [[mismatch; 28]; 28];
    // For nucleotide, we index by ASCII value mod 28 — this is a hack.
    // Instead, we use the matrix purely through score() which compares indices.
    // We'll build a separate approach: treat a/c/g/t as 0/1/2/3 etc.
    // Actually we already have a nucleotide-aware extend function, but gapped_extend
    // uses ScoringMatrix. Let's build an ASCII-aware matrix.
    // We'll use the actual byte values as indices with a larger array.
    // Since 28 is too small for ASCII, we'll use a trick: map through nt_to_2bit.
    // For simplicity, just use the ScoringMatrix with ncbi2bit encoding in query/subject.
    // Set diagonal entries for 2-bit codes 0,1,2,3 (A,C,G,T):
    for (i, row) in scores.iter_mut().enumerate().take(4) {
        for (j, cell) in row.iter_mut().enumerate().take(4) {
            *cell = if i == j { match_score } else { mismatch };
        }
    }
    ScoringMatrix {
        name: MatrixType::Blosum62, // doesn't matter for nt
        scores,
        min_score: mismatch,
    }
}

// ─── Helper: ascii amino acid → Ncbistdaa ────────────────────────────────────

pub fn aa_to_ncbistdaa(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&c| match c.to_ascii_uppercase() {
        b'A' => 1,  b'B' => 2,  b'C' => 3,  b'D' => 4,  b'E' => 5,
        b'F' => 6,  b'G' => 7,  b'H' => 8,  b'I' => 9,  b'K' => 10,
        b'L' => 11, b'M' => 12, b'N' => 13, b'P' => 14, b'Q' => 15,
        b'R' => 16, b'S' => 17, b'T' => 18, b'V' => 19, b'W' => 20,
        b'X' => 21, b'Y' => 22, b'Z' => 23, _ => 21,
    }).collect()
}

// ─── BLASTX ──────────────────────────────────────────────────────────────────

/// Run a BLASTX search (translate nucleotide query in 6 frames, search protein database).
///
/// Query coordinates in `Hsp` are nucleotide-based (0-indexed, forward strand).
/// `hsp.query_frame` holds the frame (+1/+2/+3/-1/-2/-3).
pub fn blastx_search(
    db: &BlastDb,
    nt_query: &[u8], // ASCII nucleotide
    params: &SearchParams,
) -> Vec<SearchResult> {
    // Translate in 6 frames using the specified query genetic code
    let frames = six_frame_translate_with_code(nt_query, params.query_gencode);

    // We must merge results across frames: same subject OID may appear multiple times.
    use std::collections::HashMap;
    let mut merged: HashMap<u32, SearchResult> = HashMap::new();

    for frame in &frames {
        // Filter frames based on strand selection
        let skip = match params.strand.as_str() {
            "plus"  => frame.frame < 0,
            "minus" => frame.frame > 0,
            _ => false,
        };
        if skip { continue; }

        let prot_ascii = strip_stops(&frame.protein);
        if prot_ascii.len() < params.word_size { continue; }

        // Apply SEG masking
        let mut prot_masked = prot_ascii.clone();
        if params.filter_low_complexity {
            crate::mask::apply_seg(&mut prot_masked);
        }
        let query_ncbi = aa_to_ncbistdaa(&prot_masked);

        let mut frame_params = params.clone();
        frame_params.filter_low_complexity = false; // already applied
        frame_params.comp_adjust = false; // apply below with actual coords

        let frame_results = blast_search(db, &query_ncbi, &frame_params);

        for mut result in frame_results {
            // Convert query amino acid coordinates → nucleotide coordinates
            for hsp in &mut result.hsps {
                let (nt_start, nt_end) = frame.aa_to_nt(hsp.query_start, hsp.query_end);
                hsp.query_start = nt_start;
                hsp.query_end   = nt_end;
                hsp.query_frame = frame.frame;
            }

            merged.entry(result.subject_oid)
                .and_modify(|e| e.hsps.extend(result.hsps.iter().cloned()))
                .or_insert(result);
        }
    }

    let mut results: Vec<SearchResult> = merged.into_values().map(|mut r| {
        r.hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
        if let Some(max) = params.max_hsps {
            r.hsps.truncate(max);
        }
        r
    }).collect();

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    results.truncate(params.max_target_seqs);
    results
}

// ─── TBLASTN ─────────────────────────────────────────────────────────────────

/// Run a TBLASTN search (protein query against translated nucleotide database).
///
/// Subject coordinates in `Hsp` are nucleotide-based (0-indexed, forward strand).
/// `hsp.subject_frame` holds the subject frame (+1/+2/+3/-1/-2/-3).
pub fn tblastn_search(
    db: &BlastDb,
    query: &[u8], // ASCII amino acid query
    params: &SearchParams,
) -> Vec<SearchResult> {
    // Apply SEG masking to protein query
    let mut query_work = query.to_vec();
    if params.filter_low_complexity {
        crate::mask::apply_seg(&mut query_work);
    }
    let query_ncbi = aa_to_ncbistdaa(&query_work);

    let matrix = ScoringMatrix::from_type(params.matrix);
    let gap = GapPenalty::new(params.gap_open, params.gap_extend);
    let ka = match lookup_ka_params(params.matrix, gap) {
        Some(k) => k,
        None => lookup_ka_params(MatrixType::Blosum62, GapPenalty::blosum62_default()).unwrap(),
    };

    // DB length in amino acid units (≈ nt_len / 3 × 6 frames)
    let db_aa_len = db.volume_length() * 2; // 6 frames × ~1/3 = ×2
    let num_seqs = db.num_sequences() as u64;
    let (eff_q, eff_db) = ka.effective_lengths(query_ncbi.len(), db_aa_len, num_seqs);

    let threshold = neighbor_threshold(params.matrix, params.word_size);
    let lookup = ProteinLookup::build(&query_ncbi, params.word_size, &matrix, threshold);

    let oids: Vec<u32> = (0..db.num_sequences()).collect();

    let mut results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
        let nt_subject = match db.get_sequence_nucleotide(oid) {
            Ok(s) => s,
            Err(_) => return None,
        };
        if nt_subject.len() < 3 { return None; }

        let frames = six_frame_translate_with_code(&nt_subject, params.db_gencode);
        let mut all_hsps: Vec<Hsp> = Vec::new();

        for frame in &frames {
            let subj_prot = strip_stops(&frame.protein);
            if subj_prot.len() < params.word_size { continue; }

            let mut subj_masked = subj_prot.clone();
            if params.filter_low_complexity {
                crate::mask::apply_seg(&mut subj_masked);
            }
            let subj_ncbi = aa_to_ncbistdaa(&subj_masked);

            let mut hsps = search_one_protein(
                &query_ncbi, &subj_ncbi, &lookup, &matrix, &ka,
                params, eff_q, eff_db,
            );

            // Convert subject amino acid coords → nucleotide coords
            for hsp in &mut hsps {
                let (nt_start, nt_end) = frame.aa_to_nt(hsp.subject_start, hsp.subject_end);
                hsp.subject_start = nt_start;
                hsp.subject_end   = nt_end;
                hsp.subject_frame = frame.frame;
            }
            all_hsps.extend(hsps);
        }

        if all_hsps.is_empty() { return None; }
        all_hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
        if let Some(max) = params.max_hsps {
            all_hsps.truncate(max);
        }

        let header = db.get_header(oid).unwrap_or_default();
        Some(SearchResult {
            subject_oid: oid,
            subject_title: header.title,
            subject_accession: header.accession,
            subject_len: nt_subject.len(), // report nucleotide length
            hsps: all_hsps,
        })
    }).collect();

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    results.truncate(params.max_target_seqs);
    results
}

// ─── TBLASTX ─────────────────────────────────────────────────────────────────

/// Run a TBLASTX search (translate both query and database, ungapped).
///
/// Both `hsp.query_frame` and `hsp.subject_frame` are set.
/// Coordinates are in nucleotide space (0-based, forward strand).
/// Note: TBLASTX is inherently expensive (6 × 6 frame combinations per OID).
pub fn tblastx_search(
    db: &BlastDb,
    nt_query: &[u8], // ASCII nucleotide query
    params: &SearchParams,
) -> Vec<SearchResult> {
    let query_frames = six_frame_translate_with_code(nt_query, params.query_gencode);

    let matrix = ScoringMatrix::from_type(params.matrix);
    let gap = GapPenalty::new(params.gap_open, params.gap_extend);
    let ka = match lookup_ka_params(params.matrix, gap) {
        Some(k) => k,
        None => lookup_ka_params(MatrixType::Blosum62, GapPenalty::blosum62_default()).unwrap(),
    };

    let db_aa_len = db.volume_length() * 2;
    let num_seqs = db.num_sequences() as u64;

    // Build lookup tables for query frames (filtered by strand selection)
    let q_frames_ncbi: Vec<(TranslatedFrame, Vec<u8>, ProteinLookup)> = query_frames.iter()
        .filter_map(|frame| {
            // Filter frames based on strand selection
            let skip = match params.strand.as_str() {
                "plus"  => frame.frame < 0,
                "minus" => frame.frame > 0,
                _ => false,
            };
            if skip { return None; }

            let prot = strip_stops(&frame.protein);
            if prot.len() < params.word_size { return None; }
            let mut masked = prot.clone();
            if params.filter_low_complexity { crate::mask::apply_seg(&mut masked); }
            let ncbi = aa_to_ncbistdaa(&masked);
            let threshold = neighbor_threshold(params.matrix, params.word_size);
            let lookup = ProteinLookup::build(&ncbi, params.word_size, &matrix, threshold);
            Some((frame.clone(), ncbi, lookup))
        })
        .collect();

    let oids: Vec<u32> = (0..db.num_sequences()).collect();

    let mut results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
        let nt_subject = match db.get_sequence_nucleotide(oid) {
            Ok(s) => s,
            Err(_) => return None,
        };
        if nt_subject.len() < 3 { return None; }

        let subj_frames = six_frame_translate_with_code(&nt_subject, params.db_gencode);
        let mut all_hsps: Vec<Hsp> = Vec::new();

        for (q_frame, q_ncbi, q_lookup) in &q_frames_ncbi {
            let (eff_q, eff_db) = ka.effective_lengths(q_ncbi.len(), db_aa_len, num_seqs);

            for s_frame in &subj_frames {
                let subj_prot = strip_stops(&s_frame.protein);
                if subj_prot.len() < params.word_size { continue; }
                let mut subj_masked = subj_prot.clone();
                if params.filter_low_complexity { crate::mask::apply_seg(&mut subj_masked); }
                let subj_ncbi = aa_to_ncbistdaa(&subj_masked);

                let mut hsps = search_one_protein(
                    q_ncbi, &subj_ncbi, q_lookup, &matrix, &ka, params, eff_q, eff_db,
                );

                for hsp in &mut hsps {
                    let (qs, qe) = q_frame.aa_to_nt(hsp.query_start, hsp.query_end);
                    let (ss, se) = s_frame.aa_to_nt(hsp.subject_start, hsp.subject_end);
                    hsp.query_start   = qs; hsp.query_end   = qe;
                    hsp.subject_start = ss; hsp.subject_end = se;
                    hsp.query_frame   = q_frame.frame;
                    hsp.subject_frame = s_frame.frame;
                }
                all_hsps.extend(hsps);
            }
        }

        if all_hsps.is_empty() { return None; }
        all_hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
        if let Some(max) = params.max_hsps {
            all_hsps.truncate(max);
        }

        let header = db.get_header(oid).unwrap_or_default();
        Some(SearchResult {
            subject_oid: oid,
            subject_title: header.title,
            subject_accession: header.accession,
            subject_len: nt_subject.len(),
            hsps: all_hsps,
        })
    }).collect();

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    results.truncate(params.max_target_seqs);
    results
}

/// Run a discontiguous megablast search.
pub fn dc_megablast_search(
    db: &BlastDb,
    query: &[u8],
    params: &SearchParams,
    template_type: u8,
    template_length: usize,
) -> Vec<SearchResult> {
    let ka = blastn_ka_params(params.match_score, params.mismatch, params.gap_open, params.gap_extend);
    let db_len = db.volume_length();
    let num_seqs = db.num_sequences() as u64;
    let (eff_query_len, eff_db_len) = ka.effective_lengths(query.len(), db_len, num_seqs);

    let dc_lookup = DiscontiguousLookup::build(query, template_type, template_length);

    let oids: Vec<u32> = (0..db.num_sequences()).collect();

    let mut results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
        let subject = match db.get_sequence_nucleotide(oid) {
            Ok(s) => s,
            Err(_) => return None,
        };
        if subject.is_empty() { return None; }

        let seed_hits = dc_lookup.scan_subject(&subject);
        if seed_hits.is_empty() { return None; }

        let nt_matrix = build_nt_matrix(params.match_score, params.mismatch);
        let query_2bit: Vec<u8> = query.iter().map(|&b| nt_to_2bit(b)).collect();
        let subject_2bit: Vec<u8> = subject.iter().map(|&b| nt_to_2bit(b)).collect();
        let mut hsps = Vec::new();
        let diag_offset = query.len();
        let mut diag_hit = vec![false; query.len() + subject.len() + 1];

        for (q_pos, s_pos) in seed_hits {
            let q_pos = q_pos as usize;
            let s_pos = s_pos as usize;
            let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
            if diag < diag_hit.len() && diag_hit[diag] { continue; }

            let hit = ungapped_extend_nucleotide(query, &subject, q_pos, s_pos,
                params.match_score, params.mismatch, params.x_drop_ungapped);
            if hit.score > 0 {
                if diag < diag_hit.len() { diag_hit[diag] = true; }

                let center_q = (hit.q_start + hit.q_end) / 2;
                let center_s = (hit.s_start + hit.s_end) / 2;
                let gh = gapped_extend(&query_2bit, &subject_2bit, center_q, center_s, &nt_matrix,
                    params.gap_open, params.gap_extend, params.x_drop_gapped);
                if gh.score <= 0 { continue; }

                let evalue = ka.evalue(gh.score, eff_query_len, eff_db_len);
                if evalue > params.evalue_threshold { continue; }

                let query_aln = gh.query_aln.iter().map(|&b| bit2_to_ascii(b)).collect();
                let subject_aln = gh.subject_aln.iter().map(|&b| bit2_to_ascii(b)).collect();

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
        }

        if hsps.is_empty() { return None; }
        hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));

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

impl SearchParams {
    /// Save search parameters to a JSON-like text format.
    pub fn save_strategy<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<()> {
        writeln!(out, "blast_search_strategy_v1")?;
        writeln!(out, "word_size={}", self.word_size)?;
        writeln!(out, "matrix={:?}", self.matrix)?;
        writeln!(out, "gap_open={}", self.gap_open)?;
        writeln!(out, "gap_extend={}", self.gap_extend)?;
        writeln!(out, "evalue_threshold={}", self.evalue_threshold)?;
        writeln!(out, "max_target_seqs={}", self.max_target_seqs)?;
        writeln!(out, "x_drop_ungapped={}", self.x_drop_ungapped)?;
        writeln!(out, "x_drop_gapped={}", self.x_drop_gapped)?;
        writeln!(out, "x_drop_final={}", self.x_drop_final)?;
        writeln!(out, "match_score={}", self.match_score)?;
        writeln!(out, "mismatch={}", self.mismatch)?;
        writeln!(out, "filter_low_complexity={}", self.filter_low_complexity)?;
        writeln!(out, "comp_adjust={}", self.comp_adjust)?;
        writeln!(out, "two_hit={}", self.two_hit)?;
        writeln!(out, "two_hit_window={}", self.two_hit_window)?;
        writeln!(out, "strand={}", self.strand)?;
        writeln!(out, "query_gencode={}", self.query_gencode)?;
        writeln!(out, "db_gencode={}", self.db_gencode)?;
        writeln!(out, "soft_masking={}", self.soft_masking)?;
        Ok(())
    }

    /// Load search parameters from a saved strategy file.
    pub fn load_strategy<R: std::io::BufRead>(reader: &mut R) -> std::io::Result<Self> {
        let mut params = SearchParams::blastp_defaults();
        let mut first = true;

        for line in reader.lines() {
            let line = line?;
            let line = line.trim().to_string();
            if line.is_empty() { continue; }

            if first {
                if line != "blast_search_strategy_v1" {
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                        "Not a valid search strategy file"));
                }
                first = false;
                continue;
            }

            let mut parts = line.splitn(2, '=');
            let key = parts.next().unwrap_or("");
            let val = parts.next().unwrap_or("");

            match key {
                "word_size" => { params.word_size = val.parse().unwrap_or(params.word_size); }
                "matrix" => {
                    params.matrix = match val {
                        "Blosum45" => MatrixType::Blosum45,
                        "Blosum50" => MatrixType::Blosum50,
                        "Blosum62" => MatrixType::Blosum62,
                        "Blosum80" => MatrixType::Blosum80,
                        "Blosum90" => MatrixType::Blosum90,
                        "Pam30" => MatrixType::Pam30,
                        "Pam70" => MatrixType::Pam70,
                        "Pam250" => MatrixType::Pam250,
                        _ => params.matrix,
                    };
                }
                "gap_open" => { params.gap_open = val.parse().unwrap_or(params.gap_open); }
                "gap_extend" => { params.gap_extend = val.parse().unwrap_or(params.gap_extend); }
                "evalue_threshold" => { params.evalue_threshold = val.parse().unwrap_or(params.evalue_threshold); }
                "max_target_seqs" => { params.max_target_seqs = val.parse().unwrap_or(params.max_target_seqs); }
                "x_drop_ungapped" => { params.x_drop_ungapped = val.parse().unwrap_or(params.x_drop_ungapped); }
                "x_drop_gapped" => { params.x_drop_gapped = val.parse().unwrap_or(params.x_drop_gapped); }
                "x_drop_final" => { params.x_drop_final = val.parse().unwrap_or(params.x_drop_final); }
                "match_score" => { params.match_score = val.parse().unwrap_or(params.match_score); }
                "mismatch" => { params.mismatch = val.parse().unwrap_or(params.mismatch); }
                "filter_low_complexity" => { params.filter_low_complexity = val.parse().unwrap_or(params.filter_low_complexity); }
                "comp_adjust" => { params.comp_adjust = val.parse().unwrap_or(params.comp_adjust); }
                "two_hit" => { params.two_hit = val.parse().unwrap_or(params.two_hit); }
                "two_hit_window" => { params.two_hit_window = val.parse().unwrap_or(params.two_hit_window); }
                "strand" => { params.strand = val.to_string(); }
                "query_gencode" => { params.query_gencode = val.parse().unwrap_or(params.query_gencode); }
                "db_gencode" => { params.db_gencode = val.parse().unwrap_or(params.db_gencode); }
                "soft_masking" => { params.soft_masking = val.parse().unwrap_or(params.soft_masking); }
                _ => {} // ignore unknown keys
            }
        }
        Ok(params)
    }
}

