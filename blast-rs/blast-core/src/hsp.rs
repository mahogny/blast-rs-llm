//! HSP and SearchResult data structures.

/// A high-scoring segment pair (one alignment).
#[derive(Debug, Clone)]
pub struct Hsp {
    pub score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_start: usize,   // 0-based, inclusive
    pub query_end: usize,     // 0-based, exclusive
    pub subject_start: usize, // 0-based, inclusive
    pub subject_end: usize,   // 0-based, exclusive
    pub num_identities: usize,
    pub num_gaps: usize,
    pub alignment_length: usize,
    /// Query alignment string (with gaps)
    pub query_aln: Vec<u8>,
    /// Match-line string (| for identical, space for mismatch/gap)
    pub midline: Vec<u8>,
    /// Subject alignment string (with gaps)
    pub subject_aln: Vec<u8>,
    /// Reading frame of the query translation.
    /// +1/+2/+3 for forward strand, -1/-2/-3 for reverse complement.
    /// 0 for blastp/blastn (no translation).
    pub query_frame: i32,
    /// Reading frame of the subject translation.
    /// +1/+2/+3 for forward strand, -1/-2/-3 for reverse complement.
    /// 0 for blastp/blastn (no translation).
    pub subject_frame: i32,
}

impl Hsp {
    pub fn percent_identity(&self) -> f64 {
        if self.alignment_length == 0 { return 0.0; }
        100.0 * self.num_identities as f64 / self.alignment_length as f64
    }
}

/// All HSPs for one database sequence.
#[derive(Debug, Clone)]
pub struct SearchResult {
    pub subject_oid: u32,
    pub subject_title: String,
    pub subject_accession: String,
    pub subject_len: usize,
    pub hsps: Vec<Hsp>,
}

impl SearchResult {
    /// Best E-value across all HSPs.
    pub fn best_evalue(&self) -> f64 {
        self.hsps.iter().map(|h| h.evalue).fold(f64::INFINITY, f64::min)
    }
}
