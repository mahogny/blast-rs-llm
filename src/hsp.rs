//! HSP (High-Scoring Segment Pair) and SearchResult data structures.
//!
//! These are the core result types returned by all BLAST search functions.
//! A [`SearchResult`] contains one or more [`Hsp`]s for a single database subject.

/// A high-scoring segment pair — one local alignment between query and subject.
///
/// Coordinates are 0-based with exclusive end (Rust convention).
/// For translated searches, `query_frame` and `subject_frame` indicate the
/// reading frame used (+1/+2/+3 forward, -1/-2/-3 reverse, 0 = no translation).
#[derive(Debug, Clone)]
pub struct Hsp {
    /// Raw alignment score (matrix-dependent).
    pub score: i32,
    /// Normalized bit score: `(lambda * score - ln(K)) / ln(2)`.
    pub bit_score: f64,
    /// Expected number of alignments with this score by chance.
    pub evalue: f64,
    /// Query start position (0-based, inclusive).
    pub query_start: usize,
    /// Query end position (0-based, exclusive).
    pub query_end: usize,
    /// Subject start position (0-based, inclusive).
    pub subject_start: usize,
    /// Subject end position (0-based, exclusive).
    pub subject_end: usize,
    /// Number of identical residue pairs in the alignment.
    pub num_identities: usize,
    /// Total number of gap characters in query + subject alignment strings.
    pub num_gaps: usize,
    /// Length of the alignment including gaps.
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
    /// Percentage of identical residues: `100 * identities / alignment_length`.
    pub fn percent_identity(&self) -> f64 {
        if self.alignment_length == 0 { return 0.0; }
        100.0 * self.num_identities as f64 / self.alignment_length as f64
    }
}

/// All alignments (HSPs) for one database subject sequence.
///
/// Returned as part of a `Vec<SearchResult>` from search functions like [`crate::blastp`].
/// Results are sorted by best E-value. Each result may contain multiple HSPs
/// if the query aligns to different regions of the same subject.
#[derive(Debug, Clone)]
pub struct SearchResult {
    /// Ordinal ID in the database (0-based).
    pub subject_oid: u32,
    /// Full title/description from the database header.
    pub subject_title: String,
    /// Accession identifier (first word of the FASTA header).
    pub subject_accession: String,
    /// Length of the subject sequence in residues/bases.
    pub subject_len: usize,
    /// High-scoring segment pairs, sorted by E-value (best first).
    pub hsps: Vec<Hsp>,
    /// Taxonomy IDs for this subject (from v5 database .pot/.not files).
    /// Empty if database doesn't have taxonomy data.
    pub taxids: Vec<i32>,
}

impl SearchResult {
    /// Best E-value across all HSPs.
    pub fn best_evalue(&self) -> f64 {
        self.hsps.iter().map(|h| h.evalue).fold(f64::INFINITY, f64::min)
    }
}
