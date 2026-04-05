//! # blast-rs
//!
//! A pure-Rust implementation of NCBI BLAST (Basic Local Alignment Search Tool).
//!
//! Provides protein and nucleotide sequence search against BLAST databases,
//! matching NCBI BLAST+ in correctness and performance. No dependency on the
//! NCBI C++ toolkit.
//!
//! ## Quick start
//!
//! ```no_run
//! use blast_rs::{BlastDb, SearchParams, blastp, parse_fasta};
//!
//! let db = BlastDb::open("mydb".as_ref()).unwrap();
//! let fasta = std::fs::read("query.faa").unwrap();
//! let sequences = parse_fasta(&fasta);
//!
//! let params = SearchParams::blastp()
//!     .evalue(1e-5)
//!     .max_target_seqs(50)
//!     .num_threads(8);
//!
//! for (title, seq) in &sequences {
//!     let results = blastp(&db, seq, &params);
//!     for r in &results {
//!         println!("{}\t{}\t{:.2e}", title, r.subject_accession, r.best_evalue());
//!     }
//! }
//! ```
//!
//! ## Search modes
//!
//! | Function | Query | Database | Description |
//! |----------|-------|----------|-------------|
//! | [`blastp`] | Protein | Protein | Protein-protein search |
//! | [`blastn`] | Nucleotide | Nucleotide | Nucleotide-nucleotide search |
//! | [`blastx`] | Nucleotide | Protein | Translated query (6 frames) |
//! | [`tblastn`] | Protein | Nucleotide | Translated database (6 frames) |
//! | [`tblastx`] | Nucleotide | Nucleotide | Both sides translated |
//! | [`psiblast`] | Protein | Protein | Iterative search with PSSM |
//!
//! ## Database I/O
//!
//! Read existing BLAST databases (v4 and v5) or create new ones:
//!
//! ```no_run
//! use blast_rs::{BlastDbBuilder, SequenceEntry};
//! use blast_rs::db::index::SeqType;
//! use std::path::Path;
//!
//! let mut builder = BlastDbBuilder::new(SeqType::Protein, "My database");
//! builder.add(SequenceEntry {
//!     title: "Human insulin".into(),
//!     accession: "P01308".into(),
//!     sequence: b"MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT".to_vec(),
//!     taxid: Some(9606),
//! });
//! builder.write(Path::new("mydb")).unwrap();
//! ```

pub mod db;
pub mod matrix;
pub mod stats;
mod tests;
pub mod lookup;
pub mod extend;
pub mod align;
pub mod search;
pub mod hsp;
pub mod translate;
pub mod mask;
pub mod pssm;
pub mod compo;
pub mod api;

// ── Primary search API ──────────────────────────────────────────────────────

pub use api::{blastp, blastn, blastx, tblastn, tblastx, psiblast, PsiblastParams, parse_fasta};

// ── Search configuration ────────────────────────────────────────────────────

pub use search::SearchParams;
pub use matrix::{ScoringMatrix, MatrixType};
pub use stats::KarlinAltschul;

// ── Results ─────────────────────────────────────────────────────────────────

pub use hsp::{Hsp, SearchResult};

// ── Database I/O ────────────────────────────────────────────────────────────

pub use crate::db::{BlastDb, BlastDefLine, BlastDbBuilder, SequenceEntry, TaxDb};

// ── Sequence utilities ──────────────────────────────────────────────────────

pub use translate::{six_frame_translate, six_frame_translate_with_code, reverse_complement, TranslatedFrame, get_codon_table};
pub use mask::{apply_dust, apply_seg, apply_seg_ncbistdaa, apply_repeat_mask, repeat_mask,
               apply_lowercase_mask_protein, apply_lowercase_mask_nucleotide, lowercase_mask};

// ── Advanced / low-level ────────────────────────────────────────────────────

pub use search::{blast_search, blastn_search, blastx_search, tblastn_search, tblastx_search, dc_megablast_search};
pub use pssm::{Pssm, build_pssm, psiblast_search, search_with_pssm};
pub use compo::{composition_ncbistdaa, adjust_evalue, adjust_evalue_with_mode, BACKGROUND_FREQ};
