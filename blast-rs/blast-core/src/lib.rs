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

pub use matrix::{ScoringMatrix, MatrixType};
pub use stats::KarlinAltschul;
pub use hsp::{Hsp, SearchResult};
pub use search::{SearchParams, blast_search};
pub use translate::{six_frame_translate, six_frame_translate_with_code, reverse_complement, TranslatedFrame, get_codon_table};
pub use mask::{apply_dust, apply_seg, apply_seg_ncbistdaa};
pub use pssm::{Pssm, build_pssm, psiblast_search, search_with_pssm};
pub use compo::{composition_ncbistdaa, adjust_evalue, adjust_evalue_with_mode, BACKGROUND_FREQ};

// Top-level API functions
pub use api::{blastp, blastn, blastx, tblastn, tblastx, psiblast, PsiblastParams, parse_fasta};

// Low-level search functions for advanced users
pub use search::{blastn_search, blastx_search, tblastn_search, tblastx_search};

// Re-export BlastDb types so users need only one dependency
pub use blast_db::{BlastDb, BlastDefLine, BlastDbBuilder, SequenceEntry};
