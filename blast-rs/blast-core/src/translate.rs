//! Genetic code, reverse complement, and six-frame translation.

/// 2-bit encoding: A=0, C=1, G=2, T=3, other=4 (ambiguous).
#[inline]
pub fn nt_to_2bit(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' | b'U' => 3, _ => 4,
    }
}

/// Complement of a 2-bit-encoded base.
#[inline]
fn complement_2bit(b: u8) -> u8 {
    match b { 0 => 3, 1 => 2, 2 => 1, 3 => 0, _ => 4 }
}

/// Complement of an ASCII nucleotide character.
pub fn complement_ascii(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        b'R' => b'Y', b'Y' => b'R', b'M' => b'K', b'K' => b'M',
        b'S' => b'S', b'W' => b'W', b'H' => b'D', b'D' => b'H',
        b'B' => b'V', b'V' => b'B', b'N' => b'N', _ => b'N',
    }
}

/// Reverse complement of an ASCII nucleotide sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| complement_ascii(c)).collect()
}

/// Standard genetic code indexed by 2-bit codon: A=0, C=1, G=2, T=3.
/// index = b1*16 + b2*4 + b3.
///
/// ```text
/// Rows (b1 b2): AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
/// ```
static CODON_TABLE: [u8; 64] = [
    // A__ (0-15)
    b'K', b'N', b'K', b'N',  // AA{A,C,G,T}
    b'T', b'T', b'T', b'T',  // AC{A,C,G,T}
    b'R', b'S', b'R', b'S',  // AG{A,C,G,T}
    b'I', b'I', b'M', b'I',  // AT{A,C,G,T}
    // C__ (16-31)
    b'Q', b'H', b'Q', b'H',  // CA{A,C,G,T}
    b'P', b'P', b'P', b'P',  // CC{A,C,G,T}
    b'R', b'R', b'R', b'R',  // CG{A,C,G,T}
    b'L', b'L', b'L', b'L',  // CT{A,C,G,T}
    // G__ (32-47)
    b'E', b'D', b'E', b'D',  // GA{A,C,G,T}
    b'A', b'A', b'A', b'A',  // GC{A,C,G,T}
    b'G', b'G', b'G', b'G',  // GG{A,C,G,T}
    b'V', b'V', b'V', b'V',  // GT{A,C,G,T}
    // T__ (48-63)
    b'*', b'Y', b'*', b'Y',  // TA{A,C,G,T}
    b'S', b'S', b'S', b'S',  // TC{A,C,G,T}
    b'*', b'C', b'W', b'C',  // TG{A,C,G,T}
    b'L', b'F', b'L', b'F',  // TT{A,C,G,T}
];

/// Translate a single codon from 2-bit encoded bases.
/// Returns `b'X'` for ambiguous/unknown bases, `b'*'` for stop codons.
#[inline]
fn translate_codon(b1: u8, b2: u8, b3: u8) -> u8 {
    if b1 >= 4 || b2 >= 4 || b3 >= 4 { return b'X'; }
    CODON_TABLE[(b1 as usize) * 16 + (b2 as usize) * 4 + (b3 as usize)]
}

/// Translate a nucleotide sequence starting at `offset`, yielding ASCII amino acids.
/// Stop codons (`*`) are included in the output.
pub fn translate_from(seq: &[u8], offset: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity((seq.len().saturating_sub(offset)) / 3 + 1);
    let mut i = offset;
    while i + 2 < seq.len() {
        out.push(translate_codon(nt_to_2bit(seq[i]), nt_to_2bit(seq[i+1]), nt_to_2bit(seq[i+2])));
        i += 3;
    }
    out
}

/// A single reading frame's translated protein.
#[derive(Debug, Clone)]
pub struct TranslatedFrame {
    /// Frame number: +1, +2, +3 (forward); -1, -2, -3 (reverse complement).
    pub frame: i32,
    /// ASCII amino acid sequence (may contain `*` for stop codons).
    pub protein: Vec<u8>,
    /// Length of the source nucleotide sequence.
    pub nt_len: usize,
}

impl TranslatedFrame {
    /// Convert amino acid coordinates (0-based, exclusive end) back to
    /// nucleotide coordinates in the **original** (forward-strand) sequence.
    ///
    /// Returns `(nt_start, nt_end)` both 0-based, start inclusive, end exclusive.
    /// For negative frames `nt_start < nt_end` still holds (we return the lower bound first).
    pub fn aa_to_nt(&self, aa_start: usize, aa_end: usize) -> (usize, usize) {
        let offset = (self.frame.unsigned_abs() as usize) - 1;
        if self.frame > 0 {
            let s = offset + aa_start * 3;
            let e = offset + aa_end   * 3;
            (s, e)
        } else {
            // Coords in rev-comp space → map back to original strand.
            // revcomp position p → original position nt_len - 1 - p
            // revcomp range [aa_start*3+offset .. aa_end*3+offset)
            //   → original [nt_len - aa_end*3 - offset .. nt_len - aa_start*3 - offset)
            let s = self.nt_len.saturating_sub(offset + aa_end   * 3);
            let e = self.nt_len.saturating_sub(offset + aa_start * 3);
            (s, e)
        }
    }
}

/// Translate all 6 reading frames of a nucleotide sequence.
pub fn six_frame_translate(nt_seq: &[u8]) -> [TranslatedFrame; 6] {
    let nt_len = nt_seq.len();
    let rc = reverse_complement(nt_seq);
    [
        TranslatedFrame { frame:  1, protein: translate_from(nt_seq, 0), nt_len },
        TranslatedFrame { frame:  2, protein: translate_from(nt_seq, 1), nt_len },
        TranslatedFrame { frame:  3, protein: translate_from(nt_seq, 2), nt_len },
        TranslatedFrame { frame: -1, protein: translate_from(&rc,    0), nt_len },
        TranslatedFrame { frame: -2, protein: translate_from(&rc,    1), nt_len },
        TranslatedFrame { frame: -3, protein: translate_from(&rc,    2), nt_len },
    ]
}

/// Remove stop codons from a translated protein for use as search query/subject.
/// Splits on `*` and keeps the longest segment (or truncates at the first stop).
pub fn strip_stops(protein: &[u8]) -> Vec<u8> {
    // Return everything up to the first stop (mirrors BLAST behaviour for queries).
    match protein.iter().position(|&b| b == b'*') {
        Some(pos) => protein[..pos].to_vec(),
        None => protein.to_vec(),
    }
}
