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

/// Standard genetic code (code 1) indexed by 2-bit codon: A=0, C=1, G=2, T=3.
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

// ---------------------------------------------------------------------------
// NCBI genetic code tables (1-33)
// ---------------------------------------------------------------------------
//
// All tables are derived from the standard code (code 1) by applying the
// documented codon reassignments from NCBI.
//
// Codon index: b1*16 + b2*4 + b3 where A=0 C=1 G=2 T=3.
// Key indices:
//   AAA=0   AGA=8   AGG=10  ATA=12
//   CTN=28..31  CTG=30
//   TAA=48  TAG=50  TCA=52  TGA=56  TTA=60

/// Helper: build a codon table by starting from standard code and applying changes.
const fn make_table(changes: &[(usize, u8)]) -> [u8; 64] {
    let mut t = CODON_TABLE;
    let mut i = 0;
    while i < changes.len() {
        t[changes[i].0] = changes[i].1;
        i += 1;
    }
    t
}

/// Code 2: Vertebrate Mitochondrial
/// AGA=* AGG=* ATA=M TGA=W
static CODON_TABLE_2: [u8; 64] = make_table(&[
    (8, b'*'), (10, b'*'), (12, b'M'), (56, b'W'),
]);

/// Code 3: Yeast Mitochondrial
/// ATA=M CTN=T (28..31) TGA=W
static CODON_TABLE_3: [u8; 64] = make_table(&[
    (12, b'M'), (28, b'T'), (29, b'T'), (30, b'T'), (31, b'T'), (56, b'W'),
]);

/// Code 4: Mold, Protozoan, Coelenterate Mitochondrial
/// TGA=W
static CODON_TABLE_4: [u8; 64] = make_table(&[(56, b'W')]);

/// Code 5: Invertebrate Mitochondrial
/// AGA=S AGG=S ATA=M TGA=W
static CODON_TABLE_5: [u8; 64] = make_table(&[
    (8, b'S'), (10, b'S'), (12, b'M'), (56, b'W'),
]);

/// Code 6: Ciliate, Dasycladacean, Hexamita Nuclear
/// TAA=Q TAG=Q
static CODON_TABLE_6: [u8; 64] = make_table(&[(48, b'Q'), (50, b'Q')]);

/// Code 9: Echinoderm and Flatworm Mitochondrial
/// AAA=N AGA=S AGG=S TGA=W
static CODON_TABLE_9: [u8; 64] = make_table(&[
    (0, b'N'), (8, b'S'), (10, b'S'), (56, b'W'),
]);

/// Code 10: Euplotid Nuclear
/// TGA=C
static CODON_TABLE_10: [u8; 64] = make_table(&[(56, b'C')]);

// Code 11: Bacterial, Archaeal and Plant Plastid - same as standard code 1.

/// Code 12: Alternative Yeast Nuclear
/// CTG=S
static CODON_TABLE_12: [u8; 64] = make_table(&[(30, b'S')]);

/// Code 13: Ascidian Mitochondrial
/// AGA=G AGG=G ATA=M TGA=W
static CODON_TABLE_13: [u8; 64] = make_table(&[
    (8, b'G'), (10, b'G'), (12, b'M'), (56, b'W'),
]);

/// Code 14: Alternative Flatworm Mitochondrial
/// AAA=N AGA=S AGG=S TAA=Y TGA=W
static CODON_TABLE_14: [u8; 64] = make_table(&[
    (0, b'N'), (8, b'S'), (10, b'S'), (48, b'Y'), (56, b'W'),
]);

/// Code 15: Blepharisma Nuclear
/// TAG=Q
static CODON_TABLE_15: [u8; 64] = make_table(&[(50, b'Q')]);

/// Code 16: Chlorophycean Mitochondrial
/// TAG=L
static CODON_TABLE_16: [u8; 64] = make_table(&[(50, b'L')]);

/// Code 21: Trematode Mitochondrial
/// AAA=N AGA=S AGG=S ATA=M TGA=W
static CODON_TABLE_21: [u8; 64] = make_table(&[
    (0, b'N'), (8, b'S'), (10, b'S'), (12, b'M'), (56, b'W'),
]);

/// Code 22: Scenedesmus obliquus Mitochondrial
/// TCA=* TAG=L
static CODON_TABLE_22: [u8; 64] = make_table(&[(52, b'*'), (50, b'L')]);

/// Code 23: Thraustochytrium Mitochondrial
/// TTA=*
static CODON_TABLE_23: [u8; 64] = make_table(&[(60, b'*')]);

/// Code 24: Rhabdopleuridae Mitochondrial
/// AGA=S AGG=K TGA=W
static CODON_TABLE_24: [u8; 64] = make_table(&[
    (8, b'S'), (10, b'K'), (56, b'W'),
]);

/// Code 25: Candidate Division SR1 and Gracilibacteria
/// TGA=G
static CODON_TABLE_25: [u8; 64] = make_table(&[(56, b'G')]);

/// Code 26: Pachysolen tannophilus Nuclear
/// CTG=A
static CODON_TABLE_26: [u8; 64] = make_table(&[(30, b'A')]);

/// Code 27: Karyorelictea Nuclear
/// TAA=Q TAG=Q TGA=W
static CODON_TABLE_27: [u8; 64] = make_table(&[
    (48, b'Q'), (50, b'Q'), (56, b'W'),
]);

/// Code 29: Mesodinium Nuclear
/// TAA=Y TAG=Y
static CODON_TABLE_29: [u8; 64] = make_table(&[(48, b'Y'), (50, b'Y')]);

/// Code 30: Peritrich Nuclear
/// TAA=E TAG=E
static CODON_TABLE_30: [u8; 64] = make_table(&[(48, b'E'), (50, b'E')]);

/// Code 31: Blastocrithidia Nuclear
/// TGA=W TAG=E TAA=E
static CODON_TABLE_31: [u8; 64] = make_table(&[
    (48, b'E'), (50, b'E'), (56, b'W'),
]);

/// Code 33: Cephalodiscidae Mitochondrial
/// AGA=S AGG=K TAA=Y TGA=W
static CODON_TABLE_33: [u8; 64] = make_table(&[
    (8, b'S'), (10, b'K'), (48, b'Y'), (56, b'W'),
]);

/// Return the codon table for a given NCBI genetic code number (1-33).
///
/// Panics if the code is not a recognized NCBI genetic code.
pub fn get_codon_table(code: u8) -> &'static [u8; 64] {
    match code {
        1 | 11 => &CODON_TABLE,
        2  => &CODON_TABLE_2,
        3  => &CODON_TABLE_3,
        4  => &CODON_TABLE_4,
        5  => &CODON_TABLE_5,
        6  => &CODON_TABLE_6,
        9  => &CODON_TABLE_9,
        10 => &CODON_TABLE_10,
        12 => &CODON_TABLE_12,
        13 => &CODON_TABLE_13,
        14 => &CODON_TABLE_14,
        15 => &CODON_TABLE_15,
        16 => &CODON_TABLE_16,
        21 => &CODON_TABLE_21,
        22 => &CODON_TABLE_22,
        23 => &CODON_TABLE_23,
        24 => &CODON_TABLE_24,
        25 => &CODON_TABLE_25,
        26 => &CODON_TABLE_26,
        27 => &CODON_TABLE_27,
        29 => &CODON_TABLE_29,
        30 => &CODON_TABLE_30,
        31 => &CODON_TABLE_31,
        33 => &CODON_TABLE_33,
        _  => panic!("Unsupported NCBI genetic code: {}", code),
    }
}

/// Translate a single codon from 2-bit encoded bases using the standard code.
/// Returns `b'X'` for ambiguous/unknown bases, `b'*'` for stop codons.
#[inline]
fn translate_codon(b1: u8, b2: u8, b3: u8) -> u8 {
    if b1 >= 4 || b2 >= 4 || b3 >= 4 { return b'X'; }
    CODON_TABLE[(b1 as usize) * 16 + (b2 as usize) * 4 + (b3 as usize)]
}

/// Translate a single codon from 2-bit encoded bases using the given table.
/// Returns `b'X'` for ambiguous/unknown bases.
#[inline]
fn translate_codon_with_table(b1: u8, b2: u8, b3: u8, table: &[u8; 64]) -> u8 {
    if b1 >= 4 || b2 >= 4 || b3 >= 4 { return b'X'; }
    table[(b1 as usize) * 16 + (b2 as usize) * 4 + (b3 as usize)]
}

/// Translate a nucleotide sequence starting at `offset`, yielding ASCII amino acids.
/// Uses the standard genetic code (code 1).
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

/// Translate a nucleotide sequence starting at `offset` using the specified
/// NCBI genetic code. Stop codons (`*`) are included in the output.
pub fn translate_from_with_code(seq: &[u8], offset: usize, genetic_code: u8) -> Vec<u8> {
    let table = get_codon_table(genetic_code);
    let mut out = Vec::with_capacity((seq.len().saturating_sub(offset)) / 3 + 1);
    let mut i = offset;
    while i + 2 < seq.len() {
        out.push(translate_codon_with_table(
            nt_to_2bit(seq[i]), nt_to_2bit(seq[i+1]), nt_to_2bit(seq[i+2]), table,
        ));
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

/// Translate all 6 reading frames of a nucleotide sequence using the standard
/// genetic code (code 1).
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

/// Translate all 6 reading frames of a nucleotide sequence using the specified
/// NCBI genetic code.
pub fn six_frame_translate_with_code(nt_seq: &[u8], genetic_code: u8) -> [TranslatedFrame; 6] {
    let nt_len = nt_seq.len();
    let rc = reverse_complement(nt_seq);
    [
        TranslatedFrame { frame:  1, protein: translate_from_with_code(nt_seq, 0, genetic_code), nt_len },
        TranslatedFrame { frame:  2, protein: translate_from_with_code(nt_seq, 1, genetic_code), nt_len },
        TranslatedFrame { frame:  3, protein: translate_from_with_code(nt_seq, 2, genetic_code), nt_len },
        TranslatedFrame { frame: -1, protein: translate_from_with_code(&rc,    0, genetic_code), nt_len },
        TranslatedFrame { frame: -2, protein: translate_from_with_code(&rc,    1, genetic_code), nt_len },
        TranslatedFrame { frame: -3, protein: translate_from_with_code(&rc,    2, genetic_code), nt_len },
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
