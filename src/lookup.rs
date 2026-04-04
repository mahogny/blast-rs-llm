//! Lookup tables for k-mer seeding.

use crate::matrix::ScoringMatrix;

/// NCBI-style backbone cell: stores up to 3 hits inline to avoid
/// pointer dereference for the common case (single cache line access).
#[derive(Clone, Copy)]
struct BackboneCell {
    /// Number of hits: 0 = empty, 1-3 = inline, >3 = overflow
    num_used: u32,
    /// Inline storage for up to 3 query positions, or overflow cursor
    data: [u32; 3],
}

/// Protein lookup table using neighboring words.
///
/// Uses NCBI-style architecture:
/// - Shift-based indexing (charsize=5) for fast rolling hash
/// - Backbone cells with up to 3 inline hits (no pointer dereference)
/// - Overflow array for words with >3 hits
/// - Presence bitfield for quick empty-cell rejection
pub struct ProteinLookup {
    pub word_size: usize,
    pub charsize: u32,
    pub mask: u32,
    /// Presence bitfield for quick rejection (~4 KB, L1-resident)
    presence: Vec<u64>,
    /// Backbone: one cell per possible word code. Most accesses are single cache line.
    backbone: Vec<BackboneCell>,
    /// Overflow array for words with >3 hits
    overflow: Vec<u32>,
    #[allow(dead_code)]
    capacity: usize,
}

/// Encode a protein word using shift-based encoding (NCBI style).
/// charsize bits per residue, matching ComputeTableIndex.
#[inline]
pub fn encode_protein_word(residues: &[u8]) -> u32 {
    let mut code = 0u32;
    for &r in residues {
        code = (code << 5) | (r as u32 & 0x1F);
    }
    code
}

/// Constant: max hits stored inline per backbone cell (matches NCBI AA_HITS_PER_CELL)
const HITS_PER_CELL: usize = 3;

impl ProteinLookup {
    /// Build lookup table from query (Ncbistdaa encoded).
    pub fn build(query: &[u8], word_size: usize, matrix: &ScoringMatrix, threshold: i32) -> Self {
        let charsize = 5u32;
        let capacity = 1usize << (word_size as u32 * charsize);
        let mask = (capacity - 1) as u32;
        let qlen = query.len();

        let empty_cell = BackboneCell { num_used: 0, data: [0; 3] };

        if qlen < word_size {
            return ProteinLookup {
                word_size, charsize, mask,
                presence: vec![0u64; capacity.div_ceil(64)],
                backbone: vec![empty_cell; capacity],
                overflow: Vec::new(),
                capacity,
            };
        }

        // Phase 1: collect hits into temporary Vecs
        let mut tmp = vec![Vec::new(); capacity];
        for q_pos in 0..=(qlen - word_size) {
            let query_word = &query[q_pos..q_pos + word_size];
            enumerate_neighbors_shift(query_word, word_size, charsize, matrix, threshold, &mut |code| {
                tmp[code as usize].push(q_pos as u32);
            });
        }

        // Phase 2: pack into backbone + overflow (NCBI style)
        let mut presence = vec![0u64; capacity.div_ceil(64)];
        let mut backbone = vec![empty_cell; capacity];
        let mut overflow = Vec::new();

        for code in 0..capacity {
            let hits = &tmp[code];
            if hits.is_empty() { continue; }

            presence[code / 64] |= 1u64 << (code % 64);
            let n = hits.len();
            backbone[code].num_used = n as u32;

            if n <= HITS_PER_CELL {
                // Store inline — no pointer dereference needed at lookup time
                for (i, &pos) in hits.iter().enumerate() {
                    backbone[code].data[i] = pos;
                }
            } else {
                // Overflow: store cursor into overflow array
                backbone[code].data[0] = overflow.len() as u32;
                overflow.extend_from_slice(hits);
            }
        }

        ProteinLookup { word_size, charsize, mask, presence, backbone, overflow, capacity }
    }

    /// Get hits for a word code. Inlined for the hot scanning loop.
    /// For ≤3 hits (common case): returns pointer into the backbone cell (single cache line).
    /// For >3 hits: returns slice into overflow array.
    #[inline]
    pub fn get_hits(&self, code: u32) -> &[u32] {
        let code = code as usize;
        // Quick presence check (L1 cache — 4 KB bitfield)
        let word = code >> 6;
        let bit = code & 63;
        if unsafe { *self.presence.get_unchecked(word) } & (1u64 << bit) == 0 {
            return &[];
        }
        let cell = unsafe { self.backbone.get_unchecked(code) };
        let n = cell.num_used as usize;
        if n <= HITS_PER_CELL {
            // Inline hits — same cache line as num_used, no pointer chase
            unsafe { cell.data.get_unchecked(..n) }
        } else {
            // Overflow
            let start = cell.data[0] as usize;
            unsafe { self.overflow.get_unchecked(start..start + n) }
        }
    }

    /// Look up query positions for a subject word (Ncbistdaa encoded).
    #[inline]
    pub fn lookup(&self, subject_word: &[u8]) -> Option<&[u32]> {
        let hits = self.get_hits(encode_protein_word(subject_word));
        if hits.is_empty() { None } else { Some(hits) }
    }
}

/// Enumerate neighbors using shift-based (charsize=5) encoding.
fn enumerate_neighbors_shift(
    query_word: &[u8],
    word_size: usize,
    charsize: u32,
    matrix: &ScoringMatrix,
    threshold: i32,
    callback: &mut impl FnMut(u32),
) {
    let mut max_suffix = vec![0i32; word_size + 1];
    for i in (0..word_size).rev() {
        let q = query_word[i];
        let best = (0u8..28).map(|r| matrix.score(q, r)).max().unwrap_or(0);
        max_suffix[i] = max_suffix[i + 1] + best;
    }
    let mut score_rows = vec![[0i32; 28]; word_size];
    for (i, row) in score_rows.iter_mut().enumerate() {
        let q = query_word[i];
        for r in 0u8..28 {
            row[r as usize] = matrix.score(q, r);
        }
    }
    enumerate_rec_shift(&score_rows, &max_suffix, word_size, charsize, threshold, 0, 0, 0, callback);
}

/// Recursive neighbor enumeration with shift-based encoding.
#[allow(clippy::too_many_arguments)]
fn enumerate_rec_shift(
    score_rows: &[[i32; 28]],
    max_suffix: &[i32],
    word_size: usize,
    charsize: u32,
    threshold: i32,
    pos: usize,
    current_score: i32,
    current_code: u32,
    callback: &mut impl FnMut(u32),
) {
    if pos == word_size {
        if current_score >= threshold {
            callback(current_code);
        }
        return;
    }
    if current_score + max_suffix[pos] < threshold { return; }
    let scores = &score_rows[pos];
    for r in 0u8..28 {
        let new_score = current_score + scores[r as usize];
        if new_score + max_suffix[pos + 1] < threshold { continue; }
        // Shift-based encoding: (code << charsize) | residue
        let new_code = (current_code << charsize) | r as u32;
        enumerate_rec_shift(score_rows, max_suffix, word_size, charsize, threshold, pos + 1, new_score, new_code, callback);
    }
}

/// Nucleotide lookup table.
/// Word size typically 11. Uses a flat array indexed by 2-bit packed word.
pub struct NucleotideLookup {
    pub word_size: usize,
    /// positions[word_code] = list of query positions
    pub table: Vec<Vec<u32>>,
    pub capacity: usize,
}

impl NucleotideLookup {
    pub fn build(query: &[u8], word_size: usize) -> Self {
        assert!(word_size <= 16, "word_size must be <= 16 for 32-bit index");
        let capacity = 1usize << (2 * word_size);
        let mut table = vec![Vec::new(); capacity];

        if query.len() < word_size {
            return NucleotideLookup { word_size, table, capacity };
        }

        // Convert query to 2-bit encoded
        let encoded: Vec<u8> = query.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();

        // Build initial word from first (word_size-1) bases
        let mut word_code: u32 = 0;
        let mask = (capacity - 1) as u32;

        for &b in encoded.iter().take(word_size - 1) {
            if b > 3 {
                // Ambiguous base — reset
                word_code = 0;
                continue;
            }
            word_code = ((word_code << 2) | b as u32) & mask;
        }

        let mut valid = true;
        // Check if first word_size-1 bases are valid
        for &b in encoded.iter().take(word_size - 1) {
            if b > 3 { valid = false; break; }
        }

        for i in word_size - 1..query.len() {
            let b = encoded[i];
            if b > 3 {
                valid = false;
                word_code = 0;
                continue;
            }
            if !valid {
                // Try to rebuild
                let start = i + 1 - word_size;
                valid = true;
                word_code = 0;
                for &eb in encoded.iter().take(i + 1).skip(start) {
                    if eb > 3 { valid = false; break; }
                    word_code = ((word_code << 2) | eb as u32) & mask;
                }
                if !valid { continue; }
            } else {
                word_code = ((word_code << 2) | b as u32) & mask;
            }
            let q_pos = (i + 1 - word_size) as u32;
            table[word_code as usize].push(q_pos);
        }

        NucleotideLookup { word_size, table, capacity }
    }

    /// Returns query positions for a given 2-bit encoded word.
    pub fn lookup(&self, word_code: u32) -> &[u32] {
        &self.table[word_code as usize]
    }

    /// Scan subject sequence for hits against this lookup table.
    /// Returns iterator of (query_pos, subject_pos) hits.
    pub fn scan_subject<'a>(&'a self, subject: &'a [u8]) -> impl Iterator<Item = (u32, u32)> + 'a {
        let mask = (self.capacity - 1) as u32;
        let ws = self.word_size;
        let encoded: Vec<u8> = subject.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();

        let mut hits = Vec::new();
        if encoded.len() < ws { return hits.into_iter(); }

        let mut word_code: u32 = 0;
        let mut valid = true;

        for &b in encoded.iter().take(ws - 1) {
            if b > 3 { valid = false; word_code = 0; continue; }
            word_code = ((word_code << 2) | b as u32) & mask;
        }
        // Recheck first window validity
        for &b in encoded.iter().take(ws - 1) {
            if b > 3 { valid = false; break; }
        }

        for i in ws - 1..encoded.len() {
            let b = encoded[i];
            if b > 3 {
                valid = false;
                word_code = 0;
                continue;
            }
            if !valid {
                let start = i + 1 - ws;
                valid = true;
                word_code = 0;
                for &eb in encoded.iter().take(i + 1).skip(start) {
                    if eb > 3 { valid = false; break; }
                    word_code = ((word_code << 2) | eb as u32) & mask;
                }
                if !valid { continue; }
            } else {
                word_code = ((word_code << 2) | b as u32) & mask;
            }
            let s_pos = (i + 1 - ws) as u32;
            for &q_pos in &self.table[word_code as usize] {
                hits.push((q_pos, s_pos));
            }
        }
        hits.into_iter()
    }
}

/// NCBI discontiguous megablast templates.
/// Template type 0: coding (optimized for coding sequences)
/// Template type 1: optimal (maximizes sensitivity)
/// Template type 2: two simultaneous templates
pub fn get_discontiguous_template(template_type: u8, template_length: usize) -> Vec<bool> {
    match (template_type, template_length) {
        // Coding template, length 21, 11 care positions
        (0, 21) => vec![
            true, true, true, false, true, true, false, true, false, true, false,
            false, true, true, false, true, false, true, true, true, true,
        ],
        // Optimal template, length 21, 11 care positions
        (1, 21) => vec![
            true, true, false, true, false, true, true, false, false, true, false,
            true, false, false, true, true, false, true, true, true, true,
        ],
        // Coding template, length 18, 11 care positions
        (0, 18) => vec![
            true, true, true, false, true, true, false, true, false, true,
            false, true, true, false, true, true, true, true,
        ],
        // Optimal template, length 18, 11 care positions
        (1, 18) => vec![
            true, true, false, true, false, true, true, false, true, false,
            true, false, true, true, true, true, true, true,
        ],
        // Default: contiguous (all care)
        _ => vec![true; template_length],
    }
}

/// Discontiguous megablast lookup table.
/// Uses spaced seed templates for more sensitive nucleotide searching.
pub struct DiscontiguousLookup {
    pub template_length: usize,
    pub word_size: usize, // number of care positions
    /// Template mask: true = care position, false = don't care
    pub template: Vec<bool>,
    /// table[code] = list of query positions
    pub table: Vec<Vec<u32>>,
    pub capacity: usize,
}

impl DiscontiguousLookup {
    pub fn build(query: &[u8], template_type: u8, template_length: usize) -> Self {
        let template = get_discontiguous_template(template_type, template_length);
        let word_size = template.iter().filter(|&&b| b).count();
        assert!(word_size <= 16, "word_size must be <= 16 for 32-bit index");
        let capacity = 1usize << (2 * word_size);
        let mut table = vec![Vec::new(); capacity];

        let encoded: Vec<u8> = query.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();

        if encoded.len() < template_length {
            return DiscontiguousLookup { template_length, word_size, template, table, capacity };
        }

        for pos in 0..=(encoded.len() - template_length) {
            let mut code = 0u32;
            let mut valid = true;
            for (i, &care) in template.iter().enumerate() {
                if care {
                    let b = encoded[pos + i];
                    if b > 3 { valid = false; break; }
                    code = (code << 2) | b as u32;
                }
            }
            if valid {
                table[code as usize].push(pos as u32);
            }
        }

        DiscontiguousLookup { template_length, word_size, template, table, capacity }
    }

    pub fn scan_subject<'a>(&'a self, subject: &'a [u8]) -> Vec<(u32, u32)> {
        let encoded: Vec<u8> = subject.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();
        let mut hits = Vec::new();

        if encoded.len() < self.template_length { return hits; }

        for pos in 0..=(encoded.len() - self.template_length) {
            let mut code = 0u32;
            let mut valid = true;
            for (i, &care) in self.template.iter().enumerate() {
                if care {
                    let b = encoded[pos + i];
                    if b > 3 { valid = false; break; }
                    code = (code << 2) | b as u32;
                }
            }
            if valid {
                for &q_pos in &self.table[code as usize] {
                    hits.push((q_pos, pos as u32));
                }
            }
        }
        hits
    }
}
