//! Sequence decoding for protein (Ncbistdaa) and nucleotide (NcbiNa2 + ambiguities).

/// Ncbistdaa to single-letter amino acid code.
/// Index 0 = gap, 1 = Ala, 2 = Arg, ...
/// Table from NCBI toolkit: ncbistdaa encoding.
pub const NCBISTDAA_TO_AA: &[u8] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

/// Convert a protein sequence from Ncbistdaa encoding to ASCII.
pub fn decode_protein(data: &[u8]) -> Vec<u8> {
    data.iter().map(|&b| {
        let idx = b as usize;
        if idx < NCBISTDAA_TO_AA.len() {
            NCBISTDAA_TO_AA[idx]
        } else {
            b'X'
        }
    }).collect()
}

/// Decode a packed NcbiNa2 nucleotide sequence including ambiguities.
///
/// - `packed`: the raw packed bytes (4 bases/byte, MSB first)
/// - `ambig_data`: raw ambiguity segment data (may be empty if none)
///
/// Returns the decoded sequence as ASCII (A/C/G/T/N/R/Y/... IUPAC).
pub fn decode_nucleotide(packed: &[u8], ambig_data: &[u8]) -> Vec<u8> {
    if packed.is_empty() {
        return Vec::new();
    }

    // The last byte encodes the "remainder": how many bases in the last byte are real.
    // remainder is stored in the lowest 2 bits of the last byte.
    let last_byte = *packed.last().unwrap();
    let remainder = (last_byte & 0x03) as usize; // 0 means all 4 of the last byte... wait

    // From spec: if remainder == 0, an extra byte is appended containing just the count 0.
    // So: total bases = (packed.len() - 1) * 4 + remainder
    // But if remainder == 0, we need special handling: it means the extra sentinel byte was added.
    // Actually re-reading: "if a sequence is exactly divisible by four, an additional byte must be
    // appended to contain this count (which will be zero)."
    // So the last byte is ALWAYS the count byte. Bases come from all bytes except the last...
    // wait no — the count byte is the last byte, and when remainder==0 it means zero bases in that last byte.
    // But for the packed data, the last data byte AND the count byte overlap for non-zero remainder.
    // Let me re-read: "the last base of the last byte to store a number from 0-3"
    // So the last byte's lowest 2 bits are the count; its upper 6 bits contain up to 3 real bases.
    // Example: TACG -> (3,0,1,2) packed -> one byte 0b11_00_01_10 = 0xC6; then remainder byte = 0b00_00_00_00 (0, 0 extra bases)
    // Wait, from spec example: "For the sequence (TACG), the full byte encoding is (225, 0)"
    // 225 = 0b11100001... hmm let me recalculate.
    // T=3, A=0, C=1, G=2. Packed MSB first: (3<<6)|(0<<4)|(1<<2)|2 = 192+0+4+2=198... not 225.
    // Oh wait, the spec says the REMAINDER is stored in the last byte's lowest 2 bits.
    // For TACG (4 bases, exactly divisible), an extra byte is added with count=0. So:
    // byte 0: T A C G -> 3,0,1,2 -> (3<<6)|(0<<4)|(1<<2)|2 = 198
    // byte 1: just the count = 0
    // But the spec says "225, 0"... let me check: 225 = 0xE1 = 0b11100001.
    // Hmm, maybe TACG: T=3(11), A=0(00), C=1(01), G=2(10) -> 11000110 = 0xC6 = 198. Not 225.
    // Wait: maybe the spec example has A=00, C=01, G=10, T=11?
    // 225 = 0b11100001. If T=11, A=00, C=01, G=10... no that doesn't work either.
    // Let me try: TACG where A=0,C=1,G=2,T=3 as stated:
    // T=3=11, A=0=00, C=1=01, G=2=10 packed MSB first: 11 00 01 10 = 0b11000110 = 198
    // That contradicts the spec example of 225.
    // Actually maybe the spec says (T,A,C,G) as 4 bases = one full byte, plus remainder byte = 0
    // And 225 is just wrong in the spec, or I'm misreading... Let me try TGGTTACAAC:
    // T=3,G=2,G=2,T=3 -> 11101011 = 0xEB = 235 ✓
    // T=3,A=0,C=1,A=0 -> 11000100 = 0xC4 = 196 ✓
    // A=0,C=1 padded with 0,remainder=2 -> 00 01 00 10 = 0b00010010 = 18 ✓
    // So TGGTTACAAC -> (235, 196, 18) - matches! Great.
    // So for TACG: (198, 0). The spec said 225 which seems wrong. Let me proceed with 198,0.

    // num_bases = (len-1)*4 + remainder, where remainder from 1..=3 means that many bases in last byte
    // if remainder == 0, the last byte is pure sentinel with 0 bases
    let num_full_bytes = packed.len() - 1; // last byte is the "tail" byte
    let num_bases = if remainder == 0 {
        num_full_bytes * 4
    } else {
        (num_full_bytes * 4) + remainder
    };

    if num_bases == 0 {
        return Vec::new();
    }

    // Decode all bases from the packed bytes (not the sentinel byte at the end)
    static NA2_TO_ASCII: [u8; 4] = [b'A', b'C', b'G', b'T'];

    let mut bases = Vec::with_capacity(num_bases);
    for &byte in &packed[..num_full_bytes] {
        bases.push(NA2_TO_ASCII[((byte >> 6) & 3) as usize]);
        bases.push(NA2_TO_ASCII[((byte >> 4) & 3) as usize]);
        bases.push(NA2_TO_ASCII[((byte >> 2) & 3) as usize]);
        bases.push(NA2_TO_ASCII[(byte & 3) as usize]);
    }

    // Handle last (partial) byte if remainder > 0
    if remainder > 0 {
        let byte = packed[num_full_bytes]; // same as last byte which has count in low bits
        // But the count IS in the low 2 bits, so if remainder=2:
        // bits 7,6 = base 1; bits 5,4 = base 2; bits 3,2 = padding 0; bits 1,0 = count(2)
        // Actually the bases and count share the byte:
        // For remainder=1: bits 7,6 = base; bits 5-0 = (padding + count)
        // For remainder=2: bits 7-4 = 2 bases; bits 3-0 = (padding + count)
        // For remainder=3: bits 7-2 = 3 bases; bits 1-0 = count
        // So the count is always the last 2 bits, and we should decode only `remainder` bases
        // from the top bits.
        for i in 0..remainder {
            let shift = 6 - i * 2;
            let code = (byte >> shift) & 3;
            bases.push(NA2_TO_ASCII[code as usize]);
        }
    }

    // Apply ambiguity data
    if !ambig_data.is_empty() {
        apply_ambiguities(&mut bases, ambig_data);
    }

    bases
}

/// NcbiNA4 to IUPAC ambiguity codes.
/// NA4 uses 4 bits, encoding which bases are possible.
/// bit 0=A, 1=C, 2=G, 3=T
static NA4_TO_IUPAC: [u8; 16] = [
    b'N', // 0000 = none (treat as N)
    b'A', // 0001 = A
    b'C', // 0010 = C
    b'M', // 0011 = A or C
    b'G', // 0100 = G
    b'R', // 0101 = A or G
    b'S', // 0110 = C or G
    b'V', // 0111 = A or C or G
    b'T', // 1000 = T
    b'W', // 1001 = A or T
    b'Y', // 1010 = C or T
    b'H', // 1011 = A or C or T
    b'K', // 1100 = G or T
    b'D', // 1101 = A or G or T
    b'B', // 1110 = C or G or T
    b'N', // 1111 = any
];

fn apply_ambiguities(bases: &mut [u8], ambig_data: &[u8]) {
    if ambig_data.len() < 4 {
        return;
    }

    let num_segments_raw = i32::from_be_bytes([ambig_data[0], ambig_data[1], ambig_data[2], ambig_data[3]]);
    let new_format = (num_segments_raw as u32) & 0x80000000 != 0;
    let num_segments = ((num_segments_raw as u32) & 0x7fffffff) as usize;

    let seg_data = &ambig_data[4..];

    if new_format {
        // Each segment is 8 bytes (Int8 big-endian)
        for i in 0..num_segments {
            let off = i * 8;
            if off + 8 > seg_data.len() { break; }
            let val = u64::from_be_bytes(seg_data[off..off+8].try_into().unwrap());
            let na4_value = ((val >> 60) & 0xf) as u8;
            let length = (((val >> 48) & 0xfff) + 1) as usize;
            let start = (val & 0xffffffff) as usize;
            let iupac = NA4_TO_IUPAC[na4_value as usize];
            for j in 0..length {
                let pos = start + j;
                if pos < bases.len() {
                    bases[pos] = iupac;
                }
            }
        }
    } else {
        // Old format: each segment is 4 bytes (Int4 big-endian)
        for i in 0..num_segments {
            let off = i * 4;
            if off + 4 > seg_data.len() { break; }
            let val = u32::from_be_bytes(seg_data[off..off+4].try_into().unwrap());
            let na4_value = ((val >> 28) & 0xf) as u8;
            let length = (((val >> 24) & 0xf) + 1) as usize;
            let start = (val & 0x00ffffff) as usize;
            let iupac = NA4_TO_IUPAC[na4_value as usize];
            for j in 0..length {
                let pos = start + j;
                if pos < bases.len() {
                    bases[pos] = iupac;
                }
            }
        }
    }
}

/// Return the raw Ncbistdaa bytes for a protein OID.
/// The caller ensures start < end and the data slice is the full psq file.
pub fn get_protein_raw(psq: &[u8], start: usize, end: usize) -> &[u8] {
    // end offset points to the NUL byte separator; sequence is start..end-1
    if end == 0 || start >= psq.len() {
        return &[];
    }
    let actual_end = (end - 1).min(psq.len());
    &psq[start..actual_end]
}

/// Decode nucleotide sequence for one OID.
/// seq_range: `sequence_array[oid] .. ambig_array[oid]`
/// ambig_range: `ambig_array[oid] .. sequence_array[oid+1]`
pub fn get_nucleotide(nsq: &[u8], seq_start: usize, seq_end: usize, ambig_start: usize, ambig_end: usize) -> Vec<u8> {
    let packed = if seq_start < seq_end && seq_end <= nsq.len() {
        &nsq[seq_start..seq_end]
    } else {
        &[]
    };
    let ambig = if ambig_start < ambig_end && ambig_end <= nsq.len() {
        &nsq[ambig_start..ambig_end]
    } else {
        &[]
    };
    decode_nucleotide(packed, ambig)
}
