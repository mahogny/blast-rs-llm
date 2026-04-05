//! BlastDB v4/v5 writer – creates protein (.pin/.psq/.phr) and nucleotide (.nin/.nsq/.nhr) databases.
//! V5 additionally writes LMDB accession index (.pdb/.ndb), OID→SeqIDs (.pos/.nos),
//! and OID→TaxIDs (.pot/.not) files.

use std::fs;
use std::path::{Path, PathBuf};
use byteorder::{BigEndian, LittleEndian, WriteBytesExt};
use lmdb::{self, DatabaseFlags, Environment, EnvironmentFlags, Transaction, WriteFlags};
use crate::db::error::{DbError, Result};
use crate::db::index::SeqType;

/// One sequence to be added to a BLAST database.
///
/// The `sequence` field should contain raw ASCII bytes — amino acids (A-Z) for protein,
/// nucleotides (A/C/G/T/N) for nucleotide. The builder handles encoding internally.
pub struct SequenceEntry {
    /// Full description (e.g., "Human insulin precursor").
    pub title: String,
    /// Accession identifier (e.g., "P01308").
    pub accession: String,
    /// Raw ASCII sequence bytes (uppercase preferred).
    pub sequence: Vec<u8>,
    /// NCBI taxonomy ID (e.g., 9606 for human). Optional.
    pub taxid: Option<u32>,
}

/// Builder for creating BLAST databases from sequences.
///
/// Supports v4 and v5 formats, single and multi-volume output.
///
/// ```no_run
/// use blast_rs::{BlastDbBuilder, SequenceEntry};
/// use blast_rs::db::index::SeqType;
/// use std::path::Path;
///
/// let mut builder = BlastDbBuilder::new(SeqType::Protein, "My DB");
/// builder.add(SequenceEntry {
///     title: "test".into(), accession: "P001".into(),
///     sequence: b"MKFLILLF".to_vec(), taxid: None,
/// });
/// builder.write(Path::new("mydb")).unwrap();
/// ```
pub struct BlastDbBuilder {
    pub seq_type: SeqType,
    pub db_title: String,
    pub entries: Vec<SequenceEntry>,
}

impl BlastDbBuilder {
    pub fn new(seq_type: SeqType, db_title: impl Into<String>) -> Self {
        BlastDbBuilder { seq_type, db_title: db_title.into(), entries: Vec::new() }
    }

    pub fn add(&mut self, entry: SequenceEntry) {
        self.entries.push(entry);
    }

    /// Write all database files to `base_path` (without extension) using v4 format.
    pub fn write(&self, base_path: &Path) -> Result<()> {
        match self.seq_type {
            SeqType::Protein    => self.write_protein(base_path, 4),
            SeqType::Nucleotide => self.write_nucleotide(base_path, 4),
        }
    }

    /// Write a multi-volume database, splitting when sequence data exceeds `max_file_size` bytes.
    /// Creates `.00`, `.01`, ... volume files plus a `.pal`/`.nal` alias file.
    /// If all data fits in one volume, writes a single-volume database (no alias).
    pub fn write_multivolume(&self, base_path: &Path, format_version: i32, max_file_size: u64) -> Result<()> {
        if self.entries.is_empty() {
            return if format_version == 5 { self.write_v5(base_path) } else { self.write(base_path) };
        }

        // Split entries into volumes based on cumulative sequence size
        let mut volumes: Vec<(usize, usize)> = Vec::new(); // (start_idx, end_idx) inclusive
        let mut vol_start = 0;
        let mut vol_bytes = 0u64;

        for (i, entry) in self.entries.iter().enumerate() {
            let entry_size = entry.sequence.len() as u64 + 1; // +1 for NUL terminator
            if vol_bytes + entry_size > max_file_size && i > vol_start {
                volumes.push((vol_start, i));
                vol_start = i;
                vol_bytes = 0;
            }
            vol_bytes += entry_size;
        }
        volumes.push((vol_start, self.entries.len()));

        // Single volume — no alias needed
        if volumes.len() == 1 {
            return if format_version == 5 { self.write_v5(base_path) } else { self.write(base_path) };
        }

        // Write each volume with .00, .01, ... suffix.
        // NCBI uses mydb.00.pin, mydb.00.psq etc. Since Rust's with_extension replaces
        // everything after the last dot, we write to a temp name then rename the files.
        let mut vol_names = Vec::new();
        let base_str = base_path.to_string_lossy().to_string();
        let base_dir = base_path.parent().unwrap_or(Path::new("."));
        for (vi, &(start, end)) in volumes.iter().enumerate() {
            let vol_suffix = format!("{:02}", vi);

            // Create sub-builder for this volume
            let mut vol_builder = BlastDbBuilder::new(self.seq_type, &self.db_title);
            for entry in &self.entries[start..end] {
                vol_builder.add(SequenceEntry {
                    title: entry.title.clone(),
                    accession: entry.accession.clone(),
                    sequence: entry.sequence.clone(),
                    taxid: entry.taxid,
                });
            }

            // Write to a temp path without dots in the stem, then rename
            let temp_path = base_dir.join(format!("__vol_tmp_{}", vi));
            if format_version == 5 {
                vol_builder.write_v5(&temp_path)?;
            } else {
                vol_builder.write(&temp_path)?;
            }

            // Rename files: __vol_tmp_0.pin → mydb.00.pin etc.
            let exts: &[&str] = match self.seq_type {
                SeqType::Protein    => &["pin", "psq", "phr", "pdb", "pos", "pot"],
                SeqType::Nucleotide => &["nin", "nsq", "nhr", "ndb", "nos", "not"],
            };
            for ext in exts {
                let src = temp_path.with_extension(ext);
                if src.exists() {
                    let dst = PathBuf::from(format!("{}.{}.{}", base_str, vol_suffix, ext));
                    fs::rename(&src, &dst)?;
                }
            }

            vol_names.push(format!("{}.{}",
                base_path.file_name().unwrap_or_default().to_string_lossy(), vol_suffix));
        }

        // Write alias file (.pal or .nal)
        let alias_ext = match self.seq_type {
            SeqType::Protein    => "pal",
            SeqType::Nucleotide => "nal",
        };
        let alias_content = format!(
            "#\n# Alias file created by blast-rs\n#\nTITLE {}\nDBLIST {}\n",
            self.db_title,
            vol_names.join(" "),
        );
        fs::write(base_path.with_extension(alias_ext), alias_content)?;

        Ok(())
    }

    /// Write a v5 database: sequence/header/index files (format_version=5) plus
    /// LMDB accession index (.pdb/.ndb), OID→SeqIDs (.pos/.nos), and OID→TaxIDs (.pot/.not).
    pub fn write_v5(&self, base_path: &Path) -> Result<()> {
        match self.seq_type {
            SeqType::Protein    => self.write_protein(base_path, 5)?,
            SeqType::Nucleotide => self.write_nucleotide(base_path, 5)?,
        }
        self.write_lmdb(base_path)?;
        self.write_oid_seqids(base_path)?;
        self.write_oid_taxids(base_path)?;
        Ok(())
    }

    /// Write the LMDB accession index (.pdb/.ndb).
    fn write_lmdb(&self, base: &Path) -> Result<()> {
        let ext = match self.seq_type {
            SeqType::Protein    => "pdb",
            SeqType::Nucleotide => "ndb",
        };
        let lmdb_path = base.with_extension(ext);

        let env = Environment::new()
            .set_flags(EnvironmentFlags::NO_SUB_DIR)
            .set_max_dbs(4)
            .set_map_size(1 << 30) // 1 GiB virtual; grows as needed
            .open(&lmdb_path)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB create '{}': {}", lmdb_path.display(), e)))?;

        // Create named databases
        let db_acc2oid = env.create_db(Some("acc2oid"), DatabaseFlags::DUP_SORT | DatabaseFlags::DUP_FIXED)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB create 'acc2oid': {}", e)))?;
        let db_volname = env.create_db(Some("volname"), DatabaseFlags::INTEGER_KEY)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB create 'volname': {}", e)))?;
        let db_volinfo = env.create_db(Some("volinfo"), DatabaseFlags::INTEGER_KEY)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB create 'volinfo': {}", e)))?;

        let mut txn = env.begin_rw_txn()
            .map_err(|e| DbError::InvalidFormat(format!("LMDB write txn: {}", e)))?;

        // acc2oid: accession bytes → LE u32 OID
        for (oid, entry) in self.entries.iter().enumerate() {
            let oid_bytes = (oid as u32).to_le_bytes();
            txn.put(db_acc2oid, &entry.accession.as_bytes(), &oid_bytes, WriteFlags::empty())
                .map_err(|e| DbError::InvalidFormat(format!("LMDB put acc2oid: {}", e)))?;
        }

        // volname: key=native u32 0 → volume base name
        let vol_key: u32 = 0;
        let vol_name = base.file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        txn.put(db_volname, &vol_key.to_ne_bytes(), &vol_name.as_bytes(), WriteFlags::empty())
            .map_err(|e| DbError::InvalidFormat(format!("LMDB put volname: {}", e)))?;

        // volinfo: key=native u32 0 → LE u32 num_oids
        let num_oids = (self.entries.len() as u32).to_le_bytes();
        txn.put(db_volinfo, &vol_key.to_ne_bytes(), &num_oids, WriteFlags::empty())
            .map_err(|e| DbError::InvalidFormat(format!("LMDB put volinfo: {}", e)))?;

        txn.commit()
            .map_err(|e| DbError::InvalidFormat(format!("LMDB commit: {}", e)))?;

        Ok(())
    }

    /// Write the OID→SeqIDs file (.pos/.nos).
    fn write_oid_seqids(&self, base: &Path) -> Result<()> {
        let ext = match self.seq_type {
            SeqType::Protein    => "pos",
            SeqType::Nucleotide => "nos",
        };
        let num_oids = self.entries.len() as u64;

        // Build data section: for each OID, encode the accession as a length-prefixed string
        let mut data_section: Vec<u8> = Vec::new();
        let mut end_offsets: Vec<u64> = Vec::new();

        for entry in &self.entries {
            let acc = entry.accession.as_bytes();
            if acc.len() < 0xFF {
                data_section.push(acc.len() as u8);
            } else {
                data_section.push(0xFF);
                data_section.extend_from_slice(&(acc.len() as u32).to_le_bytes());
            }
            data_section.extend_from_slice(acc);
            end_offsets.push(data_section.len() as u64);
        }

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(&num_oids.to_le_bytes());
        for &off in &end_offsets {
            buf.extend_from_slice(&off.to_le_bytes());
        }
        buf.extend_from_slice(&data_section);

        fs::write(base.with_extension(ext), &buf)?;
        Ok(())
    }

    /// Write the OID→TaxIDs file (.pot/.not).
    fn write_oid_taxids(&self, base: &Path) -> Result<()> {
        let ext = match self.seq_type {
            SeqType::Protein    => "pot",
            SeqType::Nucleotide => "not",
        };
        let num_oids = self.entries.len() as u64;

        // Build data section: for each OID, write its taxid(s) as i32 LE
        let mut data_section: Vec<u8> = Vec::new();
        let mut end_offsets: Vec<u64> = Vec::new(); // in units of i32 elements

        let mut elem_count: u64 = 0;
        for entry in &self.entries {
            if let Some(taxid) = entry.taxid {
                data_section.extend_from_slice(&(taxid as i32).to_le_bytes());
                elem_count += 1;
            }
            end_offsets.push(elem_count);
        }

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(&num_oids.to_le_bytes());
        for &off in &end_offsets {
            buf.extend_from_slice(&off.to_le_bytes());
        }
        buf.extend_from_slice(&data_section);

        fs::write(base.with_extension(ext), &buf)?;
        Ok(())
    }

    // ── Protein ─────────────────────────────────────────────────────────────

    fn write_protein(&self, base: &Path, format_version: i32) -> Result<()> {
        let mut seq_data: Vec<u8> = Vec::new();
        let mut hdr_data: Vec<u8> = Vec::new();
        let mut sequence_array: Vec<u32> = Vec::new();
        let mut header_array: Vec<u32> = Vec::new();
        let mut max_seq_length = 0u32;
        let mut volume_length  = 0u64;

        for entry in &self.entries {
            header_array.push(hdr_data.len() as u32);
            hdr_data.extend(encode_defline_ber(&entry.title, &entry.accession, entry.taxid));

            sequence_array.push(seq_data.len() as u32);
            let enc = encode_protein_seq(&entry.sequence);
            let slen = enc.len() as u32;
            seq_data.extend_from_slice(&enc);
            seq_data.push(0x00); // NUL terminator
            volume_length  += slen as u64;
            if slen > max_seq_length { max_seq_length = slen; }
        }
        // Sentinel entries (num_oids+1 each)
        header_array.push(hdr_data.len() as u32);
        sequence_array.push(seq_data.len() as u32);

        fs::write(base.with_extension("psq"), &seq_data)?;
        fs::write(base.with_extension("phr"), &hdr_data)?;
        fs::write(base.with_extension("pin"), build_index(
            SeqType::Protein, &self.db_title,
            format_version,
            self.entries.len() as u32, volume_length, max_seq_length,
            &header_array, &sequence_array, None,
        ))?;
        Ok(())
    }

    // ── Nucleotide ───────────────────────────────────────────────────────────

    fn write_nucleotide(&self, base: &Path, format_version: i32) -> Result<()> {
        let mut seq_data: Vec<u8> = Vec::new();
        let mut hdr_data: Vec<u8> = Vec::new();
        let mut sequence_array: Vec<u32> = Vec::new();
        let mut ambig_array:    Vec<u32> = Vec::new();
        let mut header_array:   Vec<u32> = Vec::new();
        let mut max_seq_length = 0u32;
        let mut volume_length  = 0u64;

        for entry in &self.entries {
            header_array.push(hdr_data.len() as u32);
            hdr_data.extend(encode_defline_ber(&entry.title, &entry.accession, entry.taxid));

            sequence_array.push(seq_data.len() as u32);
            let (packed, ambig) = encode_nucleotide_seq(&entry.sequence);
            seq_data.extend_from_slice(&packed);
            ambig_array.push(seq_data.len() as u32); // end of packed = start of ambig
            seq_data.extend_from_slice(&ambig);

            let slen = entry.sequence.len() as u32;
            volume_length  += slen as u64;
            if slen > max_seq_length { max_seq_length = slen; }
        }
        // Sentinel entries (num_oids+1 each)
        header_array.push(hdr_data.len() as u32);
        sequence_array.push(seq_data.len() as u32);
        ambig_array.push(seq_data.len() as u32); // unused sentinel

        fs::write(base.with_extension("nsq"), &seq_data)?;
        fs::write(base.with_extension("nhr"), &hdr_data)?;
        fs::write(base.with_extension("nin"), build_index(
            SeqType::Nucleotide, &self.db_title,
            format_version,
            self.entries.len() as u32, volume_length, max_seq_length,
            &header_array, &sequence_array, Some(&ambig_array),
        ))?;
        Ok(())
    }
}

// ─── Index file serialisation ────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn build_index(
    seq_type: SeqType,
    title: &str,
    format_version: i32,
    num_oids: u32,
    volume_length: u64,
    max_seq_length: u32,
    header_array: &[u32],
    sequence_array: &[u32],
    ambig_array: Option<&[u32]>,
) -> Vec<u8> {
    let mut buf: Vec<u8> = Vec::new();

    // format_version (BigEndian)
    buf.write_i32::<BigEndian>(format_version).unwrap();
    // seq_type: 1=protein, 0=nucleotide
    buf.write_i32::<BigEndian>(match seq_type {
        SeqType::Protein    => 1,
        SeqType::Nucleotide => 0,
    }).unwrap();

    write_be_string(&mut buf, title);
    write_be_string(&mut buf, ""); // create_date (empty)

    buf.write_u32::<BigEndian>(num_oids).unwrap();
    buf.write_u64::<LittleEndian>(volume_length).unwrap(); // LITTLE-endian!
    buf.write_u32::<BigEndian>(max_seq_length).unwrap();

    for &v in header_array   { buf.write_u32::<BigEndian>(v).unwrap(); }
    for &v in sequence_array { buf.write_u32::<BigEndian>(v).unwrap(); }
    if let Some(arr) = ambig_array {
        for &v in arr { buf.write_u32::<BigEndian>(v).unwrap(); }
    }

    buf
}

/// Write a length-prefixed string as used in the index (i32 BE length + bytes).
fn write_be_string(buf: &mut Vec<u8>, s: &str) {
    buf.write_i32::<BigEndian>(s.len() as i32).unwrap();
    buf.extend_from_slice(s.as_bytes());
}

// ─── Protein encoding ────────────────────────────────────────────────────────

fn encode_protein_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&c| ascii_to_ncbistdaa(c)).collect()
}

fn ascii_to_ncbistdaa(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 1,  b'B' => 2,  b'C' => 3,  b'D' => 4,  b'E' => 5,
        b'F' => 6,  b'G' => 7,  b'H' => 8,  b'I' => 9,  b'K' => 10,
        b'L' => 11, b'M' => 12, b'N' => 13, b'P' => 14, b'Q' => 15,
        b'R' => 16, b'S' => 17, b'T' => 18, b'V' => 19, b'W' => 20,
        b'X' => 21, b'Y' => 22, b'Z' => 23, b'U' => 24, b'*' => 25,
        b'O' => 26, b'J' => 27,
        _ => 21, // unknown → X
    }
}

// ─── Nucleotide encoding ─────────────────────────────────────────────────────

/// Returns `(packed_na2_bytes, ambig_data_bytes)`.
///
/// Packed NcbiNa2: 4 bases/byte MSB-first (A=0, C=1, G=2, T=3).
/// Last byte stores the remainder count in its low 2 bits (0 = sentinel, bases in prev bytes only).
/// IUPAC ambiguous bases are encoded as A in the packed stream and listed in the ambig_data section.
fn encode_nucleotide_seq(seq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    // Collect ambiguity runs: (start_pos, length_minus_1, na4_code)
    let mut segs: Vec<(u32, u32, u8)> = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        let na4 = iupac_to_na4(seq[i].to_ascii_uppercase());
        if !matches!(na4, 1 | 2 | 4 | 8) {
            // ambiguous base — collect run of identical na4
            let start = i;
            let code  = na4;
            let mut j = i + 1;
            while j < seq.len() && iupac_to_na4(seq[j].to_ascii_uppercase()) == code {
                j += 1;
            }
            segs.push((start as u32, (j - start - 1) as u32, code));
            i = j;
        } else {
            i += 1;
        }
    }

    let packed = pack_na2(seq);

    let ambig = if segs.is_empty() {
        Vec::new()
    } else {
        let mut data: Vec<u8> = Vec::new();
        // num_segments with high bit set → new (8-byte) format
        let n = (segs.len() as u32) | 0x8000_0000;
        data.write_u32::<BigEndian>(n).unwrap();
        for (start, len_m1, na4) in &segs {
            // bits 63-60: na4 code, bits 59-48: length-1, bits 47-32: 0, bits 31-0: start
            let val: u64 = ((*na4 as u64) << 60)
                | ((*len_m1 as u64 & 0xfff) << 48)
                | (*start as u64 & 0xffff_ffff);
            data.write_u64::<BigEndian>(val).unwrap();
        }
        data
    };

    (packed, ambig)
}

fn pack_na2(seq: &[u8]) -> Vec<u8> {
    let n = seq.len();
    if n == 0 {
        return vec![0x00]; // sentinel
    }
    let num_full = n / 4;
    let remainder = n % 4;
    let mut out = Vec::with_capacity(num_full + 1);

    for i in 0..num_full {
        let mut byte = 0u8;
        for j in 0..4usize {
            byte |= base_to_2bit(seq[i * 4 + j].to_ascii_uppercase()) << ((3 - j) * 2);
        }
        out.push(byte);
    }

    if remainder == 0 {
        out.push(0x00); // sentinel: 0 extra bases
    } else {
        // Pack `remainder` bases in high bits; store count in low 2 bits.
        let mut byte = remainder as u8;
        for j in 0..remainder {
            byte |= base_to_2bit(seq[num_full * 4 + j].to_ascii_uppercase()) << ((3 - j) * 2);
        }
        out.push(byte);
    }

    out
}

fn base_to_2bit(c: u8) -> u8 {
    match c { b'C' => 1, b'G' => 2, b'T' | b'U' => 3, _ => 0 }
}

/// Map an uppercase IUPAC character to a 4-bit NA4 code (bit per base: A=1,C=2,G=4,T=8).
fn iupac_to_na4(c: u8) -> u8 {
    match c {
        b'A' => 0x1, b'C' => 0x2, b'G' => 0x4, b'T' | b'U' => 0x8,
        b'M' => 0x3, b'R' => 0x5, b'S' => 0x6, b'V' => 0x7,
        b'W' => 0x9, b'Y' => 0xa, b'H' => 0xb, b'K' => 0xc,
        b'D' => 0xd, b'B' => 0xe, b'N' => 0xf,
        _    => 0xf, // treat unknown as N
    }
}

// ─── BER encoding ────────────────────────────────────────────────────────────

fn ber_length(len: usize) -> Vec<u8> {
    if len < 128 {
        vec![len as u8]
    } else {
        let mut bytes: Vec<u8> = Vec::new();
        let mut n = len;
        while n > 0 { bytes.push((n & 0xff) as u8); n >>= 8; }
        bytes.reverse();
        let mut out = vec![0x80 | bytes.len() as u8];
        out.extend_from_slice(&bytes);
        out
    }
}

fn ber_tlv(tag: u8, content: &[u8]) -> Vec<u8> {
    let mut out = vec![tag];
    out.extend(ber_length(content.len()));
    out.extend_from_slice(content);
    out
}

fn ber_visible_string(s: &str) -> Vec<u8> {
    ber_tlv(0x1a, s.as_bytes())
}

/// Minimal two's-complement BER INTEGER encoding.
fn ber_integer(n: i64) -> Vec<u8> {
    if n == 0 {
        return ber_tlv(0x02, &[0x00]);
    }
    let raw = n.to_be_bytes();
    // Strip redundant leading bytes
    let mut start = 0;
    while start < 7 {
        let b    = raw[start];
        let next = raw[start + 1];
        if (n >= 0 && b == 0x00 && next & 0x80 == 0)
        || (n <  0 && b == 0xff && next & 0x80 != 0)
        {
            start += 1;
        } else {
            break;
        }
    }
    ber_tlv(0x02, &raw[start..])
}

/// Encode one sequence as a complete BER Blast-def-line-set (SEQUENCE OF one Blast-def-line).
///
/// Structure written:
/// ```text
/// SEQUENCE (0x30)                    ← Blast-def-line-set
///   SEQUENCE (0x30)                  ← Blast-def-line
///     [0] EXPLICIT (0xa0)            ← title
///       VisibleString (0x1a) title
///     [1] EXPLICIT (0xa1)            ← seqid (Seq-id TLVs directly, no SEQUENCE OF wrapper)
///       [0] context (0xa0)           ← local Object-id
///         VisibleString (0x1a) accession
///     [2] EXPLICIT (0xa2)  [opt]     ← taxid
///       INTEGER (0x02) taxid
/// ```
pub fn encode_defline_ber(title: &str, accession: &str, taxid: Option<u32>) -> Vec<u8> {
    let title_field = ber_tlv(0xa0, &ber_visible_string(title));

    // The parser iterates Seq-id TLVs directly inside [1] (no SEQUENCE OF wrapper).
    let seq_id      = ber_tlv(0xa0, &ber_visible_string(accession)); // local Object-id (str)
    let seqid_field = ber_tlv(0xa1, &seq_id);

    let taxid_field = taxid
        .map(|t| ber_tlv(0xa2, &ber_integer(t as i64)))
        .unwrap_or_default();

    let mut defline_content = Vec::new();
    defline_content.extend_from_slice(&title_field);
    defline_content.extend_from_slice(&seqid_field);
    defline_content.extend_from_slice(&taxid_field);

    let defline = ber_tlv(0x30, &defline_content);
    ber_tlv(0x30, &defline) // Blast-def-line-set
}
