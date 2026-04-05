use std::path::{Path, PathBuf};
use std::fs;
use memmap2::Mmap;
use crate::db::error::{DbError, Result};
use crate::db::index::{IndexFile, SeqType};
use crate::db::header::{BlastDefLine, parse_def_line_set};
use crate::db::sequence::{get_protein_raw, get_nucleotide, decode_protein};
use crate::db::lmdb_v5::LmdbV5;
use crate::db::oid_seqids::OidSeqIds;
use crate::db::oid_taxids::OidTaxIds;

/// A single BLAST database volume.
struct Volume {
    index: IndexFile,
    seq_mmap: Mmap,
    hdr_mmap: Mmap,
    lmdb: Option<LmdbV5>,
    oid_seqids: Option<OidSeqIds>,
    oid_taxids: Option<OidTaxIds>,
}

/// Handle to an open BLAST database (one or more volumes).
///
/// Supports v4 and v5 formats, single and multi-volume databases, and alias files.
/// Database files are memory-mapped for efficient access.
///
/// ```no_run
/// use blast_rs::BlastDb;
/// use std::path::Path;
///
/// let db = BlastDb::open(Path::new("nr")).unwrap();
/// println!("{} sequences, {} total residues", db.num_sequences(), db.volume_length());
/// ```
pub struct BlastDb {
    volumes: Vec<Volume>,
    /// Cumulative OID count per volume: volume i covers OIDs [cum_oids[i], cum_oids[i+1]).
    cum_oids: Vec<u32>,
    /// Cached metadata from first volume / alias
    title_str: String,
    total_seqs: u32,
    total_length: u64,
    db_seq_type: SeqType,
    db_format_version: i32,
}

/// Append ".ext" to a path string, handling dotted stems correctly.
/// Unlike Path::with_extension, this always appends rather than replacing.
fn path_with_ext(base: &Path, ext: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", base.display(), ext))
}

impl BlastDb {
    /// Open a single volume from a base path (without extension).
    fn open_single_volume(path: &Path) -> Result<Volume> {
        let (is_protein, index_ext, seq_ext, hdr_ext) = if path_with_ext(path, "pin").exists() {
            (true, "pin", "psq", "phr")
        } else if path_with_ext(path, "nin").exists() {
            (false, "nin", "nsq", "nhr")
        } else {
            return Err(DbError::InvalidFormat(format!(
                "No .pin or .nin index file found at {}",
                path.display()
            )));
        };

        let index_data = fs::read(path_with_ext(path, index_ext))?;
        let index = IndexFile::parse(&index_data)?;

        let seq_file = fs::File::open(path_with_ext(path, seq_ext))?;
        let seq_mmap = unsafe { Mmap::map(&seq_file)? };

        let hdr_file = fs::File::open(path_with_ext(path, hdr_ext))?;
        let hdr_mmap = unsafe { Mmap::map(&hdr_file)? };

        let (lmdb, oid_seqids, oid_taxids) = if index.format_version == 5 {
            let (lmdb_ext, seqids_ext, taxids_ext) =
                if is_protein { ("pdb", "pos", "pot") }
                else          { ("ndb", "nos", "not") };

            let lmdb_path = path_with_ext(path, lmdb_ext);
            let lmdb = if lmdb_path.exists() {
                Some(LmdbV5::open(&lmdb_path)?)
            } else { None };

            let oid_seqids_path = path_with_ext(path, seqids_ext);
            let oid_seqids = if oid_seqids_path.exists() {
                Some(OidSeqIds::open(&oid_seqids_path)?)
            } else { None };

            let oid_taxids_path = path_with_ext(path, taxids_ext);
            let oid_taxids = if oid_taxids_path.exists() {
                Some(OidTaxIds::open(&oid_taxids_path)?)
            } else { None };

            (lmdb, oid_seqids, oid_taxids)
        } else {
            (None, None, None)
        };

        Ok(Volume { index, seq_mmap, hdr_mmap, lmdb, oid_seqids, oid_taxids })
    }

    /// Parse a BLAST alias file (.pal or .nal) and return volume paths.
    fn parse_alias_file(alias_path: &Path) -> Result<(Option<String>, Vec<PathBuf>)> {
        let content = fs::read_to_string(alias_path)?;
        let parent = alias_path.parent().unwrap_or(Path::new("."));
        let mut title = None;
        let mut db_list = Vec::new();

        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') { continue; }
            if let Some(val) = line.strip_prefix("TITLE ") {
                title = Some(val.trim().to_string());
            } else if let Some(val) = line.strip_prefix("DBLIST ") {
                for name in val.split_whitespace() {
                    let name = name.trim_matches('"');
                    let vol_path = if Path::new(name).is_absolute() {
                        PathBuf::from(name)
                    } else {
                        parent.join(name)
                    };
                    db_list.push(vol_path);
                }
            }
        }

        if db_list.is_empty() {
            return Err(DbError::InvalidFormat("Alias file has no DBLIST".into()));
        }
        Ok((title, db_list))
    }

    /// Open a BLAST database by base path (without extension).
    /// Supports single volumes, multi-volume databases, and alias files (.pal/.nal).
    pub fn open(path: &Path) -> Result<Self> {
        // Check for alias files first
        let pal = path.with_extension("pal");
        let nal = path.with_extension("nal");

        if pal.exists() || nal.exists() {
            let alias_path = if pal.exists() { pal } else { nal };
            let (alias_title, vol_paths) = Self::parse_alias_file(&alias_path)?;

            let mut volumes = Vec::new();
            for vp in &vol_paths {
                volumes.push(Self::open_single_volume(vp)?);
            }
            return Self::from_volumes(volumes, alias_title);
        }

        // Check for multi-volume: path.00.pin, path.01.pin, etc.
        // Note: we must check for the full filename (e.g. "db.00.pin") by
        // constructing the path directly, because Path::with_extension would
        // replace ".00" (treating it as an extension) instead of appending ".pin".
        let mut vol_paths = Vec::new();
        for i in 0..1000u32 {
            let vol_base = PathBuf::from(format!("{}.{:02}", path.display(), i));
            let pin_path = PathBuf::from(format!("{}.{:02}.pin", path.display(), i));
            let nin_path = PathBuf::from(format!("{}.{:02}.nin", path.display(), i));
            if pin_path.exists() || nin_path.exists() {
                vol_paths.push(vol_base);
            } else {
                break;
            }
        }

        if vol_paths.len() > 1 {
            let mut volumes = Vec::new();
            for vp in &vol_paths {
                volumes.push(Self::open_single_volume(vp)?);
            }
            return Self::from_volumes(volumes, None);
        }

        // Single volume
        let vol = Self::open_single_volume(path)?;
        Self::from_volumes(vec![vol], None)
    }

    fn from_volumes(volumes: Vec<Volume>, alias_title: Option<String>) -> Result<Self> {
        if volumes.is_empty() {
            return Err(DbError::InvalidFormat("No volumes found".into()));
        }

        let db_seq_type = volumes[0].index.seq_type;
        let db_format_version = volumes[0].index.format_version;

        let mut cum_oids = vec![0u32];
        let mut total_length = 0u64;
        let mut total_seqs = 0u32;

        for vol in &volumes {
            total_seqs += vol.index.num_oids;
            total_length += vol.index.volume_length;
            cum_oids.push(total_seqs);
        }

        let title_str = alias_title.unwrap_or_else(|| {
            volumes.iter().map(|v| v.index.title.as_str()).collect::<Vec<_>>().join("; ")
        });

        Ok(BlastDb {
            volumes,
            cum_oids,
            title_str,
            total_seqs,
            total_length,
            db_seq_type,
            db_format_version,
        })
    }

    /// Resolve a global OID to (volume_index, local_oid).
    #[inline]
    fn resolve_oid(&self, oid: u32) -> Result<(usize, u32)> {
        if oid >= self.total_seqs {
            return Err(DbError::OidOutOfRange(oid));
        }
        // Binary search for the right volume
        let mut lo = 0;
        let mut hi = self.volumes.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if oid < self.cum_oids[mid + 1] {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        let local_oid = oid - self.cum_oids[lo];
        Ok((lo, local_oid))
    }

    // ---- Metadata ----

    pub fn num_sequences(&self) -> u32 {
        self.total_seqs
    }

    pub fn seq_type(&self) -> SeqType {
        self.db_seq_type
    }

    pub fn volume_length(&self) -> u64 {
        self.total_length
    }

    pub fn title(&self) -> &str {
        &self.title_str
    }

    pub fn format_version(&self) -> i32 {
        self.db_format_version
    }

    pub fn is_v5(&self) -> bool {
        self.db_format_version == 5
    }

    /// Number of volumes in this database.
    pub fn num_volumes(&self) -> usize {
        self.volumes.len()
    }

    // ---- Sequence access ----

    pub fn get_sequence_protein_raw(&self, oid: u32) -> Result<&[u8]> {
        let (vi, local) = self.resolve_oid(oid)?;
        let vol = &self.volumes[vi];
        let start = vol.index.sequence_array[local as usize] as usize;
        let end = vol.index.sequence_array[local as usize + 1] as usize;
        Ok(get_protein_raw(&vol.seq_mmap, start, end))
    }

    pub fn get_sequence_protein(&self, oid: u32) -> Result<Vec<u8>> {
        Ok(decode_protein(self.get_sequence_protein_raw(oid)?))
    }

    pub fn get_sequence_nucleotide(&self, oid: u32) -> Result<Vec<u8>> {
        let (vi, local) = self.resolve_oid(oid)?;
        let vol = &self.volumes[vi];
        let ambig = vol.index.ambig_array.as_ref()
            .ok_or_else(|| DbError::InvalidFormat("No ambig_array for nucleotide db".into()))?;
        let seq_start = vol.index.sequence_array[local as usize] as usize;
        let seq_end = ambig[local as usize] as usize;
        let ambig_start = seq_end;
        let ambig_end = vol.index.sequence_array[local as usize + 1] as usize;
        Ok(get_nucleotide(&vol.seq_mmap, seq_start, seq_end, ambig_start, ambig_end))
    }

    // ---- Header access ----

    pub fn get_header(&self, oid: u32) -> Result<BlastDefLine> {
        let (vi, local) = self.resolve_oid(oid)?;
        let vol = &self.volumes[vi];
        let start = vol.index.header_array[local as usize] as usize;
        let end = vol.index.header_array[local as usize + 1] as usize;
        let data = &vol.hdr_mmap[start..end];
        let deflines = parse_def_line_set(data)?;
        Ok(deflines.into_iter().next().unwrap_or_default())
    }

    pub fn get_headers(&self, oid: u32) -> Result<Vec<BlastDefLine>> {
        let (vi, local) = self.resolve_oid(oid)?;
        let vol = &self.volumes[vi];
        let start = vol.index.header_array[local as usize] as usize;
        let end = vol.index.header_array[local as usize + 1] as usize;
        let data = &vol.hdr_mmap[start..end];
        parse_def_line_set(data)
    }

    // ---- V5 accession index ----

    pub fn lookup_accession(&self, accession: &str) -> Option<Result<Vec<u32>>> {
        // Search across all volumes
        let mut all_oids = Vec::new();
        let mut found = false;
        for (vi, vol) in self.volumes.iter().enumerate() {
            if let Some(ref lmdb) = vol.lmdb {
                found = true;
                match lmdb.get_oids_for_accession(accession) {
                    Ok(oids) => {
                        let offset = self.cum_oids[vi];
                        all_oids.extend(oids.iter().map(|&o| o + offset));
                    }
                    Err(e) => return Some(Err(e)),
                }
            }
        }
        if found { Some(Ok(all_oids)) } else { None }
    }

    pub fn iter_accessions<F>(&self, mut f: F) -> Option<Result<()>>
    where
        F: FnMut(&str, u32),
    {
        let mut found = false;
        for (vi, vol) in self.volumes.iter().enumerate() {
            if let Some(ref lmdb) = vol.lmdb {
                found = true;
                let offset = self.cum_oids[vi];
                let result = lmdb.iter_accessions(|acc, oid| {
                    f(acc, oid + offset);
                });
                if let Err(e) = result { return Some(Err(e)); }
            }
        }
        if found { Some(Ok(())) } else { None }
    }

    pub fn get_volumes_info(&self) -> Option<Result<Vec<(String, u32)>>> {
        // Try first volume's LMDB
        let vol = &self.volumes[0];
        let lmdb = vol.lmdb.as_ref()?;
        Some(lmdb.get_volumes_info())
    }

    // ---- V5 OID→SeqIDs ----

    pub fn get_seqids(&self, oid: u32) -> Option<Result<Vec<String>>> {
        let (vi, local) = self.resolve_oid(oid).ok()?;
        let r = self.volumes[vi].oid_seqids.as_ref()?;
        Some(r.get_seqids(local))
    }

    // ---- V5 OID→TaxIDs ----

    pub fn get_taxids(&self, oid: u32) -> Option<Result<Vec<i32>>> {
        let (vi, local) = self.resolve_oid(oid).ok()?;
        let r = self.volumes[vi].oid_taxids.as_ref()?;
        Some(r.get_taxids(local))
    }
}
