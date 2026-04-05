//! Reader for the OID→TaxIDs file (`.pot` / `.not`) used in BlastDB v5.
//!
//! File layout (all integers little-endian):
//!   [0..7]              Uint8   num_oids
//!   `[8..8+8*num_oids-1]` end-offset array: `end_offset(i)` = offset in the data
//!                                section (in units of i32, i.e. 4-byte elements)
//!                                where OID i's tax-id data ends
//!   [8+8*num_oids..]    i32[]   concatenated tax-id lists (little-endian i32)

use std::fs;
use std::path::Path;
use memmap2::Mmap;
use crate::db::error::{DbError, Result};

pub struct OidTaxIds {
    mmap: Mmap,
    num_oids: u64,
    data_offset: usize,
}

impl OidTaxIds {
    pub fn open(path: &Path) -> Result<Self> {
        let file = fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.len() < 8 {
            return Err(DbError::InvalidFormat("OID-to-taxids file too short".into()));
        }
        let num_oids = u64::from_le_bytes(mmap[0..8].try_into().unwrap());
        let data_offset = 8 + 8 * num_oids as usize;
        if data_offset > mmap.len() {
            return Err(DbError::InvalidFormat(
                "OID-to-taxids file truncated".into()
            ));
        }
        Ok(OidTaxIds { mmap, num_oids, data_offset })
    }

    pub fn num_oids(&self) -> u64 {
        self.num_oids
    }

    /// Return all taxonomy IDs for the given OID.
    pub fn get_taxids(&self, oid: u32) -> Result<Vec<i32>> {
        if oid as u64 >= self.num_oids {
            return Err(DbError::OidOutOfRange(oid));
        }
        let (begin_elem, end_elem) = self.elem_range(oid);
        let begin = self.data_offset + begin_elem * 4;
        let end = self.data_offset + end_elem * 4;
        if end > self.mmap.len() {
            return Ok(vec![]);
        }
        let data = &self.mmap[begin..end];
        let taxids = data.chunks_exact(4)
            .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect();
        Ok(taxids)
    }

    fn elem_range(&self, oid: u32) -> (usize, usize) {
        let index_base = 8usize;
        let end_elem = {
            let off = index_base + oid as usize * 8;
            u64::from_le_bytes(self.mmap[off..off+8].try_into().unwrap()) as usize
        };
        let begin_elem = if oid == 0 {
            0
        } else {
            let off = index_base + (oid as usize - 1) * 8;
            u64::from_le_bytes(self.mmap[off..off+8].try_into().unwrap()) as usize
        };
        (begin_elem, end_elem)
    }
}
