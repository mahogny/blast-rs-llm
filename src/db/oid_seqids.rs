//! Reader for the OID→SeqIDs file (`.pos` / `.nos`) used in BlastDB v5.
//!
//! File layout (all integers little-endian):
//!   [0..7]              Uint8   num_oids
//!   `[8..8+8*num_oids-1]` end-offset array: `end_offset(i)` = byte offset in the data
//!                                section where OID i's seq-id data ends
//!   [8+8*num_oids..]    data    concatenated seq-id records for each OID
//!
//! Each seq-id record:
//!   if first byte < 0xFF:  1 byte length + `length` ASCII bytes
//!   if first byte == 0xFF: 0xFF + 4 bytes Uint4 length + `length` ASCII bytes

use std::fs;
use std::path::Path;
use memmap2::Mmap;
use crate::db::error::{DbError, Result};

pub struct OidSeqIds {
    mmap: Mmap,
    num_oids: u64,
    /// Byte offset within mmap where the data section begins.
    data_offset: usize,
}

impl OidSeqIds {
    pub fn open(path: &Path) -> Result<Self> {
        let file = fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.len() < 8 {
            return Err(DbError::InvalidFormat("OID-to-seqids file too short".into()));
        }
        let num_oids = u64::from_le_bytes(mmap[0..8].try_into().unwrap());
        // Index array: num_oids × 8 bytes starting at offset 8.
        // Data section starts immediately after.
        let data_offset = 8 + 8 * num_oids as usize;
        if data_offset > mmap.len() {
            return Err(DbError::InvalidFormat(
                format!("OID-to-seqids file truncated: need {} bytes for index", data_offset)
            ));
        }
        Ok(OidSeqIds { mmap, num_oids, data_offset })
    }

    pub fn num_oids(&self) -> u64 {
        self.num_oids
    }

    /// Return all seq-id strings for the given OID.
    pub fn get_seqids(&self, oid: u32) -> Result<Vec<String>> {
        if oid as u64 >= self.num_oids {
            return Err(DbError::OidOutOfRange(oid));
        }
        let (begin, end) = self.data_range(oid);
        let data = &self.mmap[begin..end];
        Ok(parse_seqids(data))
    }

    fn data_range(&self, oid: u32) -> (usize, usize) {
        // end_offset[oid] is at index_base + oid * 8
        let index_base = 8usize;
        let end_off = {
            let off = index_base + oid as usize * 8;
            u64::from_le_bytes(self.mmap[off..off+8].try_into().unwrap()) as usize
        };
        let begin_off = if oid == 0 {
            0
        } else {
            let off = index_base + (oid as usize - 1) * 8;
            u64::from_le_bytes(self.mmap[off..off+8].try_into().unwrap()) as usize
        };
        (self.data_offset + begin_off, self.data_offset + end_off)
    }
}

fn parse_seqids(mut data: &[u8]) -> Vec<String> {
    let mut ids = Vec::new();
    while !data.is_empty() {
        let len_byte = data[0];
        data = &data[1..];
        let len = if len_byte == 0xFF {
            if data.len() < 4 { break; }
            let l = u32::from_le_bytes([data[0], data[1], data[2], data[3]]) as usize;
            data = &data[4..];
            l
        } else {
            len_byte as usize
        };
        if len > data.len() { break; }
        ids.push(String::from_utf8_lossy(&data[..len]).into_owned());
        data = &data[len..];
    }
    ids
}
