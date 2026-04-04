//! Performance benchmarks for blast-rs.
//!
//! Run with: cargo bench
//!
//! Benchmarks cover the main pipeline stages:
//! - Lookup table construction
//! - Ungapped extension
//! - Gapped extension
//! - Full blastp search (varying DB sizes)
//! - Full blastn search
//! - SEG/DUST masking
//! - Six-frame translation

use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use tempfile::TempDir;

use blast_rs::db::index::SeqType;
use blast_rs::{
    BlastDb, BlastDbBuilder, SequenceEntry, SearchParams, ScoringMatrix, MatrixType,
    blastp, blastn,
    apply_seg, apply_dust,
    six_frame_translate,
};
use blast_rs::lookup::{ProteinLookup, NucleotideLookup};
use blast_rs::extend::{ungapped_extend, gapped_extend, ungapped_extend_nucleotide};
use blast_rs::search::aa_to_ncbistdaa;

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Generate a deterministic pseudo-random protein sequence of given length.
fn random_protein(len: usize, seed: u64) -> Vec<u8> {
    let aa = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            aa[((state >> 33) as usize) % aa.len()]
        })
        .collect()
}

/// Generate a deterministic pseudo-random nucleotide sequence of given length.
fn random_nucleotide(len: usize, seed: u64) -> Vec<u8> {
    let nt = b"ACGT";
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            nt[((state >> 33) as usize) % nt.len()]
        })
        .collect()
}

/// Build a protein database with `n` sequences of ~300 aa each.
fn build_protein_db(n: usize) -> (TempDir, BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("benchdb");
    let mut builder = BlastDbBuilder::new(SeqType::Protein, "bench protein db");
    for i in 0..n {
        let seq = random_protein(300, i as u64);
        builder.add(SequenceEntry {
            title: format!("protein_{}", i),
            accession: format!("P{:06}", i),
            sequence: seq,
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = BlastDb::open(&base).unwrap();
    (tmp, db)
}

/// Build a nucleotide database with `n` sequences of ~1000 bp each.
fn build_nucleotide_db(n: usize) -> (TempDir, BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("benchdb");
    let mut builder = BlastDbBuilder::new(SeqType::Nucleotide, "bench nt db");
    for i in 0..n {
        let seq = random_nucleotide(1000, i as u64 + 10000);
        builder.add(SequenceEntry {
            title: format!("nt_seq_{}", i),
            accession: format!("N{:06}", i),
            sequence: seq,
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = BlastDb::open(&base).unwrap();
    (tmp, db)
}

// ── Lookup table benchmarks ─────────────────────────────────────────────────

fn bench_protein_lookup_build(c: &mut Criterion) {
    let query = random_protein(300, 42);
    let query_ncbi = aa_to_ncbistdaa(&query);
    let matrix = ScoringMatrix::from_type(MatrixType::Blosum62);

    c.bench_function("protein_lookup_build_300aa", |b| {
        b.iter(|| ProteinLookup::build(&query_ncbi, 3, &matrix, 11));
    });
}

fn bench_nucleotide_lookup_build(c: &mut Criterion) {
    let query = random_nucleotide(1000, 42);

    c.bench_function("nucleotide_lookup_build_1000bp", |b| {
        b.iter(|| NucleotideLookup::build(&query, 11));
    });
}

// ── Extension benchmarks ────────────────────────────────────────────────────

fn bench_ungapped_extend(c: &mut Criterion) {
    let query = aa_to_ncbistdaa(&random_protein(300, 42));
    let subject = aa_to_ncbistdaa(&random_protein(300, 99));
    let matrix = ScoringMatrix::from_type(MatrixType::Blosum62);

    c.bench_function("ungapped_extend_300aa", |b| {
        b.iter(|| ungapped_extend(&query, &subject, 150, 150, &matrix, 7));
    });
}

fn bench_gapped_extend(c: &mut Criterion) {
    // Use identical sequences to force a full-length alignment
    let seq = aa_to_ncbistdaa(&random_protein(300, 42));
    let matrix = ScoringMatrix::from_type(MatrixType::Blosum62);

    c.bench_function("gapped_extend_300aa_identical", |b| {
        b.iter(|| gapped_extend(&seq, &seq, 150, 150, &matrix, 11, 1, 25));
    });
}

fn bench_ungapped_extend_nucleotide(c: &mut Criterion) {
    let query = random_nucleotide(1000, 42);
    let subject = random_nucleotide(1000, 99);

    c.bench_function("ungapped_extend_nt_1000bp", |b| {
        b.iter(|| ungapped_extend_nucleotide(&query, &subject, 500, 500, 2, -3, 20));
    });
}

// ── Full search benchmarks ──────────────────────────────────────────────────

fn bench_blastp(c: &mut Criterion) {
    let mut group = c.benchmark_group("blastp");
    group.sample_size(10);

    let query = random_protein(300, 42);

    for db_size in [100, 1000] {
        let (_tmp, db) = build_protein_db(db_size);
        let params = SearchParams::blastp()
            .evalue(10.0)
            .num_threads(1)
            .filter_low_complexity(false)
            .comp_adjust(false);

        group.bench_with_input(
            BenchmarkId::new("single_thread", db_size),
            &db_size,
            |b, _| {
                b.iter(|| blastp(&db, &query, &params));
            },
        );
    }

    // Multi-threaded with 1000 seqs
    {
        let (_tmp, db) = build_protein_db(1000);
        let params = SearchParams::blastp()
            .evalue(10.0)
            .num_threads(4)
            .filter_low_complexity(false)
            .comp_adjust(false);

        group.bench_function("4_threads_1000", |b| {
            b.iter(|| blastp(&db, &query, &params));
        });
    }

    group.finish();
}

fn bench_blastn(c: &mut Criterion) {
    let mut group = c.benchmark_group("blastn");
    group.sample_size(10);

    let query = random_nucleotide(500, 42);

    for db_size in [100, 1000] {
        let (_tmp, db) = build_nucleotide_db(db_size);
        let params = SearchParams::blastn()
            .evalue(10.0)
            .num_threads(1)
            .filter_low_complexity(false);

        group.bench_with_input(
            BenchmarkId::new("single_thread", db_size),
            &db_size,
            |b, _| {
                b.iter(|| blastn(&db, &query, &params));
            },
        );
    }

    group.finish();
}

// ── Masking benchmarks ──────────────────────────────────────────────────────

fn bench_seg_masking(c: &mut Criterion) {
    let seq = random_protein(1000, 42);

    c.bench_function("seg_mask_1000aa", |b| {
        b.iter(|| {
            let mut s = seq.clone();
            apply_seg(&mut s);
        });
    });
}

fn bench_dust_masking(c: &mut Criterion) {
    let seq = random_nucleotide(10000, 42);

    c.bench_function("dust_mask_10000bp", |b| {
        b.iter(|| {
            let mut s = seq.clone();
            apply_dust(&mut s);
        });
    });
}

// ── Translation benchmark ───────────────────────────────────────────────────

fn bench_six_frame_translate(c: &mut Criterion) {
    let seq = random_nucleotide(3000, 42);

    c.bench_function("six_frame_translate_3000bp", |b| {
        b.iter(|| six_frame_translate(&seq));
    });
}

// ── Database I/O benchmark ──────────────────────────────────────────────────

fn bench_db_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("db_build");
    group.sample_size(10);

    let entries: Vec<SequenceEntry> = (0..1000)
        .map(|i| SequenceEntry {
            title: format!("protein_{}", i),
            accession: format!("P{:06}", i),
            sequence: random_protein(300, i as u64),
            taxid: None,
        })
        .collect();

    group.bench_function("protein_1000_seqs", |b| {
        b.iter(|| {
            let tmp = TempDir::new().unwrap();
            let base = tmp.path().join("benchdb");
            let mut builder = BlastDbBuilder::new(SeqType::Protein, "bench");
            for e in &entries {
                builder.add(SequenceEntry {
                    title: e.title.clone(),
                    accession: e.accession.clone(),
                    sequence: e.sequence.clone(),
                    taxid: None,
                });
            }
            builder.write(&base).unwrap();
        });
    });

    group.finish();
}

// ── Criterion groups ────────────────────────────────────────────────────────

criterion_group!(
    lookup,
    bench_protein_lookup_build,
    bench_nucleotide_lookup_build,
);

criterion_group!(
    extension,
    bench_ungapped_extend,
    bench_gapped_extend,
    bench_ungapped_extend_nucleotide,
);

criterion_group!(
    search,
    bench_blastp,
    bench_blastn,
);

criterion_group!(
    masking,
    bench_seg_masking,
    bench_dust_masking,
);

criterion_group!(
    misc,
    bench_six_frame_translate,
    bench_db_build,
);

criterion_main!(lookup, extension, search, masking, misc);
