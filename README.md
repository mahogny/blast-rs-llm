# blast-rs

A pure-Rust implementation of BLAST (Basic Local Alignment Search Tool). Reads and writes databases created by NCBI `makeblastdb`, implements the core BLAST algorithm for protein and nucleotide search, and produces output in all standard BLAST formats. No dependency on the NCBI C++ toolkit.

This implementation passes many tests, showing similiar behavior to the original BLAST. Further randomized testing would however be beneficial

The aim of this code is to offer a stable version that can be called as a library. It is however currently up to 20% slower than regular BLAST. Further optimization is needed to fully replace BLAST


## Features

- Read BLAST databases v4 and v5 (protein and nucleotide)
- Build BLAST databases from FASTA input (`makeblastdb`)
- `blastp` — protein-protein search with BLOSUM/PAM scoring matrices
- `blastn` — nucleotide-nucleotide search
- Parallel search via Rayon
- All standard BLAST output formats (0–18)
- Custom tabular column selection

## Using blast-rs as a library

`blast-rs` exposes a high-level Rust API so you can embed BLAST searches in your own program without shelling out to the CLI.

Add to `Cargo.toml`:

```toml
[dependencies]
blast-rs = { git = "https://github.com/henriksson-lab/blast-rs" }
```

### Protein search (blastp)

```rust
use blast_rs::{BlastDb, SearchParams, blastp, parse_fasta};

fn main() -> anyhow::Result<()> {
    let db = BlastDb::open("mydb".as_ref())?;

    let fasta = std::fs::read("query.faa")?;
    let sequences = parse_fasta(&fasta);

    let params = SearchParams::blastp()
        .evalue(1e-5)
        .max_target_seqs(50)
        .num_threads(8);

    for (title, seq) in &sequences {
        let results = blastp(&db, seq, &params);
        for r in &results {
            println!("{}\t{}\t{:.2e}", title, r.subject_accession, r.best_evalue());
        }
    }
    Ok(())
}
```

### Nucleotide search (blastn)

```rust
use blast_rs::{BlastDb, SearchParams, blastn};

let db = BlastDb::open("nt".as_ref())?;
let params = SearchParams::blastn()
    .evalue(1e-10)
    .word_size(15)
    .match_score(2)
    .mismatch(-3);

let results = blastn(&db, b"ATGCGTACGTAGCTAGC", &params);
```

### Translated search (blastx / tblastn / tblastx)

```rust
use blast_rs::{BlastDb, SearchParams, blastx, tblastn, tblastx};

// Nucleotide query vs protein database (6-frame translation of query)
let results = blastx(&db, nt_query, &SearchParams::blastx().evalue(1e-5));

// Protein query vs nucleotide database (6-frame translation of subjects)
let results = tblastn(&db, aa_query, &SearchParams::tblastn().evalue(1e-5));

// Both sides translated (expensive)
let results = tblastx(&db, nt_query, &SearchParams::tblastx());
```

### Iterative search with PSSM (psiblast)

```rust
use blast_rs::{BlastDb, SearchParams, psiblast, PsiblastParams};

let db = BlastDb::open("mydb".as_ref())?;
let search = SearchParams::blastp().evalue(10.0).matrix(blast_rs::MatrixType::Blosum62);

let params = PsiblastParams::new(search)
    .num_iterations(3)
    .inclusion_evalue(0.001);

let (results, pssm) = psiblast(&db, b"MKTLLLTLVV...", &params);
// `pssm` can be used for subsequent custom searches via `blast_rs::search_with_pssm`
```

### Parsing FASTA in memory

```rust
use blast_rs::parse_fasta;

let input = std::fs::read("sequences.faa")?;
for (title, seq) in parse_fasta(&input) {
    println!("{}: {} residues", title, seq.len());
}
```

### SearchParams builder reference

| Method | Type | Description |
|--------|------|-------------|
| `.evalue(f64)` | all | E-value threshold (default 10.0) |
| `.max_target_seqs(usize)` | all | Max hits returned (default 500) |
| `.num_threads(usize)` | all | Worker threads; 0 = all (default 0) |
| `.word_size(usize)` | all | Seed word length |
| `.gap_open(i32)` | protein | Gap open penalty |
| `.gap_extend(i32)` | protein | Gap extend penalty |
| `.matrix(MatrixType)` | protein | Scoring matrix |
| `.filter_low_complexity(bool)` | all | SEG/DUST masking (default true) |
| `.comp_adjust(bool)` | protein | Composition-based statistics (default true) |
| `.match_score(i32)` | blastn | Match reward (default 2) |
| `.mismatch(i32)` | blastn | Mismatch penalty (default -3) |

## Building

Requires Rust 1.70+ and a system LMDB library (for v5 database support).

```sh
cargo build --release
```

The binary is at `target/release/blast-cli`.

## Usage

### Build a database

```sh
# Protein database
blast-cli makeblastdb -i sequences.faa --dbtype prot -o mydb

# Nucleotide database
blast-cli makeblastdb -i sequences.fna --dbtype nucl -o mydb --title "My genome DB"
```

FASTA headers are parsed as `>accession description` (the first whitespace-delimited word becomes the accession). IUPAC ambiguity codes (N, R, Y, M, K, S, W, H, B, V, D) are preserved in nucleotide databases.

### Protein search (blastp)

```sh
blast-cli blastp -q query.faa -d mydb
```

### Nucleotide search (blastn)

```sh
blast-cli blastn -q query.fna -d mydb
```

### Common options

| Flag | Default | Description |
|------|---------|-------------|
| `-q / --query` | — | Query FASTA file |
| `-d / --db` | — | Database base path (no extension) |
| `-o / --out` | stdout | Output file |
| `--evalue` | 10 | E-value threshold |
| `--outfmt` | 0 | Output format (see below) |
| `--max-target-seqs` | 500 | Maximum hits returned |
| `--num_threads` | 0 (all) | Worker threads |
| `--matrix` | BLOSUM62 | Scoring matrix (`blastp` only) |
| `--gapopen` | matrix default | Gap open penalty |
| `--gapextend` | matrix default | Gap extend penalty |
| `--word-size` | 3 (prot) / 11 (nucl) | Seed word size |
| `--reward` | 2 | Match reward (`blastn` only) |
| `--penalty` | -3 | Mismatch penalty (`blastn` only) |

### Output formats

`--outfmt` accepts a format number, optionally followed by column names for tabular formats:

```sh
blast-cli blastp -q query.faa -d mydb --outfmt 6
blast-cli blastp -q query.faa -d mydb --outfmt "6 qseqid sseqid pident evalue bitscore"
blast-cli blastp -q query.faa -d mydb --outfmt 15   # JSON
blast-cli blastp -q query.faa -d mydb --outfmt 5    # BLAST XML
```

| Format | Description |
|--------|-------------|
| 0 | Pairwise (default) |
| 1 | Query-anchored with identities |
| 2 | Query-anchored, no identities |
| 3 | Flat query-anchored with identities |
| 4 | Flat query-anchored, no identities |
| 5 | BLAST XML |
| 6 | Tabular |
| 7 | Tabular with comment lines |
| 8 | Seqalign (text ASN.1) |
| 9 | Seqalign (binary ASN.1) |
| 10 | Comma-separated (CSV) |
| 11 | BLAST archive (ASN.1) |
| 12 | Seqalign (JSON) |
| 13 | Multiple-file BLAST JSON |
| 14 | Multiple-file BLAST XML2 |
| 15 | Single-file BLAST JSON |
| 16 | Single-file BLAST XML2 |
| 18 | Organism Report |

Tabular columns available: `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident positive gaps ppos qseq sseq btop staxid salltitles qcovs qcovhsp score`

### Inspect a database

```sh
# Dump all sequences as FASTA
blast-cli dumpdb -d mydb

# Show headers only
blast-cli dumpdb -d mydb --headers-only

# v5 databases: list all accessions
blast-cli dumpdb -d mydb --list-accessions

# v5 databases: look up OIDs for an accession
blast-cli dumpdb -d mydb --lookup NP_001234

# v5 databases: show volume info
blast-cli dumpdb -d mydb --volumes
```

## Database format support

| Feature | Supported |
|---------|-----------|
| v4 protein (`.pin / .psq / .phr`) | read + write |
| v4 nucleotide (`.nin / .nsq / .nhr`) | read + write |
| v5 LMDB accession index (`.pdb / .ndb`) | read + write |
| v5 OID→SeqIds (`.pos / .nos`) | read + write |
| v5 OID→TaxIds (`.pot / .not`) | read + write |
| Multi-volume databases | not yet |

## Scoring matrices

BLOSUM45, BLOSUM62, BLOSUM80, PAM30, PAM70, PAM250.
Karlin-Altschul parameters are pre-computed for standard gap penalty combinations.

## Dependencies

| Crate | Purpose |
|-------|---------|
| `byteorder` | Big/little-endian binary I/O |
| `memmap2` | Memory-mapped database files |
| `lmdb` | v5 LMDB accession index |
| `rayon` | Parallel search |
| `clap` | CLI argument parsing |
| `thiserror` | Error types |

## Performance benchmarks

Measured with `cargo bench` (Criterion). Database sequences are ~300aa (protein) or ~1000bp (nucleotide), randomly generated.

| Benchmark | Time | Description |
|-----------|------|-------------|
| **Lookup table construction** | | |
| protein_lookup_build (300aa) | 1.2 ms | Neighbor word enumeration, BLOSUM62, word_size=3 |
| nucleotide_lookup_build (1000bp) | 75 ms | Exact word hashing, word_size=11 |
| **Extension** | | |
| ungapped_extend (300aa) | 49 ns | X-drop ungapped protein extension |
| gapped_extend (300aa, identical) | 72 µs | Banded Smith-Waterman DP with traceback |
| ungapped_extend_nt (1000bp) | 77 ns | X-drop ungapped nucleotide extension |
| **End-to-end search** | | |
| blastp, 100 seqs, 1 thread | 3.4 ms | 300aa query vs 100×300aa DB |
| blastp, 1000 seqs, 1 thread | 3.8 ms | 300aa query vs 1000×300aa DB |
| blastp, 1000 seqs, 4 threads | 3.7 ms | Same, with Rayon parallelism |
| blastn, 100 seqs, 1 thread | 176 ms | 500bp query vs 100×1000bp DB |
| blastn, 1000 seqs, 1 thread | 176 ms | 500bp query vs 1000×1000bp DB |
| **Masking** | | |
| SEG (1000aa) | 205 µs | Low-complexity protein masking |
| DUST (10000bp) | 2.5 ms | Low-complexity nucleotide masking |
| **Misc** | | |
| six_frame_translate (3000bp) | 226 µs | All 6 reading frames |
| db_build (1000 protein seqs) | 2.1 ms | Write v4 database (1000×300aa) |

Run benchmarks yourself:

```sh
cargo bench
```

HTML reports are generated in `target/criterion/`.

### Comparison with NCBI BLAST+ (blastp)

Measured on synthetic protein databases (sequences ~300aa). Both tools find the same number of hits. NCBI BLAST+ 2.17.0.

**Single-threaded:**

| Query | DB size | NCBI BLAST+ | blast-rs | Ratio | Hits |
|-------|---------|-------------|----------|-------|------|
| 50aa  | 100     | 0.06 s      | 0.01 s   | **0.2x (5x faster)** | 5 |
| 50aa  | 1,000   | 0.07 s      | 0.04 s   | **0.6x** | 5 |
| 50aa  | 10,000  | 0.11 s      | 0.11 s   | 1.0x (parity)  | 5 |
| 300aa | 100     | 0.07 s      | 0.02 s   | **0.3x (3x faster)** | 5 |
| 300aa | 1,000   | 0.08 s      | 0.08 s   | 1.0x (parity) | 5 |
| 300aa | 10,000  | 0.24 s      | 0.33 s   | 1.4x   | 5 |
| 1000aa| 100     | 0.08 s      | 0.06 s   | **0.8x** | 20 |
| 1000aa| 1,000   | 0.12 s      | 0.14 s   | 1.2x   | 20 |
| 1000aa| 10,000  | 0.64 s      | 0.74 s   | 1.2x   | 20 |

**Multi-threaded scaling (1000aa query, 10,000 seq DB):**

| Threads | NCBI BLAST+ | blast-rs | Ratio |
|---------|-------------|----------|-------|
| 1       | 0.64 s      | 0.74 s   | 1.2x  |
| 2       | 0.32 s      | 0.42 s   | 1.3x  |
| 4       | 0.21 s      | 0.26 s   | 1.2x  |
| 8       | 0.16 s      | 0.16 s   | **1.0x (parity)** |

**Memory usage (blast-rs uses 4-8x less):**

| | NCBI BLAST+ | blast-rs |
|--|-------------|----------|
| Typical | 36-39 MB | 5-10 MB |

blast-rs is faster than NCBI BLAST+ for small/medium databases and within 1.2x for large databases (10,000 sequences). The gap closes to parity with multi-threading. Memory usage is 4-8x lower. Both tools produce identical hit counts

## License

Dual-licensed under MIT or Public Domain (Unlicense), at your option. See [LICENSE](LICENSE).
