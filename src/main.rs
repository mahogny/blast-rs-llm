mod output;

use std::path::{Path, PathBuf};
use std::fs;
use std::io::{self, BufRead, Write};
use clap::{Parser, Subcommand};
use blast_rs::{BlastDb, BlastDbBuilder, SequenceEntry};
use blast_rs::db::index::SeqType;
use blast_rs::{
    SearchParams,
    matrix::MatrixType,
    api::{blastp as api_blastp, blastn as api_blastn, blastx as api_blastx,
          tblastn as api_tblastn, tblastx as api_tblastx,
          psiblast as api_psiblast, PsiblastParams},
};

#[derive(Parser)]
#[command(name = "blast-cli", about = "Pure-Rust BLAST implementation")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Protein-protein BLAST search
    Blastp(BlastArgs),
    /// Nucleotide-nucleotide BLAST search
    Blastn(BlastArgs),
    /// Translate nucleotide query, search protein database
    Blastx(BlastArgs),
    /// Protein query against translated nucleotide database
    Tblastn(BlastArgs),
    /// Translated nucleotide query vs translated nucleotide database
    Tblastx(BlastArgs),
    /// Iterative protein search using position-specific scoring matrix (PSI-BLAST)
    Psiblast(PsiblastArgs),
    /// Dump sequences from a BLAST database (for testing)
    Dumpdb(DumpArgs),
    /// Build a BLAST database from a FASTA file
    Makeblastdb(MakeDbArgs),
    /// Retrieve sequences or information from a BLAST database
    Blastdbcmd(BlastdbcmdArgs),
    /// Create a BLAST database alias
    BlastdbAliastool(AliasArgs),
    /// Re-format BLAST archive (format 11) output into other formats
    BlastFormatter(FormatterArgs),
}

#[derive(clap::Args)]
struct BlastArgs {
    /// Query FASTA file
    #[arg(short = 'q', long)]
    query: PathBuf,
    /// Database path (without extension)
    #[arg(short, long)]
    db: PathBuf,
    /// Output file (default: stdout)
    #[arg(short, long)]
    out: Option<PathBuf>,
    /// E-value threshold
    #[arg(long, default_value = "10")]
    evalue: f64,
    /// Output format: 0=pairwise, 5=XML, 6=tabular, 7=tabular+comments, 10=CSV, 15=JSON, 16=XML2, 17=FASTA, etc.
    /// For tabular formats, columns can be specified: e.g. "6 qseqid sseqid pident evalue"
    #[arg(long = "outfmt", default_value = "0")]
    outfmt: String,
    /// Maximum target sequences
    #[arg(long, default_value = "500")]
    max_target_seqs: usize,
    /// Scoring matrix (blastp only)
    #[arg(long, default_value = "BLOSUM62")]
    matrix: String,
    /// Gap open penalty
    #[arg(long = "gapopen")]
    gap_open: Option<i32>,
    /// Gap extend penalty
    #[arg(long = "gapextend")]
    gap_extend: Option<i32>,
    /// Word size
    #[arg(long)]
    word_size: Option<usize>,
    /// Number of threads (0 = all available)
    #[arg(long = "num_threads", default_value = "0")]
    num_threads: usize,
    /// Match score (blastn/blastx/tblastn only)
    #[arg(long = "reward", default_value = "2")]
    match_score: i32,
    /// Mismatch penalty (blastn only, positive value will be negated)
    #[arg(long = "penalty", default_value = "-3")]
    mismatch: i32,
    /// Disable low-complexity (SEG/DUST) filtering
    #[arg(long = "no-lc-filter")]
    no_lc_filter: bool,
    /// Disable composition-based statistics adjustment
    /// Composition-based statistics mode: 0=off, 1=unconditional (default), 2=conditional, 3=forced
    #[arg(long = "comp-based-stats", default_value = "1")]
    comp_based_stats: u8,

    /// Disable composition-based statistics (alias for --comp-based-stats 0)
    #[arg(long = "no-comp-adjust")]
    no_comp_adjust: bool,
    /// Query strand: both, plus, minus (nucleotide searches only)
    #[arg(long, default_value = "both")]
    strand: String,
    /// Genetic code for query (blastx/tblastx)
    #[arg(long = "query_gencode", default_value = "1")]
    query_gencode: u8,
    /// Genetic code for database (tblastn/tblastx)
    #[arg(long = "db_gencode", default_value = "1")]
    db_gencode: u8,
    /// Maximum HSPs per subject
    #[arg(long = "max_hsps")]
    max_hsps: Option<usize>,
    /// Culling limit
    #[arg(long = "culling_limit")]
    culling_limit: Option<usize>,
    /// Search task preset (e.g. megablast, blastn-short, blastp-fast)
    #[arg(long)]
    task: Option<String>,
    /// Restrict search to sequences with these taxonomy IDs (comma-separated)
    #[arg(long)]
    taxids: Option<String>,
    /// File of taxonomy IDs to restrict search (one per line)
    #[arg(long)]
    taxidlist: Option<PathBuf>,
    /// File of sequence IDs to restrict search (one per line)
    #[arg(long)]
    seqidlist: Option<PathBuf>,
    /// File of sequence IDs to EXCLUDE from search (one per line)
    #[arg(long)]
    negative_seqidlist: Option<PathBuf>,
    /// X-dropoff for ungapped extensions
    #[arg(long = "xdrop_ungap")]
    xdrop_ungap: Option<i32>,
    /// X-dropoff for preliminary gapped extensions
    #[arg(long = "xdrop_gap")]
    xdrop_gap: Option<i32>,
    /// X-dropoff for final gapped extensions
    #[arg(long = "xdrop_gap_final")]
    xdrop_gap_final: Option<i32>,
    /// Soft masking (use masking for lookup table only, not extensions)
    #[arg(long = "soft_masking")]
    soft_masking: bool,
    /// Number of one-line descriptions to show
    #[arg(long = "num_descriptions")]
    num_descriptions: Option<usize>,
    /// Number of alignments to show
    #[arg(long = "num_alignments")]
    num_alignments: Option<usize>,
    /// Parse deflines in query input
    #[arg(long = "lcase_masking")]
    lcase_masking: bool,
    /// Word score threshold for protein lookup table
    #[arg(long = "threshold")]
    threshold: Option<i32>,
    /// Window size for 2-hit algorithm
    #[arg(long = "window_size")]
    window_size: Option<usize>,
    /// GI list file to restrict search
    #[arg(long)]
    gilist: Option<PathBuf>,
    /// Negative GI list file to exclude from search
    #[arg(long)]
    negative_gilist: Option<PathBuf>,
    /// Produce HTML output
    #[arg(long)]
    html: bool,
    /// Export search strategy to file
    #[arg(long = "export_search_strategy")]
    export_strategy: Option<PathBuf>,
    /// Import search strategy from file
    #[arg(long = "import_search_strategy")]
    import_strategy: Option<PathBuf>,
    /// Discontiguous megablast template type (0=coding, 1=optimal)
    #[arg(long = "template_type")]
    template_type: Option<u8>,
    /// Discontiguous megablast template length (18 or 21)
    #[arg(long = "template_length")]
    template_length: Option<usize>,
}

#[derive(clap::Args)]
struct DumpArgs {
    /// Database path (without extension)
    #[arg(short, long)]
    db: PathBuf,
    /// Max sequences to dump (0 = all)
    #[arg(long, default_value = "0")]
    max: usize,
    /// Show headers only (no sequences)
    #[arg(long)]
    headers_only: bool,
    /// List all accessions from the v5 LMDB index
    #[arg(long)]
    list_accessions: bool,
    /// Look up OIDs for a specific accession (v5 only)
    #[arg(long)]
    lookup: Option<String>,
    /// Show v5 volume info from LMDB
    #[arg(long)]
    volumes: bool,
}

#[derive(clap::Args)]
struct PsiblastArgs {
    /// Query FASTA file (protein)
    #[arg(short = 'q', long)]
    query: PathBuf,
    /// Database path (without extension)
    #[arg(short, long)]
    db: PathBuf,
    /// Output file (default: stdout)
    #[arg(short, long)]
    out: Option<PathBuf>,
    /// E-value threshold for reporting hits
    #[arg(long, default_value = "10")]
    evalue: f64,
    /// E-value threshold for including hits in PSSM construction
    #[arg(long, default_value = "0.001")]
    inclusion_ethresh: f64,
    /// Number of PSI-BLAST iterations
    #[arg(long = "num_iterations", default_value = "3")]
    num_iterations: u32,
    /// Output format
    #[arg(long = "outfmt", default_value = "0")]
    outfmt: String,
    /// Maximum target sequences
    #[arg(long, default_value = "500")]
    max_target_seqs: usize,
    /// Scoring matrix
    #[arg(long, default_value = "BLOSUM62")]
    matrix: String,
    /// Gap open penalty
    #[arg(long = "gapopen")]
    gap_open: Option<i32>,
    /// Gap extend penalty
    #[arg(long = "gapextend")]
    gap_extend: Option<i32>,
    /// Number of threads (0 = all)
    #[arg(long = "num_threads", default_value = "0")]
    num_threads: usize,
    /// Save PSSM checkpoint file (binary)
    #[arg(long = "out_pssm")]
    out_pssm: Option<PathBuf>,
    /// Save ASCII PSSM file
    #[arg(long = "out_ascii_pssm")]
    out_ascii_pssm: Option<PathBuf>,
    /// Load PSSM checkpoint file to use instead of query
    #[arg(long = "in_pssm")]
    in_pssm: Option<PathBuf>,
}

#[derive(clap::Args)]
struct BlastdbcmdArgs {
    /// Database path (without extension)
    #[arg(short, long)]
    db: PathBuf,
    /// Accession or sequence identifier to retrieve
    #[arg(long)]
    entry: Option<String>,
    /// File of accessions to retrieve (one per line)
    #[arg(long)]
    entry_batch: Option<PathBuf>,
    /// Output file (default: stdout)
    #[arg(short, long)]
    out: Option<PathBuf>,
    /// Show database info
    #[arg(long)]
    info: bool,
    /// Output format: fasta (default), accession
    #[arg(long = "outfmt", default_value = "fasta")]
    outfmt: String,
    /// Range of sequence to extract (e.g. "10-50")
    #[arg(long)]
    range: Option<String>,
}

#[derive(clap::Args)]
struct MakeDbArgs {
    /// Input FASTA file
    #[arg(short = 'i', long = "in")]
    input: PathBuf,
    /// Sequence type: prot or nucl
    #[arg(long, default_value = "prot")]
    dbtype: String,
    /// Output database base path (without extension)
    #[arg(short, long)]
    out: PathBuf,
    /// Database title (defaults to output path)
    #[arg(long)]
    title: Option<String>,
    /// Parse FASTA headers as "accession description" (first word = accession)
    #[arg(long, default_value = "true")]
    parse_seqids: bool,
    /// Database format version (4 or 5)
    #[arg(long, default_value = "4")]
    dbversion: u32,
    /// Maximum file size in bytes before splitting into volumes (0 = no splitting)
    #[arg(long = "max-file-size", default_value = "0")]
    max_file_size: u64,
}

#[derive(clap::Args)]
struct AliasArgs {
    /// List of database paths to combine (space-separated)
    #[arg(long)]
    dblist: String,
    /// Output alias file path (without extension)
    #[arg(short, long)]
    out: PathBuf,
    /// Database type: prot or nucl
    #[arg(long, default_value = "prot")]
    dbtype: String,
    /// Title for the alias database
    #[arg(long)]
    title: Option<String>,
}

#[derive(clap::Args)]
struct FormatterArgs {
    /// Input BLAST archive file (format 11 output)
    #[arg(long = "archive", short = 'a')]
    archive: PathBuf,
    /// Output format (same as --outfmt in search commands)
    #[arg(long = "outfmt", default_value = "0")]
    outfmt: String,
    /// Output file (default: stdout)
    #[arg(short, long)]
    out: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Blastp(args)      => run_blastp(args),
        Commands::Blastn(args)      => run_blastn(args),
        Commands::Blastx(args)      => run_blastx(args),
        Commands::Tblastn(args)     => run_tblastn(args),
        Commands::Tblastx(args)     => run_tblastx(args),
        Commands::Psiblast(args)    => run_psiblast(args),
        Commands::Dumpdb(args)      => run_dumpdb(args),
        Commands::Makeblastdb(args) => run_makeblastdb(args),
        Commands::Blastdbcmd(args)  => run_blastdbcmd(args),
        Commands::BlastdbAliastool(args) => run_aliastool(args),
        Commands::BlastFormatter(args) => run_blast_formatter(args),
    }
}

fn apply_task_preset(params: &mut SearchParams, task: &str, program: &str) {
    match (program, task) {
        // blastn presets
        ("blastn", "megablast") => {
            params.word_size = 28;
            params.two_hit = false;
        }
        ("blastn", "dc-megablast") => {
            params.word_size = 11;
        }
        ("blastn", "blastn-short") => {
            params.word_size = 7;
            params.match_score = 1;
            params.mismatch = -3;
            params.evalue_threshold = 1000.0;
        }
        // blastp presets
        ("blastp", "blastp-short") => {
            params.word_size = 2;
            params.evalue_threshold = 20000.0;
            params.comp_adjust = 0;
        }
        ("blastp", "blastp-fast") => {
            params.word_size = 6;
        }
        // blastx presets
        ("blastx", "blastx-fast") => {
            params.word_size = 6;
        }
        // tblastn presets
        ("tblastn", "tblastn-fast") => {
            params.word_size = 6;
        }
        _ => {
            eprintln!("Warning: unknown task '{}' for program '{}', using defaults", task, program);
        }
    }
}

fn build_oid_filter(db: &BlastDb, args: &BlastArgs) -> Option<std::collections::HashSet<u32>> {
    use std::collections::HashSet;

    // Collect GI list filter
    let gi_filter: Option<HashSet<i64>> = if let Some(ref gi_file) = args.gilist {
        let content = fs::read_to_string(gi_file).unwrap_or_else(|e| {
            eprintln!("Error reading gilist: {}", e);
            std::process::exit(1);
        });
        Some(content.lines().filter_map(|l| l.trim().parse().ok()).collect())
    } else {
        None
    };

    let neg_gi_filter: Option<HashSet<i64>> = if let Some(ref gi_file) = args.negative_gilist {
        let content = fs::read_to_string(gi_file).unwrap_or_else(|e| {
            eprintln!("Error reading negative_gilist: {}", e);
            std::process::exit(1);
        });
        Some(content.lines().filter_map(|l| l.trim().parse().ok()).collect())
    } else {
        None
    };

    // Collect taxid filter
    let tax_filter: Option<HashSet<u32>> = if let Some(ref taxids_str) = args.taxids {
        Some(taxids_str.split(',').filter_map(|s| s.trim().parse().ok()).collect())
    } else if let Some(ref taxid_file) = args.taxidlist {
        let content = fs::read_to_string(taxid_file).unwrap_or_else(|e| {
            eprintln!("Error reading taxidlist: {}", e);
            std::process::exit(1);
        });
        Some(content.lines().filter_map(|l| l.trim().parse().ok()).collect())
    } else {
        None
    };

    // Collect seqid filter
    let seqid_filter: Option<HashSet<String>> = if let Some(ref seqid_file) = args.seqidlist {
        let content = fs::read_to_string(seqid_file).unwrap_or_else(|e| {
            eprintln!("Error reading seqidlist: {}", e);
            std::process::exit(1);
        });
        Some(content.lines().filter(|l| !l.trim().is_empty()).map(|l| l.trim().to_string()).collect())
    } else {
        None
    };

    // Collect negative seqid filter
    let neg_seqid_filter: Option<HashSet<String>> = if let Some(ref neg_file) = args.negative_seqidlist {
        let content = fs::read_to_string(neg_file).unwrap_or_else(|e| {
            eprintln!("Error reading negative_seqidlist: {}", e);
            std::process::exit(1);
        });
        Some(content.lines().filter(|l| !l.trim().is_empty()).map(|l| l.trim().to_string()).collect())
    } else {
        None
    };

    if tax_filter.is_none() && seqid_filter.is_none() && neg_seqid_filter.is_none()
        && gi_filter.is_none() && neg_gi_filter.is_none()
    {
        return None;
    }

    let mut allowed_oids = HashSet::new();
    let total = db.num_sequences();

    for oid in 0..total {
        let mut include = true;

        // Tax filter
        if let Some(ref taxids) = tax_filter {
            if let Some(Ok(oid_taxids)) = db.get_taxids(oid) {
                if !oid_taxids.iter().any(|t| *t >= 0 && taxids.contains(&(*t as u32))) {
                    include = false;
                }
            } else if let Ok(header) = db.get_header(oid) {
                if !taxids.contains(&header.taxid) {
                    include = false;
                }
            } else {
                include = false;
            }
        }

        // Seqid filter (positive)
        if include {
            if let Some(ref seqids) = seqid_filter {
                if let Ok(header) = db.get_header(oid) {
                    if !seqids.contains(&header.accession) {
                        include = false;
                    }
                } else {
                    include = false;
                }
            }
        }

        // Negative seqid filter
        if include {
            if let Some(ref neg_seqids) = neg_seqid_filter {
                if let Ok(header) = db.get_header(oid) {
                    if neg_seqids.contains(&header.accession) {
                        include = false;
                    }
                }
            }
        }

        // GI filter (positive)
        if include {
            if let Some(ref gis) = gi_filter {
                if let Ok(header) = db.get_header(oid) {
                    let gi: Option<i64> = header.accession.parse().ok();
                    if let Some(gi) = gi {
                        if !gis.contains(&gi) { include = false; }
                    } else {
                        include = false; // Can't determine GI
                    }
                }
            }
        }

        // Negative GI filter
        if include {
            if let Some(ref neg_gis) = neg_gi_filter {
                if let Ok(header) = db.get_header(oid) {
                    let gi: Option<i64> = header.accession.parse().ok();
                    if let Some(gi) = gi {
                        if neg_gis.contains(&gi) { include = false; }
                    }
                }
            }
        }

        if include {
            allowed_oids.insert(oid);
        }
    }

    Some(allowed_oids)
}

fn run_blastp(args: &BlastArgs) {
    let db = BlastDb::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Error opening database: {}", e);
        std::process::exit(1);
    });

    if args.num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.num_threads)
            .build_global()
            .ok();
    }

    let matrix: MatrixType = args.matrix.parse().unwrap_or_else(|e: String| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });

    let fmt = output::OutputFormat::parse(&args.outfmt).unwrap_or_else(|e| {
        eprintln!("Error parsing --outfmt: {}", e);
        std::process::exit(1);
    });

    let mut params = SearchParams::blastp_defaults();
    if let Some(ref task) = args.task {
        apply_task_preset(&mut params, task, "blastp");
    }
    params.matrix = matrix;
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    params.filter_low_complexity = !args.no_lc_filter;
    params.comp_adjust = if args.no_comp_adjust { 0 } else { args.comp_based_stats };
    if let Some(go) = args.gap_open { params.gap_open = go; }
    if let Some(ge) = args.gap_extend { params.gap_extend = ge; }
    if let Some(ws) = args.word_size { params.word_size = ws; }
    params.strand = args.strand.clone();
    params.query_gencode = args.query_gencode;
    params.db_gencode = args.db_gencode;
    params.max_hsps = args.max_hsps;
    params.culling_limit = args.culling_limit;
    if let Some(xu) = args.xdrop_ungap { params.x_drop_ungapped = xu; }
    if let Some(xg) = args.xdrop_gap { params.x_drop_gapped = xg; }
    if let Some(xf) = args.xdrop_gap_final { params.x_drop_final = xf; }
    if let Some(ws) = args.window_size { params.two_hit_window = ws; }
    // params.soft_masking = args.soft_masking;  // TODO: add to SearchParams
    // params.lcase_masking = args.lcase_masking;  // TODO: add to SearchParams

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let out: Box<dyn Write> = match &args.out {
        Some(p) => Box::new(fs::File::create(p).unwrap()),
        None => Box::new(io::stdout()),
    };
    let mut out = io::BufWriter::new(out);

    let db_path = args.db.to_string_lossy();
    let db_title = db.title().to_string();
    let db_num_seqs = db.num_sequences() as u64;
    let db_len = db.volume_length();
    let matrix_name = args.matrix.as_str();

    if fmt.fmt_id == 5 {
        output::write_xml_header(&mut out, "blastp", &db_path, "blast-cli 0.1.0").unwrap();
    } else if fmt.fmt_id == 13 || fmt.fmt_id == 15 {
        output::write_json_header(&mut out).unwrap();
    }

    if args.html {
        writeln!(out, "<!DOCTYPE html><html><head><title>BLAST Results</title>\
            <style>body {{ font-family: monospace; white-space: pre; }}</style>\
            </head><body>").unwrap();
    }

    let filter = build_oid_filter(&db, args);

    for (iter_num, (query_title, query_seq)) in queries.iter().enumerate() {
        let mut results = api_blastp(&db, query_seq, &params);
        if let Some(ref f) = filter {
            results.retain(|r| f.contains(&r.subject_oid));
        }

        // JSON comma separator between multi-query entries
        if (fmt.fmt_id == 13 || fmt.fmt_id == 15) && iter_num > 0 {
            writeln!(out, ",").unwrap();
        }

        let ctx = output::SearchContext {
            program: "blastp",
            db_path: &db_path,
            db_title: &db_title,
            db_num_seqs,
            db_len,
            query_title,
            query_seq,
            query_len: query_seq.len(),
            matrix: matrix_name,
            gap_open: params.gap_open,
            gap_extend: params.gap_extend,
            evalue_threshold: args.evalue,
            iter_num: iter_num + 1,
            num_descriptions: None,
            num_alignments: None, taxdb: None,
        };

        output::write_results(&mut out, &fmt, &ctx, &results, Some(&db)).unwrap();
    }

    if args.html {
        writeln!(out, "</body></html>").unwrap();
    }

    if fmt.fmt_id == 5 {
        output::write_xml_footer(&mut out).unwrap();
    } else if fmt.fmt_id == 13 || fmt.fmt_id == 15 {
        output::write_json_footer(&mut out).unwrap();
    }
}

fn run_blastn(args: &BlastArgs) {
    let db = BlastDb::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Error opening database: {}", e);
        std::process::exit(1);
    });

    if args.num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.num_threads)
            .build_global()
            .ok();
    }

    let fmt = output::OutputFormat::parse(&args.outfmt).unwrap_or_else(|e| {
        eprintln!("Error parsing --outfmt: {}", e);
        std::process::exit(1);
    });

    let mut params = if let Some(ref strat_path) = args.import_strategy {
        let file = fs::File::open(strat_path).unwrap_or_else(|e| {
            eprintln!("Error reading strategy: {}", e);
            std::process::exit(1);
        });
        let mut reader = io::BufReader::new(file);
        SearchParams::load_strategy(&mut reader).unwrap_or_else(|e| {
            eprintln!("Error parsing strategy: {}", e);
            std::process::exit(1);
        })
    } else {
        SearchParams::blastn_defaults()
    };
    if let Some(ref task) = args.task {
        apply_task_preset(&mut params, task, "blastn");
    }
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    params.match_score = args.match_score;
    params.mismatch = args.mismatch;
    params.filter_low_complexity = !args.no_lc_filter;
    params.comp_adjust = if args.no_comp_adjust { 0 } else { args.comp_based_stats };
    if let Some(go) = args.gap_open { params.gap_open = go; }
    if let Some(ge) = args.gap_extend { params.gap_extend = ge; }
    if let Some(ws) = args.word_size { params.word_size = ws; }
    params.strand = args.strand.clone();
    params.query_gencode = args.query_gencode;
    params.db_gencode = args.db_gencode;
    params.max_hsps = args.max_hsps;
    params.culling_limit = args.culling_limit;
    if let Some(xu) = args.xdrop_ungap { params.x_drop_ungapped = xu; }
    if let Some(xg) = args.xdrop_gap { params.x_drop_gapped = xg; }
    if let Some(xf) = args.xdrop_gap_final { params.x_drop_final = xf; }
    if let Some(ws) = args.window_size { params.two_hit_window = ws; }
    // params.soft_masking = args.soft_masking;  // TODO: add to SearchParams
    // params.lcase_masking = args.lcase_masking;  // TODO: add to SearchParams

    let use_dc = args.task.as_deref() == Some("dc-megablast");

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let out: Box<dyn Write> = match &args.out {
        Some(p) => Box::new(fs::File::create(p).unwrap()),
        None => Box::new(io::stdout()),
    };
    let mut out = io::BufWriter::new(out);

    let db_path = args.db.to_string_lossy();
    let db_title = db.title().to_string();
    let db_num_seqs = db.num_sequences() as u64;
    let db_len = db.volume_length();

    if fmt.fmt_id == 5 {
        output::write_xml_header(&mut out, "blastn", &db_path, "blast-cli 0.1.0").unwrap();
    } else if fmt.fmt_id == 13 || fmt.fmt_id == 15 {
        output::write_json_header(&mut out).unwrap();
    }

    if args.html {
        writeln!(out, "<!DOCTYPE html><html><head><title>BLAST Results</title>\
            <style>body {{ font-family: monospace; white-space: pre; }}</style>\
            </head><body>").unwrap();
    }

    let filter = build_oid_filter(&db, args);

    for (iter_num, (query_title, query_seq)) in queries.iter().enumerate() {
        if (fmt.fmt_id == 13 || fmt.fmt_id == 15) && iter_num > 0 { writeln!(out, ",").unwrap(); }
        let mut results = if use_dc {
            let tt = args.template_type.unwrap_or(0);
            let tl = args.template_length.unwrap_or(21);
            blast_rs::search::dc_megablast_search(&db, query_seq, &params, tt, tl)
        } else {
            api_blastn(&db, query_seq, &params)
        };
        if let Some(ref f) = filter {
            results.retain(|r| f.contains(&r.subject_oid));
        }

        let ctx = output::SearchContext {
            program: "blastn",
            db_path: &db_path,
            db_title: &db_title,
            db_num_seqs,
            db_len,
            query_title,
            query_seq,
            query_len: query_seq.len(),
            matrix: "N/A",
            gap_open: params.gap_open,
            gap_extend: params.gap_extend,
            evalue_threshold: args.evalue,
            iter_num: iter_num + 1,
            num_descriptions: None,
            num_alignments: None, taxdb: None,
        };

        output::write_results(&mut out, &fmt, &ctx, &results, Some(&db)).unwrap();
    }

    if args.html {
        writeln!(out, "</body></html>").unwrap();
    }

    if fmt.fmt_id == 5 {
        output::write_xml_footer(&mut out).unwrap();
    } else if fmt.fmt_id == 13 || fmt.fmt_id == 15 {
        output::write_json_footer(&mut out).unwrap();
    }

    if let Some(ref strat_path) = args.export_strategy {
        let mut f = fs::File::create(strat_path).unwrap();
        params.save_strategy(&mut f).unwrap();
        eprintln!("Search strategy saved to {}", strat_path.display());
    }
}

// ── Shared helpers ────────────────────────────────────────────────────────────

fn make_output(args_out: &Option<PathBuf>) -> io::BufWriter<Box<dyn Write>> {
    let out: Box<dyn Write> = match args_out {
        Some(p) => Box::new(fs::File::create(p).unwrap()),
        None    => Box::new(io::stdout()),
    };
    io::BufWriter::new(out)
}

fn parse_fmt(outfmt: &str) -> output::OutputFormat {
    output::OutputFormat::parse(outfmt).unwrap_or_else(|e| {
        eprintln!("Error parsing --outfmt: {}", e);
        std::process::exit(1);
    })
}

fn parse_matrix(matrix: &str) -> MatrixType {
    matrix.parse().unwrap_or_else(|e: String| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    })
}

fn open_db(db: &Path) -> BlastDb {
    BlastDb::open(db).unwrap_or_else(|e| {
        eprintln!("Error opening database: {}", e);
        std::process::exit(1);
    })
}

fn set_threads(n: usize) {
    if n > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }
}

// ── BLASTX ───────────────────────────────────────────────────────────────────

fn run_blastx(args: &BlastArgs) {
    let db = open_db(&args.db);
    set_threads(args.num_threads);

    let fmt = parse_fmt(&args.outfmt);
    let mut params = SearchParams::blastx_defaults();
    if let Some(ref task) = args.task {
        apply_task_preset(&mut params, task, "blastx");
    }
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    params.filter_low_complexity = !args.no_lc_filter;
    params.comp_adjust = if args.no_comp_adjust { 0 } else { args.comp_based_stats };
    if let Some(go) = args.gap_open  { params.gap_open = go; }
    if let Some(ge) = args.gap_extend { params.gap_extend = ge; }
    if let Some(ws) = args.word_size  { params.word_size = ws; }
    params.strand = args.strand.clone();
    params.query_gencode = args.query_gencode;
    params.db_gencode = args.db_gencode;
    params.max_hsps = args.max_hsps;
    params.culling_limit = args.culling_limit;
    if let Some(xu) = args.xdrop_ungap { params.x_drop_ungapped = xu; }
    if let Some(xg) = args.xdrop_gap { params.x_drop_gapped = xg; }
    if let Some(xf) = args.xdrop_gap_final { params.x_drop_final = xf; }
    if let Some(ws) = args.window_size { params.two_hit_window = ws; }
    // params.soft_masking = args.soft_masking;  // TODO: add to SearchParams
    // params.lcase_masking = args.lcase_masking;  // TODO: add to SearchParams

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let mut out = make_output(&args.out);
    let db_path  = args.db.to_string_lossy();
    let db_title = db.title().to_string();
    let db_num_seqs = db.num_sequences() as u64;
    let db_len   = db.volume_length();

    if fmt.fmt_id == 5 { output::write_xml_header(&mut out, "blastx", &db_path, "blast-cli 0.1.0").unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_header(&mut out).unwrap(); }

    if args.html {
        writeln!(out, "<!DOCTYPE html><html><head><title>BLAST Results</title>\
            <style>body {{ font-family: monospace; white-space: pre; }}</style>\
            </head><body>").unwrap();
    }

    let filter = build_oid_filter(&db, args);

    for (iter_num, (query_title, query_seq)) in queries.iter().enumerate() {
        if (fmt.fmt_id == 13 || fmt.fmt_id == 15) && iter_num > 0 { writeln!(out, ",").unwrap(); }
        let mut results = api_blastx(&db, query_seq, &params);
        if let Some(ref f) = filter {
            results.retain(|r| f.contains(&r.subject_oid));
        }
        let ctx = output::SearchContext {
            program: "blastx", db_path: &db_path, db_title: &db_title,
            db_num_seqs, db_len, query_title,
            query_seq, query_len: query_seq.len(),
            matrix: "N/A", gap_open: params.gap_open, gap_extend: params.gap_extend,
            evalue_threshold: args.evalue, iter_num: iter_num + 1,
            num_descriptions: None, num_alignments: None, taxdb: None,
        };
        output::write_results(&mut out, &fmt, &ctx, &results, Some(&db)).unwrap();
    }

    if args.html { writeln!(out, "</body></html>").unwrap(); }

    if fmt.fmt_id == 5 { output::write_xml_footer(&mut out).unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_footer(&mut out).unwrap(); }
}

// ── TBLASTN ──────────────────────────────────────────────────────────────────

fn run_tblastn(args: &BlastArgs) {
    let db = open_db(&args.db);
    set_threads(args.num_threads);

    let matrix = parse_matrix(&args.matrix);
    let fmt = parse_fmt(&args.outfmt);
    let mut params = SearchParams::tblastn_defaults();
    if let Some(ref task) = args.task {
        apply_task_preset(&mut params, task, "tblastn");
    }
    params.matrix = matrix;
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    params.filter_low_complexity = !args.no_lc_filter;
    params.comp_adjust = if args.no_comp_adjust { 0 } else { args.comp_based_stats };
    if let Some(go) = args.gap_open   { params.gap_open = go; }
    if let Some(ge) = args.gap_extend { params.gap_extend = ge; }
    if let Some(ws) = args.word_size  { params.word_size = ws; }
    params.strand = args.strand.clone();
    params.query_gencode = args.query_gencode;
    params.db_gencode = args.db_gencode;
    params.max_hsps = args.max_hsps;
    params.culling_limit = args.culling_limit;
    if let Some(xu) = args.xdrop_ungap { params.x_drop_ungapped = xu; }
    if let Some(xg) = args.xdrop_gap { params.x_drop_gapped = xg; }
    if let Some(xf) = args.xdrop_gap_final { params.x_drop_final = xf; }
    if let Some(ws) = args.window_size { params.two_hit_window = ws; }
    // params.soft_masking = args.soft_masking;  // TODO: add to SearchParams
    // params.lcase_masking = args.lcase_masking;  // TODO: add to SearchParams

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let mut out = make_output(&args.out);
    let db_path     = args.db.to_string_lossy();
    let db_title    = db.title().to_string();
    let db_num_seqs = db.num_sequences() as u64;
    let db_len      = db.volume_length();
    let matrix_name = args.matrix.as_str();

    if fmt.fmt_id == 5 { output::write_xml_header(&mut out, "tblastn", &db_path, "blast-cli 0.1.0").unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_header(&mut out).unwrap(); }

    if args.html {
        writeln!(out, "<!DOCTYPE html><html><head><title>BLAST Results</title>\
            <style>body {{ font-family: monospace; white-space: pre; }}</style>\
            </head><body>").unwrap();
    }

    let filter = build_oid_filter(&db, args);

    for (iter_num, (query_title, query_seq)) in queries.iter().enumerate() {
        if (fmt.fmt_id == 13 || fmt.fmt_id == 15) && iter_num > 0 { writeln!(out, ",").unwrap(); }
        let mut results = api_tblastn(&db, query_seq, &params);
        if let Some(ref f) = filter {
            results.retain(|r| f.contains(&r.subject_oid));
        }
        let ctx = output::SearchContext {
            program: "tblastn", db_path: &db_path, db_title: &db_title,
            db_num_seqs, db_len, query_title,
            query_seq, query_len: query_seq.len(),
            matrix: matrix_name, gap_open: params.gap_open, gap_extend: params.gap_extend,
            evalue_threshold: args.evalue, iter_num: iter_num + 1,
            num_descriptions: None, num_alignments: None, taxdb: None,
        };
        output::write_results(&mut out, &fmt, &ctx, &results, Some(&db)).unwrap();
    }

    if args.html { writeln!(out, "</body></html>").unwrap(); }

    if fmt.fmt_id == 5 { output::write_xml_footer(&mut out).unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_footer(&mut out).unwrap(); }
}

// ── TBLASTX ──────────────────────────────────────────────────────────────────

fn run_tblastx(args: &BlastArgs) {
    let db = open_db(&args.db);
    set_threads(args.num_threads);

    let fmt = parse_fmt(&args.outfmt);
    let mut params = SearchParams::tblastx_defaults();
    if let Some(ref task) = args.task {
        apply_task_preset(&mut params, task, "tblastx");
    }
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    params.filter_low_complexity = !args.no_lc_filter;
    params.comp_adjust = if args.no_comp_adjust { 0 } else { args.comp_based_stats };
    if let Some(ws) = args.word_size { params.word_size = ws; }
    params.strand = args.strand.clone();
    params.query_gencode = args.query_gencode;
    params.db_gencode = args.db_gencode;
    params.max_hsps = args.max_hsps;
    params.culling_limit = args.culling_limit;
    if let Some(xu) = args.xdrop_ungap { params.x_drop_ungapped = xu; }
    if let Some(xg) = args.xdrop_gap { params.x_drop_gapped = xg; }
    if let Some(xf) = args.xdrop_gap_final { params.x_drop_final = xf; }
    if let Some(ws) = args.window_size { params.two_hit_window = ws; }
    // params.soft_masking = args.soft_masking;  // TODO: add to SearchParams
    // params.lcase_masking = args.lcase_masking;  // TODO: add to SearchParams

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let mut out = make_output(&args.out);
    let db_path     = args.db.to_string_lossy();
    let db_title    = db.title().to_string();
    let db_num_seqs = db.num_sequences() as u64;
    let db_len      = db.volume_length();

    if fmt.fmt_id == 5 { output::write_xml_header(&mut out, "tblastx", &db_path, "blast-cli 0.1.0").unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_header(&mut out).unwrap(); }

    if args.html {
        writeln!(out, "<!DOCTYPE html><html><head><title>BLAST Results</title>\
            <style>body {{ font-family: monospace; white-space: pre; }}</style>\
            </head><body>").unwrap();
    }

    let filter = build_oid_filter(&db, args);

    for (iter_num, (query_title, query_seq)) in queries.iter().enumerate() {
        if (fmt.fmt_id == 13 || fmt.fmt_id == 15) && iter_num > 0 { writeln!(out, ",").unwrap(); }
        let mut results = api_tblastx(&db, query_seq, &params);
        if let Some(ref f) = filter {
            results.retain(|r| f.contains(&r.subject_oid));
        }
        let ctx = output::SearchContext {
            program: "tblastx", db_path: &db_path, db_title: &db_title,
            db_num_seqs, db_len, query_title,
            query_seq, query_len: query_seq.len(),
            matrix: "N/A", gap_open: params.gap_open, gap_extend: params.gap_extend,
            evalue_threshold: args.evalue, iter_num: iter_num + 1,
            num_descriptions: None, num_alignments: None, taxdb: None,
        };
        output::write_results(&mut out, &fmt, &ctx, &results, Some(&db)).unwrap();
    }

    if args.html { writeln!(out, "</body></html>").unwrap(); }

    if fmt.fmt_id == 5 { output::write_xml_footer(&mut out).unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_footer(&mut out).unwrap(); }
}

// ── PSI-BLAST ────────────────────────────────────────────────────────────────

fn run_psiblast(args: &PsiblastArgs) {
    let db = open_db(&args.db);
    set_threads(args.num_threads);

    let matrix = parse_matrix(&args.matrix);
    let fmt = parse_fmt(&args.outfmt);
    let mut params = SearchParams::blastp_defaults();
    params.matrix = matrix;
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    if let Some(go) = args.gap_open   { params.gap_open = go; }
    if let Some(ge) = args.gap_extend { params.gap_extend = ge; }

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let mut out = make_output(&args.out);
    let db_path     = args.db.to_string_lossy();
    let db_title    = db.title().to_string();
    let db_num_seqs = db.num_sequences() as u64;
    let db_len      = db.volume_length();
    let matrix_name = args.matrix.as_str();

    if fmt.fmt_id == 5 { output::write_xml_header(&mut out, "psiblast", &db_path, "blast-cli 0.1.0").unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_header(&mut out).unwrap(); }

    let psi_params = PsiblastParams::new(params.clone())
        .num_iterations(args.num_iterations)
        .inclusion_evalue(args.inclusion_ethresh);

    for (iter_num, (query_title, query_seq)) in queries.iter().enumerate() {
        if (fmt.fmt_id == 13 || fmt.fmt_id == 15) && iter_num > 0 { writeln!(out, ",").unwrap(); }
        let (results, pssm) = api_psiblast(&db, query_seq, &psi_params);
        let ctx = output::SearchContext {
            program: "psiblast", db_path: &db_path, db_title: &db_title,
            db_num_seqs, db_len, query_title,
            query_seq, query_len: query_seq.len(),
            matrix: matrix_name, gap_open: params.gap_open, gap_extend: params.gap_extend,
            evalue_threshold: args.evalue, iter_num: iter_num + 1,
            num_descriptions: None, num_alignments: None, taxdb: None,
        };
        output::write_results(&mut out, &fmt, &ctx, &results, Some(&db)).unwrap();

        // Save PSSM checkpoint if requested
        if let Some(ref pssm_path) = args.out_pssm {
            let mut f = fs::File::create(pssm_path).unwrap_or_else(|e| {
                eprintln!("Error creating PSSM file: {}", e);
                std::process::exit(1);
            });
            pssm.write_checkpoint(&mut f).unwrap();
            eprintln!("PSSM checkpoint saved to {}", pssm_path.display());
        }

        // Save ASCII PSSM if requested
        if let Some(ref ascii_path) = args.out_ascii_pssm {
            let ncbi_query = blast_rs::search::aa_to_ncbistdaa(query_seq);
            let mut f = fs::File::create(ascii_path).unwrap_or_else(|e| {
                eprintln!("Error creating ASCII PSSM file: {}", e);
                std::process::exit(1);
            });
            pssm.write_ascii(&mut f, &ncbi_query).unwrap();
            eprintln!("ASCII PSSM saved to {}", ascii_path.display());
        }
    }

    if fmt.fmt_id == 5 { output::write_xml_footer(&mut out).unwrap(); }
    else if fmt.fmt_id == 13 || fmt.fmt_id == 15 { output::write_json_footer(&mut out).unwrap(); }
}

fn run_aliastool(args: &AliasArgs) {
    let ext = match args.dbtype.to_ascii_lowercase().as_str() {
        "prot" | "protein" => "pal",
        "nucl" | "nucleotide" => "nal",
        other => {
            eprintln!("Unknown dbtype '{}'. Use 'prot' or 'nucl'.", other);
            std::process::exit(1);
        }
    };

    let title = args.title.clone().unwrap_or_else(|| {
        args.out.to_string_lossy().into_owned()
    });

    let dblist = args.dblist.clone();

    let alias_path = args.out.with_extension(ext);
    let mut f = fs::File::create(&alias_path).unwrap_or_else(|e| {
        eprintln!("Error creating alias file: {}", e);
        std::process::exit(1);
    });

    writeln!(f, "#").unwrap();
    writeln!(f, "# Alias file created by blast-cli blastdb-aliastool").unwrap();
    writeln!(f, "#").unwrap();
    writeln!(f, "TITLE {}", title).unwrap();
    writeln!(f, "DBLIST {}", dblist).unwrap();

    eprintln!("Alias file written: {}", alias_path.display());
}

fn run_makeblastdb(args: &MakeDbArgs) {
    let seq_type = match args.dbtype.to_ascii_lowercase().as_str() {
        "prot" | "protein" => SeqType::Protein,
        "nucl" | "nucleotide" => SeqType::Nucleotide,
        other => {
            eprintln!("Unknown dbtype '{}'. Use 'prot' or 'nucl'.", other);
            std::process::exit(1);
        }
    };

    let db_title = args.title.clone().unwrap_or_else(|| {
        args.out.to_string_lossy().into_owned()
    });

    let sequences = read_fasta(&args.input).unwrap_or_else(|e| {
        eprintln!("Error reading input FASTA: {}", e);
        std::process::exit(1);
    });

    if sequences.is_empty() {
        eprintln!("No sequences found in input file.");
        std::process::exit(1);
    }

    let mut builder = BlastDbBuilder::new(seq_type, &db_title);
    for (header, seq) in sequences {
        let (accession, title) = if args.parse_seqids {
            let mut words = header.splitn(2, ' ');
            let acc = words.next().unwrap_or(&header).to_string();
            let ttl = words.next().unwrap_or("").to_string();
            (acc, ttl)
        } else {
            (header.clone(), header.clone())
        };
        builder.add(SequenceEntry { title, accession, sequence: seq, taxid: None });
    }

    if args.max_file_size > 0 {
        builder.write_multivolume(&args.out, args.dbversion as i32, args.max_file_size).unwrap_or_else(|e| {
            eprintln!("Error writing database: {}", e);
            std::process::exit(1);
        });
    } else if args.dbversion == 5 {
        builder.write_v5(&args.out).unwrap_or_else(|e| {
            eprintln!("Error writing database: {}", e);
            std::process::exit(1);
        });
    } else {
        builder.write(&args.out).unwrap_or_else(|e| {
            eprintln!("Error writing database: {}", e);
            std::process::exit(1);
        });
    }

    let n = builder.entries.len();
    let ext = match seq_type { SeqType::Protein => "pin", SeqType::Nucleotide => "nin" };
    eprintln!(
        "Database written (v{}): {} sequences → {}.{}",
        args.dbversion,
        n,
        args.out.display(),
        ext
    );
}

fn run_dumpdb(args: &DumpArgs) {
    let db = BlastDb::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Error opening database: {}", e);
        std::process::exit(1);
    });

    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    // v5: accession lookup
    if let Some(acc) = &args.lookup {
        match db.lookup_accession(acc) {
            None => eprintln!("Database is v4 (no LMDB accession index)."),
            Some(Err(e)) => eprintln!("Lookup error: {}", e),
            Some(Ok(oids)) => {
                if oids.is_empty() {
                    writeln!(out, "Accession '{}' not found.", acc).unwrap();
                } else {
                    for oid in oids {
                        writeln!(out, "{}\t{}", acc, oid).unwrap();
                    }
                }
            }
        }
        return;
    }

    // v5: list volume info
    if args.volumes {
        match db.get_volumes_info() {
            None => eprintln!("Database is v4 (no LMDB volume info)."),
            Some(Err(e)) => eprintln!("Volume info error: {}", e),
            Some(Ok(vols)) => {
                writeln!(out, "DB version: {}", db.format_version()).unwrap();
                writeln!(out, "Volumes ({}):", vols.len()).unwrap();
                for (name, n_oids) in &vols {
                    writeln!(out, "  {} ({} OIDs)", name, n_oids).unwrap();
                }
            }
        }
        return;
    }

    // v5: list all accessions
    if args.list_accessions {
        match db.iter_accessions(|acc, oid| {
            writeln!(out, "{}\t{}", acc, oid).ok();
        }) {
            None => eprintln!("Database is v4 (no LMDB accession index)."),
            Some(Err(e)) => eprintln!("Iteration error: {}", e),
            Some(Ok(())) => {}
        }
        return;
    }

    // Standard sequence dump
    let n = if args.max == 0 { db.num_sequences() } else { args.max as u32 };

    for oid in 0..n {
        let header = db.get_header(oid).unwrap_or_default();

        // For v5, prefer the .pos/.nos seqid list for the defline.
        let accession = if let Some(Ok(ids)) = db.get_seqids(oid) {
            ids.into_iter().next().unwrap_or_else(|| header.accession.clone())
        } else {
            header.accession.clone()
        };

        let title = if header.title.is_empty() {
            format!("OID={}", oid)
        } else {
            header.title.clone()
        };

        if args.headers_only {
            writeln!(out, ">{} {} [taxid={}]", accession, title, header.taxid).unwrap();
            continue;
        }

        match db.seq_type() {
            blast_rs::db::index::SeqType::Protein => {
                let seq = db.get_sequence_protein(oid).unwrap_or_default();
                writeln!(out, ">{} {}", accession, title).unwrap();
                for chunk in seq.chunks(60) {
                    writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or("?")).unwrap();
                }
            }
            blast_rs::db::index::SeqType::Nucleotide => {
                let seq = db.get_sequence_nucleotide(oid).unwrap_or_default();
                writeln!(out, ">{} {}", accession, title).unwrap();
                for chunk in seq.chunks(60) {
                    writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or("?")).unwrap();
                }
            }
        }
    }
}

// ---- FASTA reader ----

fn run_blastdbcmd(args: &BlastdbcmdArgs) {
    let db = open_db(&args.db);
    let stdout = io::stdout();
    let out: Box<dyn Write> = match &args.out {
        Some(p) => Box::new(fs::File::create(p).unwrap()),
        None => Box::new(stdout.lock()),
    };
    let mut out = io::BufWriter::new(out);

    if args.info {
        writeln!(out, "Database: {}", args.db.display()).unwrap();
        writeln!(out, "  {} sequences; {} total letters", db.num_sequences(), db.volume_length()).unwrap();
        writeln!(out, "  Title: {}", db.title()).unwrap();
        writeln!(out, "  Type: {}", if db.seq_type() == SeqType::Protein { "Protein" } else { "Nucleotide" }).unwrap();
        writeln!(out, "  Format version: {}", db.format_version()).unwrap();
        writeln!(out, "  Volumes: {}", db.num_volumes()).unwrap();
        return;
    }

    let accessions: Vec<String> = if let Some(ref entry) = args.entry {
        entry.split(',').map(|s| s.trim().to_string()).collect()
    } else if let Some(ref batch_path) = args.entry_batch {
        let content = fs::read_to_string(batch_path).unwrap_or_else(|e| {
            eprintln!("Error reading entry batch file: {}", e);
            std::process::exit(1);
        });
        content.lines().filter(|l| !l.trim().is_empty()).map(|l| l.trim().to_string()).collect()
    } else {
        eprintln!("Error: specify --entry or --entry_batch or --info");
        std::process::exit(1);
    };

    let range = args.range.as_ref().and_then(|r| {
        let parts: Vec<&str> = r.split('-').collect();
        if parts.len() == 2 {
            let start: usize = parts[0].parse().ok()?;
            let end: usize = parts[1].parse().ok()?;
            Some((start.saturating_sub(1), end)) // convert 1-based to 0-based
        } else {
            None
        }
    });

    for acc in &accessions {
        // Try accession lookup (v5) first, then scan headers
        let oid = if let Some(result) = db.lookup_accession(acc) {
            match result {
                Ok(oids) if !oids.is_empty() => Some(oids[0]),
                _ => None,
            }
        } else {
            // Linear scan for v4 databases
            let mut found = None;
            for oid in 0..db.num_sequences() {
                if let Ok(header) = db.get_header(oid) {
                    if header.accession == *acc {
                        found = Some(oid);
                        break;
                    }
                }
            }
            found
        };

        let oid = match oid {
            Some(o) => o,
            None => {
                eprintln!("Accession '{}' not found.", acc);
                continue;
            }
        };

        if args.outfmt == "accession" {
            let header = db.get_header(oid).unwrap_or_default();
            writeln!(out, "{}", header.accession).unwrap();
            continue;
        }

        // FASTA output (default)
        let header = db.get_header(oid).unwrap_or_default();
        let title = if header.title.is_empty() { acc.clone() } else { header.title };

        match db.seq_type() {
            SeqType::Protein => {
                let seq = db.get_sequence_protein(oid).unwrap_or_default();
                let seq = if let Some((start, end)) = range {
                    seq[start..end.min(seq.len())].to_vec()
                } else { seq };
                writeln!(out, ">{} {}", header.accession, title).unwrap();
                for chunk in seq.chunks(60) {
                    out.write_all(chunk).unwrap();
                    writeln!(out).unwrap();
                }
            }
            SeqType::Nucleotide => {
                let seq = db.get_sequence_nucleotide(oid).unwrap_or_default();
                let seq = if let Some((start, end)) = range {
                    seq[start..end.min(seq.len())].to_vec()
                } else { seq };
                writeln!(out, ">{} {}", header.accession, title).unwrap();
                for chunk in seq.chunks(60) {
                    out.write_all(chunk).unwrap();
                    writeln!(out).unwrap();
                }
            }
        }
    }
}

fn read_fasta(path: &Path) -> io::Result<Vec<(String, Vec<u8>)>> {
    let file = fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_title = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() { continue; }
        if let Some(stripped) = line.strip_prefix('>') {
            if !current_title.is_empty() {
                sequences.push((current_title.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_title = stripped.to_string();
        } else {
            current_seq.extend_from_slice(line.as_bytes());
        }
    }
    if !current_title.is_empty() {
        sequences.push((current_title, current_seq));
    }
    Ok(sequences)
}

// ── BLAST Formatter ─────────────────────────────────────────────────────────

fn run_blast_formatter(args: &FormatterArgs) {
    let archive = std::fs::read_to_string(&args.archive).unwrap_or_else(|e| {
        eprintln!("Error reading archive: {}", e);
        std::process::exit(1);
    });

    let fmt = output::OutputFormat::parse(&args.outfmt).unwrap_or_else(|e| {
        eprintln!("Error parsing --outfmt: {}", e);
        std::process::exit(1);
    });

    let mut out: Box<dyn Write> = match &args.out {
        Some(p) => Box::new(std::io::BufWriter::new(
            fs::File::create(p).unwrap_or_else(|e| { eprintln!("Error: {}", e); std::process::exit(1); })
        )),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };

    // Parse the text ASN.1 archive to extract search metadata and results
    let (program, db_path, query_title, query_len, results) = parse_archive_asn1(&archive);

    let ctx = output::SearchContext {
        program: &program,
        db_path: &db_path,
        db_title: &db_path,
        db_num_seqs: 0,
        db_len: 0,
        query_title: &query_title,
        query_seq: &[],
        query_len,
        matrix: "BLOSUM62",
        gap_open: 11,
        gap_extend: 1,
        evalue_threshold: 10.0,
        iter_num: 1,
        num_descriptions: None,
        num_alignments: None, taxdb: None,
    };

    output::write_results(&mut out, &fmt, &ctx, &results, None).unwrap();
}

/// Parse a text ASN.1 BLAST archive (format 11) into SearchResults.
/// Splits by "type partial" alignment blocks for robust extraction.
fn parse_archive_asn1(text: &str) -> (String, String, String, usize, Vec<blast_rs::SearchResult>) {
    use blast_rs::{SearchResult, Hsp};

    let mut program = String::from("blastp");
    let mut db_path = String::new();
    let mut query_title = String::new();
    let mut query_len = 0usize;
    let mut results: Vec<SearchResult> = Vec::new();

    // Extract metadata from the request section
    for line in text.lines() {
        let t = line.trim();
        if t.starts_with("program \"") {
            program = extract_quoted(t, "program ");
        } else if t.starts_with("subject database \"") {
            db_path = extract_quoted(t, "subject database ");
        } else if t.starts_with("title \"") && query_title.is_empty() {
            query_title = extract_quoted(t, "title ");
        } else if t.starts_with("length ") && query_len == 0 {
            query_len = t[7..].trim().parse().unwrap_or(0);
        }
    }

    // Split text into alignment blocks (each starts with "type partial")
    let blocks: Vec<&str> = text.split("type partial").skip(1).collect();

    for block in blocks {
        let mut score = 0i32;
        let mut evalue = 0.0f64;
        let mut bitscore = 0.0f64;
        let mut nident = 0usize;
        let mut _qid = String::new();
        let mut sid = String::new();
        let mut starts = Vec::new();
        let mut lens = Vec::new();

        let mut last_id_name = String::new();
        let mut in_starts = false;
        let mut in_lens = false;
        let mut in_ids = false;
        let mut id_count = 0;

        for line in block.lines() {
            let t = line.trim();

            // Track which score field we're reading
            if t.contains("id str \"score\"") { last_id_name = "score".into(); }
            else if t.contains("id str \"e_value\"") { last_id_name = "e_value".into(); }
            else if t.contains("id str \"bit_score\"") { last_id_name = "bit_score".into(); }
            else if t.contains("id str \"num_ident\"") { last_id_name = "num_ident".into(); }

            // Extract integer values
            if let Some(rest) = t.strip_prefix("value int ") {
                let v: i32 = rest.trim().parse().unwrap_or(0);
                match last_id_name.as_str() {
                    "score" => score = v,
                    "num_ident" => nident = v as usize,
                    _ => {}
                }
                last_id_name.clear();
            }

            // Extract real values: value real { mantissa, 10, exponent }
            if t.starts_with("value real {") {
                let inner = t.trim_start_matches("value real {").trim_end_matches('}').trim();
                let parts: Vec<&str> = inner.split(',').map(|s| s.trim()).collect();
                if parts.len() == 3 {
                    let m: f64 = parts[0].parse().unwrap_or(0.0);
                    let e: i32 = parts[2].parse().unwrap_or(0);
                    let val = m * 10.0f64.powi(e);
                    match last_id_name.as_str() {
                        "e_value" => evalue = val,
                        "bit_score" => bitscore = val,
                        _ => {}
                    }
                }
                last_id_name.clear();
            }

            // IDs
            if t == "ids {" { in_ids = true; id_count = 0; }
            if in_ids && t.starts_with("local str \"") {
                let id = extract_quoted(t, "local str ");
                if id_count == 0 { _qid = id; } else { sid = id; }
                id_count += 1;
                if id_count >= 2 { in_ids = false; }
            }

            // Starts array
            if t == "starts {" { in_starts = true; starts.clear(); continue; }
            if in_starts {
                if t.starts_with('}') { in_starts = false; continue; }
                if let Ok(v) = t.trim_end_matches(',').trim().parse::<usize>() {
                    starts.push(v);
                }
                continue;
            }

            // Lens array
            if t.starts_with("lens {") { in_lens = true; lens.clear(); continue; }
            if in_lens {
                if t.starts_with('}') { in_lens = false; continue; }
                if let Ok(v) = t.trim_end_matches(',').trim().parse::<usize>() {
                    lens.push(v);
                }
                continue;
            }
        }

        // Build HSP if we have valid data
        if !sid.is_empty() && !starts.is_empty() && !lens.is_empty() {
            let q_start = starts.first().copied().unwrap_or(0);
            let s_start = starts.get(1).copied().unwrap_or(0);
            let aln_len = lens.first().copied().unwrap_or(0);

            let hsp = Hsp {
                score, bit_score: bitscore, evalue,
                query_start: q_start, query_end: q_start + aln_len,
                subject_start: s_start, subject_end: s_start + aln_len,
                num_identities: nident, num_gaps: 0,
                alignment_length: aln_len,
                query_aln: vec![], midline: vec![], subject_aln: vec![],
                query_frame: 0, subject_frame: 0,
            };

            if let Some(existing) = results.iter_mut().find(|r| r.subject_accession == sid) {
                existing.hsps.push(hsp);
            } else {
                results.push(SearchResult {
                    subject_oid: 0,
                    subject_title: sid.clone(),
                    subject_accession: sid,
                    subject_len: s_start + aln_len,
                    hsps: vec![hsp],
                    taxids: vec![],
                });
            }
        }
    }

    (program, db_path, query_title, query_len, results)
}

fn extract_quoted(line: &str, prefix: &str) -> String {
    line.strip_prefix(prefix)
        .unwrap_or(line)
        .trim()
        .trim_start_matches('"')
        .trim_end_matches(['"', ',', ' '])
        .to_string()
}

