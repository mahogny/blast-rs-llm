//! All BLAST output format writers.
//!
//! Formats implemented:
//!   0  – Pairwise
//!   1  – Query-anchored with identities
//!   2  – Query-anchored no identities (positive marks)
//!   3  – Flat query-anchored with identities
//!   4  – Flat query-anchored no identities
//!   5  – BLAST XML (version 1)
//!   6  – Tabular (12 standard columns, or user-specified)
//!   7  – Tabular with comment lines
//!   8  – Text ASN.1 (stub)
//!   9  – Binary ASN.1 (stub)
//!  10  – Comma-separated values (CSV)
//!  11  – BLAST archive (stub)
//!  12  – Seqalign JSON (stub)
//!  13  – Multiple-file JSON (stub)
//!  14  – Multiple-file XML2 (stub)
//!  15  – Single-file BLAST JSON
//!  16  – Single-file BLAST XML2
//!  17  – Subject sequences (FASTA)
//!  18  – Organism report (stub)

use std::io::{self, Write};
use blast_core::hsp::{Hsp, SearchResult};
use blast_db::BlastDb;

// ─── Public API ────────────────────────────────────────────────────────────

/// Parsed `--outfmt` argument.
#[derive(Debug, Clone)]
pub struct OutputFormat {
    pub fmt_id: u32,
    /// For tabular formats 6/7/10: selected columns (empty = default 12).
    pub columns: Vec<TabularColumn>,
}

impl OutputFormat {
    /// Parse `"6"` or `"6 qseqid sseqid pident"` style strings.
    pub fn parse(s: &str) -> Result<Self, String> {
        let mut tokens = s.split_whitespace();
        let fmt_id: u32 = tokens
            .next()
            .unwrap_or("0")
            .parse()
            .map_err(|_| format!("invalid outfmt number: '{}'", s))?;

        let mut columns: Vec<TabularColumn> = tokens
            .map(|tok| TabularColumn::parse(tok).ok_or_else(|| format!("unknown column: '{}'", tok)))
            .collect::<Result<Vec<_>, _>>()?;

        if columns.is_empty() && matches!(fmt_id, 6 | 7 | 10) {
            columns = TabularColumn::default_columns();
        }

        Ok(OutputFormat { fmt_id, columns })
    }
}

/// Context passed to every formatter.
pub struct SearchContext<'a> {
    pub program: &'a str,
    pub db_path: &'a str,
    pub db_title: &'a str,
    pub db_num_seqs: u64,
    pub db_len: u64,
    pub query_title: &'a str,
    pub query_seq: &'a [u8],
    pub query_len: usize,
    pub matrix: &'a str,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue_threshold: f64,
    pub iter_num: usize, // 1-based, for multi-query XML/JSON
}

/// Write all results for one query in the requested format.
pub fn write_results(
    out: &mut impl Write,
    fmt: &OutputFormat,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    db: Option<&BlastDb>,
) -> io::Result<()> {
    match fmt.fmt_id {
        0 => fmt0_pairwise(out, ctx, results, false, false),
        1 => fmt0_pairwise(out, ctx, results, true, false),
        2 => fmt0_pairwise(out, ctx, results, true, true),
        3 => fmt34_flat(out, ctx, results, false),
        4 => fmt34_flat(out, ctx, results, true),
        5 => fmt5_xml(out, ctx, results),
        6 | 7 | 10 => fmt6_tabular(out, ctx, results, fmt),
        8 | 9 | 11 | 12 | 13 | 14 | 18 => {
            writeln!(out, "# Format {} is not implemented.", fmt.fmt_id)
        }
        15 => fmt15_json(out, ctx, results),
        16 => fmt16_xml2(out, ctx, results),
        17 => fmt17_fasta(out, ctx, results, db),
        _ => writeln!(out, "# Unknown output format {}.", fmt.fmt_id),
    }
}

/// Write the XML document-open boilerplate (call once before first query).
pub fn write_xml_header(out: &mut impl Write, program: &str, db: &str, version: &str) -> io::Result<()> {
    writeln!(out, r#"<?xml version="1.0"?>"#)?;
    writeln!(out, r#"<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">"#)?;
    writeln!(out, "<BlastOutput>")?;
    writeln!(out, "  <BlastOutput_program>{}</BlastOutput_program>", xml_escape(program))?;
    writeln!(out, "  <BlastOutput_version>{}</BlastOutput_version>", xml_escape(version))?;
    writeln!(out, "  <BlastOutput_db>{}</BlastOutput_db>", xml_escape(db))?;
    writeln!(out, "  <BlastOutput_iterations>")
}

pub fn write_xml_footer(out: &mut impl Write) -> io::Result<()> {
    writeln!(out, "  </BlastOutput_iterations>")?;
    writeln!(out, "</BlastOutput>")
}

/// Write the JSON array open `[` (call once).
pub fn write_json_header(out: &mut impl Write) -> io::Result<()> {
    writeln!(out, "[")
}
pub fn write_json_footer(out: &mut impl Write) -> io::Result<()> {
    writeln!(out, "]")
}

// ─── Format 0/1/2: Pairwise ────────────────────────────────────────────────

/// Format 0: classic pairwise.
/// Format 1: midline shows '|' for identical, ' ' for mismatch/conservative (no '+').
/// Format 2: midline shows '+' for identical OR conservative, ' ' for mismatch.
fn fmt0_pairwise(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    query_anchored: bool,    // true → formats 1/2
    positive_marks: bool,    // true → format 2 (use '+' for all conserved, no '|')
) -> io::Result<()> {
    let qid = first_word(ctx.query_title);

    writeln!(out, "BLAST{} {}  Reference: Rust BLAST implementation", if query_anchored {"X"} else {""}, ctx.program.to_uppercase())?;
    writeln!(out)?;
    writeln!(out, "Database: {}", ctx.db_title)?;
    writeln!(out, "           {} sequences; {} total letters", ctx.db_num_seqs, ctx.db_len)?;
    writeln!(out)?;
    writeln!(out, "Query= {}", ctx.query_title)?;
    writeln!(out, "         ({} letters)", ctx.query_len)?;
    writeln!(out)?;

    if results.is_empty() {
        writeln!(out, "                                                                 Score     E")?;
        writeln!(out, "Sequences producing significant alignments:                      (Bits)  Value")?;
        writeln!(out)?;
        writeln!(out, "***** No significant similarity found. *****")?;
        writeln!(out)?;
        return Ok(());
    }

    writeln!(out, "                                                                 Score     E")?;
    writeln!(out, "Sequences producing significant alignments:                      (Bits)  Value")?;
    writeln!(out)?;
    for r in results {
        let title = subject_title(r);
        let best = &r.hsps[0];
        writeln!(out, "{:<67} {:.0}  {:.2e}", &title[..title.len().min(67)], best.bit_score, best.evalue)?;
    }
    writeln!(out)?;

    writeln!(out, "ALIGNMENTS")?;
    for r in results {
        let title = subject_title(r);
        writeln!(out, "> {}", title)?;
        writeln!(out, "   Length = {}", r.subject_len)?;
        writeln!(out)?;

        for (hi, hsp) in r.hsps.iter().enumerate() {
            write_hsp_header(out, hsp, hi + 1)?;
            write_alignment_blocks(out, hsp, query_anchored, positive_marks, 60)?;
        }
    }
    Ok(())
}

// ─── Format 3/4: Flat query-anchored ───────────────────────────────────────

fn fmt34_flat(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    positive_marks: bool,
) -> io::Result<()> {
    let qid = first_word(ctx.query_title);
    writeln!(out, "Query= {} ({} letters)", ctx.query_title, ctx.query_len)?;
    writeln!(out)?;

    if results.is_empty() {
        writeln!(out, "***** No significant similarity found. *****")?;
        return Ok(());
    }

    for r in results {
        let title = subject_title(r);
        writeln!(out, "> {}", title)?;
        writeln!(out, "   Length = {}", r.subject_len)?;
        for (hi, hsp) in r.hsps.iter().enumerate() {
            write_hsp_header(out, hsp, hi + 1)?;
            write_alignment_blocks(out, hsp, true, positive_marks, 60)?;
        }
    }
    Ok(())
}

// ─── Format 5: BLAST XML 1 ─────────────────────────────────────────────────

fn fmt5_xml(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, "    <Iteration>")?;
    writeln!(out, "      <Iteration_iter-num>{}</Iteration_iter-num>", ctx.iter_num)?;
    writeln!(out, "      <Iteration_query-ID>Query_{}</Iteration_query-ID>", ctx.iter_num)?;
    writeln!(out, "      <Iteration_query-def>{}</Iteration_query-def>", xml_escape(ctx.query_title))?;
    writeln!(out, "      <Iteration_query-len>{}</Iteration_query-len>", ctx.query_len)?;
    writeln!(out, "      <Iteration_hits>")?;

    for (hi, r) in results.iter().enumerate() {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        writeln!(out, "        <Hit>")?;
        writeln!(out, "          <Hit_num>{}</Hit_num>", hi + 1)?;
        writeln!(out, "          <Hit_id>{}</Hit_id>", xml_escape(sid))?;
        writeln!(out, "          <Hit_def>{}</Hit_def>", xml_escape(&r.subject_title))?;
        writeln!(out, "          <Hit_len>{}</Hit_len>", r.subject_len)?;
        writeln!(out, "          <Hit_hsps>")?;

        for (hj, hsp) in r.hsps.iter().enumerate() {
            let positives = count_positive_chars(&hsp.midline);
            writeln!(out, "            <Hsp>")?;
            writeln!(out, "              <Hsp_num>{}</Hsp_num>", hj + 1)?;
            writeln!(out, "              <Hsp_bit-score>{:.4}</Hsp_bit-score>", hsp.bit_score)?;
            writeln!(out, "              <Hsp_score>{}</Hsp_score>", hsp.score)?;
            writeln!(out, "              <Hsp_evalue>{:.6e}</Hsp_evalue>", hsp.evalue)?;
            writeln!(out, "              <Hsp_query-from>{}</Hsp_query-from>", hsp.query_start + 1)?;
            writeln!(out, "              <Hsp_query-to>{}</Hsp_query-to>", hsp.query_end)?;
            writeln!(out, "              <Hsp_hit-from>{}</Hsp_hit-from>", hsp.subject_start + 1)?;
            writeln!(out, "              <Hsp_hit-to>{}</Hsp_hit-to>", hsp.subject_end)?;
            writeln!(out, "              <Hsp_query-frame>0</Hsp_query-frame>")?;
            writeln!(out, "              <Hsp_hit-frame>0</Hsp_hit-frame>")?;
            writeln!(out, "              <Hsp_identity>{}</Hsp_identity>", hsp.num_identities)?;
            writeln!(out, "              <Hsp_positive>{}</Hsp_positive>", positives)?;
            writeln!(out, "              <Hsp_gaps>{}</Hsp_gaps>", hsp.num_gaps)?;
            writeln!(out, "              <Hsp_align-len>{}</Hsp_align-len>", hsp.alignment_length)?;
            writeln!(out, "              <Hsp_qseq>{}</Hsp_qseq>", xml_escape(str_from_aln(&hsp.query_aln)))?;
            writeln!(out, "              <Hsp_hseq>{}</Hsp_hseq>", xml_escape(str_from_aln(&hsp.subject_aln)))?;
            writeln!(out, "              <Hsp_midline>{}</Hsp_midline>", xml_escape(str_from_aln(&hsp.midline)))?;
            writeln!(out, "            </Hsp>")?;
        }

        writeln!(out, "          </Hit_hsps>")?;
        writeln!(out, "        </Hit>")?;
    }

    writeln!(out, "      </Iteration_hits>")?;
    writeln!(out, "      <Iteration_stat>")?;
    writeln!(out, "        <Statistics>")?;
    writeln!(out, "          <Statistics_db-num>{}</Statistics_db-num>", ctx.db_num_seqs)?;
    writeln!(out, "          <Statistics_db-len>{}</Statistics_db-len>", ctx.db_len)?;
    writeln!(out, "          <Statistics_hsp-len>0</Statistics_hsp-len>")?;
    writeln!(out, "          <Statistics_eff-space>0</Statistics_eff-space>")?;
    writeln!(out, "          <Statistics_kappa>0</Statistics_kappa>")?;
    writeln!(out, "          <Statistics_lambda>0</Statistics_lambda>")?;
    writeln!(out, "          <Statistics_entropy>0</Statistics_entropy>")?;
    writeln!(out, "        </Statistics>")?;
    writeln!(out, "      </Iteration_stat>")?;
    writeln!(out, "    </Iteration>")
}

// ─── Formats 6/7/10: Tabular / CSV ─────────────────────────────────────────

/// All recognised tabular column identifiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TabularColumn {
    Qseqid,
    Sseqid,
    Pident,
    Length,
    Mismatch,
    Gapopen,
    Qstart,
    Qend,
    Sstart,
    Send,
    Evalue,
    Bitscore,
    // Extended columns
    Qlen,
    Slen,
    Nident,
    Positive,
    Gaps,
    Ppos,
    Qseq,
    Sseq,
    Btop,
    Staxid,
    Salltitles,
    Qcovs,
    QcovHsp,
    Score,
}

impl TabularColumn {
    fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "qseqid"     => Self::Qseqid,
            "sseqid"     => Self::Sseqid,
            "pident"     => Self::Pident,
            "length"     => Self::Length,
            "mismatch"   => Self::Mismatch,
            "gapopen"    => Self::Gapopen,
            "qstart"     => Self::Qstart,
            "qend"       => Self::Qend,
            "sstart"     => Self::Sstart,
            "send"       => Self::Send,
            "evalue"     => Self::Evalue,
            "bitscore"   => Self::Bitscore,
            "qlen"       => Self::Qlen,
            "slen"       => Self::Slen,
            "nident"     => Self::Nident,
            "positive"   => Self::Positive,
            "gaps"       => Self::Gaps,
            "ppos"       => Self::Ppos,
            "qseq"       => Self::Qseq,
            "sseq"       => Self::Sseq,
            "btop"       => Self::Btop,
            "staxid"     => Self::Staxid,
            "salltitles" => Self::Salltitles,
            "qcovs"      => Self::Qcovs,
            "qcovhsp"    => Self::QcovHsp,
            "score"      => Self::Score,
            _ => return None,
        })
    }

    fn default_columns() -> Vec<Self> {
        vec![
            Self::Qseqid, Self::Sseqid, Self::Pident, Self::Length,
            Self::Mismatch, Self::Gapopen,
            Self::Qstart, Self::Qend, Self::Sstart, Self::Send,
            Self::Evalue, Self::Bitscore,
        ]
    }

    fn header_name(&self) -> &'static str {
        match self {
            Self::Qseqid     => "qseqid",
            Self::Sseqid     => "sseqid",
            Self::Pident     => "pident",
            Self::Length     => "length",
            Self::Mismatch   => "mismatch",
            Self::Gapopen    => "gapopen",
            Self::Qstart     => "qstart",
            Self::Qend       => "qend",
            Self::Sstart     => "sstart",
            Self::Send       => "send",
            Self::Evalue     => "evalue",
            Self::Bitscore   => "bitscore",
            Self::Qlen       => "qlen",
            Self::Slen       => "slen",
            Self::Nident     => "nident",
            Self::Positive   => "positive",
            Self::Gaps       => "gaps",
            Self::Ppos       => "ppos",
            Self::Qseq       => "qseq",
            Self::Sseq       => "sseq",
            Self::Btop       => "btop",
            Self::Staxid     => "staxid",
            Self::Salltitles => "salltitles",
            Self::Qcovs      => "qcovs",
            Self::QcovHsp    => "qcovhsp",
            Self::Score      => "score",
        }
    }
}

fn fmt6_tabular(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    fmt: &OutputFormat,
) -> io::Result<()> {
    let sep = if fmt.fmt_id == 10 { "," } else { "\t" };
    let with_comments = fmt.fmt_id == 7;

    if with_comments {
        writeln!(out, "# BLAST{}", ctx.program.to_uppercase())?;
        writeln!(out, "# Query: {}", ctx.query_title)?;
        writeln!(out, "# Database: {}", ctx.db_path)?;
        let col_names: Vec<&str> = fmt.columns.iter().map(|c| c.header_name()).collect();
        writeln!(out, "# Fields: {}", col_names.join(", "))?;
        writeln!(out, "# {} hits found", results.iter().map(|r| r.hsps.len()).sum::<usize>())?;
    }

    let qid = first_word(ctx.query_title);

    for r in results {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        for hsp in &r.hsps {
            let fields: Vec<String> = fmt.columns.iter().map(|col| {
                tabular_field(col, qid, sid, hsp, r, ctx)
            }).collect();
            writeln!(out, "{}", fields.join(sep))?;
        }
    }

    if with_comments {
        writeln!(out, "# BLAST processed {} queries", 1)?;
    }

    Ok(())
}

fn tabular_field(
    col: &TabularColumn,
    qid: &str,
    sid: &str,
    hsp: &Hsp,
    r: &SearchResult,
    ctx: &SearchContext<'_>,
) -> String {
    let aln_len = hsp.alignment_length;
    let mismatches = aln_len.saturating_sub(hsp.num_identities + hsp.num_gaps);
    let gap_opens = count_gap_opens(&hsp.query_aln) + count_gap_opens(&hsp.subject_aln);
    let positives = count_positive_chars(&hsp.midline);

    match col {
        TabularColumn::Qseqid     => qid.to_string(),
        TabularColumn::Sseqid     => sid.to_string(),
        TabularColumn::Pident     => format!("{:.2}", hsp.percent_identity()),
        TabularColumn::Length     => aln_len.to_string(),
        TabularColumn::Mismatch   => mismatches.to_string(),
        TabularColumn::Gapopen    => gap_opens.to_string(),
        TabularColumn::Qstart     => (hsp.query_start + 1).to_string(),
        TabularColumn::Qend       => hsp.query_end.to_string(),
        TabularColumn::Sstart     => (hsp.subject_start + 1).to_string(),
        TabularColumn::Send       => hsp.subject_end.to_string(),
        TabularColumn::Evalue     => format!("{:.2e}", hsp.evalue),
        TabularColumn::Bitscore   => format!("{:.1}", hsp.bit_score),
        TabularColumn::Score      => hsp.score.to_string(),
        TabularColumn::Qlen       => ctx.query_len.to_string(),
        TabularColumn::Slen       => r.subject_len.to_string(),
        TabularColumn::Nident     => hsp.num_identities.to_string(),
        TabularColumn::Positive   => positives.to_string(),
        TabularColumn::Gaps       => hsp.num_gaps.to_string(),
        TabularColumn::Ppos       => {
            if aln_len == 0 { "0.00".to_string() }
            else { format!("{:.2}", 100.0 * positives as f64 / aln_len as f64) }
        }
        TabularColumn::Qseq       => str_from_aln(&hsp.query_aln).to_string(),
        TabularColumn::Sseq       => str_from_aln(&hsp.subject_aln).to_string(),
        TabularColumn::Btop       => build_btop(&hsp.query_aln, &hsp.subject_aln),
        TabularColumn::Staxid     => "N/A".to_string(),
        TabularColumn::Salltitles => r.subject_title.clone(),
        TabularColumn::Qcovs      => {
            if ctx.query_len == 0 { "0".to_string() }
            else {
                let cov_len = hsp.query_end - hsp.query_start;
                format!("{}", (100 * cov_len / ctx.query_len).min(100))
            }
        }
        TabularColumn::QcovHsp    => {
            if ctx.query_len == 0 { "0".to_string() }
            else {
                let cov_len = hsp.query_end - hsp.query_start;
                format!("{}", (100 * cov_len / ctx.query_len).min(100))
            }
        }
    }
}

// ─── Format 15: Single-file BLAST JSON ─────────────────────────────────────

fn fmt15_json(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, "  {{")?;
    writeln!(out, "    \"BlastOutput2\": [{{")?;
    writeln!(out, "      \"report\": {{")?;
    writeln!(out, "        \"program\": {},", json_str(ctx.program))?;
    writeln!(out, "        \"version\": \"blast-rs 0.1.0\",")?;
    writeln!(out, "        \"reference\": \"Rust BLAST implementation\",")?;
    writeln!(out, "        \"search_target\": {{ \"db\": {} }},", json_str(ctx.db_path))?;
    writeln!(out, "        \"params\": {{")?;
    writeln!(out, "          \"matrix\": {},", json_str(ctx.matrix))?;
    writeln!(out, "          \"expect\": {},", ctx.evalue_threshold)?;
    writeln!(out, "          \"gap_open\": {},", ctx.gap_open)?;
    writeln!(out, "          \"gap_extend\": {}", ctx.gap_extend)?;
    writeln!(out, "        }},")?;
    writeln!(out, "        \"results\": {{")?;
    writeln!(out, "          \"search\": {{")?;
    writeln!(out, "            \"query_id\": {},", json_str(&format!("Query_{}", ctx.iter_num)))?;
    writeln!(out, "            \"query_title\": {},", json_str(ctx.query_title))?;
    writeln!(out, "            \"query_len\": {},", ctx.query_len)?;
    writeln!(out, "            \"hits\": [")?;

    for (hi, r) in results.iter().enumerate() {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        let comma_r = if hi + 1 < results.len() { "," } else { "" };
        writeln!(out, "              {{")?;
        writeln!(out, "                \"num\": {},", hi + 1)?;
        writeln!(out, "                \"description\": [{{")?;
        writeln!(out, "                  \"id\": {},", json_str(sid))?;
        writeln!(out, "                  \"title\": {}", json_str(&r.subject_title))?;
        writeln!(out, "                }}],")?;
        writeln!(out, "                \"len\": {},", r.subject_len)?;
        writeln!(out, "                \"hsps\": [")?;

        for (hj, hsp) in r.hsps.iter().enumerate() {
            let positives = count_positive_chars(&hsp.midline);
            let mismatches = hsp.alignment_length.saturating_sub(hsp.num_identities + hsp.num_gaps);
            let comma_h = if hj + 1 < r.hsps.len() { "," } else { "" };
            writeln!(out, "                  {{")?;
            writeln!(out, "                    \"num\": {},", hj + 1)?;
            writeln!(out, "                    \"bit_score\": {:.4},", hsp.bit_score)?;
            writeln!(out, "                    \"score\": {},", hsp.score)?;
            writeln!(out, "                    \"evalue\": {:.6e},", hsp.evalue)?;
            writeln!(out, "                    \"identity\": {},", hsp.num_identities)?;
            writeln!(out, "                    \"positive\": {},", positives)?;
            writeln!(out, "                    \"gaps\": {},", hsp.num_gaps)?;
            writeln!(out, "                    \"align_len\": {},", hsp.alignment_length)?;
            writeln!(out, "                    \"query_from\": {},", hsp.query_start + 1)?;
            writeln!(out, "                    \"query_to\": {},", hsp.query_end)?;
            writeln!(out, "                    \"hit_from\": {},", hsp.subject_start + 1)?;
            writeln!(out, "                    \"hit_to\": {},", hsp.subject_end)?;
            writeln!(out, "                    \"qseq\": {},", json_str(str_from_aln(&hsp.query_aln)))?;
            writeln!(out, "                    \"hseq\": {},", json_str(str_from_aln(&hsp.subject_aln)))?;
            writeln!(out, "                    \"midline\": {}", json_str(str_from_aln(&hsp.midline)))?;
            writeln!(out, "                  }}{}", comma_h)?;
        }

        writeln!(out, "                ]")?;
        writeln!(out, "              }}{}", comma_r)?;
    }

    writeln!(out, "            ],")?;
    writeln!(out, "            \"stat\": {{")?;
    writeln!(out, "              \"db_num\": {},", ctx.db_num_seqs)?;
    writeln!(out, "              \"db_len\": {}", ctx.db_len)?;
    writeln!(out, "            }}")?;
    writeln!(out, "          }}")?;
    writeln!(out, "        }}")?;
    writeln!(out, "      }}")?;
    writeln!(out, "    }}]")?;
    writeln!(out, "  }}")
}

// ─── Format 16: Single-file BLAST XML2 ─────────────────────────────────────

fn fmt16_xml2(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, r#"<?xml version="1.0"?>"#)?;
    writeln!(out, r#"<BlastXML2>"#)?;
    writeln!(out, "  <BlastXML2_report>")?;
    writeln!(out, "    <Report>")?;
    writeln!(out, "      <Report_program>{}</Report_program>", xml_escape(ctx.program))?;
    writeln!(out, "      <Report_version>blast-rs 0.1.0</Report_version>")?;
    writeln!(out, "      <Report_search-target>")?;
    writeln!(out, "        <Target>")?;
    writeln!(out, "          <Target_db>{}</Target_db>", xml_escape(ctx.db_path))?;
    writeln!(out, "        </Target>")?;
    writeln!(out, "      </Report_search-target>")?;
    writeln!(out, "      <Report_params>")?;
    writeln!(out, "        <Parameters>")?;
    writeln!(out, "          <Parameters_matrix>{}</Parameters_matrix>", xml_escape(ctx.matrix))?;
    writeln!(out, "          <Parameters_expect>{}</Parameters_expect>", ctx.evalue_threshold)?;
    writeln!(out, "          <Parameters_gap-open>{}</Parameters_gap-open>", ctx.gap_open)?;
    writeln!(out, "          <Parameters_gap-extend>{}</Parameters_gap-extend>", ctx.gap_extend)?;
    writeln!(out, "        </Parameters>")?;
    writeln!(out, "      </Report_params>")?;
    writeln!(out, "      <Report_results>")?;
    writeln!(out, "        <Results>")?;
    writeln!(out, "          <Results_search>")?;
    writeln!(out, "            <Search>")?;
    writeln!(out, "              <Search_query-id>Query_{}</Search_query-id>", ctx.iter_num)?;
    writeln!(out, "              <Search_query-title>{}</Search_query-title>", xml_escape(ctx.query_title))?;
    writeln!(out, "              <Search_query-len>{}</Search_query-len>", ctx.query_len)?;
    writeln!(out, "              <Search_hits>")?;

    for (hi, r) in results.iter().enumerate() {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        writeln!(out, "                <Hit>")?;
        writeln!(out, "                  <Hit_num>{}</Hit_num>", hi + 1)?;
        writeln!(out, "                  <Hit_description>")?;
        writeln!(out, "                    <HitDescr>")?;
        writeln!(out, "                      <HitDescr_id>{}</HitDescr_id>", xml_escape(sid))?;
        writeln!(out, "                      <HitDescr_title>{}</HitDescr_title>", xml_escape(&r.subject_title))?;
        writeln!(out, "                    </HitDescr>")?;
        writeln!(out, "                  </Hit_description>")?;
        writeln!(out, "                  <Hit_len>{}</Hit_len>", r.subject_len)?;
        writeln!(out, "                  <Hit_hsps>")?;

        for (hj, hsp) in r.hsps.iter().enumerate() {
            let positives = count_positive_chars(&hsp.midline);
            writeln!(out, "                    <Hsp>")?;
            writeln!(out, "                      <Hsp_num>{}</Hsp_num>", hj + 1)?;
            writeln!(out, "                      <Hsp_bit-score>{:.4}</Hsp_bit-score>", hsp.bit_score)?;
            writeln!(out, "                      <Hsp_score>{}</Hsp_score>", hsp.score)?;
            writeln!(out, "                      <Hsp_evalue>{:.6e}</Hsp_evalue>", hsp.evalue)?;
            writeln!(out, "                      <Hsp_query-from>{}</Hsp_query-from>", hsp.query_start + 1)?;
            writeln!(out, "                      <Hsp_query-to>{}</Hsp_query-to>", hsp.query_end)?;
            writeln!(out, "                      <Hsp_hit-from>{}</Hsp_hit-from>", hsp.subject_start + 1)?;
            writeln!(out, "                      <Hsp_hit-to>{}</Hsp_hit-to>", hsp.subject_end)?;
            writeln!(out, "                      <Hsp_identity>{}</Hsp_identity>", hsp.num_identities)?;
            writeln!(out, "                      <Hsp_positive>{}</Hsp_positive>", positives)?;
            writeln!(out, "                      <Hsp_gaps>{}</Hsp_gaps>", hsp.num_gaps)?;
            writeln!(out, "                      <Hsp_align-len>{}</Hsp_align-len>", hsp.alignment_length)?;
            writeln!(out, "                      <Hsp_qseq>{}</Hsp_qseq>", xml_escape(str_from_aln(&hsp.query_aln)))?;
            writeln!(out, "                      <Hsp_hseq>{}</Hsp_hseq>", xml_escape(str_from_aln(&hsp.subject_aln)))?;
            writeln!(out, "                      <Hsp_midline>{}</Hsp_midline>", xml_escape(str_from_aln(&hsp.midline)))?;
            writeln!(out, "                    </Hsp>")?;
        }

        writeln!(out, "                  </Hit_hsps>")?;
        writeln!(out, "                </Hit>")?;
    }

    writeln!(out, "              </Search_hits>")?;
    writeln!(out, "              <Search_stat>")?;
    writeln!(out, "                <Statistics>")?;
    writeln!(out, "                  <Statistics_db-num>{}</Statistics_db-num>", ctx.db_num_seqs)?;
    writeln!(out, "                  <Statistics_db-len>{}</Statistics_db-len>", ctx.db_len)?;
    writeln!(out, "                </Statistics>")?;
    writeln!(out, "              </Search_stat>")?;
    writeln!(out, "            </Search>")?;
    writeln!(out, "          </Results_search>")?;
    writeln!(out, "        </Results>")?;
    writeln!(out, "      </Report_results>")?;
    writeln!(out, "    </Report>")?;
    writeln!(out, "  </BlastXML2_report>")?;
    writeln!(out, "</BlastXML2>")
}

// ─── Format 17: Subject FASTA sequences ────────────────────────────────────

fn fmt17_fasta(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    db: Option<&BlastDb>,
) -> io::Result<()> {
    let db = match db {
        Some(d) => d,
        None => {
            writeln!(out, "# Error: database not available for format 17")?;
            return Ok(());
        }
    };

    for r in results {
        let sid = subject_title(r);
        writeln!(out, ">{}", sid)?;

        let seq: Vec<u8> = match db.seq_type() {
            blast_db::index::SeqType::Protein => {
                db.get_sequence_protein(r.subject_oid).unwrap_or_default()
            }
            blast_db::index::SeqType::Nucleotide => {
                db.get_sequence_nucleotide(r.subject_oid).unwrap_or_default()
            }
        };

        for chunk in seq.chunks(60) {
            out.write_all(chunk)?;
            writeln!(out)?;
        }
    }
    Ok(())
}

// ─── Shared alignment helpers ───────────────────────────────────────────────

/// Write the Score / Expect / Identity / Gaps header for one HSP.
fn write_hsp_header(out: &mut impl Write, hsp: &Hsp, hsp_num: usize) -> io::Result<()> {
    writeln!(out, " Score = {:.1} bits ({}),  Expect = {:.2e}", hsp.bit_score, hsp.score, hsp.evalue)?;
    let positives = count_positive_chars(&hsp.midline);
    writeln!(out,
        " Identities = {}/{} ({:.0}%), Positives = {}/{} ({:.0}%), Gaps = {}/{}",
        hsp.num_identities, hsp.alignment_length, hsp.percent_identity(),
        positives, hsp.alignment_length, 100.0 * positives as f64 / hsp.alignment_length.max(1) as f64,
        hsp.num_gaps, hsp.alignment_length,
    )?;
    writeln!(out)
}

/// Print the alignment in 60-character blocks.
///
/// `query_anchored`: formats 1-4 — replace '+' midline chars with ' ' or use '+' for all.
/// `positive_marks`: formats 2/4 — use '+' for everything non-gap, ' ' for mismatch.
fn write_alignment_blocks(
    out: &mut impl Write,
    hsp: &Hsp,
    query_anchored: bool,
    positive_marks: bool,
    line_width: usize,
) -> io::Result<()> {
    let aln_len = hsp.query_aln.len();
    let mut pos = 0;
    let mut q_pos = hsp.query_start;
    let mut s_pos = hsp.subject_start;

    while pos < aln_len {
        let end = (pos + line_width).min(aln_len);
        let q_chunk = &hsp.query_aln[pos..end];
        let m_chunk = &hsp.midline[pos..end];
        let s_chunk = &hsp.subject_aln[pos..end];

        let q_chars = q_chunk.iter().filter(|&&b| b != b'-').count();
        let s_chars = s_chunk.iter().filter(|&&b| b != b'-').count();

        // Possibly remap midline for formats 1-4
        let midline_display: Vec<u8> = if query_anchored || positive_marks {
            m_chunk.iter().map(|&c| {
                if positive_marks {
                    // Format 2/4: '+' for identical, space for mismatch/gap
                    if c == b'|' { b'+' }
                    else if c == b'+' { b'+' }
                    else { b' ' }
                } else {
                    // Format 1/3: '|' for identical, space for everything else
                    if c == b'|' { b'|' } else { b' ' }
                }
            }).collect()
        } else {
            m_chunk.to_vec()
        };

        writeln!(out, "Query  {:>6}  {}  {}",
            q_pos + 1,
            std::str::from_utf8(q_chunk).unwrap_or("?"),
            q_pos + q_chars,
        )?;
        writeln!(out, "              {}",
            std::str::from_utf8(&midline_display).unwrap_or("?"),
        )?;
        writeln!(out, "Sbjct  {:>6}  {}  {}",
            s_pos + 1,
            std::str::from_utf8(s_chunk).unwrap_or("?"),
            s_pos + s_chars,
        )?;
        writeln!(out)?;

        q_pos += q_chars;
        s_pos += s_chars;
        pos = end;
    }
    Ok(())
}

/// Build BTOP (BLAST Traceback Operations) string.
/// Format: alternating run-lengths of matches and substitution pairs.
/// Matches: an integer (e.g., "5")
/// Substitutions: two characters (e.g., "TG" = query T aligned to subject G)
/// Gaps in query: "-" + subject char (e.g., "-A")
/// Gaps in subject: query char + "-" (e.g., "A-")
fn build_btop(query_aln: &[u8], subject_aln: &[u8]) -> String {
    let mut btop = String::new();
    let mut match_run = 0usize;

    let flush_matches = |s: &mut String, n: &mut usize| {
        if *n > 0 {
            s.push_str(&n.to_string());
            *n = 0;
        }
    };

    for (&q, &s) in query_aln.iter().zip(subject_aln.iter()) {
        if q == b'-' {
            flush_matches(&mut btop, &mut match_run);
            btop.push('-');
            btop.push(s as char);
        } else if s == b'-' {
            flush_matches(&mut btop, &mut match_run);
            btop.push(q as char);
            btop.push('-');
        } else if q == s {
            match_run += 1;
        } else {
            flush_matches(&mut btop, &mut match_run);
            btop.push(q as char);
            btop.push(s as char);
        }
    }
    flush_matches(&mut btop, &mut match_run);
    btop
}

// ─── Utility functions ──────────────────────────────────────────────────────

fn first_word(s: &str) -> &str {
    s.split_whitespace().next().unwrap_or(s)
}

fn subject_title(r: &SearchResult) -> String {
    if !r.subject_accession.is_empty() {
        format!("{} {}", r.subject_accession, r.subject_title)
    } else {
        r.subject_title.clone()
    }
}

fn str_from_aln(aln: &[u8]) -> &str {
    std::str::from_utf8(aln).unwrap_or("?")
}

fn count_gap_opens(aln: &[u8]) -> usize {
    let mut opens = 0;
    let mut in_gap = false;
    for &b in aln {
        if b == b'-' {
            if !in_gap { opens += 1; in_gap = true; }
        } else {
            in_gap = false;
        }
    }
    opens
}

/// Count '|' and '+' characters in the midline (identities + positives).
fn count_positive_chars(midline: &[u8]) -> usize {
    midline.iter().filter(|&&c| c == b'|' || c == b'+').count()
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
     .replace('<', "&lt;")
     .replace('>', "&gt;")
     .replace('"', "&quot;")
     .replace('\'', "&apos;")
}

fn json_str(s: &str) -> String {
    let escaped = s
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t");
    format!("\"{}\"", escaped)
}
