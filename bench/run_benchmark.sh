#!/usr/bin/env bash
# Performance benchmark: blast-rs vs NCBI BLAST (blastp)
# Generates synthetic protein data, runs timed searches, compares correctness.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BLAST_RS="$REPO_ROOT/target/release/blast-cli"
NCBI_BLASTP="$(which blastp)"
NCBI_MAKEDB="$(which makeblastdb)"
GNU_TIME="$(which time)"

RESULTS_DIR="$SCRIPT_DIR/results"
TMPDIR="$SCRIPT_DIR/tmp"
SUMMARY="$RESULTS_DIR/summary.tsv"

DB_SIZES=(100 1000 10000)
QUERY_LENGTHS=(50 300 1000)
THREAD_COUNTS=(1 2 4 8)
EVALUE="1e-5"

# ── Helpers ─────────────────────────────────────────────────────────────────

generate_fasta() {
    local num_seqs=$1 seq_len=$2 outfile=$3
    python3 - "$num_seqs" "$seq_len" "$outfile" <<'PYEOF'
import sys, random
random.seed(42)  # reproducible
num_seqs, seq_len = int(sys.argv[1]), int(sys.argv[2])
outfile = sys.argv[3]
aa = "ACDEFGHIKLMNPQRSTVWY"
# Swiss-Prot amino acid frequencies (%)
wt = [8.25,1.37,5.45,6.75,3.86,7.07,2.27,5.96,5.84,9.66,
      2.42,4.06,4.70,3.93,5.53,6.56,5.34,6.87,1.08,2.92]
with open(outfile, 'w') as f:
    for i in range(num_seqs):
        actual_len = max(20, int(seq_len * random.uniform(0.8, 1.2)))
        f.write(f">seq_{i+1} synthetic protein length~{actual_len}\n")
        seq = ''.join(random.choices(aa, weights=wt, k=actual_len))
        for j in range(0, len(seq), 70):
            f.write(seq[j:j+70] + "\n")
PYEOF
}

parse_wall_clock() {
    # Parse GNU time -v "Elapsed (wall clock)" line → seconds
    local time_log=$1
    local wall
    wall=$(grep "Elapsed (wall clock)" "$time_log" | sed 's/.*): //')
    python3 -c "
parts = '$wall'.split(':')
if len(parts) == 3:
    print(f'{float(parts[0])*3600 + float(parts[1])*60 + float(parts[2]):.3f}')
elif len(parts) == 2:
    print(f'{float(parts[0])*60 + float(parts[1]):.3f}')
else:
    print(f'{float(parts[0]):.3f}')
"
}

run_and_record() {
    local test_name=$1 tool_name=$2 threads=$3
    shift 3
    local label="${tool_name}_${test_name}"
    local time_log="$RESULTS_DIR/${label}.time"
    local output="$RESULTS_DIR/${label}.tsv"

    echo "  Running: $label ..."
    "$GNU_TIME" -v "$@" > "$output" 2> "$time_log" || true

    local wall_sec rss hits
    wall_sec=$(parse_wall_clock "$time_log")
    rss=$(grep "Maximum resident set size" "$time_log" | awk '{print $NF}')
    hits=$(wc -l < "$output" | tr -d ' ')

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$test_name" "$tool_name" "$threads" "$wall_sec" "$rss" "$hits" >> "$SUMMARY"
}

# ── Setup ───────────────────────────────────────────────────────────────────

echo "=== blast-rs vs NCBI BLAST benchmark ==="
echo ""
echo "NCBI blastp: $NCBI_BLASTP ($(blastp -version 2>&1 | head -1))"
echo "blast-rs:    $BLAST_RS"
echo "GNU time:    $GNU_TIME"
echo "CPUs:        $(nproc)"
echo ""

rm -rf "$RESULTS_DIR" "$TMPDIR"
mkdir -p "$RESULTS_DIR" "$TMPDIR"

# ── Phase 1: Generate synthetic data ────────────────────────────────────────

echo "--- Generating synthetic protein data ---"
for size in "${DB_SIZES[@]}"; do
    echo "  Database: $size sequences x ~300aa"
    generate_fasta "$size" 300 "$TMPDIR/db_${size}.faa"
done

for qlen in "${QUERY_LENGTHS[@]}"; do
    echo "  Query: 5 sequences x ~${qlen}aa"
    generate_fasta 5 "$qlen" "$TMPDIR/query_${qlen}aa.faa"
done
echo ""

# ── Phase 2: Build databases ───────────────────────────────────────────────

echo "--- Building databases ---"
for size in "${DB_SIZES[@]}"; do
    echo "  NCBI makeblastdb ($size seqs)..."
    "$NCBI_MAKEDB" -in "$TMPDIR/db_${size}.faa" -dbtype prot \
        -out "$TMPDIR/ncbi_db_${size}" > /dev/null 2>&1

    echo "  blast-rs makeblastdb ($size seqs)..."
    "$BLAST_RS" makeblastdb -i "$TMPDIR/db_${size}.faa" --dbtype prot \
        -o "$TMPDIR/rs_db_${size}" > /dev/null 2>&1
done
echo ""

# ── Phase 3: Timed searches ───────────────────────────────────────────────

echo "--- Running timed searches (72 runs) ---"
printf "test_name\ttool\tthreads\twall_time_sec\tpeak_rss_kb\tnum_hits\n" > "$SUMMARY"

for size in "${DB_SIZES[@]}"; do
  for qlen in "${QUERY_LENGTHS[@]}"; do
    for threads in "${THREAD_COUNTS[@]}"; do
        test_name="db${size}_q${qlen}aa_t${threads}"

        run_and_record "$test_name" "ncbi" "$threads" \
            "$NCBI_BLASTP" -query "$TMPDIR/query_${qlen}aa.faa" \
            -db "$TMPDIR/ncbi_db_${size}" \
            -outfmt 6 -evalue "$EVALUE" -num_threads "$threads"

        run_and_record "$test_name" "blast-rs" "$threads" \
            "$BLAST_RS" blastp -q "$TMPDIR/query_${qlen}aa.faa" \
            --db "$TMPDIR/rs_db_${size}" \
            --outfmt 6 --evalue "$EVALUE" --num_threads "$threads"
    done
  done
done
echo ""

# ── Phase 4: Correctness comparison ───────────────────────────────────────

echo "--- Correctness comparison (threads=1) ---"
echo ""

for size in "${DB_SIZES[@]}"; do
  for qlen in "${QUERY_LENGTHS[@]}"; do
    test_name="db${size}_q${qlen}aa_t1"
    ncbi_out="$RESULTS_DIR/ncbi_${test_name}.tsv"
    rs_out="$RESULTS_DIR/blast-rs_${test_name}.tsv"

    ncbi_hits=$(awk -F'\t' '{print $1"\t"$2}' "$ncbi_out" | sort -u | wc -l | tr -d ' ')
    rs_hits=$(awk -F'\t' '{print $1"\t"$2}' "$rs_out" | sort -u | wc -l | tr -d ' ')

    shared=$(comm -12 \
        <(awk -F'\t' '{print $1"\t"$2}' "$ncbi_out" | sort -u) \
        <(awk -F'\t' '{print $1"\t"$2}' "$rs_out" | sort -u) | wc -l | tr -d ' ')

    printf "  %-30s  ncbi=%4s  blast-rs=%4s  shared=%4s\n" \
        "$test_name" "$ncbi_hits" "$rs_hits" "$shared"

    diff <(sort "$ncbi_out") <(sort "$rs_out") \
        > "$RESULTS_DIR/diff_${test_name}.txt" 2>&1 || true
  done
done
echo ""

# ── Phase 5: Summary ──────────────────────────────────────────────────────

echo "=== Benchmark Summary ==="
echo ""
column -t -s $'\t' "$SUMMARY"
echo ""
echo "Results directory: $RESULTS_DIR/"
echo "Summary TSV:       $SUMMARY"
