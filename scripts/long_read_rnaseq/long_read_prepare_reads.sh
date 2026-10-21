#!/usr/bin/env bash
# Long-read RNA-seq: EndoC NanoPlot QC, chopper filtering and Pychopper.
# Edit settings, then run: bash long_read_prepare_reads.sh
# Methods versions: NanoPlot 1.42.0, chopper 0.9.2, Pychopper 2.7.10.

set -euo pipefail

ROOT="/path/to/long_read_results"
RAW_DIR="/path/to/GEO_FASTQs"
CELL="EndoC"          # EndoC or RPE
CONDITION="SCR"       # SCR, EX2, EX7 or EX2EX7
REPLICATE="1"         # 1, 2 or 3

case "$CELL" in EndoC|RPE) ;; *) echo 'Invalid CELL.' >&2; exit 1 ;; esac
case "$CONDITION" in SCR|EX2|EX7|EX2EX7) ;; *) echo 'Invalid CONDITION.' >&2; exit 1 ;; esac
case "$REPLICATE" in 1|2|3) ;; *) echo 'Invalid REPLICATE.' >&2; exit 1 ;; esac

# Names match the GEO submission metadata and count-matrix sample columns.
SAMPLE_ID="${CELL}_LongReadRNAseq_HNRNPD_${CONDITION}_rep${REPLICATE}"
INPUT_FASTQ="$RAW_DIR/${SAMPLE_ID}.fq.gz"
SAMPLE_DIR="$ROOT/$CELL/$SAMPLE_ID"
READ_DIR="$SAMPLE_DIR/reads"
FILTERED_FASTQ="$READ_DIR/filtered_q7_l500_${SAMPLE_ID}.fq.gz"
FULL_LENGTH_FASTQ="$READ_DIR/Pychopper_full_length_output_${SAMPLE_ID}.fq.gz"

[[ -s "$INPUT_FASTQ" && -r "$INPUT_FASTQ" ]] || { echo "Missing FASTQ: $INPUT_FASTQ" >&2; exit 1; }
[[ ! -e "$READ_DIR" ]] || { echo 'Read-preparation outputs exist; use a new ROOT.' >&2; exit 1; }
command -v chopper >/dev/null 2>&1 || { echo 'Activate chopper first.' >&2; exit 1; }
command -v pychopper >/dev/null 2>&1 || { echo 'Activate Pychopper first.' >&2; exit 1; }
command -v gunzip >/dev/null 2>&1 || { echo 'gunzip is required.' >&2; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo 'gzip is required.' >&2; exit 1; }
if [[ "$CELL" == EndoC ]]; then
    command -v NanoPlot >/dev/null 2>&1 || { echo 'Activate NanoPlot first.' >&2; exit 1; }
    [[ ! -e "$SAMPLE_DIR/NanoPlot" ]] || { echo 'NanoPlot output exists.' >&2; exit 1; }
fi
mkdir -p "$READ_DIR"

# 1. Raw-read QC was performed for EndoC only.
if [[ "$CELL" == EndoC ]]; then
    NanoPlot --fastq "$INPUT_FASTQ" --N50 -o "$SAMPLE_DIR/NanoPlot"
fi

# 2. Preserve the original minimum quality and length thresholds.
gunzip -c "$INPUT_FASTQ" | chopper -q 7 -l 500 | gzip > "$FILTERED_FASTQ"

# 3. Identify, orient and trim full-length cDNA reads.
pychopper -t 16 -m edlib \
    -S "$READ_DIR/statistics_${SAMPLE_ID}.tsv" \
    -r "$READ_DIR/report_${SAMPLE_ID}.pdf" \
    -u "$READ_DIR/unclassified_${SAMPLE_ID}.fq" \
    -w "$READ_DIR/rescued_${SAMPLE_ID}.fq" \
    "$FILTERED_FASTQ" "$FULL_LENGTH_FASTQ"

# Pass this full-length FASTQ to the subsequent minimap2 alignment script.
printf 'Full-length FASTQ: %s\n' "$FULL_LENGTH_FASTQ"
