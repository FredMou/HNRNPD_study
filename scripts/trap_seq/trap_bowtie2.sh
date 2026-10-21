#!/usr/bin/env bash
# Build the rRNA/tRNA index and filter one TRAP-seq sample.
# Requires Bowtie2 2.5.4. Edit settings, then: bash trap_bowtie2.sh
set -euo pipefail

ROOT="/path/to/TRAP_results"
CELL="EndoC"                  # EndoC or RPE
CONDITION="SCR"               # SCR, EX2, EX7 or EX2EX7
REPLICATE="1"                 # 1 or 2
STAGE="all"                   # index, filter or all; build the shared index once.
REFERENCE_DIR="/path/to/TRAP_reference"
REFERENCE_FASTA="$REFERENCE_DIR/rRNA_tRNA_combined.fasta"
INDEX_PREFIX="$ROOT/Bowtie2_index/rRNA_tRNA_combined_idx"

case "$CELL" in EndoC|RPE) ;; *) echo 'Invalid CELL.' >&2; exit 1 ;; esac
case "$CONDITION" in SCR|EX2|EX7|EX2EX7) ;; *) echo 'Invalid CONDITION.' >&2; exit 1 ;; esac
case "$REPLICATE" in 1|2) ;; *) echo 'Invalid REPLICATE.' >&2; exit 1 ;; esac
case "$STAGE" in index|filter|all) ;; *) echo 'Invalid STAGE.' >&2; exit 1 ;; esac
SAMPLE="${CELL}_TRAPseq_HNRNPD_${CONDITION}_rep${REPLICATE}"
READ_DIR="$ROOT/$CELL/$SAMPLE/reads"
R1="$READ_DIR/${SAMPLE}_R1.fastq.gz"
R2="$READ_DIR/${SAMPLE}_R2.fastq.gz"
OUT_DIR="$ROOT/$CELL/$SAMPLE/Bowtie2"

# Use the reference from trap_prepare_reference.sh for both cell types.
if [[ "$STAGE" == index || "$STAGE" == all ]]; then
    command -v bowtie2-build >/dev/null 2>&1 || { echo 'Activate Bowtie2 first.' >&2; exit 1; }
    [[ -s "$REFERENCE_FASTA" ]] || { echo 'Check REFERENCE_FASTA.' >&2; exit 1; }
    [[ ! -e "${INDEX_PREFIX}.1.bt2" && ! -e "${INDEX_PREFIX}.1.bt2l" ]] || { echo 'Index exists; select STAGE=filter.' >&2; exit 1; }
    mkdir -p "$(dirname "$INDEX_PREFIX")"
    bowtie2-build "$REFERENCE_FASTA" "$INDEX_PREFIX"
fi

if [[ "$STAGE" == filter || "$STAGE" == all ]]; then
    command -v bowtie2 >/dev/null 2>&1 || { echo 'Activate Bowtie2 first.' >&2; exit 1; }
    [[ -s "$R1" && -s "$R2" ]] || { echo 'Check trimmed R1/R2 paths; combine RPE lanes first.' >&2; exit 1; }
    EXT=bt2
    [[ -s "${INDEX_PREFIX}.1.bt2" ]] || EXT=bt2l
    for PART in 1 2 3 4 rev.1 rev.2; do
        [[ -s "${INDEX_PREFIX}.${PART}.${EXT}" ]] || { echo 'Bowtie2 index is missing or incomplete.' >&2; exit 1; }
    done
    [[ ! -e "$OUT_DIR" ]] || { echo 'Bowtie2 output directory already exists.' >&2; exit 1; }
    mkdir -p "$OUT_DIR"

    # Original options: retain pairs without concordant rRNA/tRNA alignments.
    bowtie2 -x "$INDEX_PREFIX" \
        -1 "$R1" -2 "$R2" \
        --un-conc-gz "$OUT_DIR/${SAMPLE}_rm_rtRNA_%.fq.gz" \
        -S /dev/null -p 10 \
        2> "$OUT_DIR/${SAMPLE}_rtRNA_stats.txt"
fi
