#!/usr/bin/env bash
# Trim one TRAP-seq sample using the sequencing service's Cutadapt 4.9 command.
# Edit settings, then: bash trap_cutadapt.sh
set -euo pipefail

RAW_DIR="/path/to/untrimmed_GEO_FASTQs"
ROOT="/path/to/TRAP_results"
CELL="EndoC"                  # EndoC or RPE
CONDITION="SCR"               # SCR, EX2, EX7 or EX2EX7
REPLICATE="1"                 # 1 or 2

case "$CELL" in EndoC) MIN_LENGTH=75 ;; RPE) MIN_LENGTH=40 ;; *) echo 'Invalid CELL.' >&2; exit 1 ;; esac
case "$CONDITION" in SCR|EX2|EX7|EX2EX7) ;; *) echo 'Invalid CONDITION.' >&2; exit 1 ;; esac
case "$REPLICATE" in 1|2) ;; *) echo 'Invalid REPLICATE.' >&2; exit 1 ;; esac
command -v cutadapt >/dev/null 2>&1 || { echo 'Activate Cutadapt 4.9 first.' >&2; exit 1; }
[[ "$(cutadapt --version)" == 4.9 ]] || { echo 'Use Cutadapt 4.9.' >&2; exit 1; }

SAMPLE="${CELL}_TRAPseq_HNRNPD_${CONDITION}_rep${REPLICATE}"
OUT_DIR="$ROOT/$CELL/$SAMPLE/reads"
[[ ! -e "$OUT_DIR" ]] || { echo 'Output directory already exists.' >&2; exit 1; }

# Keep the service's named adapter, including its leading N.
# -u 1 -U 1 removes the library-derived first base from each mate.
trim_pair() {
    cutadapt \
        -a Nextera=NCTGTCTCTTATACACATCT \
        -A Nextera=NCTGTCTCTTATACACATCT \
        --cores 8 \
        --quality-cutoff=22 \
        --minimum-length="$MIN_LENGTH" \
        --error-rate=0.11 \
        -u 1 -U 1 \
        --overlap=5 \
        --pair-filter=any \
        --output="$OUT_DIR/${1}_R1.fastq.gz" \
        --paired-output="$OUT_DIR/${1}_R2.fastq.gz" \
        "$RAW_DIR/${1}_R1.fastq.gz" \
        "$RAW_DIR/${1}_R2.fastq.gz" \
        > "$OUT_DIR/${1}_cutadapt.log"
}

if [[ "$CELL" == EndoC ]]; then
    [[ -s "$RAW_DIR/${SAMPLE}_R1.fastq.gz" && -s "$RAW_DIR/${SAMPLE}_R2.fastq.gz" ]] || { echo 'Check both EndoC FASTQs.' >&2; exit 1; }
    mkdir -p "$OUT_DIR"
    trim_pair "$SAMPLE"
else
    [[ -s "$RAW_DIR/${SAMPLE}_L3_R1.fastq.gz" && -s "$RAW_DIR/${SAMPLE}_L3_R2.fastq.gz" && -s "$RAW_DIR/${SAMPLE}_L4_R1.fastq.gz" && -s "$RAW_DIR/${SAMPLE}_L4_R2.fastq.gz" ]] || { echo 'Check all four RPE lane FASTQs.' >&2; exit 1; }
    mkdir -p "$OUT_DIR"
    trim_pair "${SAMPLE}_L3"
    trim_pair "${SAMPLE}_L4"

    # Merge after trimming, keeping matching L3/L4 order for both mates.
    cat "$OUT_DIR/${SAMPLE}_L3_R1.fastq.gz" "$OUT_DIR/${SAMPLE}_L4_R1.fastq.gz" > "$OUT_DIR/${SAMPLE}_R1.fastq.gz"
    cat "$OUT_DIR/${SAMPLE}_L3_R2.fastq.gz" "$OUT_DIR/${SAMPLE}_L4_R2.fastq.gz" > "$OUT_DIR/${SAMPLE}_R2.fastq.gz"
fi

