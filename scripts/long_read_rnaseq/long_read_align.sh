#!/usr/bin/env bash
# Minimap2 alignment of Pychopper full-length reads; minimap2 2.29 and samtools.
# Edit settings, then run: bash long_read_align.sh

set -euo pipefail

# Match these settings to long_read_prepare_reads.sh.
ROOT="/path/to/long_read_results"
CELL="EndoC"          # EndoC or RPE
CONDITION="SCR"       # SCR, EX2, EX7 or EX2EX7
REPLICATE="1"         # 1, 2 or 3
FASTA="/path/to/GRCh38.primary_assembly.genome.fa"

# Decompress the repository BED once before alignment.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
JUNCTION_BED="$SCRIPT_DIR/data/junctions.GENECODEv43.bed"

case "$CELL" in EndoC|RPE) ;; *) echo 'Invalid CELL.' >&2; exit 1 ;; esac
case "$CONDITION" in SCR|EX2|EX7|EX2EX7) ;; *) echo 'Invalid CONDITION.' >&2; exit 1 ;; esac
case "$REPLICATE" in 1|2|3) ;; *) echo 'Invalid REPLICATE.' >&2; exit 1 ;; esac

SAMPLE_ID="${CELL}_LongReadRNAseq_HNRNPD_${CONDITION}_rep${REPLICATE}"
SAMPLE_DIR="$ROOT/$CELL/$SAMPLE_ID"
INPUT_FASTQ="$SAMPLE_DIR/reads/Pychopper_full_length_output_${SAMPLE_ID}.fq.gz"
ALIGN_DIR="$SAMPLE_DIR/minimap2"
UNSORTED_SAM="$ALIGN_DIR/unsorted_PychopperFilt_${SAMPLE_ID}.sam"
SORTED_SAM="$ALIGN_DIR/sorted_PychopperFilt_${SAMPLE_ID}.sam"
SORTED_BAM="$ALIGN_DIR/sorted_PychopperFilt_${SAMPLE_ID}.bam"

command -v minimap2 >/dev/null 2>&1 || { echo 'Activate minimap2 first.' >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo 'Activate samtools first.' >&2; exit 1; }
[[ -s "$INPUT_FASTQ" && -r "$INPUT_FASTQ" ]] || { echo "Missing full-length FASTQ: $INPUT_FASTQ" >&2; exit 1; }
[[ -s "$FASTA" && -r "$FASTA" ]] || { echo 'Check FASTA path.' >&2; exit 1; }
[[ -s "$JUNCTION_BED" && -r "$JUNCTION_BED" ]] || { echo 'Check BED path; decompress data/junctions.GENECODEv43.bed.gz first.' >&2; exit 1; }
[[ ! -e "$ALIGN_DIR" ]] || { echo 'Alignment outputs exist; use a new ROOT.' >&2; exit 1; }
mkdir -p "$ALIGN_DIR"

# Preserve the original alignment options for both cell types.
minimap2 -ax splice -ub -k14 -w 4 \
    --junc-bed "$JUNCTION_BED" --secondary=no \
    "$FASTA" "$INPUT_FASTQ" > "$UNSORTED_SAM"

samtools sort -o "$SORTED_SAM" "$UNSORTED_SAM"

# The supplied RPE scripts additionally converted SAM to BAM and ran flagstat.
if [[ "$CELL" == RPE ]]; then
    samtools view -Sb "$SORTED_SAM" > "$SORTED_BAM"
    samtools flagstat "$SORTED_BAM" > "$ALIGN_DIR/sorted_PychopperFilt_${SAMPLE_ID}_flagstat"
fi

printf 'Sorted SAM: %s\n' "$SORTED_SAM"
if [[ "$CELL" == RPE ]]; then
    printf 'Sorted BAM: %s\n' "$SORTED_BAM"
fi
