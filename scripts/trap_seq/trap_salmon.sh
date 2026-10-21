#!/usr/bin/env bash
# Prepare the transcript FASTA once and quantify one STAR-aligned sample.
# Requires gffread for preparation and Salmon 1.10.3 for quantification.
# Edit settings, then: bash trap_salmon.sh
set -euo pipefail

ROOT="/path/to/TRAP_results"
REFERENCE_DIR="/path/to/GRCh38_GENCODEv43"
CELL="EndoC"                  # EndoC or RPE
CONDITION="SCR"               # SCR, EX2, EX7 or EX2EX7
REPLICATE="1"                 # 1 or 2
STAGE="all"                   # prepare, quant or all; prepare the shared FASTA once.
GENOME_FASTA="$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa"
GTF="$REFERENCE_DIR/gencode.v43.primary_assembly.annotation.gtf"
TRANSCRIPTS="$REFERENCE_DIR/gencode.v43.primary_assembly.gffread.transcripts.fa"

case "$CELL" in EndoC|RPE) ;; *) echo 'Invalid CELL.' >&2; exit 1 ;; esac
case "$CONDITION" in SCR|EX2|EX7|EX2EX7) ;; *) echo 'Invalid CONDITION.' >&2; exit 1 ;; esac
case "$REPLICATE" in 1|2) ;; *) echo 'Invalid REPLICATE.' >&2; exit 1 ;; esac
case "$STAGE" in prepare|quant|all) ;; *) echo 'Invalid STAGE.' >&2; exit 1 ;; esac
SAMPLE="${CELL}_TRAPseq_HNRNPD_${CONDITION}_rep${REPLICATE}"
BAM="$ROOT/$CELL/$SAMPLE/STAR/${SAMPLE}_Aligned.toTranscriptome.out.bam"
OUT_DIR="$ROOT/$CELL/$SAMPLE/Salmon"

# Original gffread command; only the output filename has been renamed.
if [[ "$STAGE" == prepare || "$STAGE" == all ]]; then
    command -v gffread >/dev/null 2>&1 || { echo 'Activate gffread first.' >&2; exit 1; }
    [[ -s "$GENOME_FASTA" && -s "$GTF" ]] || { echo 'Check genome FASTA and GTF.' >&2; exit 1; }
    [[ ! -e "$TRANSCRIPTS" ]] || { echo 'Transcript FASTA exists; select STAGE=quant.' >&2; exit 1; }
    gffread "$GTF" -g "$GENOME_FASTA" -w "$TRANSCRIPTS"
fi

# Alignment-based quantification: no Salmon index or further alignment required.
if [[ "$STAGE" == quant || "$STAGE" == all ]]; then
    command -v salmon >/dev/null 2>&1 || { echo 'Activate Salmon first.' >&2; exit 1; }
    [[ -s "$TRANSCRIPTS" && -s "$BAM" && -s "$GTF" ]] || { echo 'Check transcript FASTA, STAR transcriptome BAM and GTF.' >&2; exit 1; }
    [[ ! -e "$OUT_DIR" ]] || { echo 'Salmon output directory already exists.' >&2; exit 1; }
    mkdir -p "$(dirname "$OUT_DIR")"
    salmon quant -l A \
        -p 12 \
        -t "$TRANSCRIPTS" \
        -a "$BAM" \
        -g "$GTF" \
        -o "$OUT_DIR"
fi
