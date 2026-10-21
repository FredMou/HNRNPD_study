#!/usr/bin/env bash
# IsoQuant: quantify all 12 samples together for one cell type.
# Methods version: IsoQuant 3.7.0. Also requires samtools.
# Edit settings, then run: bash long_read_isoquant.sh

set -euo pipefail

ROOT="/path/to/long_read_results"   # Same ROOT as the preceding scripts.
CELL="EndoC"                      # EndoC or RPE
REFERENCE_DIR="/path/to/GRCh38_GENCODEv43"
FASTA="$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa"
GTF="$REFERENCE_DIR/gencode.v43.primary_assembly.annotation.gtf"
OUT_DIR="$ROOT/$CELL/IsoQuant"

case "$CELL" in EndoC|RPE) ;; *) echo 'Invalid CELL.' >&2; exit 1 ;; esac
command -v isoquant.py >/dev/null 2>&1 || { echo 'Activate IsoQuant first.' >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo 'Activate samtools first.' >&2; exit 1; }
[[ -s "$FASTA" && -r "$FASTA" && -s "$GTF" && -r "$GTF" ]] || { echo 'Check reference paths.' >&2; exit 1; }
[[ ! -e "$OUT_DIR" ]] || { echo 'IsoQuant outputs exist; use a new OUT_DIR.' >&2; exit 1; }

# Preserve the original sample order: SCR, EX2, EX7, EX2EX7; replicates 1-3.
BAMS=()
LABELS=()
for CONDITION in SCR EX2 EX7 EX2EX7; do
    for REPLICATE in 1 2 3; do
        SAMPLE_ID="${CELL}_LongReadRNAseq_HNRNPD_${CONDITION}_rep${REPLICATE}"
        ALIGN_DIR="$ROOT/$CELL/$SAMPLE_ID/minimap2"
        BAM="$ALIGN_DIR/sorted_PychopperFilt_${SAMPLE_ID}.bam"
        SAM="$ALIGN_DIR/sorted_PychopperFilt_${SAMPLE_ID}.sam"
        if [[ "$CELL" == EndoC && ! -s "$BAM" ]]; then
            [[ -s "$SAM" && -r "$SAM" ]] || { echo "Missing sorted SAM: $SAM" >&2; exit 1; }
        else
            [[ -s "$BAM" && -r "$BAM" ]] || { echo "Missing sorted BAM: $BAM" >&2; exit 1; }
        fi
        BAMS+=("$BAM")
        LABELS+=("$SAMPLE_ID")
    done
done

# EndoC alignment stops at sorted SAM. Convert format only; do not filter reads.
# Index both cell types' BAMs for IsoQuant's regional access.
for BAM in "${BAMS[@]}"; do
    if [[ "$CELL" == EndoC && ! -s "$BAM" ]]; then
        samtools view -Sb "${BAM%.bam}.sam" > "$BAM"
    fi
    if [[ ! -s "$BAM.bai" && ! -s "${BAM%.bam}.bai" && ! -s "$BAM.csi" ]]; then
        samtools index "$BAM"
    fi
done

mkdir -p "$OUT_DIR/tmp"
export TMPDIR="$OUT_DIR/tmp"

# Original IsoQuant analysis options; GEO sample names replace barcode labels.
# Reference counts and transcript-discovery outputs may both be generated.
isoquant.py --reference "$FASTA" \
    --complete_genedb --genedb "$GTF" \
    --bam "${BAMS[@]}" \
    --label "${LABELS[@]}" \
    --data_type nanopore \
    --prefix isoquant \
    --threads 18 \
    --sqanti_output \
    --output "$OUT_DIR/isoquant_output"
