#!/usr/bin/env bash
# STAR index generation and single-end alignment of replicates A/B.
# Requires STAR 2.7.11a and Linux-compatible zcat. Edit settings, then:
# bash iclip_star.sh

set -euo pipefail

# Use the same ROOT, CELL and ISOFORM in all scripts.
ROOT="/path/to/iclip_results"
CELL="EndoC"                 # EndoC or RPE
ISOFORM="p37"                # p37, p40, p42 or p45

case "$CELL" in
    EndoC) BAM_LABEL=EndoC ;;
    RPE) BAM_LABEL=HRPE ;;
    *) echo 'CELL must be EndoC or RPE.' >&2; exit 1 ;;
esac
case "$ISOFORM" in
    p37|p40|p42|p45) ;;
    *) echo 'Choose ISOFORM=p37, p40, p42 or p45.' >&2; exit 1 ;;
esac
SAMPLE_DIR="$ROOT/$CELL/PureCLIP_${ISOFORM}_FWD"
SAMPLE="${BAM_LABEL}_${ISOFORM}"

STAGE="all"                 # index, align or all; use align after building index once.
REFERENCE_DIR="/path/to/GRCh38_GENCODEv43"
INDEX_DIR="/path/to/STAR_index_v43"
FASTQ_DIR="/path/to/selected_cell_FASTQs"
FASTA="$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa"
GTF="$REFERENCE_DIR/gencode.v43.primary_assembly.annotation.gtf"
FASTQ_A="$FASTQ_DIR/${CELL}_iCLIP_HNRNPD_${ISOFORM}_rep1_R1.fastq.gz"
FASTQ_B="$FASTQ_DIR/${CELL}_iCLIP_HNRNPD_${ISOFORM}_rep2_R1.fastq.gz"
STAR_A="$SAMPLE_DIR/replicateA/STAR"
STAR_B="$SAMPLE_DIR/replicateB/STAR"

case "$STAGE" in index|align|all) ;; *) echo 'Invalid STAGE.' >&2; exit 1 ;; esac
command -v STAR >/dev/null 2>&1 || { echo 'Activate STAR first.' >&2; exit 1; }

# Build the shared index once; original sjdbOverhang=167 is preserved.
if [[ "$STAGE" == index || "$STAGE" == all ]]; then
    [[ -s "$FASTA" && -s "$GTF" ]] || { echo 'Check reference paths.' >&2; exit 1; }
    [[ ! -e "$INDEX_DIR" ]] || { echo 'Index exists: select STAGE=align.' >&2; exit 1; }
    mkdir -p "$INDEX_DIR"
STAR --runMode genomeGenerate \
    --runThreadN 15 \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$FASTA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 167
fi

# GEO rep1/rep2 correspond to A/B. Align post-Ultraplex R1 only, as originally.
# The deposited R2 files are not aligned; R1 names must retain UMIs.
if [[ "$STAGE" == align || "$STAGE" == all ]]; then
    command -v zcat >/dev/null 2>&1 || { echo 'Linux-compatible zcat is required.' >&2; exit 1; }
    [[ -s "$INDEX_DIR/Genome" && -s "$INDEX_DIR/SA" && -s "$INDEX_DIR/SAindex" && -s "$INDEX_DIR/genomeParameters.txt" ]] || { echo 'Check STAR index.' >&2; exit 1; }
    [[ -s "$FASTQ_A" && -s "$FASTQ_B" ]] || { echo 'Check both FASTQ paths.' >&2; exit 1; }
    [[ ! -e "$STAR_A" && ! -e "$STAR_B" ]] || { echo 'STAR outputs already exist.' >&2; exit 1; }

# Replicate A.
mkdir -p "$STAR_A"
STAR --outSAMtype BAM SortedByCoordinate \
    --runThreadN 15 \
    --genomeDir "$INDEX_DIR" \
    --readFilesIn "$FASTQ_A" \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --outFilterMultimapNmax 1 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --scoreDelOpen -1 \
    --alignIntronMin 20 \
    --outSJfilterReads Unique \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFileNamePrefix "$STAR_A/" \
    --alignEndsType Extend5pOfRead1

# Replicate B.
mkdir -p "$STAR_B"
STAR --outSAMtype BAM SortedByCoordinate \
    --runThreadN 15 \
    --genomeDir "$INDEX_DIR" \
    --readFilesIn "$FASTQ_B" \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --outFilterMultimapNmax 1 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --scoreDelOpen -1 \
    --alignIntronMin 20 \
    --outSJfilterReads Unique \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFileNamePrefix "$STAR_B/" \
    --alignEndsType Extend5pOfRead1
fi
