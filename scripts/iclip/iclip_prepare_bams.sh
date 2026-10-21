#!/usr/bin/env bash
# Chromosome filtering, per-replicate UMI deduplication, then A/B BAM merging.
# Requires samtools and umi_tools (original analysis versions). Edit settings:
# bash iclip_prepare_bams.sh

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

BAM_A="$SAMPLE_DIR/replicateA/STAR/Aligned.sortedByCoord.out.bam"
BAM_B="$SAMPLE_DIR/replicateB/STAR/Aligned.sortedByCoord.out.bam"
PREP_DIR="$SAMPLE_DIR/prepared"
FILTER_A="$PREP_DIR/${SAMPLE}_repA_aligned.f.bam"
FILTER_B="$PREP_DIR/${SAMPLE}_repB_aligned.f.bam"
DEDUP_A="$PREP_DIR/${SAMPLE}_repA_aligned.f.duprm.bam"
DEDUP_B="$PREP_DIR/${SAMPLE}_repB_aligned.f.duprm.bam"
MERGED_BAM="$PREP_DIR/${SAMPLE}_aligned.f.duprm.merged.bam"

command -v samtools >/dev/null 2>&1 || { echo 'Activate samtools first.' >&2; exit 1; }
command -v umi_tools >/dev/null 2>&1 || { echo 'Activate UMI-tools first.' >&2; exit 1; }
[[ -s "$BAM_A" && -s "$BAM_B" ]] || { echo 'Check both STAR BAMs.' >&2; exit 1; }
[[ ! -e "$PREP_DIR" ]] || { echo 'Prepared outputs already exist.' >&2; exit 1; }
mkdir -p "$PREP_DIR"

# 1. Retain chr1-22/X/Y. No -f 2: these alignments are single-end.

# Replicate A.
samtools index "$BAM_A"
samtools view -hb "$BAM_A" -o "$FILTER_A" \
    chr1:1 chr2:1 chr3:1 chr4:1 chr5:1 chr6:1 chr7:1 chr8:1 chr9:1 chr10:1 chr11:1 chr12:1 chr13:1 chr14:1 chr15:1 chr16:1 chr17:1 chr18:1 chr19:1 chr20:1 chr21:1 chr22:1 chrX:1 chrY:1
samtools index "$FILTER_A"

# Replicate B.
samtools index "$BAM_B"
samtools view -hb "$BAM_B" -o "$FILTER_B" \
    chr1:1 chr2:1 chr3:1 chr4:1 chr5:1 chr6:1 chr7:1 chr8:1 chr9:1 chr10:1 chr11:1 chr12:1 chr13:1 chr14:1 chr15:1 chr16:1 chr17:1 chr18:1 chr19:1 chr20:1 chr21:1 chr22:1 chrX:1 chrY:1
samtools index "$FILTER_B"

# 2. Remove PCR duplicates independently using original read-name UMI settings.
umi_tools dedup \
    -I "$FILTER_A" \
    -L "$PREP_DIR/${SAMPLE}_repA_aligned.f.duprm.log" \
    -S "$DEDUP_A" \
    --extract-umi-method read_id \
    --umi-separator rbc: \
    --method unique

umi_tools dedup \
    -I "$FILTER_B" \
    -L "$PREP_DIR/${SAMPLE}_repB_aligned.f.duprm.log" \
    -S "$DEDUP_B" \
    --extract-umi-method read_id \
    --umi-separator rbc: \
    --method unique

# 3. Merge coordinate-sorted A/B BAMs and index; preserve original -f.
samtools merge -f "$MERGED_BAM" "$DEDUP_A" "$DEDUP_B"
samtools index "$MERGED_BAM"

