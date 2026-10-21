#!/usr/bin/env bash
# PureCLIP crosslink-site and binding-region calling from merged A/B BAMs.
# Requires pureclip (original analysis version). Edit settings, then:
# bash iclip_pureclip.sh

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

FASTA="/path/to/GRCh38_GENCODEv43/GRCh38.primary_assembly.genome.fa"
MERGED_BAM="$SAMPLE_DIR/prepared/${SAMPLE}_aligned.f.duprm.merged.bam"
RESULT_DIR="$SAMPLE_DIR/PureCLIP"

command -v pureclip >/dev/null 2>&1 || { echo 'Activate PureCLIP first.' >&2; exit 1; }
[[ -s "$FASTA" && -s "$MERGED_BAM" && -s "$MERGED_BAM.bai" ]] || { echo 'Check FASTA, merged BAM and BAI paths.' >&2; exit 1; }
[[ ! -e "$RESULT_DIR" ]] || { echo 'PureCLIP outputs already exist.' >&2; exit 1; }
mkdir -p "$RESULT_DIR"

# Preserve long-double precision and 15 threads; output sites and regions.
pureclip -i "$MERGED_BAM" \
    -bai "$MERGED_BAM.bai" \
    -g "$FASTA" \
    -ld \
    -nt 15 \
    -o "$RESULT_DIR/PureCLIP.${SAMPLE}_crosslink_sites.bed" \
    -or "$RESULT_DIR/PureCLIP.${SAMPLE}_crosslink_regions.bed"

