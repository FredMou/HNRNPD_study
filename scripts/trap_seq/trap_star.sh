#!/usr/bin/env bash
# Build the shared STAR index and align one Bowtie2-filtered TRAP-seq sample.
# Edit settings, activate STAR, then: bash trap_star.sh
set -euo pipefail

ROOT="/path/to/TRAP_results"
REFERENCE_DIR="/path/to/GRCh38_GENCODEv43"
CELL="EndoC"                  # EndoC or RPE
CONDITION="SCR"               # SCR, EX2, EX7 or EX2EX7
REPLICATE="1"                 # 1 or 2
STAGE="all"                   # index, align or all; build the index once.
GENOME_FASTA="$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa"
GTF="$REFERENCE_DIR/gencode.v43.primary_assembly.annotation.gtf"
GENOME_DIR="$ROOT/STAR_index/genome_dir"

case "$CELL" in EndoC|RPE) ;; *) echo 'Invalid CELL.' >&2; exit 1 ;; esac
case "$CONDITION" in SCR|EX2|EX7|EX2EX7) ;; *) echo 'Invalid CONDITION.' >&2; exit 1 ;; esac
case "$REPLICATE" in 1|2) ;; *) echo 'Invalid REPLICATE.' >&2; exit 1 ;; esac
case "$STAGE" in index|align|all) ;; *) echo 'Invalid STAGE.' >&2; exit 1 ;; esac
command -v STAR >/dev/null 2>&1 || { echo 'Activate STAR first.' >&2; exit 1; }
SAMPLE="${CELL}_TRAPseq_HNRNPD_${CONDITION}_rep${REPLICATE}"
READ_DIR="$ROOT/$CELL/$SAMPLE/Bowtie2"
R1="$READ_DIR/${SAMPLE}_rm_rtRNA_1.fq.gz"
R2="$READ_DIR/${SAMPLE}_rm_rtRNA_2.fq.gz"
OUT_DIR="$ROOT/$CELL/$SAMPLE/STAR"

# Original index: genome only; supply GTF during alignment below.
if [[ "$STAGE" == index || "$STAGE" == all ]]; then
    [[ -s "$GENOME_FASTA" ]] || { echo 'Check GENOME_FASTA.' >&2; exit 1; }
    [[ ! -e "$GENOME_DIR" ]] || { echo 'Index directory exists; select STAGE=align.' >&2; exit 1; }
    mkdir -p "$GENOME_DIR"
    ( cd "$GENOME_DIR"
      STAR --runThreadN 12 \
          --runMode genomeGenerate \
          --genomeDir "$GENOME_DIR" \
          --genomeFastaFiles "$GENOME_FASTA"
    )
fi

if [[ "$STAGE" == align || "$STAGE" == all ]]; then
    [[ -s "$R1" && -s "$R2" && -s "$GTF" ]] || { echo 'Check filtered FASTQs and GTF.' >&2; exit 1; }
    [[ -s "$GENOME_DIR/Genome" && -s "$GENOME_DIR/SA" && -s "$GENOME_DIR/SAindex" ]] || { echo 'STAR index is missing or incomplete.' >&2; exit 1; }
    [[ ! -e "$OUT_DIR" ]] || { echo 'STAR output directory already exists.' >&2; exit 1; }
    mkdir -p "$OUT_DIR"
    cd "$OUT_DIR"
    STAR --runThreadN 13 \
        --outFilterType BySJout \
        --outFilterMismatchNmax 2 \
        --outFilterMultimapNmax 1 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --quantMode TranscriptomeSAM GeneCounts \
        --sjdbGTFfile "$GTF"
    # Original read-group line was disconnected from STAR; retained as a comment.
    # --outSAMattrRGline ID:sample1_ScrA LB:ribo_seq PL:ILLUMINA SM:sample1_ScrA
fi
