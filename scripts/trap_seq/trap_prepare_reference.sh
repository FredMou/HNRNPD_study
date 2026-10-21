#!/usr/bin/env bash
# Prepare the shared rRNA/tRNA reference. Requires wget, gzip, awk and sed.
# Edit REFERENCE_DIR, place the two SILVA FASTAs there, then run this script.
set -euo pipefail

REFERENCE_DIR="/path/to/TRAP_reference"
RFAM_URL="https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0"
mkdir -p "$REFERENCE_DIR"
cd "$REFERENCE_DIR"

# SILVA 138.2 inputs: use these exact release files, uncompressed.
LSU="SILVA_138.2_LSURef_NR99_tax_silva.fasta"
SSU="SILVA_138.2_SSURef_NR99_tax_silva.fasta"
[[ -s "$LSU" && -s "$SSU" ]] || { echo 'Place both SILVA 138.2 FASTAs in REFERENCE_DIR first.' >&2; exit 1; }
[[ ! -e rRNA_tRNA_combined.fasta ]] || { echo 'Combined reference already exists.' >&2; exit 1; }

# Select the original Rfam 15.0 rRNA and tRNA families (all species).
mkdir -p Rfam_15.0/fasta_files
wget -O Rfam_15.0/family.txt.gz "$RFAM_URL/database_files/family.txt.gz"
gzip -dc Rfam_15.0/family.txt.gz > Rfam_15.0/family.txt
cut -f 1,2,19 Rfam_15.0/family.txt > Rfam_15.0/family_selected_columns.tsv
awk '$0 ~ /Gene; rRNA/' Rfam_15.0/family_selected_columns.tsv > Rfam_15.0/family_rRNA.tsv
awk '$0 ~ /Gene; tRNA/' Rfam_15.0/family_selected_columns.tsv > Rfam_15.0/family_tRNA.tsv

for TYPE in rRNA tRNA; do
    [[ -s "Rfam_15.0/family_${TYPE}.tsv" ]] || { echo "No $TYPE families selected." >&2; exit 1; }
    : > "Rfam_${TYPE}.fasta"
    while IFS= read -r FAMILY; do
        wget -O "Rfam_15.0/fasta_files/${FAMILY}.fa.gz" "$RFAM_URL/fasta_files/${FAMILY}.fa.gz"
        gzip -dc "Rfam_15.0/fasta_files/${FAMILY}.fa.gz" | sed "s/^>/>${FAMILY}:/" >> "Rfam_${TYPE}.fasta"
    done < <(cut -f 1 "Rfam_15.0/family_${TYPE}.tsv")
done

# Original order and U-to-T conversion; no additional sequence filtering.
cat "$LSU" "$SSU" > silva_rRNA_NR99.fa
sed 's/U/T/g' silva_rRNA_NR99.fa > silva_rRNA_NR99_DNA.fa
sed 's/U/T/g' Rfam_rRNA.fasta > Rfam_rRNA_DNA.fasta
sed 's/U/T/g' Rfam_tRNA.fasta > Rfam_tRNA_DNA.fasta
cat silva_rRNA_NR99_DNA.fa Rfam_rRNA_DNA.fasta Rfam_tRNA_DNA.fasta > rRNA_tRNA_combined.fasta.tmp
mv rRNA_tRNA_combined.fasta.tmp rRNA_tRNA_combined.fasta
