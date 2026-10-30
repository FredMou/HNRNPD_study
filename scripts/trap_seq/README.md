# TRAP-seq analysis

Untrimmed FASTQs and processed count matrices: [GSE350756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350756), within GSE350758. Both cell types have two biological replicates per condition. For RPE, two libraries (L3 and L4) were prepared from each biological sample using the same protocol with different starting RNA concentrations. Both libraries were sequenced.

| Script | Step |
|---|---|
| `trap_cutadapt.sh` | Trim adapters and read ends; combine trimmed RPE libraries |
| `trap_prepare_reference.sh` | Combine SILVA 138.2 and Rfam 15.0 rRNA/tRNA sequences |
| `trap_bowtie2.sh` | Build the contaminant index and filter reads |
| `trap_star.sh` | Build the genome index and align filtered reads |
| `trap_salmon.sh` | Extract transcript sequences and quantify STAR transcriptome BAMs |
| `trap_gene_deseq2.R` | Gene-level differential expression from GEO counts |
| `trap_transcript_deseq2.R` | Transcript-level differential expression from GEO counts |
| `trap_functional_enrichment.R` | GO/KEGG enrichment of separate up/down sets |

## Running

Edit the settings at the top of each script. Run Bash scripts with `bash SCRIPT.sh` and R scripts with `Rscript SCRIPT.R`. Use the same absolute `ROOT` across scripts. Select `CELL` (`EndoC` or `RPE`), `CONDITION` (`SCR`, `EX2`, `EX7`, `EX2EX7`) and `REPLICATE` (`1` or `2`) where applicable. For example, `SAMPLE=EndoC_TRAPseq_HNRNPD_SCR_rep1`.

`trap_cutadapt.sh` reproduces the sequencing service’s Cutadapt 4.9 command: adapter `NCTGTCTCTTATACACATCT` for both mates, removal of the first base from each mate (`-u 1 -U 1`), 3′ quality cutoff 22, minimum length 75 nt (EndoC-βH1) or 40 nt (RPE), error rate 0.11, minimum adapter overlap 5 and `--pair-filter=any`, using 8 cores.

Set `RAW_DIR` to the untrimmed GEO FASTQs. EndoC inputs are `SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`; RPE inputs are `SAMPLE_L3_R1.fastq.gz`, `SAMPLE_L3_R2.fastq.gz` and the corresponding L4 files. RPE library files (L3 and L4) are trimmed separately, then concatenated separately for R1 and R2 within each biological sample. The original L3/L4 filename labels are retained. Outputs are `ROOT/CELL/SAMPLE/reads/SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`, ready for `trap_bowtie2.sh`; per-run Cutadapt reports are retained.

## References and shared steps

For `trap_prepare_reference.sh`, place the uncompressed `SILVA_138.2_LSURef_NR99_tax_silva.fasta` and `SILVA_138.2_SSURef_NR99_tax_silva.fasta` from the [SILVA archive](https://www.arb-silva.de/no_cache/download/archive/release_138_2/Exports/) in `REFERENCE_DIR`. The script downloads Rfam families and performs U-to-T conversion. Bowtie2 uses this reference folder.

STAR/Salmon instead share a genome/annotation `REFERENCE_DIR` containing `GRCh38.primary_assembly.genome.fa` and `gencode.v43.primary_assembly.annotation.gtf` from [GENCODE v43](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/). Salmon extracts `gencode.v43.primary_assembly.gffread.transcripts.fa`; no Salmon index is needed.

Build shared references/indexes once. Use `STAGE="all"` initially, then `"filter"` (Bowtie2), `"align"` (STAR) or `"quant"` (Salmon) for subsequent samples. Outputs connect through `ROOT/CELL/SAMPLE/Bowtie2/`, `STAR/` and `Salmon/`. Existing outputs are protected.

## R analyses

Set `GEO_DIR` to the downloaded matrices and select `CELL`. DESeq2 inputs are `CELL_TRAPseq_HNRNPD_gene_counts.txt` (`Gene_id`) or `CELL_TRAPseq_HNRNPD_transcript_counts.txt` (`Transcript_id`), followed by eight GEO sample columns. These matrices already contain the protein-coding/lncRNA selection. Original integer conversion, ≥10 counts in ≥2 samples filtering and treatment-versus-SCR tests are retained.

Workbooks under `ROOT/CELL/DESeq2_gene/` or `DESeq2_transcript/` contain all results and BH-adjusted p < 0.05 results. EndoC gene-level workbooks also include nominal p < 0.05 exploratory results, following the manuscript. Fits, sample assignments and session information are saved.

Enrichment reads these workbooks. Set `LEVEL`, `COMPARISON` and `GTF`. Input selection defaults to nominal p < 0.05 for EndoC gene-level sets and BH-adjusted p < 0.05 otherwise. The editable strict `MIN_ABS_LOG2FC=0.585` is suggested by example filenames; confirm the threshold for the relevant analysis. Transcript sets use unique parent genes. GO BP/MF/CC and KEGG retain the original default database background. Tables contain nominal p < 0.05 exploratory terms and BH-adjusted p < 0.05 terms separately, plus selected inputs and settings. Database updates can change enrichment results.

## Software

Cutadapt 4.9; Bowtie2 2.5.4; STAR 2.7.11b; Salmon 1.10.3; gffread; wget/gzip/awk/sed. R: DESeq2 (methods: 1.44.0), openxlsx; enrichment additionally uses clusterProfiler, org.Hs.eg.db and rtracklayer. No plotting packages are needed.


