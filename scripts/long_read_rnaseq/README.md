# Long-read RNA-seq analysis

Analysis scripts for EndoC-βH1 and RPE: four conditions (SCR, EX2, EX7, EX2EX7), each with three biological replicates.

## Inputs

- FASTQs and processed gene/transcript raw-count matrices: [GSE350755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350755), within [GSE350758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350758).
- GENCODE v43 primary assembly: [genome FASTA](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz) and [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz). Decompress before use.
- Exact alignment junction reference: `data/junctions.GENECODEv43.bed.gz`. Decompress once with `gunzip -k data/junctions.GENECODEv43.bed.gz`.

## Run order

| Script | Purpose |
| --- | --- |
| `long_read_prepare_reads.sh` | NanoPlot (EndoC only), chopper and Pychopper |
| `long_read_align.sh` | Minimap2 alignment and coordinate sorting |
| `long_read_isoquant.sh` | Joint quantification of all 12 samples per cell type |
| `long_read_gene_deseq2.R` | Gene-level differential expression from GEO counts |
| `long_read_transcript_deseq2.R` | Transcript-level differential expression from GEO counts |
| `long_read_functional_enrichment.R` | Separate up/down GO BP/MF/CC and KEGG enrichment |

Edit paths and use the same `ROOT` and `CELL` across scripts. Preparation/alignment process one `CONDITION` and `REPLICATE` per run; repeat for all 12 samples before IsoQuant. Input names match GEO; `sample_metadata.tsv` is a reference table, not a required script input. Outputs are under `ROOT/CELL/`; existing stage directories are not overwritten.

Run Bash files with `bash SCRIPT.sh` and R files with `Rscript SCRIPT.R`.

## Analysis settings

Read preparation uses NanoPlot for EndoC only, chopper (Q7, minimum 500 nt) and Pychopper. Minimap2 uses the included junction reference. IsoQuant quantifies all 12 samples together per cell type. Its reference gene/transcript counts were used downstream; transcript-discovery outputs were not used.

DESeq2 reads the processed `CELL_LongReadRNAseq_HNRNPD_gene_counts.txt` or `CELL_LongReadRNAseq_HNRNPD_transcript_counts.txt` from GEO. These matrices already contain the original protein-coding/lncRNA selection and duplicate handling; do not substitute unfiltered IsoQuant tables. Sample assignments come from column names. Original integer conversion, ≥10 counts in ≥3 samples, `~ passage` design, three removal-versus-SCR contrasts and DESeq2 defaults are retained. Workbooks contain all results and BH-adjusted p < 0.05 results; EndoC additionally exports nominal p < 0.05 exploratory sets. Selected sets retain `abs(log2FoldChange) > 0`. Results, fitted objects and session information are saved separately in `DESeq2_gene/` and `DESeq2_transcript/`.

Enrichment reads the preceding workbook's `all_results` sheet. Select `LEVEL` and `COMPARISON`, and supply the GTF for gene/parent-gene mapping. Input thresholds are editable; manuscript defaults are RPE `padj < 0.05`, EndoC `pvalue < 0.05`, and `abs(log2FoldChange) > 0`. Historical input-list thresholds were not recorded in the supplied enrichment examples. Parent genes are deduplicated within each direction; a gene with opposing transcript responses may enter both sets. Original enrichment settings are BH adjustment, `pvalueCutoff=1`, `qvalueCutoff=1` and the default database background. Workbooks distinguish all returned terms, nominal enrichment p < 0.05 exploratory terms and BH-adjusted p < 0.05 terms. Selected inputs, mapping status, settings and session information are saved under `enrichment/`. Database updates may change results. No plotting is included.

## Software

NanoPlot 1.42.0, chopper 0.9.2, Pychopper 2.7.10, minimap2 2.29, samtools, IsoQuant 3.7.0, Bash and gzip/gunzip. R analyses require DESeq2 1.44.0, openxlsx, clusterProfiler, org.Hs.eg.db and rtracklayer. KEGG requires internet access.
