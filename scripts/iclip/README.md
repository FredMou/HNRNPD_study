# iCLIP-seq analysis

Main analysis workflow for EndoC-βH1 and RPE cells, covering HNRNPD p37, p40, p42 and p45.

## Inputs

- Post-Ultraplex FASTQs and replicate/strand-specific BigWigs: [GSE350757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350757), part of [GSE350758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350758).
- Unmodified GENCODE v43 primary-assembly references: [genome FASTA](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz), [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz) and [GFF3](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gff3.gz). Decompress before use.

Download the inputs to local folders; references are not duplicated here.

## Run order

| Script | Analysis |
|---|---|
| `iclip_star.sh` | STAR index generation and single-end alignment |
| `iclip_prepare_bams.sh` | Main-chromosome filtering, UMI deduplication and replicate merging/indexing |
| `iclip_pureclip.sh` | Crosslink-site and binding-region calling |
| `iclip_define_binding_sites.R` | Binding-site definition and annotation |
| `iclip_functional_enrichment.R` | GO Biological Process enrichment |

Run the scripts in the order above. Edit absolute paths and use the same `ROOT`, `CELL` (`EndoC` or `RPE`) and `ISOFORM` (`p37`, `p40`, `p42` or `p45`) across scripts. Each run processes one cell/isoform, including replicates A/B. Build the STAR index once with `STAGE="all"`; use `STAGE="align"` for subsequent samples.

```bash
bash iclip_star.sh
bash iclip_prepare_bams.sh
bash iclip_pureclip.sh
Rscript iclip_define_binding_sites.R
Rscript iclip_functional_enrichment.R
```

Set `FASTQ_DIR` to the downloaded GEO FASTQs and `BW_DIR` to the four BigWigs for the selected cell/isoform. Names are read directly as deposited, for example:

- `EndoC_iCLIP_HNRNPD_p37_rep1_R1.fastq.gz`
- `EndoC_iCLIP_HNRNPD_p37_rep1_plus.bw` and `EndoC_iCLIP_HNRNPD_p37_rep1_minus.bw`

Use the matching `rep2` files for the second replicate; RPE inputs use the `RPE` prefix. GEO replicates 1/2 correspond to A/B internally. STAR aligns only post-Ultraplex R1, as in the original analysis; deposited R2 files are not used for alignment. FASTQ read names must retain UMIs. Existing output names (including `HRPE` for RPE) are unchanged, so downstream scripts remain connected. Enrichment reads the preceding binding-site table directly.

Requires STAR 2.7.11a, samtools, UMI-tools, PureCLIP, Linux-compatible zcat, R, BindingSiteFinder 2.2.0, clusterProfiler 4.12.6 and the R dependencies listed in the scripts. Excel export requires xlsx/rJava and Java. R scripts record session information.

Binding-site settings retain 5-nt sites, global cutoff 0.01, gene-wise cutoff 0.2 and default frequency-based annotation. GO BP results include separate all-returned, BH-adjusted p < 0.05 and nominal-p exploratory sheets; MF/CC/KEGG examples are commented out. Outputs are sample-specific; existing stage directories are not overwritten.
