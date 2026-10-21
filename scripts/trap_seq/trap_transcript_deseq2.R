# Transcript-level DESeq2 analysis of the processed raw-count matrices from GSE350756.
# Edit settings, then run: Rscript trap_transcript_deseq2.R
# Methods version: DESeq2 1.44.0. No plotting or gene-level analysis.

library(DESeq2)
library(openxlsx)

ROOT <- "/path/to/TRAP_results"
GEO_DIR <- "/path/to/downloaded_GEO_count_matrices"
CELL <- "EndoC"                       # EndoC or RPE
COUNT_FILE <- file.path(GEO_DIR, paste0(CELL, "_TRAPseq_HNRNPD_transcript_counts.txt"))
OUT_DIR <- file.path(ROOT, CELL, "DESeq2_transcript")

stopifnot(CELL %in% c("EndoC", "RPE"))
if (!file.exists(COUNT_FILE)) stop("Count matrix not found: ", COUNT_FILE)
if (dir.exists(OUT_DIR)) stop("Output directory exists; choose a new OUT_DIR.")

# GEO matrices already include the original protein-coding/lncRNA selection.
input <- read.delim(COUNT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
if (names(input)[1] != "Transcript_id") stop("Expected Transcript_id as the first column.")
if (anyNA(input$Transcript_id) || any(!nzchar(input$Transcript_id)) || anyDuplicated(input$Transcript_id)) {
    stop("Transcript identifiers must be non-missing and unique.")
}
rownames(input) <- input$Transcript_id
raw_counts <- as.matrix(input[, -1, drop = FALSE])
if (!is.numeric(raw_counts) || anyNA(raw_counts) ||
    any(!is.finite(raw_counts)) || any(raw_counts < 0) ||
    any(raw_counts > .Machine$integer.max)) stop("Check the raw count values.")
if (any(raw_counts != floor(raw_counts))) {
    warning("Fractional counts will be truncated, matching the original integer conversion.")
}
mode(raw_counts) <- "integer"

# Reconstruct the confirmed sample assignments from the GEO column names.
conditions <- rep(c("SCR", "EX2", "EX7", "EX2EX7"), each = 2)
replicates <- rep(1:2, times = 4)
sample_ids <- paste0(CELL, "_TRAPseq_HNRNPD_", conditions, "_rep", replicates)
if (ncol(raw_counts) != 8 || anyDuplicated(colnames(raw_counts)) ||
    !setequal(colnames(raw_counts), sample_ids)) stop("Expected the 8 GEO sample columns for CELL.")
raw_counts <- raw_counts[, sample_ids, drop = FALSE]
passage_levels <- c("Scramble", "SSS005", "SSS006", "SSS005_SSS006")
col_data <- data.frame(
    passage = factor(rep(passage_levels, each = 2), levels = passage_levels),
    row.names = sample_ids
)

# Original model and prefilter; all 8 samples are fitted together.
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = col_data,
                              design = ~ passage)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
if (nrow(dds) == 0) stop("No transcripts pass the original count filter.")
dds <- DESeq(dds)

dir.create(OUT_DIR, recursive = TRUE)
write.table(data.frame(sample_id = sample_ids, condition = conditions,
                       biological_replicate = replicates, passage = col_data$passage),
            file.path(OUT_DIR, "sample_assignments.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# Retain results() defaults, including independent filtering and BH adjustment.
contrast_names <- c(SSS005 = "EX2_vs_SCR", SSS006 = "EX7_vs_SCR",
                    SSS005_SSS006 = "EX2EX7_vs_SCR")
for (treatment in names(contrast_names)) {
    res <- as.data.frame(results(dds, contrast = c("passage", treatment, "Scramble")))
    res <- res[order(res$pvalue), , drop = FALSE]
    res <- data.frame(transcript_id = rownames(res), res, row.names = NULL)

    # Preserve the original strict abs(log2FoldChange) > 0 selection.
    bh <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0)
    tables <- list(all_results = res, BH_padj_lt_0.05 = bh)
    write.xlsx(tables, file.path(OUT_DIR, paste0(CELL, "_", contrast_names[[treatment]], ".xlsx")))
}

saveRDS(dds, file.path(OUT_DIR, "transcript_DESeq2_fit.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
