# GO/KEGG enrichment of DE genes or parent genes of DE transcripts.
# Edit settings, then run: Rscript long_read_functional_enrichment.R
library(clusterProfiler)
library(org.Hs.eg.db)
library(rtracklayer)
library(openxlsx)

ROOT <- "/path/to/long_read_results"
CELL <- "RPE"                         # EndoC or RPE
LEVEL <- "gene"                       # gene or transcript
COMPARISON <- "EX2_vs_SCR"             # EX2_vs_SCR, EX7_vs_SCR, EX2EX7_vs_SCR
GTF <- "/path/to/gencode.v43.primary_assembly.annotation.gtf"

# Manuscript defaults; historical input-list cutoffs are not recorded in examples.
DE_P_COLUMN <- if (CELL == "RPE") "padj" else "pvalue"
DE_P_CUTOFF <- 0.05
MIN_ABS_LOG2FC <- 0                    # Strict >; no extra magnitude cutoff.

RESULT_FILE <- file.path(ROOT, CELL, paste0("DESeq2_", LEVEL),
                         paste0(CELL, "_", COMPARISON, ".xlsx"))
OUT_DIR <- file.path(ROOT, CELL, "enrichment", LEVEL, COMPARISON,
                     paste0(DE_P_COLUMN, "_", DE_P_CUTOFF, "_log2FC_gt_", MIN_ABS_LOG2FC))
stopifnot(CELL %in% c("EndoC", "RPE"), LEVEL %in% c("gene", "transcript"),
          COMPARISON %in% c("EX2_vs_SCR", "EX7_vs_SCR", "EX2EX7_vs_SCR"),
          DE_P_COLUMN %in% c("pvalue", "padj"),
          is.finite(DE_P_CUTOFF), DE_P_CUTOFF > 0, DE_P_CUTOFF <= 1,
          is.finite(MIN_ABS_LOG2FC), MIN_ABS_LOG2FC >= 0)
if (!file.exists(RESULT_FILE) || !file.exists(GTF)) stop("Check result and GTF paths.")
if (dir.exists(OUT_DIR)) stop("Output directory exists; choose a new OUT_DIR.")

# Read the all-results sheet produced by the preceding DESeq2 script.
res <- read.xlsx(RESULT_FILE, sheet = "all_results")
id_column <- if (LEVEL == "gene") "gene_id" else "transcript_id"
stopifnot(all(c(id_column, DE_P_COLUMN, "log2FoldChange") %in% names(res)))
selected <- res[!is.na(res[[DE_P_COLUMN]]) & res[[DE_P_COLUMN]] < DE_P_CUTOFF &
                !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > MIN_ABS_LOG2FC, , drop = FALSE]

# Map versioned IDs using the same GENCODE release; do not strip ID versions.
annotation <- as.data.frame(import(GTF))
feature <- if (LEVEL == "gene") "gene" else "transcript"
annotation <- annotation[annotation$type == feature, , drop = FALSE]
stopifnot(all(c(id_column, "gene_id", "gene_name") %in% names(annotation)))
if (anyDuplicated(annotation[[id_column]])) stop("Annotation feature IDs are not unique.")
m <- match(selected[[id_column]], annotation[[id_column]])
selected$parent_gene_id <- annotation$gene_id[m]
selected$gene_symbol <- annotation$gene_name[m]

# Preserve the examples' SYMBOL-to-ENTREZ mapping and default database background.
symbol_to_entrez <- as.list(org.Hs.egSYMBOL2EG)
dir.create(OUT_DIR, recursive = TRUE)
write.table(selected, file.path(OUT_DIR, "selected_DE_features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
settings <- c(paste("CELL", CELL), paste("LEVEL", LEVEL), paste("COMPARISON", COMPARISON),
              paste("DE_P_COLUMN", DE_P_COLUMN), paste("DE_P_CUTOFF", DE_P_CUTOFF),
              paste("STRICT_MIN_ABS_LOG2FC", MIN_ABS_LOG2FC), paste("GTF", GTF),
              "Input selection follows manuscript defaults; historical list thresholds unverified.",
              "Enrichment: BH; pvalueCutoff=1; qvalueCutoff=1; no custom universe.",
              "Term exports: nominal pvalue<0.05 exploratory; BH p.adjust<0.05 separately.")
writeLines(settings, file.path(OUT_DIR, "settings.txt"))

for (direction in c("up", "down")) {
    features <- selected[if (direction == "up") selected$log2FoldChange > 0 else
                         selected$log2FoldChange < 0, , drop = FALSE]
    symbols <- unique(features$gene_symbol[!is.na(features$gene_symbol) & nzchar(features$gene_symbol)])
    mapped <- symbol_to_entrez[names(symbol_to_entrez) %in% symbols]
    entrez <- unique(as.character(unlist(mapped, use.names = FALSE)))
    write.table(data.frame(gene_symbol = symbols,
                            mapped_to_entrez = symbols %in% names(mapped)),
                file.path(OUT_DIR, paste0(direction, "_input_genes.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    if (!length(entrez)) {
        writeLines("No mapped input genes; enrichment was not run.",
                   file.path(OUT_DIR, paste0(direction, "_no_mapped_genes.txt")))
        next
    }
    for (category in c("BP", "MF", "CC", "KEGG")) {
        enrichment <- if (category == "KEGG") {
            enrichKEGG(entrez, organism = "hsa", pAdjustMethod = "BH", keyType = "kegg",
                       pvalueCutoff = 1, qvalueCutoff = 1)
        } else {
            enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = category,
                     pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                     keyType = "ENTREZID")
        }
        tab <- as.data.frame(enrichment)
        if (nrow(tab)) {
            tab$geneName <- vapply(strsplit(tab$geneID, "/", fixed = TRUE), function(ids) {
                converted <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
                paste(converted$SYMBOL, collapse = "/")
            }, character(1))
            nominal <- subset(tab, pvalue < 0.05)
            bh <- subset(tab, p.adjust < 0.05)
        } else {
            nominal <- bh <- tab
        }
        write.xlsx(list(all_returned = tab, nominal_p_lt_0.05_exploratory = nominal,
                        BH_padj_lt_0.05 = bh),
                   file.path(OUT_DIR, paste0(direction, "_", category, ".xlsx")))
    }
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
