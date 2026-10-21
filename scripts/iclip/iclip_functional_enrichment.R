# GO Biological Process over-representation analysis of HNRNPD target genes.
# Edit settings, then run: Rscript iclip_functional_enrichment.R
# Figure 2 reports BP; other original categories are shown as inactive examples.

ROOT <- "/path/to/iclip_results"  # Same ROOT as the preceding iCLIP scripts.
CELL <- "EndoC"                  # EndoC or RPE
ISOFORM <- "p37"                 # p37, p40, p42 or p45
file_label <- if (CELL == "RPE") "HRPE" else "EndoC"
sample <- paste(file_label, ISOFORM, sep = "_")
sample_dir <- file.path(ROOT, CELL, paste0("PureCLIP_", ISOFORM, "_FWD"))
GENE_FILE <- file.path(sample_dir, "BindingSiteFinder", paste0(sample, "_results_v43.txt"))

stopifnot(CELL %in% c("EndoC", "RPE"), ISOFORM %in% c("p37", "p40", "p42", "p45"))
output_dir <- file.path(ROOT, CELL, paste0("PureCLIP_", ISOFORM, "_FWD"), "Functional_enrichment")
# library(tidyr)  # Unused in the analysis; no explicit attachment is needed.
# library(readxl) # No longer needed: read the preceding script's TSV directly.
# Required functions use namespaces: clusterProfiler, org.Hs.eg.db and openxlsx.
stopifnot(
    requireNamespace("clusterProfiler", quietly = TRUE),
    requireNamespace("org.Hs.eg.db", quietly = TRUE),
    requireNamespace("openxlsx", quietly = TRUE),
    file.exists(GENE_FILE)
)
if (dir.exists(output_dir)) stop("Output directory exists; retain it or use a new ROOT.")

# 1. Map symbols to Entrez IDs using the original SYMBOL2EG mapping.
gene_table <- read.delim(GENE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
if (!"geneName" %in% names(gene_table)) stop("Input needs a geneName column.")
gene_list <- unique(gene_table$geneName)
gene_list <- gene_list[!is.na(gene_list) & nzchar(gene_list)]
if (!is.character(gene_list)) stop("geneName must contain gene symbols as text.")
database <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG)
gene_list_enrich <- database[names(database) %in% gene_list]
entrez_ids <- as.character(unlist(gene_list_enrich))
if (!length(entrez_ids)) stop("No gene symbols mapped to Entrez IDs.")

# 2. Original settings: BH; permissive internal cutoffs; no custom universe.
# Significance is applied separately below as p.adjust < 0.05.
org_db <- org.Hs.eg.db::org.Hs.eg.db
edo_bp <- clusterProfiler::enrichGO(
    entrez_ids, OrgDb = org_db, ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, keyType = "ENTREZID"
)
# Optional categories from the original scripts; not run for Figure 2.
# edo_mf <- clusterProfiler::enrichGO(
#     entrez_ids, OrgDb = org_db, ont = "MF", pAdjustMethod = "BH",
#     pvalueCutoff = 1, qvalueCutoff = 1, keyType = "ENTREZID"
# )
# edo_cc <- clusterProfiler::enrichGO(
#     entrez_ids, OrgDb = org_db, ont = "CC", pAdjustMethod = "BH",
#     pvalueCutoff = 1, qvalueCutoff = 1, keyType = "ENTREZID"
# )
# # KEGG uses online annotations by default; database updates may affect reruns.
# edo_kegg <- clusterProfiler::enrichKEGG(
#     entrez_ids, organism = "hsa", pAdjustMethod = "BH", keyType = "kegg",
#     pvalueCutoff = 1, qvalueCutoff = 1
# )

# 3. Export all returned terms, original BH-significant terms, and nominal terms.
# Nominal-only terms are exploratory, not BH-significant. No new tests are added.
id_to_gene <- function(ids) {
    mapped <- clusterProfiler::bitr(
        strsplit(ids, "/", fixed = TRUE)[[1]],
        fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db"
    )
    paste(mapped$SYMBOL, collapse = "/")
}
all_terms <- as.data.frame(edo_bp)
if (nrow(all_terms)) {
    all_terms$geneName <- vapply(all_terms$geneID, id_to_gene, character(1))
} else {
    all_terms$geneName <- character(0)
}
bh_terms <- if (nrow(all_terms)) subset(all_terms, p.adjust < 0.05) else all_terms
nominal_terms <- if (nrow(all_terms)) subset(all_terms, pvalue < 0.05) else all_terms

dir.create(output_dir, recursive = TRUE)
openxlsx::write.xlsx(
    list(All_returned_terms = all_terms, BH_adjusted_p_lt_0.05 = bh_terms,
         Nominal_p_lt_0.05_exploratory = nominal_terms),
    file.path(output_dir, paste0(ISOFORM, "_PureCLIP_", CELL, ".GO_BP.xlsx"))
)
write.table(data.frame(geneName = gene_list), file.path(output_dir, paste0(ISOFORM, "_GeneList_", CELL, ".tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save analysis objects, mapping and versions for traceability.
saveRDS(list(BP = edo_bp),
        file.path(output_dir, "enrichment_objects.rds"))
write.table(data.frame(SYMBOL = rep(names(gene_list_enrich), lengths(gene_list_enrich)),
                       ENTREZID = entrez_ids),
            file.path(output_dir, "symbol_entrez_mapping.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(setdiff(unique(gene_list[!is.na(gene_list)]), names(gene_list_enrich)),
           file.path(output_dir, "unmapped_symbols.txt"))
writeLines(c(paste("Run date:", Sys.Date()), paste("Input:", GENE_FILE),
             "No custom universe; BH correction separately within each analysis.",
             capture.output(sessionInfo())), file.path(output_dir, "sessionInfo.txt"))
