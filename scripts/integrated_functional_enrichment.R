# GO/KEGG enrichment of one manually selected integrated target-gene list.
# Edit settings, then: Rscript integrated_functional_enrichment.R
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

ROOT <- "/path/to/integration_results"
CELL <- "RPE"                          # EndoC or RPE
GENE_SET <- "p37_p40_p42_shared"        # Label for the supplied list, not a selection rule.
GENE_FILE <- "/path/to/RPE_p37_p40_p42_shared.xlsx"
SHEET <- 1
OUT_DIR <- file.path(ROOT, CELL, GENE_SET)

stopifnot(CELL %in% c("EndoC", "RPE"), length(GENE_SET) == 1,
          nzchar(GENE_SET), !grepl("[/\\\\]", GENE_SET),
          !GENE_SET %in% c(".", ".."))
if (!file.exists(GENE_FILE)) stop("Check GENE_FILE.")
if (dir.exists(OUT_DIR)) stop("Output directory exists; choose a new OUT_DIR.")
input <- read.xlsx(GENE_FILE, sheet = SHEET)
if (!"geneName" %in% names(input)) stop("Input must contain a geneName column of human gene symbols.")
symbols <- unique(as.character(input$geneName))
symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
if (!length(symbols)) stop("No input gene symbols.")

# Preserve the original SYMBOL-to-ENTREZ mapping and database background.
database <- as.list(org.Hs.egSYMBOL2EG)
mapped <- database[names(database) %in% symbols]
entrez <- as.character(unlist(mapped, use.names = FALSE))
if (!length(entrez)) stop("No input symbols map to Entrez IDs.")
dir.create(OUT_DIR, recursive = TRUE)
write.table(data.frame(geneName = symbols, mapped_to_entrez = symbols %in% names(mapped)),
            file.path(OUT_DIR, "input_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("CELL", CELL), paste("GENE_SET", GENE_SET),
             paste("GENE_FILE", normalizePath(GENE_FILE)), paste("SHEET", SHEET),
             "Preselected input list; no DE or fold-change selection is applied here.",
             "Enrichment: BH; pvalueCutoff=1; qvalueCutoff=1; default database background.",
             "Exports: nominal pvalue<0.05 exploratory; BH p.adjust<0.05 separately."),
           file.path(OUT_DIR, "settings.txt"))

for (category in c("BP", "MF", "CC", "KEGG")) {
    enrichment <- if (category == "KEGG") {
        enrichKEGG(entrez, organism = "hsa", pAdjustMethod = "BH", keyType = "kegg",
                   pvalueCutoff = 1, qvalueCutoff = 1)
    } else {
        enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = category, pAdjustMethod = "BH",
                 pvalueCutoff = 1, qvalueCutoff = 1, keyType = "ENTREZID")
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
                    BH_padj_lt_0.05 = bh), file.path(OUT_DIR, paste0(category, ".xlsx")))
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
