# Define HNRNPD binding sites using PureCLIP sites and replicate BigWigs.
# Edit settings, then run: Rscript iclip_define_binding_sites.R

ROOT <- "/path/to/iclip_results"  # Same ROOT as the Bash scripts.
CELL <- "EndoC"                  # EndoC or RPE
ISOFORM <- "p37"                 # p37, p40, p42 or p45
ANNOTATION <- "/path/to/GRCh38_GENCODEv43/gencode.v43.primary_assembly.annotation.gff3"
BW_DIR <- "/path/to/selected_isoform/bw_files"
RUN_WIDTH_DIAGNOSTICS <- FALSE   # Optional original chr1 width-estimation check.

stopifnot(CELL %in% c("EndoC", "RPE"), ISOFORM %in% c("p37", "p40", "p42", "p45"))
# Internal RPE output filenames and dataset names retain the original HRPE label.
file_label <- if (CELL == "RPE") "HRPE" else "EndoC"
sample <- paste(file_label, ISOFORM, sep = "_")
sample_dir <- file.path(ROOT, CELL, paste0("PureCLIP_", ISOFORM, "_FWD"))
peak_file <- file.path(sample_dir, "PureCLIP", paste0("PureCLIP.", sample, "_crosslink_sites.bed"))
output_dir <- file.path(sample_dir, "BindingSiteFinder")

# Required namespaces: rtracklayer, GenomicFeatures, GenomicRanges, GenomeInfoDb,
# S4Vectors and BindingSiteFinder. xlsx/rJava are used indirectly for Excel export.
# library(ggplot2)       # Not used directly in this script.
# library(AnnotationDbi) # No direct calls; dependencies load it as needed.
# library(dplyr)         # Not used directly; BindingSiteFinder manages its imports.
# library(reshape2)      # Unused.
# library(UpSetR)        # Unused.
# library(rJava)         # Dependency of xlsx; explicit attachment is unnecessary.
# library(xlsx)          # Called internally by exportTargetGenes().
# library(txdbmaker)     # Only needed for the constructor fallback below.
stopifnot(
    requireNamespace("BindingSiteFinder", quietly = TRUE),
    requireNamespace("rtracklayer", quietly = TRUE),
    requireNamespace("GenomicFeatures", quietly = TRUE),
    requireNamespace("GenomicRanges", quietly = TRUE),
    requireNamespace("GenomeInfoDb", quietly = TRUE),
    requireNamespace("S4Vectors", quietly = TRUE),
    requireNamespace("xlsx", quietly = TRUE)
)
stopifnot(file.exists(ANNOTATION), file.exists(peak_file))
if (dir.exists(output_dir)) stop("Output directory exists; retain it or use a new ROOT.")

# 1. Build GENCODE v43 gene and transcript-region annotations.
# Newer Bioconductor versions moved this constructor to txdbmaker.
if ("makeTxDbFromGFF" %in% getNamespaceExports("GenomicFeatures")) {
    anno_db <- GenomicFeatures::makeTxDbFromGFF(file = ANNOTATION, format = "gff3")
} else {
    stopifnot(requireNamespace("txdbmaker", quietly = TRUE))
    anno_db <- txdbmaker::makeTxDbFromGFF(file = ANNOTATION, format = "gff3")
}
anno_info <- rtracklayer::import(ANNOTATION, format = "gff3")
genes <- GenomicFeatures::genes(anno_db)
idx <- match(genes$gene_id, anno_info$gene_id)
valid <- !is.na(idx)
message("Genes without matching GFF3 metadata: ", sum(!valid))
meta <- cbind(S4Vectors::mcols(genes)[valid, ], S4Vectors::mcols(anno_info)[idx[valid], ])
meta <- meta[, !duplicated(colnames(meta))]
idx2 <- match(genes$gene_id, meta$gene_id)
valid2 <- !is.na(idx2)
genes <- genes[valid2, ]
S4Vectors::mcols(genes) <- meta[idx2[valid2], ]
regions <- GenomicRanges::GRangesList(
    CDS = GenomicFeatures::cds(anno_db),
    Intron = unlist(GenomicFeatures::intronsByTranscript(anno_db)),
    UTR3 = unlist(GenomicFeatures::threeUTRsByTranscript(anno_db)),
    UTR5 = unlist(GenomicFeatures::fiveUTRsByTranscript(anno_db))
)

# 2. Import sites and pair GEO rep1/rep2 BigWigs in the original A/B order.
# PureCLIP normally writes BED6; retain V7 when an input has a seventh column.
first_record <- readLines(peak_file, n = 1L, warn = FALSE)
n_columns <- length(strsplit(first_record, "\t", fixed = TRUE)[[1L]])
if (!n_columns %in% c(6L, 7L)) stop("Expected a nonempty PureCLIP BED6 or BED7 file.")
if (n_columns == 7L) {
    peaks <- rtracklayer::import(peak_file, format = "BED", extraCols = c(V7 = "character"))
} else {
    peaks <- rtracklayer::import(peak_file, format = "BED")
}
peaks <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")
track_stem <- paste0(CELL, "_iCLIP_HNRNPD_", ISOFORM, "_rep", 1:2, "_")
plus_tracks <- file.path(BW_DIR, paste0(track_stem, "plus.bw"))
minus_tracks <- file.path(BW_DIR, paste0(track_stem, "minus.bw"))
if (!all(file.exists(c(plus_tracks, minus_tracks)))) {
    stop("Missing GEO rep1/rep2 plus/minus BigWigs; check BW_DIR.")
}
col_data <- data.frame(
    name = paste(file_label, ISOFORM, sep = "."), id = 1:2,
    condition = factor(rep(ISOFORM, 2), levels = ISOFORM),
    clPlus = plus_tracks, clMinus = minus_tracks
)
bds_initial <- BindingSiteFinder::BSFDataSetFromBigWig(ranges = peaks, meta = col_data)

# 3. Preserve the active final settings in both supplied cell-type scripts.
# All eight final calls use these values; do not infer different cutoffs from comments.
if (RUN_WIDTH_DIAGNOSTICS) {
    width_check <- BindingSiteFinder::BSFind(
        object = bds_initial, anno.genes = genes, anno.transcriptRegionList = regions,
        est.subsetChromosome = "chr1", veryQuiet = TRUE, est.maxBsWidth = 29
    )
    print(BindingSiteFinder::processingStepsFlowChart(width_check))
    print(BindingSiteFinder::estimateBsWidthPlot(width_check))
}
bds <- BindingSiteFinder::BSFind(
    object = bds_initial, bsSize = 5,
    anno.genes = genes, anno.transcriptRegionList = regions,
    veryQuiet = TRUE, est.maxBsWidth = 29,
    cutoff.geneWiseFilter = 0.2, cutoff.globalFilter = 0.01
)

# 4. Export to a separate directory for each isoform.
dir.create(output_dir, recursive = TRUE)
results <- BindingSiteFinder::getRanges(bds)
write.table(as.data.frame(results), file.path(output_dir, paste0(sample, "_results_v43.txt")),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
BindingSiteFinder::exportToBED(bds, con = file.path(output_dir, paste0(sample, "_binding_sites.bed")))
# "xslx" is the original spelling accepted by this package's export API.
BindingSiteFinder::exportTargetGenes(bds, path = paste0(output_dir, "/"),
                                    format = "xslx", split = "none", sort = "bs")
BindingSiteFinder::exportTargetGenes(bds, path = file.path(output_dir, "SplitByTransRegion_"),
                                    format = "xslx", split = "transcriptRegion")
saveRDS(bds, file.path(output_dir, paste0(sample, "_BSFDataSet.rds")))
# Retain the original object naming convention for downstream load() workflows.
object_name <- paste0("bds", ISOFORM, "_", CELL)
assign(object_name, bds)
save(list = object_name, file = file.path(output_dir, paste0("BDSobject_", sample, ".RData")))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

