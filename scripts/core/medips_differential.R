#!/usr/bin/env Rscript
# MEDIPS Differential Analysis Script
# Compare two groups (e.g., AEG vs CTRL) for differential methylation

suppressPackageStartupMessages({
  library(MEDIPS)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(optparse)
  library(data.table)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--bam-list-1"), type="character", help="Comma-separated BAM files for group 1 (e.g., AEG)"),
  make_option(c("--bam-list-2"), type="character", help="Comma-separated BAM files for group 2 (e.g., CTRL)"),
  make_option(c("--group1-name"), type="character", default="Group1", help="Name for group 1"),
  make_option(c("--group2-name"), type="character", default="Group2", help="Name for group 2"),
  make_option(c("-o", "--output"), type="character", help="Output prefix"),
  make_option(c("-g", "--genome"), type="character", default="hg38", help="Genome assembly [default: %default]"),
  make_option(c("-w", "--window"), type="integer", default=300, help="Window size [default: %default]"),
  make_option(c("-e", "--extend"), type="integer", default=0, help="Extend reads [default: %default]"),
  make_option(c("-u", "--uniq"), type="integer", default=1, help="Uniqueness threshold [default: %default]"),
  make_option(c("-s", "--shift"), type="integer", default=0, help="Shift reads [default: %default]"),
  make_option(c("-p", "--paired"), action="store_true", default=TRUE, help="Paired-end data"),
  make_option(c("--chr_select"), type="character", default=NULL, help="Chromosomes to analyze (comma-separated)"),
  make_option(c("--p-adj"), type="character", default="BH", help="P-value adjustment method [default: %default]"),
  make_option(c("--diff-method"), type="character", default="edgeR", help="Differential method: edgeR or ttest [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$`bam-list-1`) || is.null(opt$`bam-list-2`) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("BAM lists for both groups and output prefix are required")
}

# Parse BAM file lists
bams_group1 <- strsplit(opt$`bam-list-1`, ",")[[1]]
bams_group2 <- strsplit(opt$`bam-list-2`, ",")[[1]]

cat("=== MEDIPS Differential Analysis ===\n")
cat("Group 1 (", opt$`group1-name`, "):", length(bams_group1), "samples\n")
cat("Group 2 (", opt$`group2-name`, "):", length(bams_group2), "samples\n")
cat("Output prefix:", opt$output, "\n")
cat("Genome:", opt$genome, "\n")
cat("Window size:", opt$window, "\n")

# Set up genome - MEDIPS expects genome name as string
if (opt$genome == "hg38") {
  BSgenome_name <- "BSgenome.Hsapiens.UCSC.hg38"
} else if (opt$genome == "hg19") {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  BSgenome_name <- "BSgenome.Hsapiens.UCSC.hg19"
} else {
  stop(paste("Unsupported genome:", opt$genome))
}

# Chromosome selection
if (!is.null(opt$chr_select)) {
  chr_select <- strsplit(opt$chr_select, ",")[[1]]
} else {
  chr_select <- paste0("chr", c(1:22, "X", "Y"))
}

cat("Chromosomes:", paste(chr_select, collapse=", "), "\n\n")

# Create MEDIPS SETs for Group 1
cat("Creating MEDIPS SETs for", opt$`group1-name`, "...\n")
MSet_list1 <- lapply(seq_along(bams_group1), function(i) {
  bam <- bams_group1[i]
  cat("  [", i, "/", length(bams_group1), "]", basename(bam), "\n")
  MEDIPS.createSet(
    file = bam,
    BSgenome = BSgenome_name,
    extend = opt$extend,
    shift = opt$shift,
    uniq = opt$uniq,
    window_size = opt$window,
    chr.select = chr_select,
    paired = opt$paired
  )
})

# Create MEDIPS SETs for Group 2
cat("\nCreating MEDIPS SETs for", opt$`group2-name`, "...\n")
MSet_list2 <- lapply(seq_along(bams_group2), function(i) {
  bam <- bams_group2[i]
  cat("  [", i, "/", length(bams_group2), "]", basename(bam), "\n")
  MEDIPS.createSet(
    file = bam,
    BSgenome = BSgenome_name,
    extend = opt$extend,
    shift = opt$shift,
    uniq = opt$uniq,
    window_size = opt$window,
    chr.select = chr_select,
    paired = opt$paired
  )
})

# Create coupling set (CpG density) using first sample as reference
cat("\nCalculating CpG coupling vector...\n")
CS <- MEDIPS.couplingVector(
  pattern = "CG",
  refObj = MSet_list1[[1]]
)

# Run differential MEDIPS analysis
cat("\nRunning differential methylation analysis...\n")
cat("Method:", opt$`diff-method`, "\n")
cat("P-value adjustment:", opt$`p-adj`, "\n")

mr <- MEDIPS.meth(
  MSet1 = MSet_list1,
  MSet2 = MSet_list2,
  CSet = CS,
  p.adj = opt$`p-adj`,
  diff.method = opt$`diff-method`,
  MeDIP = TRUE,
  CNV = FALSE,
  minRowSum = 10
)

# Convert results to data.table
cat("\nProcessing results...\n")
results_dt <- as.data.table(mr)

# Save full results
full_output <- paste0(opt$output, ".differential.tsv")
fwrite(results_dt, full_output, sep="\t")
cat("Full results saved to:", full_output, "\n")

# Extract significant DMRs (various thresholds)
cat("\nExtracting significant DMRs...\n")

# Find the log2 fold change column
logFC_col <- grep("logFC$|log2FC$|edgeR\\.logFC$", names(results_dt), value=TRUE)
pval_col <- grep("\\.p\\.value$|pvalue$|PValue$", names(results_dt), value=TRUE)
padj_col <- grep("adj\\.p\\.value$|padj$|FDR$", names(results_dt), value=TRUE)

if (length(logFC_col) > 0 && length(padj_col) > 0) {
  logFC_col <- logFC_col[1]
  padj_col <- padj_col[1]
  
  cat("Using columns: logFC=", logFC_col, ", padj=", padj_col, "\n")
  
  # DMRs at padj < 0.05
  dmr_005 <- results_dt[get(padj_col) < 0.05]
  dmr_005_output <- paste0(opt$output, ".DMRs.padj005.tsv")
  fwrite(dmr_005, dmr_005_output, sep="\t")
  cat("DMRs (padj < 0.05):", nrow(dmr_005), "-> saved to", dmr_005_output, "\n")
  
  # DMRs at padj < 0.01
  dmr_001 <- results_dt[get(padj_col) < 0.01]
  dmr_001_output <- paste0(opt$output, ".DMRs.padj001.tsv")
  fwrite(dmr_001, dmr_001_output, sep="\t")
  cat("DMRs (padj < 0.01):", nrow(dmr_001), "-> saved to", dmr_001_output, "\n")
  
  # Hyper-methylated in Group1 (positive logFC)
  hyper <- results_dt[get(padj_col) < 0.05 & get(logFC_col) > 0]
  hyper_output <- paste0(opt$output, ".DMRs.hyper.", opt$`group1-name`, ".tsv")
  fwrite(hyper, hyper_output, sep="\t")
  cat("Hyper-methylated in", opt$`group1-name`, ":", nrow(hyper), "-> saved to", hyper_output, "\n")
  
  # Hypo-methylated in Group1 (negative logFC)
  hypo <- results_dt[get(padj_col) < 0.05 & get(logFC_col) < 0]
  hypo_output <- paste0(opt$output, ".DMRs.hypo.", opt$`group1-name`, ".tsv")
  fwrite(hypo, hypo_output, sep="\t")
  cat("Hypo-methylated in", opt$`group1-name`, ":", nrow(hypo), "-> saved to", hypo_output, "\n")
  
  # Create BED files for downstream analysis
  if (nrow(hyper) > 0) {
    hyper_bed <- hyper[, .(chr, start, stop)]
    hyper_bed_output <- paste0(opt$output, ".DMRs.hyper.", opt$`group1-name`, ".bed")
    fwrite(hyper_bed, hyper_bed_output, sep="\t", col.names=FALSE)
    cat("Hyper BED:", hyper_bed_output, "\n")
  }
  
  if (nrow(hypo) > 0) {
    hypo_bed <- hypo[, .(chr, start, stop)]
    hypo_bed_output <- paste0(opt$output, ".DMRs.hypo.", opt$`group1-name`, ".bed")
    fwrite(hypo_bed, hypo_bed_output, sep="\t", col.names=FALSE)
    cat("Hypo BED:", hypo_bed_output, "\n")
  }
  
} else {
  cat("Warning: Could not identify logFC or adjusted p-value columns\n")
  cat("Available columns:", paste(names(results_dt), collapse=", "), "\n")
}

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total windows analyzed:", nrow(results_dt), "\n")
cat("Group 1 (", opt$`group1-name`, ") samples:", length(bams_group1), "\n")
cat("Group 2 (", opt$`group2-name`, ") samples:", length(bams_group2), "\n")

# Save summary
summary_output <- paste0(opt$output, ".summary.txt")
summary_data <- c(
  paste("Comparison:", opt$`group1-name`, "vs", opt$`group2-name`),
  paste("Group 1 samples:", length(bams_group1)),
  paste("Group 2 samples:", length(bams_group2)),
  paste("Total windows:", nrow(results_dt)),
  paste("Window size:", opt$window),
  paste("Genome:", opt$genome),
  paste("Differential method:", opt$`diff-method`),
  paste("P-value adjustment:", opt$`p-adj`)
)
if (exists("dmr_005")) {
  summary_data <- c(summary_data,
    paste("DMRs (padj < 0.05):", nrow(dmr_005)),
    paste("  Hyper in", opt$`group1-name`, ":", nrow(hyper)),
    paste("  Hypo in", opt$`group1-name`, ":", nrow(hypo))
  )
}
writeLines(summary_data, summary_output)
cat("Summary saved to:", summary_output, "\n")

cat("\n=== Differential Analysis Complete ===\n")
