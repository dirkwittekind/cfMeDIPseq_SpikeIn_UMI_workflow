#!/usr/bin/env Rscript
#' DMR Analysis using limma for cfMeDIP-seq Windows
#' Computes differential methylation between groups using limma/eBayes.
#' Based on standard workflow for AEG vs CTRL discrimination.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(argparse)
})

# Parse arguments
parser <- ArgumentParser(description = "DMR analysis with limma for cfMeDIP-seq")
parser$add_argument("--ann", required = TRUE, help = "Sample annotation TSV")
parser$add_argument("--matrix", required = TRUE, help = "Feature matrix TSV")
parser$add_argument("--out-dir", required = TRUE, help = "Output directory")
parser$add_argument("--prefix", default = "dmr", help = "Output file prefix")
parser$add_argument("--case-group", default = "AEG", help = "Case group name")
parser$add_argument("--ctrl-group", default = "CTRL", help = "Control group name")
parser$add_argument("--covariates", nargs = "+", default = NULL, help = "Covariate columns to adjust for")
parser$add_argument("--top-n", type = "integer", default = 1000, help = "Number of top DMRs to output separately")
parser$add_argument("--fdr-threshold", type = "double", default = 0.05, help = "FDR threshold for significance")

args <- parser$parse_args()

cat("=" , rep("=", 59), "\n", sep = "")
cat("DMR Analysis with limma\n")
cat("=" , rep("=", 59), "\n", sep = "")

# Load annotation
cat("\nLoading annotation:", args$ann, "\n")
ann <- fread(args$ann, sep = "\t", data.table = FALSE)
colnames(ann) <- trimws(colnames(ann))
ann$sample_id <- trimws(as.character(ann$sample_id))
ann$group <- trimws(as.character(ann$group))

# Map common control names
ann$group[ann$group == "Kontrolle"] <- "CTRL"
ann$group[ann$group == "Control"] <- "CTRL"
ann$group[ann$group == "control"] <- "CTRL"

keep <- ann$sample_id
cat("  Samples:", length(keep), "\n")
cat("  Groups:", paste(table(ann$group), names(table(ann$group)), collapse = ", "), "\n")

# Load matrix
cat("\nLoading matrix:", args$matrix, "\n")
x <- fread(args$matrix, sep = "\t", data.table = FALSE)

# Clean column names: remove quotes, #, and .bw suffix from BigWig output
clean_col <- function(c) {
  c <- gsub("^[#'\"]+|['\"]+$", "", c)
  c <- trimws(c)
  if (grepl("\\.bw$", c)) {
    c <- sub("\\.bw$", "", c)
  }
  return(c)
}
colnames(x) <- sapply(colnames(x), clean_col)

feat_col <- colnames(x)[1]
x[[feat_col]] <- trimws(as.character(x[[feat_col]]))

# Detect matrix orientation
sample_cols <- intersect(colnames(x), keep)
cat("  Samples found in matrix:", length(sample_cols), "\n")

if (length(sample_cols) >= 2) {
  # Samples in columns
  X <- as.matrix(x[, c(feat_col, sample_cols)])
  rownames(X) <- X[, 1]
  X <- X[, -1, drop = FALSE]
} else {
  stop("Matrix orientation unclear: Samples not found in columns. Please check your matrix format.")
}

# Reorder columns to match annotation
X <- X[, keep, drop = FALSE]
mode(X) <- "numeric"
X[is.na(X)] <- 0

cat("  Features:", nrow(X), "\n")
cat("  Samples:", ncol(X), "\n")

# Log-transform (if signal is CPM/normalized)
cat("\nLog-transforming data (log2(x + 1))...\n")
Y <- log2(X + 1)

# Create design matrix
cat("\nBuilding design matrix...\n")
cat("  Case group:", args$case_group, "\n")
cat("  Control group:", args$ctrl_group, "\n")

# Ensure factor levels: control first, then case
ann_ordered <- ann[match(keep, ann$sample_id), ]
group_factor <- factor(ann_ordered$group, levels = c(args$ctrl_group, args$case_group))

if (!is.null(args$covariates) && length(args$covariates) > 0) {
  # Build design with covariates
  cat("  Covariates:", paste(args$covariates, collapse = ", "), "\n")
  
  cov_data <- data.frame(group = group_factor)
  for (cov in args$covariates) {
    if (cov %in% colnames(ann_ordered)) {
      cov_vals <- ann_ordered[[cov]]
      # Try numeric first, then factor
      cov_num <- suppressWarnings(as.numeric(cov_vals))
      if (all(!is.na(cov_num))) {
        cov_data[[cov]] <- cov_num
      } else {
        cov_data[[cov]] <- factor(cov_vals)
      }
    }
  }
  
  formula_str <- paste("~ group +", paste(args$covariates, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = cov_data)
} else {
  design <- model.matrix(~ group_factor)
}

cat("  Design matrix columns:", paste(colnames(design), collapse = ", "), "\n")

# Run limma
cat("\nFitting linear models...\n")
fit <- lmFit(Y, design)

cat("Running empirical Bayes moderation...\n")
fit <- eBayes(fit, trend = TRUE)

# Get coefficient name for case group
coef_name <- grep(args$case_group, colnames(design), value = TRUE)
if (length(coef_name) == 0) {
  coef_name <- colnames(design)[2]  # Default to second column
}
cat("  Testing coefficient:", coef_name, "\n")

# Extract results
cat("\nExtracting results...\n")
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
tt$feature <- rownames(tt)

# Add direction interpretation
tt$direction <- ifelse(tt$logFC > 0,
                       paste0("hyper_in_", args$case_group),
                       paste0("hypo_in_", args$case_group))

# Summary statistics
n_sig <- sum(tt$adj.P.Val < args$fdr_threshold)
n_hyper <- sum(tt$logFC > 0 & tt$adj.P.Val < args$fdr_threshold)
n_hypo <- sum(tt$logFC < 0 & tt$adj.P.Val < args$fdr_threshold)

cat("\n" , rep("=", 60), "\n", sep = "")
cat("Results Summary\n")
cat(rep("=", 60), "\n", sep = "")
cat("Total features:", nrow(tt), "\n")
cat("Significant (FDR <", args$fdr_threshold, "):", n_sig, "\n")
cat("  Hyper in", args$case_group, ":", n_hyper, "\n")
cat("  Hypo in", args$case_group, ":", n_hypo, "\n")

# Create output directory
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# Save full results
out_file <- file.path(args$out_dir, paste0(args$prefix, ".tsv"))
fwrite(tt, out_file, sep = "\t")
cat("\nFull results:", out_file, "\n")

# Compress
system(paste("gzip -f", shQuote(out_file)))
cat("Compressed to:", paste0(out_file, ".gz"), "\n")

# Save top N
top_file <- file.path(args$out_dir, paste0(args$prefix, ".top", args$top_n, ".tsv"))
fwrite(head(tt, args$top_n), top_file, sep = "\t")
cat("Top", args$top_n, "DMRs:", top_file, "\n")

# Save significant features
if (n_sig > 0) {
  sig_file <- file.path(args$out_dir, paste0(args$prefix, ".significant.tsv"))
  fwrite(tt[tt$adj.P.Val < args$fdr_threshold, ], sig_file, sep = "\t")
  cat("Significant features:", sig_file, "\n")
}

# Generate volcano plot
cat("\nGenerating volcano plot...\n")
volcano_file <- file.path(args$out_dir, paste0(args$prefix, ".volcano.png"))

png(volcano_file, width = 800, height = 600, res = 100)
par(mar = c(5, 5, 4, 2))

# Color by significance and direction
cols <- rep("gray60", nrow(tt))
cols[tt$adj.P.Val < args$fdr_threshold & tt$logFC > 0] <- "red"
cols[tt$adj.P.Val < args$fdr_threshold & tt$logFC < 0] <- "blue"

plot(tt$logFC, -log10(tt$P.Value),
     pch = 16, cex = 0.5, col = cols,
     xlab = paste("log2 FC (", args$ctrl_group, " -> ", args$case_group, ")", sep = ""),
     ylab = "-log10(P-value)",
     main = paste("DMR Volcano Plot:", args$case_group, "vs", args$ctrl_group))

# Add threshold lines
abline(h = -log10(0.05), lty = 2, col = "gray40")
abline(v = c(-0.5, 0.5), lty = 2, col = "gray40")

# Legend
legend("topright",
       legend = c(paste("Hyper in", args$case_group, "(n=", n_hyper, ")"),
                  paste("Hypo in", args$case_group, "(n=", n_hypo, ")"),
                  "Not significant"),
       pch = 16, col = c("red", "blue", "gray60"),
       bty = "n")

dev.off()
cat("Volcano plot:", volcano_file, "\n")

# Generate QQ plot
qq_file <- file.path(args$out_dir, paste0(args$prefix, ".qq.png"))

png(qq_file, width = 600, height = 600, res = 100)
par(mar = c(5, 5, 4, 2))

# Expected vs observed p-values
n <- nrow(tt)
expected <- -log10((1:n) / (n + 1))
observed <- sort(-log10(tt$P.Value))

plot(expected, observed,
     pch = 16, cex = 0.5,
     xlab = "Expected -log10(P-value)",
     ylab = "Observed -log10(P-value)",
     main = "QQ Plot")
abline(0, 1, col = "red", lwd = 2)

# Calculate genomic inflation factor
lambda <- median(qchisq(1 - tt$P.Value, df = 1)) / qchisq(0.5, df = 1)
legend("topleft", 
       legend = paste("λ =", round(lambda, 3)),
       bty = "n")

dev.off()
cat("QQ plot:", qq_file, "\n")

# Save summary statistics
summary_df <- data.frame(
  metric = c("n_features", "n_significant", "n_hyper", "n_hypo", 
             "fdr_threshold", "lambda_gc"),
  value = c(nrow(tt), n_sig, n_hyper, n_hypo, 
            args$fdr_threshold, round(lambda, 4))
)
summary_file <- file.path(args$out_dir, paste0(args$prefix, ".summary.tsv"))
fwrite(summary_df, summary_file, sep = "\t")
cat("Summary:", summary_file, "\n")

cat("\nDMR analysis complete!\n")
