#!/usr/bin/env Rscript
# Fragment Profile Analysis Script
# Analyzes cfDNA fragment size distribution and ratios

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(data.table)
  library(optparse)
  library(ggplot2)
})

# Parse arguments
option_list <- list(
  make_option(c("-b", "--bam"), type="character", help="Input BAM file"),
  make_option(c("-o", "--output"), type="character", help="Output prefix"),
  make_option(c("-r", "--regions"), type="character", default=NULL, help="BED file with regions of interest"),
  make_option(c("-w", "--window"), type="integer", default=5000000, help="Window size for genome-wide analysis [default: %default]"),
  make_option(c("--min_size"), type="integer", default=50, help="Minimum fragment size [default: %default]"),
  make_option(c("--max_size"), type="integer", default=500, help="Maximum fragment size [default: %default]"),
  make_option(c("--short_min"), type="integer", default=100, help="Short fragment minimum [default: %default]"),
  make_option(c("--short_max"), type="integer", default=150, help="Short fragment maximum [default: %default]"),
  make_option(c("--long_min"), type="integer", default=151, help="Long fragment minimum [default: %default]"),
  make_option(c("--long_max"), type="integer", default=220, help="Long fragment maximum [default: %default]"),
  make_option(c("--chrom_sizes"), type="character", default=NULL, help="Chromosome sizes file"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$bam) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("BAM file and output prefix are required")
}

cat("=== Fragment Profile Analysis ===\n")
cat("BAM:", opt$bam, "\n")
cat("Output:", opt$output, "\n")

# Function to extract insert sizes from BAM
extract_insert_sizes <- function(bam_file, regions = NULL, min_size = 50, max_size = 500) {
  param <- ScanBamParam(
    flag = scanBamFlag(
      isPaired = TRUE,
      isProperPair = TRUE,
      isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE,
      isDuplicate = FALSE
    ),
    what = c("isize", "rname", "pos")
  )
  
  if (!is.null(regions)) {
    param <- ScanBamParam(param, which = regions)
  }
  
  bam <- scanBam(bam_file, param = param)[[1]]
  
  # Get absolute insert sizes and filter
  isizes <- abs(bam$isize)
  isizes <- isizes[!is.na(isizes) & isizes >= min_size & isizes <= max_size]
  
  return(isizes)
}

# Function to calculate fragment ratio
calc_fragment_ratio <- function(isizes, short_min, short_max, long_min, long_max) {
  short_frags <- sum(isizes >= short_min & isizes <= short_max)
  long_frags <- sum(isizes >= long_min & isizes <= long_max)
  
  if (long_frags == 0) {
    return(NA)
  }
  
  return(short_frags / long_frags)
}

# Extract genome-wide insert sizes
cat("\nExtracting insert sizes...\n")
insert_sizes <- extract_insert_sizes(
  opt$bam, 
  min_size = opt$min_size, 
  max_size = opt$max_size
)

cat("Total fragments:", length(insert_sizes), "\n")

# Calculate basic statistics
size_stats <- data.table(
  metric = c("mean", "median", "sd", "min", "max", "n"),
  value = c(
    mean(insert_sizes),
    median(insert_sizes),
    sd(insert_sizes),
    min(insert_sizes),
    max(insert_sizes),
    length(insert_sizes)
  )
)

# Calculate fragment ratio
frag_ratio <- calc_fragment_ratio(
  insert_sizes,
  opt$short_min, opt$short_max,
  opt$long_min, opt$long_max
)

short_count <- sum(insert_sizes >= opt$short_min & insert_sizes <= opt$short_max)
long_count <- sum(insert_sizes >= opt$long_min & insert_sizes <= opt$long_max)

cat("Short fragments (", opt$short_min, "-", opt$short_max, "bp):", short_count, "\n")
cat("Long fragments (", opt$long_min, "-", opt$long_max, "bp):", long_count, "\n")
cat("Fragment ratio:", round(frag_ratio, 4), "\n")

# Create fragment size distribution
cat("\nCreating fragment size distribution...\n")
size_dist <- data.table(size = insert_sizes)[, .N, by = size][order(size)]
setnames(size_dist, c("fragment_size", "count"))
size_dist[, fraction := count / sum(count)]

# Save distribution
dist_output <- paste0(opt$output, ".fragment_distribution.tsv")
fwrite(size_dist, dist_output, sep="\t")
cat("Distribution saved to:", dist_output, "\n")

# Create plot
cat("Creating distribution plot...\n")
plot_output <- paste0(opt$output, ".fragment_distribution.pdf")

p <- ggplot(size_dist, aes(x = fragment_size, y = fraction)) +
  geom_line(color = "steelblue", size = 0.8) +
  geom_vline(xintercept = c(opt$short_min, opt$short_max), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = c(opt$long_min, opt$long_max), linetype = "dashed", color = "blue", alpha = 0.5) +
  annotate("rect", xmin = opt$short_min, xmax = opt$short_max, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "red") +
  annotate("rect", xmin = opt$long_min, xmax = opt$long_max, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "blue") +
  labs(
    title = paste("Fragment Size Distribution -", basename(opt$bam)),
    subtitle = paste("Ratio (short/long):", round(frag_ratio, 3)),
    x = "Fragment Size (bp)",
    y = "Fraction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

ggsave(plot_output, p, width = 8, height = 5)
cat("Plot saved to:", plot_output, "\n")

# Genome-wide fragment ratio analysis
if (!is.null(opt$chrom_sizes)) {
  cat("\nPerforming genome-wide fragment ratio analysis...\n")
  
  # Read chromosome sizes
  chrom_sizes <- fread(opt$chrom_sizes, header = FALSE, col.names = c("chr", "size"))
  
  # Filter to standard chromosomes
  std_chroms <- paste0("chr", c(1:22, "X", "Y"))
  chrom_sizes <- chrom_sizes[chr %in% std_chroms]
  
  # Create windows
  windows <- lapply(1:nrow(chrom_sizes), function(i) {
    chr <- chrom_sizes$chr[i]
    size <- chrom_sizes$size[i]
    starts <- seq(1, size, by = opt$window)
    ends <- pmin(starts + opt$window - 1, size)
    data.table(chr = chr, start = starts, end = ends)
  })
  windows <- rbindlist(windows)
  
  # Convert to GRanges
  windows_gr <- GRanges(
    seqnames = windows$chr,
    ranges = IRanges(start = windows$start, end = windows$end)
  )
  
  # Read BAM with position info
  cat("Reading BAM file for genome-wide analysis...\n")
  param <- ScanBamParam(
    flag = scanBamFlag(
      isPaired = TRUE,
      isProperPair = TRUE,
      isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE,
      isDuplicate = FALSE
    ),
    what = c("isize", "rname", "pos")
  )
  
  bam_data <- scanBam(opt$bam, param = param)[[1]]
  
  # Create fragment GRanges
  valid_idx <- !is.na(bam_data$isize) & abs(bam_data$isize) >= opt$min_size & abs(bam_data$isize) <= opt$max_size
  frag_gr <- GRanges(
    seqnames = bam_data$rname[valid_idx],
    ranges = IRanges(start = bam_data$pos[valid_idx], width = 1)
  )
  frag_gr$isize <- abs(bam_data$isize[valid_idx])
  
  # Calculate ratio per window
  cat("Calculating ratios per window...\n")
  overlaps <- findOverlaps(frag_gr, windows_gr)
  
  frag_dt <- data.table(
    window_idx = subjectHits(overlaps),
    isize = frag_gr$isize[queryHits(overlaps)]
  )
  
  window_ratios <- frag_dt[, {
    short <- sum(isize >= opt$short_min & isize <= opt$short_max)
    long <- sum(isize >= opt$long_min & isize <= opt$long_max)
    ratio <- ifelse(long > 0, short / long, NA)
    list(short_count = short, long_count = long, ratio = ratio, total = .N)
  }, by = window_idx]
  
  # Merge with window coordinates
  windows$window_idx <- 1:nrow(windows)
  window_results <- merge(windows, window_ratios, by = "window_idx", all.x = TRUE)
  window_results[, window_idx := NULL]
  
  # Save genome-wide results
  gw_output <- paste0(opt$output, ".genome_wide_ratios.tsv")
  fwrite(window_results, gw_output, sep="\t")
  cat("Genome-wide ratios saved to:", gw_output, "\n")
}

# Save summary
summary_output <- paste0(opt$output, ".fragment_summary.tsv")
summary_dt <- data.table(
  sample = basename(opt$bam),
  total_fragments = length(insert_sizes),
  mean_size = mean(insert_sizes),
  median_size = median(insert_sizes),
  sd_size = sd(insert_sizes),
  short_count = short_count,
  long_count = long_count,
  fragment_ratio = frag_ratio
)
fwrite(summary_dt, summary_output, sep="\t")
cat("Summary saved to:", summary_output, "\n")

cat("\nFragment profile analysis complete!\n")
