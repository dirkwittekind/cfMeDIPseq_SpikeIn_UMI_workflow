#!/usr/bin/env Rscript
# MEDIPS Analysis Script
# Methylation quantification using MEDIPS, MeDEStrand, and QSEA

suppressPackageStartupMessages({
  library(MEDIPS)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(optparse)
  library(data.table)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-b", "--bam"), type="character", help="Input BAM file"),
  make_option(c("-o", "--output"), type="character", help="Output prefix"),
  make_option(c("-g", "--genome"), type="character", default="hg38", help="Genome assembly [default: %default]"),
  make_option(c("-w", "--window"), type="integer", default=300, help="Window size [default: %default]"),
  make_option(c("-e", "--extend"), type="integer", default=0, help="Extend reads [default: %default]"),
  make_option(c("-u", "--uniq"), type="integer", default=1, help="Uniqueness threshold [default: %default]"),
  make_option(c("-s", "--shift"), type="integer", default=0, help="Shift reads [default: %default]"),
  make_option(c("-p", "--paired"), action="store_true", default=TRUE, help="Paired-end data"),
  make_option(c("--chr_select"), type="character", default=NULL, help="Chromosomes to analyze (comma-separated)"),
  make_option(c("--method"), type="character", default="medips", help="Method: medips, medestrand, qsea, or all [default: medips only]"),
  make_option(c("--skip_qsea"), action="store_true", default=TRUE, help="Skip QSEA analysis (very slow) [default: TRUE]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$bam) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("BAM file and output prefix are required")
}

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

cat("=== MEDIPS Analysis ===\n")
cat("BAM:", opt$bam, "\n")
cat("Output:", opt$output, "\n")
cat("Genome:", opt$genome, "\n")
cat("Window:", opt$window, "\n")
cat("Chromosomes:", paste(chr_select, collapse=", "), "\n")

# Create MEDIPS SET
cat("\nCreating MEDIPS SET...\n")
MSet <- MEDIPS.createSet(
  file = opt$bam,
  BSgenome = BSgenome_name,
  extend = opt$extend,
  shift = opt$shift,
  uniq = opt$uniq,
  window_size = opt$window,
  chr.select = chr_select,
  paired = opt$paired
)

# Calculate coupling set (CpG density)
cat("Calculating coupling set...\n")
CS <- MEDIPS.couplingVector(
  pattern = "CG",
  refObj = MSet
)

# Run MEDIPS analysis
if (opt$method %in% c("medips", "all")) {
  cat("\nRunning MEDIPS methylation calling...\n")
  
  tryCatch({
    # Check for potential integer overflow before running
    # Use length of genome_count slot to get number of windows
    n_regions <- length(MSet@genome_count)
    potential_overflow <- as.numeric(opt$window) * as.numeric(n_regions)
    if (potential_overflow > 2^31 - 1) {
      warning(sprintf("Potential integer overflow detected: window_size(%d) * regions(%d) = %.0f > 2^31",
                      opt$window, n_regions, potential_overflow))
      cat("Consider using a larger window size (>= 1000bp) for whole-genome analysis\n")
    }
    
    # Calculate methylation values
    mr <- MEDIPS.meth(
      MSet1 = MSet,
      CSet = CS,
      p.adj = "bonferroni",
      diff.method = "edgeR",
      MeDIP = TRUE,
      CNV = FALSE,
      minRowSum = 10
    )
    
    # Convert to data.table
    medips_dt <- as.data.table(mr)
    
    # Save results
    medips_output <- paste0(opt$output, ".medips.tsv")
    fwrite(medips_dt, medips_output, sep="\t")
    cat("MEDIPS results saved to:", medips_output, "\n")
    
    # Save RMS (relative methylation score) as bedGraph
    rms_output <- paste0(opt$output, ".medips.rms.bedGraph")
    # Column names may vary - find the RMS column dynamically
    rms_col <- grep("rms$|RMS$", names(medips_dt), value=TRUE)
    if (length(rms_col) > 0) {
      rms_col <- rms_col[1]
      rms_dt <- medips_dt[, c("chr", "start", "stop", rms_col), with=FALSE]
      setnames(rms_dt, c("chr", "start", "end", "rms"))
      fwrite(rms_dt, rms_output, sep="\t", col.names=FALSE)
      cat("RMS bedGraph saved to:", rms_output, "\n")
    } else {
      cat("Warning: No RMS column found in MEDIPS output\n")
      cat("Available columns:", paste(names(medips_dt), collapse=", "), "\n")
    }
    
  }, error = function(e) {
    cat("\nERROR in MEDIPS analysis:", conditionMessage(e), "\n")
    cat("This may be caused by:\n")
    cat("  - Integer overflow: Use larger window size (meth_qc_window_size >= 1000)\n")
    cat("  - Memory issues: Reduce parallel jobs or increase RAM\n")
    cat("  - BAM file issues: Check BAM is properly indexed and paired-end\n")
    stop(e)
  })
}

# Run MeDEStrand analysis
if (opt$method %in% c("medestrand", "all")) {
  cat("\nRunning MeDEStrand analysis...\n")
  
  tryCatch({
    suppressPackageStartupMessages(library(MeDEStrand))
    
    MeDESet <- MeDEStrand.createSet(
      file = opt$bam,
      BSgenome = BSgenome_name,
      extend = opt$extend,
      shift = opt$shift,
      uniq = opt$uniq,
      window_size = opt$window,
      chr.select = chr_select,
      paired = opt$paired
    )
    
    # Calculate strand-specific methylation
    mede_result <- MeDEStrand.meth(
      MSet1 = MeDESet,
      CSet = CS,
      minRowSum = 10
    )
    
    mede_dt <- as.data.table(mede_result)
    mede_output <- paste0(opt$output, ".medestrand.tsv")
    fwrite(mede_dt, mede_output, sep="\t")
    cat("MeDEStrand results saved to:", mede_output, "\n")
    
  }, error = function(e) {
    cat("Warning: MeDEStrand analysis failed:", conditionMessage(e), "\n")
  })
}

# Run QSEA analysis (optional - very slow for whole genome)
if (opt$method %in% c("qsea", "all") && !opt$skip_qsea) {
  cat("\nRunning QSEA analysis...\n")
  
  tryCatch({
    suppressPackageStartupMessages(library(qsea))
    
    # Create sample table
    sample_table <- data.frame(
      sample_name = basename(opt$bam),
      file_name = opt$bam,
      group = "sample",
      stringsAsFactors = FALSE
    )
    
    # Create QSEA set
    qs <- createQseaSet(
      sampleTable = sample_table,
      BSgenome = BSgenome_name,
      window_size = opt$window,
      chr.select = chr_select
    )
    
    # Add coverage
    qs <- addCoverage(qs, uniquePos = TRUE, paired = opt$paired)
    
    # Add CpG enrichment
    qs <- addPatternDensity(qs, "CG", name = "CpG")
    
    # Normalize
    qs <- addLibraryFactors(qs)
    
    # Calculate beta values
    qs <- addEnrichmentParameters(qs, enrichmentPattern = "CpG")
    
    # Get results
    qsea_result <- makeTable(qs, samples = 1, norm_method = "beta")
    qsea_dt <- as.data.table(qsea_result)
    
    qsea_output <- paste0(opt$output, ".qsea.tsv")
    fwrite(qsea_dt, qsea_output, sep="\t")
    cat("QSEA results saved to:", qsea_output, "\n")
    
  }, error = function(e) {
    cat("Warning: QSEA analysis failed:", conditionMessage(e), "\n")
  })
}

# Calculate summary statistics
cat("\n=== Summary Statistics ===\n")

# Get slot names safely
slot_names <- slotNames(MSet)
cat("Available slots:", paste(slot_names, collapse=", "), "\n")

# Total windows
n_windows <- length(MSet@genome_count)
cat("Total windows:", n_windows, "\n")

# Total reads - try different slot names
total_reads <- tryCatch({
  if ("number_reads" %in% slot_names) {
    MSet@number_reads
  } else if ("reads" %in% slot_names) {
    MSet@reads
  } else {
    sum(MSet@genome_count)
  }
}, error = function(e) sum(MSet@genome_count))
cat("Total reads:", total_reads, "\n")

# Genome coverage - try different slot names
genome_cov <- tryCatch({
  if ("genome_coverage" %in% slot_names) {
    MSet@genome_coverage
  } else {
    sum(MSet@genome_count > 0) / n_windows
  }
}, error = function(e) sum(MSet@genome_count > 0) / n_windows)
cat("Genome-wide coverage:", round(genome_cov, 4), "\n")

# Save summary
summary_output <- paste0(opt$output, ".summary.txt")
summary_data <- c(
  paste("Sample:", basename(opt$bam)),
  paste("Total reads:", total_reads),
  paste("Total windows:", n_windows),
  paste("Genome-wide coverage:", round(genome_cov, 4)),
  paste("Window size:", opt$window),
  paste("Genome:", opt$genome)
)
writeLines(summary_data, summary_output)

cat("\nAnalysis complete!\n")
