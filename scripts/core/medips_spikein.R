args <- commandArgs(trailingOnly = TRUE)

sample_id    <- args[1]
sample_bam   <- args[2]
ispaired     <- as.logical(args[3])
ws           <- as.numeric(args[4])
bsgenome     <- args[5]
bsgenome_pkg <- args[6]

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(MEDIPS))

# Safety checks
if (!file.exists(sample_bam)) stop("BAM not found: ", sample_bam)
if (is.na(ws) || ws <= 0) stop("Invalid window_size: ", ws)
if (!nzchar(bsgenome)) stop("BSgenome package name is empty")

# If bsgenome not installed, try installing from tarball (should normally be handled by Snakemake install rule)
if (!(bsgenome %in% rownames(installed.packages()))) {
  if (nzchar(bsgenome_pkg) && file.exists(bsgenome_pkg) && grepl("\\.tar\\.gz$", bsgenome_pkg)) {
    message("Installing BSgenome from tarball: ", bsgenome_pkg)
    install.packages(bsgenome_pkg, repos = NULL, type = "source")
  } else {
    stop(
      "BSgenome package not installed in this R env: ", bsgenome, "\n",
      "Installation must be done BEFORE running this script (handled by Snakemake install rule)."
    )
  }
}

suppressPackageStartupMessages(library(bsgenome, character.only = TRUE))

# Determine spike-in chromosomes from BAM header
bam_hdr <- Rsamtools::scanBamHeader(sample_bam)[[1]]$targets
chr <- names(bam_hdr)[grepl("^(160b|320b)_", names(bam_hdr))]
if (length(chr) == 0) {
  stop("No spike-in contigs found (expected names starting with 160b_ or 320b_) in BAM header: ", sample_bam)
}

dir.create("meth_qc_quant_spikein", showWarnings = FALSE, recursive = TRUE)

###################################
###         MEDIPS QC           ###
###################################

saturation <- MEDIPS.saturation(
  file = sample_bam,
  BSgenome = bsgenome,
  window_size = ws,
  paired = ispaired,
  uniq = 0,
  chr.select = chr
)

png(file = paste0("meth_qc_quant_spikein/", sample_id, "_saturation.png"),
    res = 300, width = 5, height = 5, units = "in")
MEDIPS.plotSaturation(saturationObj = saturation)
dev.off()

coverage <- MEDIPS.seqCoverage(
  file = sample_bam,
  BSgenome = bsgenome,
  paired = ispaired,
  uniq = 0,
  chr.select = chr
)

png(file = paste0("meth_qc_quant_spikein/", sample_id, "_seqCoverage.png"),
    res = 300, width = 5, height = 5, units = "in")
MEDIPS.plotSeqCoverage(
  seqCoverageObj = coverage,
  type = "pie",
  cov.level = c(0, 1, 2, 3, 4, 5)
)
dev.off()

# Revised enrichment function (kept as in your pipeline)
MEDIPS.CpGenrichNew <- function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=FALSE) {

  if (is.null(BSgenome)) stop("Must specify a BSgenome library.")

  fileName <- basename(file)
  path <- dirname(file)
  if (path == "") path <- getwd()
  if (!fileName %in% dir(path)) stop(paste("File", fileName, "not found in", path))

  dataset <- get(ls(paste("package:", BSgenome, sep=""))[1])

  if (!paired) {
    GRange.Reads <- getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)
  } else {
    GRange.Reads <- getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)
  }

  suppressPackageStartupMessages(library(gtools))
  if (length(unique(seqlevels(GRange.Reads))) > 1) chromosomes <- mixedsort(unique(seqlevels(GRange.Reads)))
  if (length(unique(seqlevels(GRange.Reads))) == 1) chromosomes <- unique(seqlevels(GRange.Reads))

  readsChars <- unlist(getSeq(dataset, GRange.Reads, as.character=TRUE))

  regions.CG <- sum(vcountPattern("CG", readsChars))
  regions.C  <- sum(vcountPattern("C", readsChars))
  regions.G  <- sum(vcountPattern("G", readsChars))
  all.genomic <- sum(width(readsChars))

  regions.relH <- as.numeric(regions.CG) / as.numeric(all.genomic) * 100
  regions.GoGe <- (as.numeric(regions.CG) * as.numeric(all.genomic)) / (as.numeric(regions.C) * as.numeric(regions.G))

  CG <- DNAStringSet("CG")
  pdict0 <- PDict(CG)
  params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
  genome.CG <- sum(bsapply(params, pdict = pdict0))
  params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify=TRUE)
  alphabet <- bsapply(params)
  genome.l <- sum(as.numeric(alphabet))
  genome.C <- as.numeric(sum(alphabet[2,]))
  genome.G <- as.numeric(sum(alphabet[3,]))
  genome.relH <- genome.CG / genome.l * 100
  genome.GoGe <- (genome.CG * genome.l) / (genome.C * genome.G)

  enrichment.score.relH <- regions.relH / genome.relH
  enrichment.score.GoGe <- regions.GoGe / genome.GoGe

  return(list(
    genome = BSgenome,
    regions.CG = regions.CG, regions.C = regions.C, regions.G = regions.G,
    regions.relH = regions.relH, regions.GoGe = regions.GoGe,
    genome.C = genome.C, genome.G = genome.G, genome.CG = genome.CG,
    genome.relH = genome.relH, genome.GoGe = genome.GoGe,
    enrichment.score.relH = enrichment.score.relH,
    enrichment.score.GoGe = enrichment.score.GoGe
  ))
}

cpg_enrich <- MEDIPS.CpGenrichNew(
  file = sample_bam,
  BSgenome = bsgenome,
  paired = ispaired,
  uniq = 0,
  chr.select = chr
)

cov_sum <- table(coverage$cov.res)
numberCpG <- length(coverage$cov.res)

idx_6 <- (names(cov_sum) != "0" & names(cov_sum) != "1" & names(cov_sum) != "2" &
          names(cov_sum) != "3" & names(cov_sum) != "4" & names(cov_sum) != "5")

coverage.fracCpGwoRead <- if (length(cov_sum[names(cov_sum) == "0"]) == 0) 0 else cov_sum[names(cov_sum) == "0"] / numberCpG
coverage.fracCpGw1Read <- if (length(cov_sum[names(cov_sum) == "1"]) == 0) 0 else cov_sum[names(cov_sum) == "1"] / numberCpG
coverage.fracCpGw2Read <- if (length(cov_sum[names(cov_sum) == "2"]) == 0) 0 else cov_sum[names(cov_sum) == "2"] / numberCpG
coverage.fracCpGw3Read <- if (length(cov_sum[names(cov_sum) == "3"]) == 0) 0 else cov_sum[names(cov_sum) == "3"] / numberCpG
coverage.fracCpGw4Read <- if (length(cov_sum[names(cov_sum) == "4"]) == 0) 0 else cov_sum[names(cov_sum) == "4"] / numberCpG
coverage.fracCpGw5Read <- if (length(cov_sum[names(cov_sum) == "5"]) == 0) 0 else cov_sum[names(cov_sum) == "5"] / numberCpG

medips_qc <- data.frame(
  Sample = sample_id,
  saturation.numberReads = saturation$numberReads,
  saturation.maxEstCor = saturation$maxEstCor[2],
  saturation.maxTruCor = saturation$maxTruCor[2],
  coverage.fracReadsWoCpG = coverage$numberReadsWO / coverage$numberReads,
  coverage.fracCpGwoRead = coverage.fracCpGwoRead,
  coverage.fracCpGw1Read = coverage.fracCpGw1Read,
  coverage.fracCpGw2Read = coverage.fracCpGw2Read,
  coverage.fracCpGw3Read = coverage.fracCpGw3Read,
  coverage.fracCpGw4Read = coverage.fracCpGw4Read,
  coverage.fracCpGw5Read = coverage.fracCpGw5Read,
  coverage.fracCpGgt5Reads = sum(cov_sum[idx_6]) / numberCpG,
  enrichment.regions.CG = cpg_enrich$regions.CG,
  enrichment.regions.C = cpg_enrich$regions.C,
  enrichment.regions.G = cpg_enrich$regions.G,
  enrichment.regions.relH = cpg_enrich$regions.relH,
  enrichment.regions.GoGe = cpg_enrich$regions.GoGe,
  enrichment.genome.CG = cpg_enrich$genome.CG,
  enrichment.genome.C = cpg_enrich$genome.C,
  enrichment.genome.G = cpg_enrich$genome.G,
  enrichment.genome.relH = cpg_enrich$genome.relH,
  enrichment.genome.GoGe = cpg_enrich$genome.GoGe,
  enrichment.relH = cpg_enrich$enrichment.score.relH,
  enrichment.GoGe = cpg_enrich$enrichment.score.GoGe
)

write.table(
  medips_qc,
  file = paste0("meth_qc_quant_spikein/", sample_id, "_meth_qc.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

###################################
###   MEDIPS quantification     ###
###################################

mset <- MEDIPS.createSet(
  file = sample_bam,
  BSgenome = bsgenome,
  paired = ispaired,
  window_size = ws,
  uniq = 0,
  chr.select = chr
)

cset <- MEDIPS.couplingVector(pattern = "CG", refObj = mset)
CpG <- cset@genome_CF

meth <- MEDIPS.meth(
  MSet1 = mset,
  CSet = cset,
  CNV = FALSE,
  MeDIP = FALSE
)

meth <- data.frame(meth[, 1:5], CpG)
colnames(meth)[c(4, 5)] <- c("count", "rpkm")

bin_id <- 1:nrow(meth)
strand <- rep(".", nrow(meth))
grange_cpg <- data.frame(meth[, 1:3], bin_id, CpG, strand)
colnames(grange_cpg)[1] <- "#chr"

write.table(
  grange_cpg,
  file = paste0("meth_qc_quant_spikein/", sample_id, "_Granges_CpGs.bed"),
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)

for (i in 4:5) {
  out <- meth[, i]
  feature <- colnames(meth)[i]
  write.table(
    out,
    file = paste0("meth_qc_quant_spikein/", sample_id, "_", feature, ".txt"),
    row.names = FALSE, col.names = sample_id, quote = FALSE
  )
}

save(meth, file = paste0("meth_qc_quant_spikein/", sample_id, "_meth_quant.RData"))

