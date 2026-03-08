# =============================================================================
# MEDIPIPE Core Rules - Spike-in Analysis
# =============================================================================
# For cfMeDIP-seq experiments with spike-in controls for normalization.
# Supports custom BSgenome packages for spike-in sequences.
# =============================================================================

import os

# -----------------------------------------------------------------------------
# Spike-in Configuration
# -----------------------------------------------------------------------------
# Check for spike-in mode (backward compatible with multiple config keys)
USE_SPIKEIN = bool(config.get("use_spikein", False) or config.get("spike_in", False))

# Spike-in specific paths
SPIKEIN_RLIB = os.path.join(WORK_DIR, "spikein_bsgenome/Rlib")
SPIKEIN_OK = os.path.join(WORK_DIR, "spikein_bsgenome/.installed.ok")
SPIKEIN_BSG = config.get("spike_in_bsgenome", "")
SPIKEIN_PKG = config.get("spike_in_bsgenome_pkg", "")
SPIKEIN_CHR = config.get("spike_in_chr", "")

# -----------------------------------------------------------------------------
# Combined BWA Index (spike-in aware)
# -----------------------------------------------------------------------------
def get_bwa_index_spikein():
    """Get BWA index path - returns combined index if spike-in enabled"""
    if USE_SPIKEIN:
        return config.get("spike_idx", get_bwa_index())
    return get_bwa_index()


# =============================================================================
# Rules (only active when USE_SPIKEIN is True)
# =============================================================================

if USE_SPIKEIN:
    
    # -------------------------------------------------------------------------
    # Install Custom BSgenome for Spike-ins
    # -------------------------------------------------------------------------
    localrules: install_spikein_bsgenome
    
    rule install_spikein_bsgenome:
        """Install custom BSgenome tarball into a persistent library under the workdir.
        
        This avoids conda-env library locks and survives env recreation.
        The spike-in BSgenome package (e.g., BSgenome.Hsapiens.UCSC.hg38plusspikeins)
        is installed from a local tarball specified in config['spike_in_bsgenome_pkg'].
        """
        input:
            tarball=SPIKEIN_PKG
        output:
            ok=SPIKEIN_OK
        conda:
            os.path.join(ENV_DIR, "medipipe_r.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p "{SPIKEIN_RLIB}"
            mkdir -p "$(dirname "{output.ok}")"

            if [ ! -s "{input.tarball}" ]; then
              echo "[ERROR] spike_in_bsgenome_pkg tarball missing: {input.tarball}" >&2
              exit 1
            fi

            export R_LIBS_USER="$(pwd)/{SPIKEIN_RLIB}"
            # Clean up any lock files from previous failed installs
            find "$R_LIBS_USER" -maxdepth 1 -type d -name "00LOCK-*" -exec rm -rf {{}} + 2>/dev/null || true

            Rscript --vanilla - <<'RS'
            tb <- "{input.tarball}"
            pk <- "{SPIKEIN_BSG}"

            if (!file.exists(tb)) stop(paste("Tarball not found:", tb))
            if (!nzchar(pk)) stop("Config key 'spike_in_bsgenome' is empty")

            message("R_LIBS_USER = ", Sys.getenv("R_LIBS_USER"))
            message("Installing custom BSgenome from: ", tb)

            install.packages(tb, repos = NULL, type = "source")
            suppressPackageStartupMessages(library(pk, character.only = TRUE))
            message("OK: loaded ", pk)
RS

            echo "OK" > "{output.ok}"
            """


    # -------------------------------------------------------------------------
    # Extract Spike-in Reads from Deduplicated BAM
    # -------------------------------------------------------------------------
    rule extract_spikein_bam:
        """Extract spike-in reads from deduplicated BAM based on contig names.
        
        Requires config['spike_in_chr'] to be set with space-separated
        spike-in contig names (e.g., "spikein_1 spikein_2 spikein_3").
        """
        input:
            bam=get_dedup_bam
        output:
            bam=os.path.join(OUTPUT_DIR, "spikein_bam/{sample}_spikein.bam"),
            bai=os.path.join(OUTPUT_DIR, "spikein_bam/{sample}_spikein.bam.bai"),
            stats=os.path.join(OUTPUT_DIR, "stats/{sample}_spikein.stats.txt")
        params:
            spikein_chr=SPIKEIN_CHR
        threads: config.get("threads", 24)
        log:
            os.path.join(WORK_DIR, "logs/spikein/{sample}_extract.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_core.yaml")
        shell:
            r"""
            set -euo pipefail
            
            if [ -z "{params.spikein_chr}" ] || [ "{params.spikein_chr}" = "" ]; then
              echo "[ERROR] config['spike_in_chr'] is empty. Set spike-in contig names (e.g., 'spike1 spike2 ...')." >&2
              exit 1
            fi

            mkdir -p $(dirname {output.bam})

            samtools view -@ {threads} -hb "{input.bam}" {params.spikein_chr} | \
            samtools sort -@ {threads} -o "{output.bam}" - 2> {log}
            
            samtools index -@ {threads} "{output.bam}"
            samtools stats -@ {threads} "{output.bam}" > "{output.stats}"
            """


    # -------------------------------------------------------------------------
    # Insert Size Metrics for Spike-in
    # -------------------------------------------------------------------------
    rule insert_size_spikein:
        """Collect insert size metrics for spike-in BAM using Picard"""
        input:
            bam=os.path.join(OUTPUT_DIR, "spikein_bam/{sample}_spikein.bam")
        output:
            metrics=os.path.join(OUTPUT_DIR, "insert_size/{sample}_spikein_insert_size_metrics.txt"),
            hist=os.path.join(OUTPUT_DIR, "insert_size/{sample}_spikein_insert_size_histogram.pdf")
        log:
            os.path.join(WORK_DIR, "logs/picard/{sample}_insert_size_spikein.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_core.yaml")
        shell:
            """
            mkdir -p $(dirname {output.metrics})
            
            picard CollectInsertSizeMetrics \\
                M=0.05 \\
                I={input.bam} \\
                O={output.metrics} \\
                H={output.hist} \\
                2> {log}
            """


    # -------------------------------------------------------------------------
    # MEDIPS QC/Quantification for Spike-ins
    # -------------------------------------------------------------------------
    rule meth_qc_quant_spikein:
        """Run MEDIPS QC/quant on spike-in-only BAM using the custom spike-in BSgenome.
        
        Requires install_spikein_bsgenome rule to complete first.
        Uses the persistent R library installed in the workdir.
        """
        input:
            ok=SPIKEIN_OK,
            bam=os.path.join(OUTPUT_DIR, "spikein_bam/{sample}_spikein.bam")
        output:
            qc=os.path.join(OUTPUT_DIR, "meth_quant_spikein/{sample}_meth_qc.txt"),
            rdata=os.path.join(OUTPUT_DIR, "meth_quant_spikein/{sample}_meth_quant.RData"),
            granges=os.path.join(OUTPUT_DIR, "meth_quant_spikein/{sample}_Granges_CpGs.bed"),
            counts=os.path.join(OUTPUT_DIR, "meth_quant_spikein/{sample}_count.txt"),
            rpkm=os.path.join(OUTPUT_DIR, "meth_quant_spikein/{sample}_rpkm.txt")
        params:
            paired="TRUE" if PAIRED_END else "FALSE",
            window_size=str(WINDOW_SIZE),
            bsgenome=SPIKEIN_BSG,
            bsgenome_pkg=SPIKEIN_PKG,
            rlib=SPIKEIN_RLIB,
            script=os.path.join(PIPE_DIR, "scripts/core/medips_spikein.R"),
            outdir=os.path.join(OUTPUT_DIR, "meth_quant_spikein")
        resources:
            mem_mb=60000
        threads: 8
        log:
            os.path.join(WORK_DIR, "logs/meth_quant_spikein/{sample}.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_r.yaml")
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.qc})

            export R_LIBS_USER="$(pwd)/{params.rlib}"

            Rscript --vanilla "{params.script}" \
              "{wildcards.sample}" \
              "{input.bam}" \
              "{params.paired}" \
              "{params.window_size}" \
              "{params.bsgenome}" \
              "{params.bsgenome_pkg}" \
              "$(dirname {output.qc})" \
              2> "{log}"
            """


    # -------------------------------------------------------------------------
    # Export Spike-in Methylation Array to TSV
    # -------------------------------------------------------------------------
    rule export_spikein_meth_tsv:
        """Export spike-in methylation array RData to TSV format"""
        input:
            rdata=os.path.join(OUTPUT_DIR, "meth_quant_spikein/{sample}_meth_quant.RData")
        output:
            tsv=os.path.join(OUTPUT_DIR, "meth_quant_spikein/{sample}_meth_quant.tsv")
        params:
            rlib=SPIKEIN_RLIB
        conda:
            os.path.join(ENV_DIR, "medipipe_r.yaml")
        shell:
            r"""
            set -euo pipefail
            export R_LIBS_USER="$(pwd)/{params.rlib}"

            Rscript --vanilla - <<'RS'
            in_rdata <- "{input.rdata}"
            out_tsv  <- "{output.tsv}"

            load(in_rdata)  # expects object named 'meth'
            if (!exists("meth")) stop("Object 'meth' not found in RData: ", in_rdata)

            write.table(
              meth,
              file = out_tsv,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE
            )
RS
            """
