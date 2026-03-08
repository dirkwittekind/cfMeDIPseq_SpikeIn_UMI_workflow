# =============================================================================
# MEDIPIPE Core Workflow - Stages 1-5
# =============================================================================
# This Snakefile runs the core cfMeDIP-seq analysis pipeline:
#   Stage 1: Preprocessing (FASTQ → BAM)
#   Stage 2: Quality Control
#   Stage 3: Methylation Quantification
#   Stage 4: Fragment Profile Analysis
#   Stage 5: Tissue-of-Origin Analysis
# =============================================================================

import os
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Determine rule paths
PIPE_DIR = config.get("pipe_dir", "/home/dirk/medipipe_warp")
RULES_DIR = os.path.join(PIPE_DIR, "rules")

# Include common functions and variables
include: os.path.join(RULES_DIR, "core/common.smk")

# Include core workflow rules
include: os.path.join(RULES_DIR, "core/reads_qc.smk")
include: os.path.join(RULES_DIR, "core/mapping.smk")
include: os.path.join(RULES_DIR, "core/deduplication.smk")
include: os.path.join(RULES_DIR, "core/meth_qc_quant.smk")
include: os.path.join(RULES_DIR, "core/fragment_profile.smk")
include: os.path.join(RULES_DIR, "core/tissue_of_origin.smk")

# Include spike-in rules (only active when USE_SPIKEIN=True)
include: os.path.join(RULES_DIR, "core/spikein.smk")

# Include QC rules
include: os.path.join(RULES_DIR, "qc/overamp_qc.smk")
include: os.path.join(RULES_DIR, "qc/qc_report.smk")

# Include differential methylation analysis
include: os.path.join(RULES_DIR, "core/differential_methylation.smk")

# Include ML discrimination workflow (advanced)
include: os.path.join(RULES_DIR, "advanced/ml_discrimination.smk")


# =============================================================================
# Main Target Rule
# =============================================================================
rule all:
    """Main target: run all enabled core workflow stages"""
    input:
        get_core_workflow_targets()


# =============================================================================
# Stage-specific Target Rules (for running individual stages)
# =============================================================================
rule stage1_preprocessing:
    """Stage 1: Generate deduplicated BAM files from FASTQ"""
    input:
        expand(os.path.join(OUTPUT_DIR, "dedup_bam/{sample}_dedup.bam"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "dedup_bam/{sample}_dedup.bam.bai"), sample=SAMPLE_IDS)


rule stage2_qc:
    """Stage 2: Run QC analyses"""
    input:
        expand(os.path.join(OUTPUT_DIR, "fastqc/{sample}_R1_fastqc.zip"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "fastqc/{sample}_R2_fastqc.zip"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "stats/{sample}_dedup.stats.txt"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "insert_size/{sample}_insert_size_metrics.txt"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "qc/gc_bias/{sample}_gc_bias.txt"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "qc/preseq/{sample}_lc_extrap.txt"), sample=SAMPLE_IDS)


rule stage3_methylation:
    """Stage 3: Methylation quantification"""
    input:
        expand(os.path.join(OUTPUT_DIR, "meth_quant/{sample}_meth_qc.txt"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "meth_quant/{sample}_count.txt"), sample=SAMPLE_IDS)


rule stage4_fragment:
    """Stage 4: Fragment profile analysis"""
    input:
        expand(os.path.join(OUTPUT_DIR, "fragment_profile/{sample}_fragment_profile.txt"), sample=SAMPLE_IDS)


rule stage5_too:
    """Stage 5: Tissue-of-origin analysis"""
    input:
        os.path.join(OUTPUT_DIR, "tissue_of_origin/too_fractions.tsv")


rule differential_analysis:
    """Run differential methylation analysis between groups"""
    input:
        get_differential_targets() if 'get_differential_targets' in dir() else []


# =============================================================================
# Utility Rules
# =============================================================================
rule clean_temp:
    """Remove temporary files from work directory"""
    shell:
        """
        echo "Cleaning temporary files in {WORK_DIR}..."
        rm -rf {WORK_DIR}/aligned
        rm -rf {WORK_DIR}/sorted
        rm -rf {WORK_DIR}/trimmed
        rm -rf {WORK_DIR}/umi_extracted
        echo "Temporary files cleaned."
        """


rule multiqc:
    """Generate MultiQC report"""
    input:
        expand(os.path.join(OUTPUT_DIR, "fastqc/{sample}_R1_fastqc.zip"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "stats/{sample}_dedup.stats.txt"), sample=SAMPLE_IDS),
        expand(os.path.join(OUTPUT_DIR, "insert_size/{sample}_insert_size_metrics.txt"), sample=SAMPLE_IDS)
    output:
        report=os.path.join(OUTPUT_DIR, "multiqc/multiqc_report.html")
    params:
        outdir=os.path.join(OUTPUT_DIR, "multiqc"),
        indir=OUTPUT_DIR
    log:
        os.path.join(WORK_DIR, "logs/multiqc/multiqc.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p {params.outdir}
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
