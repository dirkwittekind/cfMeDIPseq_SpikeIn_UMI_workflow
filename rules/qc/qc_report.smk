# =============================================================================
# MEDIPIPE QC Rules - Comprehensive QC Report Generation
# =============================================================================
# Automatically generates per-sample and summary QC reports after deduplication

import os
from pathlib import Path

# -----------------------------------------------------------------------------
# Comprehensive QC Report Generator
# Triggered after all samples have dedup stats
# -----------------------------------------------------------------------------
rule generate_qc_report:
    """Generate comprehensive QC report for all samples
    
    This rule collects metrics from:
    - Deduplicated BAM stats (read counts, insert size, quality)
    - Methylation QC (relH, GoGe)
    - Spike-in analysis (if available)
    - UMI stats (if available)
    
    Outputs:
    - qc_report.txt (human-readable)
    - qc_report.tsv (machine-readable)
    - qc_report.json (programmatic access)
    """
    input:
        stats=expand(os.path.join(OUTPUT_DIR, "dedup_bam_pe/{sample}_dedup.bam.stats.txt"), sample=SAMPLE_IDS)
    output:
        txt=os.path.join(OUTPUT_DIR, "qc_reports/qc_report.txt"),
        tsv=os.path.join(OUTPUT_DIR, "qc_reports/qc_report.tsv"),
        json=os.path.join(OUTPUT_DIR, "qc_reports/qc_report.json")
    params:
        output_dir=OUTPUT_DIR,
        spikein_dir=lambda wildcards: os.path.join(Path(OUTPUT_DIR).parent, "output_spikein") if USE_SPIKEIN else "",
        pipe_dir=PIPE_DIR
    log:
        os.path.join(WORK_DIR, "logs/qc/qc_report.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.txt})
        
        # Run QC report generator
        python3 {params.pipe_dir}/scripts/qc/generate_qc_report.py \
            --output-dir {params.output_dir} \
            $([ -n "{params.spikein_dir}" ] && [ -d "{params.spikein_dir}" ] && echo "--spikein-dir {params.spikein_dir}") \
            --report qc_reports/qc_report.txt \
            2>&1 | tee {log}
        """


# -----------------------------------------------------------------------------
# Per-Sample QC Report
# Generated individually for each sample after its dedup completes
# -----------------------------------------------------------------------------
rule sample_qc_report:
    """Generate individual QC report for a single sample
    
    This allows QC reports to be generated as soon as each sample
    finishes processing, without waiting for all samples.
    """
    input:
        stats=os.path.join(OUTPUT_DIR, "dedup_bam_pe/{sample}_dedup.bam.stats.txt")
    output:
        txt=os.path.join(OUTPUT_DIR, "qc_reports/per_sample/{sample}_qc.txt"),
        json=os.path.join(OUTPUT_DIR, "qc_reports/per_sample/{sample}_qc.json")
    params:
        output_dir=OUTPUT_DIR,
        spikein_dir=lambda wildcards: os.path.join(Path(OUTPUT_DIR).parent, "output_spikein") if USE_SPIKEIN else "",
        pipe_dir=PIPE_DIR,
        sample="{sample}"
    log:
        os.path.join(WORK_DIR, "logs/qc/qc_report_{sample}.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.txt})
        
        # Run QC report generator for single sample
        python3 {params.pipe_dir}/scripts/qc/generate_qc_report.py \
            --output-dir {params.output_dir} \
            $([ -n "{params.spikein_dir}" ] && [ -d "{params.spikein_dir}" ] && echo "--spikein-dir {params.spikein_dir}") \
            --samples {params.sample} \
            --report qc_reports/per_sample/{params.sample}_qc.txt \
            2>&1 | tee {log}
        """


# -----------------------------------------------------------------------------
# Helper function to include QC reports in workflow targets
# -----------------------------------------------------------------------------
def get_qc_report_targets():
    """Get QC report output targets"""
    targets = []
    
    # Summary report (generated after all samples)
    targets.append(os.path.join(OUTPUT_DIR, "qc_reports/qc_report.txt"))
    
    # Per-sample reports (generated incrementally)
    for sample in SAMPLE_IDS:
        targets.append(os.path.join(OUTPUT_DIR, f"qc_reports/per_sample/{sample}_qc.txt"))
    
    return targets
