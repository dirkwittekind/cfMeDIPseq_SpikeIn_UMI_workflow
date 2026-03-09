# =============================================================================
# MEDIPIPE Core Rules - Methylation QC and Quantification
# =============================================================================

# -----------------------------------------------------------------------------
# MEDIPS + MeDEStrand + QSEA Analysis
# -----------------------------------------------------------------------------
rule meth_qc_quant:
    """Run MEDIPS, MeDEStrand, and QSEA for methylation QC and quantification
    
    Outputs:
      - {sample}.medips.tsv: Methylation counts per window
      - {sample}.medips.rms.bedGraph: Relative methylation scores (bedGraph)
      - {sample}.summary.txt: QC summary with read counts, coverage, window size
    """
    input:
        bam=get_dedup_bam
    output:
        medips=os.path.join(OUTPUT_DIR, "meth_quant/{sample}.medips.tsv"),
        summary=os.path.join(OUTPUT_DIR, "meth_quant/{sample}.summary.txt"),
        bedgraph=os.path.join(OUTPUT_DIR, "meth_quant/{sample}.medips.rms.bedGraph")
    params:
        genome=GENOME,
        paired=PAIRED_END,
        # Use meth_qc_window_size (default 1000) to avoid integer overflow in R
        # Small windows (300bp) cause overflow: window_size * regions > 2^31
        window_size=config.get("meth_qc_window_size", max(1000, WINDOW_SIZE)),
        script_dir=os.path.join(PIPE_DIR, "scripts/core")
    resources:
        mem_mb=60000
    threads: 8
    log:
        os.path.join(WORK_DIR, "logs/meth_quant/{sample}.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_r.yaml")
    shell:
        """
        mkdir -p $(dirname {output.summary})
        
        Rscript --vanilla {params.script_dir}/medips_analysis.R \\
            --bam {input.bam} \\
            --output $(dirname {output.summary})/{wildcards.sample} \\
            --genome {params.genome} \\
            --window {params.window_size} \\
            --paired \\
            --skip_qsea \\
            2> {log}
        """


# NOTE: Spike-in quantification rules are now in spikein.smk
# They include BSgenome installation dependency and more complete outputs
