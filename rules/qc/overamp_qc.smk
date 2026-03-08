# =============================================================================
# MEDIPIPE QC Rules - Overamplification Detection
# =============================================================================

QC_CONFIG = config.get("qc", {})

# -----------------------------------------------------------------------------
# GC Bias Metrics
# -----------------------------------------------------------------------------
rule gc_bias:
    """Collect GC bias metrics using Picard"""
    input:
        bam=get_dedup_bam
    output:
        metrics=os.path.join(OUTPUT_DIR, "qc/gc_bias/{sample}_gc_bias.txt"),
        chart=os.path.join(OUTPUT_DIR, "qc/gc_bias/{sample}_gc_bias.pdf"),
        summary=os.path.join(OUTPUT_DIR, "qc/gc_bias/{sample}_gc_summary.txt")
    params:
        ref_fasta=REF.get("fasta", "")
    threads: 4
    log:
        os.path.join(WORK_DIR, "logs/qc/{sample}_gc_bias.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.metrics})
        
        picard CollectGcBiasMetrics \\
            I={input.bam} \\
            O={output.metrics} \\
            CHART={output.chart} \\
            S={output.summary} \\
            R={params.ref_fasta} \\
            VALIDATION_STRINGENCY=LENIENT \\
            2> {log}
        """


# -----------------------------------------------------------------------------
# Preseq Library Complexity
# -----------------------------------------------------------------------------
rule preseq_lc:
    """Estimate library complexity using preseq"""
    input:
        bam=get_dedup_bam
    output:
        curve=os.path.join(OUTPUT_DIR, "qc/preseq/{sample}_lc_extrap.txt")
    threads: 2
    log:
        os.path.join(WORK_DIR, "logs/qc/{sample}_preseq.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.curve})
        preseq lc_extrap -B {input.bam} -o {output.curve} 2> {log}
        """


# -----------------------------------------------------------------------------
# Fingerprint Plot
# -----------------------------------------------------------------------------
rule fingerprint:
    """Generate fingerprint plot using deepTools"""
    input:
        bam=get_dedup_bam
    output:
        plot=os.path.join(OUTPUT_DIR, "qc/fingerprint/{sample}_fingerprint.pdf"),
        metrics=os.path.join(OUTPUT_DIR, "qc/fingerprint/{sample}_fingerprint.txt")
    threads: 8
    log:
        os.path.join(WORK_DIR, "logs/qc/{sample}_fingerprint.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.plot})
        
        plotFingerprint \\
            --bamfiles {input.bam} \\
            --plotFile {output.plot} \\
            --outRawCounts {output.metrics} \\
            --plotFileFormat pdf \\
            -p {threads} \\
            2> {log}
        """


# -----------------------------------------------------------------------------
# Coverage Analysis
# -----------------------------------------------------------------------------
rule coverage_plot:
    """Generate coverage plot using deepTools"""
    input:
        bam=get_dedup_bam
    output:
        plot=os.path.join(OUTPUT_DIR, "qc/coverage/{sample}_coverage.pdf")
    threads: 8
    log:
        os.path.join(WORK_DIR, "logs/qc/{sample}_coverage.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.plot})
        
        plotCoverage \\
            -b {input.bam} \\
            -o {output.plot} \\
            --plotFileFormat pdf \\
            -p {threads} \\
            2> {log}
        """


# -----------------------------------------------------------------------------
# Overamplification Summary
# -----------------------------------------------------------------------------
rule overamp_summary:
    """Summarize overamplification QC metrics"""
    input:
        gc_bias=expand(os.path.join(OUTPUT_DIR, "qc/gc_bias/{sample}_gc_bias.txt"), sample=SAMPLE_IDS),
        preseq=expand(os.path.join(OUTPUT_DIR, "qc/preseq/{sample}_lc_extrap.txt"), sample=SAMPLE_IDS),
        insert_sizes=expand(os.path.join(OUTPUT_DIR, "insert_size/{sample}_insert_size_metrics.txt"), sample=SAMPLE_IDS)
    output:
        summary=os.path.join(OUTPUT_DIR, "qc/overamp_summary.tsv")
    params:
        script_dir=os.path.join(PIPE_DIR, "scripts/core")
    log:
        os.path.join(WORK_DIR, "logs/qc/overamp_summary.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_ml.yaml")
    shell:
        """
        python {params.script_dir}/summarize_overamp.py \\
            --gc-dir $(dirname {input.gc_bias[0]}) \\
            --preseq-dir $(dirname {input.preseq[0]}) \\
            --insert-dir $(dirname {input.insert_sizes[0]}) \\
            --output {output.summary} \\
            2> {log}
        """
