# =============================================================================
# MEDIPIPE Core Rules - Reads QC (UMI extraction, Trimming, FastQC)
# =============================================================================

# -----------------------------------------------------------------------------
# UMI Extraction (optional)
# Twist UMI Adapter system: extracts inline UMIs from R1+R2,
# removes UMI+skip bases from reads, writes UMIs to read IDs.
# Default: 5bp UMI + 2bp skip per read
# -----------------------------------------------------------------------------
if USE_UMI:
    rule umi_extract:
        """Extract inline UMIs from paired-end reads using UMI-tools
        
        Reads the UMI pattern from config:
        - umi_regex_r1: pattern for R1 (default: 5bp UMI + 2bp skip)
        - umi_regex_r2: pattern for R2 (default: 5bp UMI + 2bp skip)
        
        Pattern syntax: (?P<umi_N>.{X})(?P<discard_N>.{Y}).*
        - umi_N: captures X bases as UMI barcode
        - discard_N: removes Y bases after UMI
        - .*: matches remaining read sequence
        
        Output FASTQs have UMIs appended to read names (e.g. @READNAME_AACGT)
        """
        input:
            r1=get_raw_fastq_r1,
            r2=get_raw_fastq_r2 if PAIRED_END else []
        output:
            r1=temp(os.path.join(WORK_DIR, "umi_extracted/{sample}_R1.fastq.gz")),
            r2=temp(os.path.join(WORK_DIR, "umi_extracted/{sample}_R2.fastq.gz")) if PAIRED_END else [],
        params:
            # Use lambda to prevent Snakemake from interpreting {X} as wildcards
            pattern_r1=lambda wildcards: UMI_REGEX_R1,
            pattern_r2=lambda wildcards: UMI_REGEX_R2,
        threads: 4
        log:
            os.path.join(WORK_DIR, "logs/umi_extract/{sample}.log")
        conda:
            UMI_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.r1})
            
            umi_tools extract \
                --extract-method=regex \
                --bc-pattern='{params.pattern_r1}' \
                --bc-pattern2='{params.pattern_r2}' \
                --stdin={input.r1} \
                --read2-in={input.r2} \
                --stdout={output.r1} \
                --read2-out={output.r2} \
                --log={log}
            """

# -----------------------------------------------------------------------------
# Adapter Trimming
# -----------------------------------------------------------------------------
rule trim_galore_pe:
    """Trim adapters and low-quality bases using Trim Galore (paired-end)"""
    input:
        unpack(get_fastq_for_trim)
    output:
        r1=temp(os.path.join(WORK_DIR, "trimmed/{sample}_R1_val_1.fq.gz")),
        r2=temp(os.path.join(WORK_DIR, "trimmed/{sample}_R2_val_2.fq.gz")),
        report1=os.path.join(OUTPUT_DIR, "trim_reports/{sample}_R1.fastq.gz_trimming_report.txt"),
        report2=os.path.join(OUTPUT_DIR, "trim_reports/{sample}_R2.fastq.gz_trimming_report.txt")
    params:
        outdir=os.path.join(WORK_DIR, "trimmed"),
        report_dir=os.path.join(OUTPUT_DIR, "trim_reports")
    threads: 8
    log:
        os.path.join(WORK_DIR, "logs/trim_galore/{sample}.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p {params.outdir} {params.report_dir}
        
        trim_galore \\
            --paired \\
            --quality 20 \\
            --stringency 3 \\
            --length 20 \\
            --cores {threads} \\
            --output_dir {params.outdir} \\
            {input.r1} {input.r2} \\
            2> {log}
        
        # Move reports to output directory
        mv {params.outdir}/{wildcards.sample}_R1.fastq.gz_trimming_report.txt {output.report1} || true
        mv {params.outdir}/{wildcards.sample}_R2.fastq.gz_trimming_report.txt {output.report2} || true
        """


rule trim_galore_se:
    """Trim adapters and low-quality bases using Trim Galore (single-end)"""
    input:
        r1=get_raw_fastq_r1
    output:
        r1=temp(os.path.join(WORK_DIR, "trimmed/{sample}_trimmed.fq.gz")),
        report=os.path.join(OUTPUT_DIR, "trim_reports/{sample}.fastq.gz_trimming_report.txt")
    params:
        outdir=os.path.join(WORK_DIR, "trimmed")
    threads: 8
    log:
        os.path.join(WORK_DIR, "logs/trim_galore/{sample}.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p {params.outdir}
        
        trim_galore \\
            --quality 20 \\
            --stringency 3 \\
            --length 20 \\
            --cores {threads} \\
            --output_dir {params.outdir} \\
            {input.r1} \\
            2> {log}
        """


# -----------------------------------------------------------------------------
# FastQC
# -----------------------------------------------------------------------------
rule fastqc_raw:
    """Run FastQC on raw reads"""
    input:
        r1=get_raw_fastq_r1,
        r2=get_raw_fastq_r2 if PAIRED_END else []
    output:
        html1=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R1_fastqc.html"),
        zip1=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R1_fastqc.zip"),
        html2=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R2_fastqc.html") if PAIRED_END else [],
        zip2=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R2_fastqc.zip") if PAIRED_END else []
    params:
        outdir=os.path.join(OUTPUT_DIR, "fastqc")
    threads: 4
    log:
        os.path.join(WORK_DIR, "logs/fastqc/{sample}_raw.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} --outdir {params.outdir} {input.r1} {input.r2} 2> {log}
        
        # Rename FastQC outputs to expected filenames (FastQC uses input basename)
        R1_BASE=$(basename {input.r1} .fastq.gz)
        R1_BASE=${{R1_BASE%.fq.gz}}
        if [ -f "{params.outdir}/${{R1_BASE}}_fastqc.html" ]; then
            mv "{params.outdir}/${{R1_BASE}}_fastqc.html" {output.html1} 2>/dev/null || true
            mv "{params.outdir}/${{R1_BASE}}_fastqc.zip" {output.zip1} 2>/dev/null || true
        fi
        
        if [ -n "{input.r2}" ]; then
            R2_BASE=$(basename {input.r2} .fastq.gz)
            R2_BASE=${{R2_BASE%.fq.gz}}
            if [ -f "{params.outdir}/${{R2_BASE}}_fastqc.html" ]; then
                mv "{params.outdir}/${{R2_BASE}}_fastqc.html" {output.html2} 2>/dev/null || true
                mv "{params.outdir}/${{R2_BASE}}_fastqc.zip" {output.zip2} 2>/dev/null || true
            fi
        fi
        """


rule fastqc_trimmed:
    """Run FastQC on trimmed reads"""
    input:
        unpack(get_trimmed_fastq)
    output:
        html1=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R1_val_1_fastqc.html"),
        zip1=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R1_val_1_fastqc.zip"),
        html2=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R2_val_2_fastqc.html") if PAIRED_END else [],
        zip2=os.path.join(OUTPUT_DIR, "fastqc/{sample}_R2_val_2_fastqc.zip") if PAIRED_END else []
    params:
        outdir=os.path.join(OUTPUT_DIR, "fastqc")
    threads: 4
    log:
        os.path.join(WORK_DIR, "logs/fastqc/{sample}_trimmed.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} --outdir {params.outdir} {input.r1} {input.r2} 2> {log}
        """
