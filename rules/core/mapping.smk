# =============================================================================
# MEDIPIPE Core Rules - Mapping (BWA alignment, sorting, indexing)
# =============================================================================

# -----------------------------------------------------------------------------
# BWA Alignment
# -----------------------------------------------------------------------------
rule bwa_align:
    """Align reads to reference genome using BWA MEM"""
    input:
        unpack(get_trimmed_fastq)
    output:
        bam=temp(os.path.join(WORK_DIR, "aligned/{sample}.bam"))
    params:
        bwa_index=get_bwa_index()
    threads: config.get("threads", 24)
    log:
        os.path.join(WORK_DIR, "logs/bwa/{sample}.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        bwa mem -M -t {threads} \\
            {params.bwa_index} \\
            {input.r1} {input.r2} \\
            2> {log} | \\
        samtools view -Sb --threads {threads} - > {output.bam}
        """


# -----------------------------------------------------------------------------
# Sort and Index
# -----------------------------------------------------------------------------
rule samtools_sort:
    """Sort and fixmate BAM file"""
    input:
        bam=os.path.join(WORK_DIR, "aligned/{sample}.bam")
    output:
        bam=temp(os.path.join(WORK_DIR, "sorted/{sample}_sorted.bam")),
        stats=os.path.join(OUTPUT_DIR, "stats/{sample}_sorted.stats.txt")
    threads: config.get("threads", 24)
    log:
        os.path.join(WORK_DIR, "logs/samtools/{sample}_sort.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {output.stats})
        
        samtools fixmate -@ {threads} -m {input.bam} - | \\
        samtools sort -@ {threads} -o {output.bam} - 2> {log}
        
        samtools index -@ {threads} {output.bam}
        samtools stats -@ {threads} {output.bam} > {output.stats}
        """
