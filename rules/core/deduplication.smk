# =============================================================================
# MEDIPIPE Core Rules - Deduplication
# =============================================================================

# -----------------------------------------------------------------------------
# Standard Deduplication (samtools markdup)
# -----------------------------------------------------------------------------
rule samtools_markdup:
    """Remove duplicates using samtools markdup"""
    input:
        bam=os.path.join(WORK_DIR, "sorted/{sample}_sorted.bam")
    output:
        bam=os.path.join(OUTPUT_DIR, "dedup_bam/{sample}_dedup.bam"),
        bai=os.path.join(OUTPUT_DIR, "dedup_bam/{sample}_dedup.bam.bai"),
        stats=os.path.join(OUTPUT_DIR, "stats/{sample}_dedup.stats.txt")
    threads: config.get("threads", 24)
    log:
        os.path.join(WORK_DIR, "logs/dedup/{sample}.log")
    conda:
        os.path.join(ENV_DIR, "medipipe_core.yaml")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        # Filter: properly paired (-f 2), exclude unmapped/secondary/qcfail/dup (-F 2828)
        samtools view -b -f 2 -F 2828 --threads {threads} {input.bam} | \\
        samtools markdup -@ {threads} -r - {output.bam} 2> {log}
        
        samtools index -@ {threads} {output.bam}
        samtools stats -@ {threads} {output.bam} > {output.stats}
        """


# -----------------------------------------------------------------------------
# UMI-based Deduplication (UMI-tools)
# Uses UMI + genomic position for duplicate detection (paired-end).
# -----------------------------------------------------------------------------
if USE_UMI:
    rule umi_dedup:
        """UMI-aware deduplication on coordinate-sorted BAM (paired-end).
        
        Uses UMI-tools dedup with --paired mode. The UMI is read from the
        read name (appended during umi_extract step, separated by '_').
        
        Produces:
        - Deduplicated BAM + index
        - Deduplication stats (TSV)
        - samtools stats on final BAM
        """
        input:
            bam=os.path.join(WORK_DIR, "sorted/{sample}_sorted.bam")
        output:
            bam=os.path.join(OUTPUT_DIR, "dedup_bam_umi/{sample}_dedup.bam"),
            bai=os.path.join(OUTPUT_DIR, "dedup_bam_umi/{sample}_dedup.bam.bai"),
            stats=os.path.join(OUTPUT_DIR, "stats/{sample}_dedup_umi.stats.txt"),
            umi_stats=os.path.join(OUTPUT_DIR, "stats/{sample}_umi_dedup_stats.tsv")
        params:
            stats_prefix=os.path.join(OUTPUT_DIR, "stats/{sample}_umi_dedup")
        threads: config.get("threads", 24)
        log:
            os.path.join(WORK_DIR, "logs/umi_dedup/{sample}.log")
        conda:
            UMI_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bam})
            
            # UMI-tools dedup (paired-end)
            # --umi-separator='_' matches the output from umi_tools extract
            umi_tools dedup \
                --paired \
                -I {input.bam} \
                -S {output.bam} \
                --umi-separator='_' \
                --output-stats={params.stats_prefix} \
                > {log} 2>&1
            
            # Normalize stats file name (umi_tools outputs *_stats.tsv or *.tsv)
            if [ -f "{params.stats_prefix}.tsv" ]; then
                mv -f "{params.stats_prefix}.tsv" "{output.umi_stats}"
            elif [ -f "{params.stats_prefix}_stats.tsv" ]; then
                mv -f "{params.stats_prefix}_stats.tsv" "{output.umi_stats}"
            else
                # Create empty stats file if none found
                : > "{output.umi_stats}"
            fi
            
            # Index and collect stats
            samtools index -@ {threads} {output.bam} {output.bai}
            samtools stats -@ {threads} {output.bam} > {output.stats}
            """


# -----------------------------------------------------------------------------
# Insert Size Metrics
# -----------------------------------------------------------------------------
rule insert_size:
    """Collect insert size metrics using Picard"""
    input:
        bam=get_dedup_bam
    output:
        metrics=os.path.join(OUTPUT_DIR, "insert_size/{sample}_insert_size_metrics.txt"),
        hist=os.path.join(OUTPUT_DIR, "insert_size/{sample}_insert_size_histogram.pdf")
    log:
        os.path.join(WORK_DIR, "logs/picard/{sample}_insert_size.log")
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
