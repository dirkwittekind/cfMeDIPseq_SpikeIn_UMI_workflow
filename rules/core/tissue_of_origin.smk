# =============================================================================
# MEDIPIPE Core Rules - Tissue-of-Origin Analysis
# =============================================================================

TOO_CONFIG = config.get("tissue_of_origin", {})
# Fixed: Coordinate-to-marker_id conversion and cfMeDIP-specific scaling implemented
# Based on Loyfer et al. 2023 (Nature) methodology
TOO_ENABLE = TOO_CONFIG.get("enable", True)

if TOO_ENABLE:
    # -------------------------------------------------------------------------
    # Generate BigWig files for ToO analysis
    # -------------------------------------------------------------------------
    rule too_bigwig:
        """Generate BigWig files for tissue-of-origin analysis"""
        input:
            bam=get_dedup_bam
        output:
            bw=os.path.join(OUTPUT_DIR, "bigwig/{sample}.bw")
        params:
            bin_size=config.get("qc", {}).get("bigwig_bin_size", 25),
            norm=config.get("qc", {}).get("bigwig_normalization", "CPM")
        threads: 8
        log:
            os.path.join(WORK_DIR, "logs/bigwig/{sample}.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_core.yaml")
        shell:
            """
            mkdir -p $(dirname {output.bw})
            
            bamCoverage \\
                -b {input.bam} \\
                -o {output.bw} \\
                --binSize {params.bin_size} \\
                --normalizeUsing {params.norm} \\
                -p {threads} \\
                2> {log}
            """

    # -------------------------------------------------------------------------
    # Extract marker region coverage
    # -------------------------------------------------------------------------
    rule too_marker_matrix:
        """Extract coverage at tissue marker regions"""
        input:
            bigwigs=expand(os.path.join(OUTPUT_DIR, "bigwig/{sample}.bw"), sample=SAMPLE_IDS),
            markers=TOO_CONFIG.get("markers_bed", os.path.join(RESOURCES_DIR, "tissue_of_origin/marker_regions.bed"))
        output:
            matrix=os.path.join(OUTPUT_DIR, "tissue_of_origin/marker_matrix.tsv"),
            npz=os.path.join(OUTPUT_DIR, "tissue_of_origin/marker_matrix.npz")
        params:
            script_dir=os.path.join(PIPE_DIR, "scripts/core")
        threads: 8
        log:
            os.path.join(WORK_DIR, "logs/too/marker_matrix.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_core.yaml")
        shell:
            """
            mkdir -p $(dirname {output.matrix})
            
            # Use multiBigwigSummary for marker extraction
            multiBigwigSummary BED-file \\
                --bwfiles {input.bigwigs} \\
                --BED {input.markers} \\
                --outFileName {output.npz} \\
                --outRawCounts {output.matrix} \\
                -p {threads} \\
                2> {log}
            """

    # -------------------------------------------------------------------------
    # NNLS Deconvolution (Loyfer et al. 2023 methodology)
    # -------------------------------------------------------------------------
    rule too_nnls_deconv:
        """Perform NNLS deconvolution to estimate tissue fractions using Loyfer et al. 2023 markers"""
        input:
            marker_matrix=os.path.join(OUTPUT_DIR, "tissue_of_origin/marker_matrix.tsv"),
            signature=TOO_CONFIG.get("signature_tsv", os.path.join(RESOURCES_DIR, "tissue_of_origin/reference_methylation_matrix.tsv")),
            markers_bed=TOO_CONFIG.get("markers_bed", os.path.join(RESOURCES_DIR, "tissue_of_origin/marker_regions.bed"))
        output:
            fractions=os.path.join(OUTPUT_DIR, "tissue_of_origin/too_fractions.tsv")
        params:
            scale_mode=TOO_CONFIG.get("scale", "cfmedip"),
            script_dir=os.path.join(PIPE_DIR, "scripts/core")
        log:
            os.path.join(WORK_DIR, "logs/too/nnls_deconv.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script_dir}/nnls_deconv.py \\
                -s {input.marker_matrix} \\
                -r {input.signature} \\
                -m {input.markers_bed} \\
                -o {output.fractions} \\
                --scale {params.scale_mode} \\
                2> {log}
            """
