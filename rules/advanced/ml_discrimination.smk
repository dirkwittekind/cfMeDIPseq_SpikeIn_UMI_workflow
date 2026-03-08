# =============================================================================
# MEDIPIPE Advanced Rules - ML Discrimination Workflow
# Full workflow for cfMeDIP-seq group discrimination
# Based on standard workflow: BigWigs -> Windows -> Matrix -> ML/DMR
# =============================================================================

ML_DISCRIM_CONFIG = config.get("ml_discrimination", {})
ADVANCED_CONFIG = config.get("advanced", {})

# Default settings
WINDOW_SIZE = ML_DISCRIM_CONFIG.get("window_size", 2000)
KBEST = ML_DISCRIM_CONFIG.get("kbest", 5000)
C_PARAM = ML_DISCRIM_CONFIG.get("C", 1.0)
RSKFOLD_SPLITS = ML_DISCRIM_CONFIG.get("rskfold_splits", 5)
RSKFOLD_REPEATS = ML_DISCRIM_CONFIG.get("rskfold_repeats", 50)
ML_SEED = ML_DISCRIM_CONFIG.get("seed", 1)
ML_MODE = ML_DISCRIM_CONFIG.get("mode", "methylOnly")
COVARIATES = ML_DISCRIM_CONFIG.get("covariates", [])
CASE_GROUP = config.get("comparison", {}).get("group1", "AEG")
CTRL_GROUP = config.get("comparison", {}).get("group2", "CTRL")
AUTOSOMAL_ONLY = ML_DISCRIM_CONFIG.get("autosomal_only", True)  # Default: use only chr1-22

# Input sources
INPUT_SOURCE = ML_DISCRIM_CONFIG.get("input_source", "bigwig")  # "bigwig" or "dedup_bam"

# Directories
ML_BASE = os.path.join(OUTPUT_DIR, "ml_discrimination")
BIGWIG_DIR = os.path.join(OUTPUT_DIR, "bigwig")
DEDUP_BAM_DIR = os.path.join(OUTPUT_DIR, "dedup_bam_pe")
WINDOWS_DIR = os.path.join(ML_BASE, "resources/windows")
MATRIX_DIR = os.path.join(ML_BASE, "matrices")
ML_OUTPUT_DIR = os.path.join(ML_BASE, f"Outputs/{CASE_GROUP}_vs_{CTRL_GROUP}")

# Genome resources
GENOME = config.get("genome", "hg38")
GENOME_FILE = config.get("ref_files", {}).get("chrom_sizes", 
    os.path.join(RESOURCES_DIR, f"genomes/{GENOME}/{GENOME}.chrom.sizes"))
PROMOTER_BED = ML_DISCRIM_CONFIG.get("promoter_bed",
    os.path.join(RESOURCES_DIR, "annotations/hg38_promoters_2kb.bed"))

if ADVANCED_CONFIG.get("ml_classification", False):
    
    # =========================================================================
    # 1. Generate BigWigs from Dedup BAMs (if needed)
    # =========================================================================
    rule ml_generate_bigwig:
        """Generate BigWig files from deduplicated BAMs using deepTools bamCoverage"""
        input:
            bam=os.path.join(DEDUP_BAM_DIR, "{sample}.bam"),
            bai=os.path.join(DEDUP_BAM_DIR, "{sample}.bam.bai")
        output:
            bw=os.path.join(BIGWIG_DIR, "{sample}.bw")
        params:
            bin_size=config.get("qc", {}).get("bigwig_bin_size", 10),
            norm=config.get("qc", {}).get("bigwig_normalization", "CPM"),
            eff_genome_size=ML_DISCRIM_CONFIG.get("effective_genome_size", 2913022398)
        threads: 12
        log:
            os.path.join(WORK_DIR, "logs/ml/bigwig/{sample}.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            bamCoverage \
                -b {input.bam} \
                -o {output.bw} \
                --binSize {params.bin_size} \
                --normalizeUsing {params.norm} \
                --effectiveGenomeSize {params.eff_genome_size} \
                --ignoreDuplicates \
                --numberOfProcessors {threads} \
                2> {log}
            """
    
    # =========================================================================
    # 2. Create Genome Windows BED
    # =========================================================================
    rule ml_create_windows:
        """Create genome-wide windows BED file (with autosomal-only option)"""
        input:
            genome=GENOME_FILE
        output:
            all_windows=os.path.join(WINDOWS_DIR, f"{GENOME}_w{WINDOW_SIZE}.all.bed"),
            autosomal=os.path.join(WINDOWS_DIR, f"{GENOME}_w{WINDOW_SIZE}.autosomal.bed"),
            nonpromoter=os.path.join(WINDOWS_DIR, f"{GENOME}_w{WINDOW_SIZE}.nonpromoter.bed")
        params:
            window_size=WINDOW_SIZE,
            promoter_bed=PROMOTER_BED
        log:
            os.path.join(WORK_DIR, "logs/ml/create_windows.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            mkdir -p {WINDOWS_DIR}
            
            # Create all windows (standard chromosomes only)
            bedtools makewindows -g {input.genome} -w {params.window_size} \
                | awk 'BEGIN{{OFS="\\t"}} $1 ~ /^chr([0-9]+|X|Y)$/ {{print $0}}' \
                > {output.all_windows} 2> {log}
            
            # Create autosomal-only windows (chr1-22, exclude chrX/chrY/chrM)
            awk 'BEGIN{{OFS="\\t"}} $1 ~ /^chr[0-9]+$/ {{print $0}}' {output.all_windows} \
                > {output.autosomal} 2>> {log}
            
            # Create nonpromoter windows (if promoter BED exists)
            if [ -f "{params.promoter_bed}" ]; then
                bedtools subtract \
                    -a {output.all_windows} \
                    -b {params.promoter_bed} \
                    > {output.nonpromoter} 2>> {log}
            else
                cp {output.all_windows} {output.nonpromoter}
            fi
            
            echo "Created windows: $(wc -l < {output.all_windows}) total, $(wc -l < {output.autosomal}) autosomal, $(wc -l < {output.nonpromoter}) nonpromoter" >> {log}
            """
    
    # =========================================================================
    # 3. Build Feature Matrix from BigWigs (autosomal-only by default)
    # =========================================================================
    # Determine which window file to use based on AUTOSOMAL_ONLY setting
    WINDOWS_BED = os.path.join(WINDOWS_DIR, f"{GENOME}_w{WINDOW_SIZE}.autosomal.bed") if AUTOSOMAL_ONLY else os.path.join(WINDOWS_DIR, f"{GENOME}_w{WINDOW_SIZE}.all.bed")
    MATRIX_SUFFIX = "autosomal" if AUTOSOMAL_ONLY else "all"
    
    rule ml_build_matrix:
        """Build feature matrix from BigWigs using multiBigwigSummary (autosomal-only by default)"""
        input:
            bigwigs=expand(os.path.join(BIGWIG_DIR, "{sample}.bw"), sample=SAMPLE_IDS),
            windows=WINDOWS_BED
        output:
            npz=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.npz"),
            tsv=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv")
        params:
            bigwig_list=lambda wildcards, input: " ".join(input.bigwigs),
            autosomal_only=AUTOSOMAL_ONLY
        threads: 12
        log:
            os.path.join(WORK_DIR, "logs/ml/build_matrix.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            mkdir -p $(dirname {output.npz})
            
            echo "Building matrix with autosomal_only={params.autosomal_only}" > {log}
            echo "Using windows: {input.windows}" >> {log}
            
            multiBigwigSummary BED-file \
                --bwfiles {params.bigwig_list} \
                --BED {input.windows} \
                --outFileName {output.npz} \
                --outRawCounts {output.tsv} \
                --numberOfProcessors {threads} \
                2>> {log}
            
            echo "Matrix shape: $(head -1 {output.tsv} | awk -F'\\t' '{{print NF}}') columns" >> {log}
            echo "Matrix rows: $(wc -l < {output.tsv})" >> {log}
            """
    
    # =========================================================================
    # 4. ML Cross-Validation Discrimination
    # =========================================================================
    rule ml_discrimination_cv:
        """Run ML CV pipeline for group discrimination (autosomal-only by default)"""
        input:
            matrix=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv"),
            ann=config.get("samples", "")
        output:
            summary=os.path.join(ML_OUTPUT_DIR, f"windows{WINDOW_SIZE}_{MATRIX_SUFFIX}__{ML_MODE}/ml_cv.summary.tsv"),
            curves=os.path.join(ML_OUTPUT_DIR, f"windows{WINDOW_SIZE}_{MATRIX_SUFFIX}__{ML_MODE}/ml_cv.curves.png"),
            preds_loo=os.path.join(ML_OUTPUT_DIR, f"windows{WINDOW_SIZE}_{MATRIX_SUFFIX}__{ML_MODE}/ml_cv.preds.LOO.tsv"),
            preds_rskf=os.path.join(ML_OUTPUT_DIR, f"windows{WINDOW_SIZE}_{MATRIX_SUFFIX}__{ML_MODE}/ml_cv.preds.RSKFold_perSampleMean.tsv"),
            importance=os.path.join(ML_OUTPUT_DIR, f"windows{WINDOW_SIZE}_{MATRIX_SUFFIX}__{ML_MODE}/ml_cv.feature_importance.tsv")
        params:
            script=os.path.join(PIPE_DIR, "scripts/ml/ml_cv_discrimination.py"),
            out_dir=os.path.join(ML_OUTPUT_DIR, f"windows{WINDOW_SIZE}_{MATRIX_SUFFIX}__{ML_MODE}"),
            kbest=KBEST,
            C=C_PARAM,
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP,
            seed=ML_SEED,
            mode=ML_MODE,
            covariates=" ".join(COVARIATES) if COVARIATES else ""
        log:
            os.path.join(WORK_DIR, "logs/ml/ml_cv.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --ann {input.ann} \
                --matrix {input.matrix} \
                --out-dir {params.out_dir} \
                --prefix ml_cv \
                --case-group {params.case_group} \
                --ctrl-group {params.ctrl_group} \
                --kbest {params.kbest} \
                --C {params.C} \
                --seed {params.seed} \
                --mode {params.mode} \
                --no-hypo-hyper \
                2> {log}
            """
    
    # =========================================================================
    # 5. Univariate Feature Driver Analysis (autosomal-only by default)
    # =========================================================================
    rule ml_univariate_drivers:
        """Rank features by univariate effect size (autosomal-only = no chrX/Y artifacts)"""
        input:
            matrix=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv"),
            ann=config.get("samples", "")
        output:
            top_features=os.path.join(ML_OUTPUT_DIR, "feature_drivers/univariate_top_windows.tsv")
        params:
            script=os.path.join(PIPE_DIR, "scripts/ml/rank_windows_univariate.py"),
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP,
            gene_bed=ML_DISCRIM_CONFIG.get("gene_bed", "")
        log:
            os.path.join(WORK_DIR, "logs/ml/univariate_drivers.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --ann {input.ann} \
                --matrix {input.matrix} \
                --out {output.top_features} \
                --case-group {params.case_group} \
                --ctrl-group {params.ctrl_group} \
                --top-n 2000 \
                --gene-bed {params.gene_bed} \
                2> {log}
            """
    
    # =========================================================================
    # 6. ML Robustness Checks
    # =========================================================================
    rule ml_robustness_checks:
        """Run robustness checks: depth confounders, permutation test, feature driver analysis"""
        input:
            matrix=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv"),
            ann=config.get("samples", ""),
            drivers=os.path.join(ML_OUTPUT_DIR, "feature_drivers/univariate_top_windows.tsv")
        output:
            robustness=os.path.join(ML_OUTPUT_DIR, "robustness/robustness_checks.tsv"),
            depth_metrics=os.path.join(ML_OUTPUT_DIR, "robustness/sample_depth_metrics.tsv"),
            robustness_json=os.path.join(ML_OUTPUT_DIR, "robustness/robustness_checks.json")
        params:
            script=os.path.join(PIPE_DIR, "scripts/ml/ml_robustness_checks.py"),
            out_dir=os.path.join(ML_OUTPUT_DIR, "robustness"),
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP,
            n_permutations=ML_DISCRIM_CONFIG.get("n_permutations", 1000)
        log:
            os.path.join(WORK_DIR, "logs/ml/robustness_checks.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --matrix {input.matrix} \
                --annotation {input.ann} \
                --drivers {input.drivers} \
                --case-group {params.case_group} \
                --control-group {params.ctrl_group} \
                --output {params.out_dir} \
                --n-permutations {params.n_permutations} \
                2> {log}
            """
    
    # =========================================================================
    # 7. Exact Permutation Test (confirmatory statistical validation)
    # =========================================================================
    EXACT_PERM_TOP_K = ML_DISCRIM_CONFIG.get("exact_perm_top_k", 100)
    EXACT_PERM_C = ML_DISCRIM_CONFIG.get("exact_perm_C", 0.1)
    
    rule ml_exact_permutation:
        """Run exact permutation test with strict leakage control (462 labelings for n=11)"""
        input:
            matrix=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv"),
            ann=config.get("samples", "")
        output:
            summary=os.path.join(ML_OUTPUT_DIR, f"exact_permutation/autosomal_windows.summary.tsv"),
            null_dist=os.path.join(ML_OUTPUT_DIR, f"exact_permutation/autosomal_windows.exact_null.tsv"),
            stability=os.path.join(ML_OUTPUT_DIR, f"exact_permutation/autosomal_windows.sample_stability.tsv"),
            p_values=os.path.join(ML_OUTPUT_DIR, f"exact_permutation/autosomal_windows.exact_p_values.json"),
            auroc_plot=os.path.join(ML_OUTPUT_DIR, f"exact_permutation/autosomal_windows.permutation_auroc.pdf")
        params:
            script=os.path.join(PIPE_DIR, "scripts/ml/ml_exact_permutation.py"),
            out_dir=os.path.join(ML_OUTPUT_DIR, "exact_permutation"),
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP,
            top_k=EXACT_PERM_TOP_K,
            C=EXACT_PERM_C
        log:
            os.path.join(WORK_DIR, "logs/ml/exact_permutation.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --matrix {input.matrix} \
                --annotation {input.ann} \
                --output {params.out_dir} \
                --feature-space autosomal_windows \
                --case-group {params.case_group} \
                --ctrl-group {params.ctrl_group} \
                --top-k {params.top_k} \
                --C {params.C} \
                2> {log}
            """
    
    # =========================================================================
    # 8. DMR Analysis with limma (autosomal-only by default)
    # =========================================================================
    rule ml_dmr_analysis:
        """Run DMR analysis using limma (autosomal-only = no sex chr confounders)"""
        input:
            matrix=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv"),
            ann=config.get("samples", "")
        output:
            dmr_gz=os.path.join(ML_OUTPUT_DIR, f"DMR/windows{WINDOW_SIZE}_{MATRIX_SUFFIX}/dmr.tsv.gz"),
            dmr_top=os.path.join(ML_OUTPUT_DIR, f"DMR/windows{WINDOW_SIZE}_{MATRIX_SUFFIX}/dmr.top1000.tsv"),
            volcano=os.path.join(ML_OUTPUT_DIR, f"DMR/windows{WINDOW_SIZE}_{MATRIX_SUFFIX}/dmr.volcano.png"),
            summary=os.path.join(ML_OUTPUT_DIR, f"DMR/windows{WINDOW_SIZE}_{MATRIX_SUFFIX}/dmr.summary.tsv")
        params:
            script=os.path.join(PIPE_DIR, "scripts/ml/dmr_windows_limma.R"),
            out_dir=os.path.join(ML_OUTPUT_DIR, f"DMR/windows{WINDOW_SIZE}_{MATRIX_SUFFIX}"),
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP,
            covariates=" ".join(COVARIATES) if COVARIATES else ""
        log:
            os.path.join(WORK_DIR, "logs/ml/dmr_analysis.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_r.yaml")
        shell:
            """
            Rscript {params.script} \
                --ann {input.ann} \
                --matrix {input.matrix} \
                --out-dir {params.out_dir} \
                --prefix dmr \
                --case-group {params.case_group} \
                --ctrl-group {params.ctrl_group} \
                --top-n 1000 \
                2> {log}
            """
    
    # =========================================================================
    # 9. Exact Permutation with Dual Model (L2 + Elastic-Net)
    # =========================================================================
    rule ml_exact_permutation_dual:
        """Run exact permutation test with BOTH L2 and Elastic-Net models (parallelized)"""
        input:
            matrix=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv"),
            ann=config.get("samples", "")
        output:
            l2_summary=os.path.join(ML_OUTPUT_DIR, f"exact_permutation/autosomal_windows.L2.summary.tsv"),
            enet_summary=os.path.join(ML_OUTPUT_DIR, f"exact_permutation/autosomal_windows.ElasticNet.summary.tsv")
        params:
            script=os.path.join(PIPE_DIR, "scripts/ml/ml_exact_permutation.py"),
            out_dir=os.path.join(ML_OUTPUT_DIR, "exact_permutation"),
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP,
            top_k=EXACT_PERM_TOP_K,
            C=EXACT_PERM_C
        threads: 40  # Use most cores for parallelized permutation
        log:
            os.path.join(WORK_DIR, "logs/ml/exact_permutation_dual.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --matrix {input.matrix} \
                --annotation {input.ann} \
                --output {params.out_dir} \
                --feature-space autosomal_windows \
                --case-group {params.case_group} \
                --ctrl-group {params.ctrl_group} \
                --top-k {params.top_k} \
                --C {params.C} \
                --model-type both \
                --n-jobs {threads} \
                2> {log}
            """
    
    # =========================================================================
    # 10. Per-Sample QC Reports
    # =========================================================================
    rule generate_sample_qc_report:
        """Generate per-sample QC report with all metrics"""
        input:
            # Ensure basic QC outputs exist
            bam_stats=os.path.join(OUTPUT_DIR, "dedup_bam_pe/{sample}.bam.stats.txt") if os.path.exists(os.path.join(OUTPUT_DIR, "dedup_bam_pe")) else []
        output:
            html=os.path.join(ML_BASE, "reports/samples/{sample}_qc_report.html"),
            tsv=os.path.join(ML_BASE, "reports/samples/{sample}_qc_report.tsv")
        params:
            script=os.path.join(PIPE_DIR, "scripts/reporting/generate_sample_report.py"),
            report_dir=os.path.join(ML_BASE, "reports/samples")
        log:
            os.path.join(WORK_DIR, "logs/reports/{sample}_qc_report.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --sample-id {wildcards.sample} \
                --output-dir {OUTPUT_DIR} \
                --report-dir {params.report_dir} \
                --format both \
                2> {log}
            """
    
    # =========================================================================
    # 11. ML Analysis Report with Interpretation
    # =========================================================================
    rule generate_ml_report:
        """Generate comprehensive ML discrimination report with interpretation"""
        input:
            exact_perm=os.path.join(ML_OUTPUT_DIR, "exact_permutation/autosomal_windows.summary.tsv"),
            robustness=os.path.join(ML_OUTPUT_DIR, "robustness/robustness_checks.tsv"),
            univariate=os.path.join(ML_OUTPUT_DIR, "feature_drivers/univariate_top_windows.tsv")
        output:
            html=os.path.join(ML_BASE, f"reports/ml_report_{CASE_GROUP}_vs_{CTRL_GROUP}.html"),
            json=os.path.join(ML_BASE, f"reports/ml_report_{CASE_GROUP}_vs_{CTRL_GROUP}.json")
        params:
            script=os.path.join(PIPE_DIR, "scripts/reporting/generate_ml_report.py"),
            ml_output_dir=ML_OUTPUT_DIR,
            report_dir=os.path.join(ML_BASE, "reports"),
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP
        log:
            os.path.join(WORK_DIR, "logs/reports/ml_report.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --ml-output-dir {params.ml_output_dir} \
                --report-dir {params.report_dir} \
                --case-group {params.case_group} \
                --ctrl-group {params.ctrl_group} \
                --format all \
                2> {log}
            """
    
    # =========================================================================
    # 12. ML Visualization (Volcano Plots, Driver Heatmaps)
    # =========================================================================
    ML_VIZ_TOP_N_HEATMAP = ML_DISCRIM_CONFIG.get("top_n_heatmap", 20)
    ML_VIZ_TOP_N_BARPLOT = ML_DISCRIM_CONFIG.get("top_n_barplot", 30)
    ELASTICNET_ALPHA = ML_DISCRIM_CONFIG.get("elasticnet_alpha", 0.01)
    
    rule ml_visualization:
        """Generate ML visualizations: volcano plots, driver heatmaps, feature importance"""
        input:
            matrix=os.path.join(MATRIX_DIR, f"windows/{GENOME}_w{WINDOW_SIZE}/{MATRIX_SUFFIX}/matrix.tsv"),
            ann=config.get("samples", "")
        output:
            l2_volcano=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.L2.volcano.pdf"),
            l2_heatmap=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.L2.driver_heatmap.pdf"),
            enet_volcano=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.ElasticNet.volcano.pdf"),
            enet_heatmap=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.ElasticNet.driver_heatmap.pdf"),
            univariate=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.univariate_stats.tsv")
        params:
            script=os.path.join(PIPE_DIR, "scripts/ml/ml_visualization.py"),
            out_dir=os.path.join(ML_OUTPUT_DIR, "visualizations"),
            case_group=CASE_GROUP,
            ctrl_group=CTRL_GROUP,
            top_k=KBEST,
            C=C_PARAM,
            alpha=ELASTICNET_ALPHA,
            top_n_heatmap=ML_VIZ_TOP_N_HEATMAP,
            top_n_barplot=ML_VIZ_TOP_N_BARPLOT
        threads: 4
        log:
            os.path.join(WORK_DIR, "logs/ml/visualization.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script} \
                --matrix {input.matrix} \
                --annotation {input.ann} \
                --output {params.out_dir} \
                --feature-space autosomal_windows \
                --case-group {params.case_group} \
                --ctrl-group {params.ctrl_group} \
                --top-k {params.top_k} \
                --C {params.C} \
                --alpha {params.alpha} \
                --model-type both \
                --top-n-heatmap {params.top_n_heatmap} \
                --top-n-barplot {params.top_n_barplot} \
                2> {log}
            """
    
    # =========================================================================
    # 13. Combined ML Discrimination Target (with reports)
    # =========================================================================
    rule ml_discrimination_all:
        """Run complete ML discrimination workflow (autosomal-only by default)"""
        input:
            ml_summary=os.path.join(ML_OUTPUT_DIR, f"windows{WINDOW_SIZE}_{MATRIX_SUFFIX}__{ML_MODE}/ml_cv.summary.tsv"),
            univariate=os.path.join(ML_OUTPUT_DIR, "feature_drivers/univariate_top_windows.tsv"),
            dmr_summary=os.path.join(ML_OUTPUT_DIR, f"DMR/windows{WINDOW_SIZE}_{MATRIX_SUFFIX}/dmr.summary.tsv"),
            robustness=os.path.join(ML_OUTPUT_DIR, "robustness/robustness_checks.tsv"),
            exact_perm=os.path.join(ML_OUTPUT_DIR, "exact_permutation/autosomal_windows.summary.tsv"),
            ml_report=os.path.join(ML_BASE, f"reports/ml_report_{CASE_GROUP}_vs_{CTRL_GROUP}.html"),
            # Visualizations (volcano plots, heatmaps)
            l2_volcano=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.L2.volcano.pdf"),
            l2_heatmap=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.L2.driver_heatmap.pdf"),
            enet_heatmap=os.path.join(ML_OUTPUT_DIR, "visualizations/autosomal_windows.ElasticNet.driver_heatmap.pdf")
        output:
            report=os.path.join(ML_OUTPUT_DIR, "ml_discrimination_report.txt")
        params:
            autosomal_only=AUTOSOMAL_ONLY
        shell:
            """
            echo "ML Discrimination Analysis Complete" > {output.report}
            echo "====================================" >> {output.report}
            echo "" >> {output.report}
            echo "Groups: {CASE_GROUP} vs {CTRL_GROUP}" >> {output.report}
            echo "Window size: {WINDOW_SIZE} bp" >> {output.report}
            echo "Autosomal-only: {params.autosomal_only}" >> {output.report}
            echo "K-best features: {KBEST}" >> {output.report}
            echo "Models: L2 Logistic Regression + Elastic-Net" >> {output.report}
            echo "" >> {output.report}
            echo "ML CV Results:" >> {output.report}
            cat {input.ml_summary} >> {output.report}
            echo "" >> {output.report}
            echo "Exact Permutation Test (confirmatory):" >> {output.report}
            cat {input.exact_perm} >> {output.report}
            echo "" >> {output.report}
            echo "Robustness Checks:" >> {output.report}
            cat {input.robustness} >> {output.report}
            echo "" >> {output.report}
            echo "DMR Results:" >> {output.report}
            cat {input.dmr_summary} >> {output.report}
            echo "" >> {output.report}
            echo "HTML Report: {input.ml_report}" >> {output.report}
            echo "" >> {output.report}
            echo "Output directory: {ML_OUTPUT_DIR}" >> {output.report}
            """
