# =============================================================================
# MEDIPIPE Advanced Rules - Machine Learning Classification
# =============================================================================

ML_CONFIG = config.get("ml", {})
ADVANCED_CONFIG = config.get("advanced", {})

if ADVANCED_CONFIG.get("ml_classification", False):
    # -------------------------------------------------------------------------
    # Prepare ML Features
    # -------------------------------------------------------------------------
    rule ml_prepare_features:
        """Prepare feature matrix for ML classification"""
        input:
            meth_counts=expand(os.path.join(OUTPUT_DIR, "meth_quant/{sample}_count.txt"), sample=SAMPLE_IDS),
            samples_file=config.get("samples", "")
        output:
            features=os.path.join(OUTPUT_DIR, "ml/feature_matrix.tsv"),
            labels=os.path.join(OUTPUT_DIR, "ml/labels.tsv")
        params:
            script_dir=os.path.join(PIPE_DIR, "scripts/ml"),
            n_features=ML_CONFIG.get("n_features", 1000),
            feature_selection=ML_CONFIG.get("feature_selection", "variance_threshold")
        log:
            os.path.join(WORK_DIR, "logs/ml/prepare_features.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            mkdir -p $(dirname {output.features})
            
            python {params.script_dir}/feature_selection.py \\
                --input-dir $(dirname {input.meth_counts[0]}) \\
                --samples {input.samples_file} \\
                --output-features {output.features} \\
                --output-labels {output.labels} \\
                --n-features {params.n_features} \\
                --method {params.feature_selection} \\
                2> {log}
            """

    # -------------------------------------------------------------------------
    # Train ML Models
    # -------------------------------------------------------------------------
    rule ml_train:
        """Train machine learning classifiers with nested cross-validation"""
        input:
            features=os.path.join(OUTPUT_DIR, "ml/feature_matrix.tsv"),
            labels=os.path.join(OUTPUT_DIR, "ml/labels.tsv")
        output:
            results=os.path.join(OUTPUT_DIR, "ml/model_results.tsv"),
            model=os.path.join(OUTPUT_DIR, "ml/best_model.pkl"),
            roc=os.path.join(OUTPUT_DIR, "ml/roc_curve.png"),
            pr=os.path.join(OUTPUT_DIR, "ml/pr_curve.png"),
            feature_importance=os.path.join(OUTPUT_DIR, "ml/feature_importance.tsv")
        params:
            script_dir=os.path.join(PIPE_DIR, "scripts/ml"),
            algorithms=ML_CONFIG.get("algorithms", ["random_forest", "gradient_boosting"]),
            cv_folds=ML_CONFIG.get("cv_folds", 5)
        threads: config.get("threads", 24)
        log:
            os.path.join(WORK_DIR, "logs/ml/train.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_ml.yaml")
        shell:
            """
            python {params.script_dir}/train_classifier.py \\
                --features {input.features} \\
                --labels {input.labels} \\
                --output-dir $(dirname {output.results}) \\
                --algorithms {params.algorithms} \\
                --cv-folds {params.cv_folds} \\
                --threads {threads} \\
                2> {log}
            """
