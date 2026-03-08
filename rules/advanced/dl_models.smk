# =============================================================================
# MEDIPIPE Advanced Rules - Deep Learning Models (PyTorch)
# =============================================================================

DL_CONFIG = config.get("dl", {})
ADVANCED_CONFIG = config.get("advanced", {})

if ADVANCED_CONFIG.get("dl_classification", False):
    # -------------------------------------------------------------------------
    # Prepare DL Dataset
    # -------------------------------------------------------------------------
    rule dl_prepare_dataset:
        """Prepare dataset for deep learning training"""
        input:
            meth_counts=expand(os.path.join(OUTPUT_DIR, "meth_quant/{sample}_count.txt"), sample=SAMPLE_IDS),
            samples_file=config.get("samples", "")
        output:
            train_data=os.path.join(OUTPUT_DIR, "dl/train_data.h5"),
            test_data=os.path.join(OUTPUT_DIR, "dl/test_data.h5"),
            metadata=os.path.join(OUTPUT_DIR, "dl/dataset_metadata.json")
        params:
            script_dir=os.path.join(PIPE_DIR, "scripts/dl"),
            test_size=0.2,
            random_state=42
        log:
            os.path.join(WORK_DIR, "logs/dl/prepare_dataset.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_dl.yaml")
        shell:
            """
            mkdir -p $(dirname {output.train_data})
            
            python {params.script_dir}/prepare_dl_dataset.py \\
                --input-dir $(dirname {input.meth_counts[0]}) \\
                --samples {input.samples_file} \\
                --output-train {output.train_data} \\
                --output-test {output.test_data} \\
                --output-metadata {output.metadata} \\
                --test-size {params.test_size} \\
                --random-state {params.random_state} \\
                2> {log}
            """

    # -------------------------------------------------------------------------
    # Train 1D CNN Model
    # -------------------------------------------------------------------------
    rule dl_train_cnn:
        """Train 1D CNN model for methylation classification"""
        input:
            train_data=os.path.join(OUTPUT_DIR, "dl/train_data.h5"),
            test_data=os.path.join(OUTPUT_DIR, "dl/test_data.h5")
        output:
            model=os.path.join(OUTPUT_DIR, "dl/cnn_model.pt"),
            results=os.path.join(OUTPUT_DIR, "dl/cnn_results.tsv"),
            history=os.path.join(OUTPUT_DIR, "dl/cnn_training_history.png")
        params:
            script_dir=os.path.join(PIPE_DIR, "scripts/dl"),
            epochs=DL_CONFIG.get("epochs", 100),
            batch_size=DL_CONFIG.get("batch_size", 32),
            lr=DL_CONFIG.get("learning_rate", 0.001),
            patience=DL_CONFIG.get("early_stopping_patience", 10)
        threads: config.get("threads", 24)
        resources:
            gpu=1 if DL_CONFIG.get("use_gpu", True) else 0
        log:
            os.path.join(WORK_DIR, "logs/dl/train_cnn.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_dl.yaml")
        shell:
            """
            python {params.script_dir}/train_dl.py \\
                --train-data {input.train_data} \\
                --test-data {input.test_data} \\
                --model-type cnn \\
                --output-model {output.model} \\
                --output-results {output.results} \\
                --output-history {output.history} \\
                --epochs {params.epochs} \\
                --batch-size {params.batch_size} \\
                --learning-rate {params.lr} \\
                --patience {params.patience} \\
                2> {log}
            """

    # -------------------------------------------------------------------------
    # Train Transformer Model
    # -------------------------------------------------------------------------
    rule dl_train_transformer:
        """Train Transformer model for methylation classification"""
        input:
            train_data=os.path.join(OUTPUT_DIR, "dl/train_data.h5"),
            test_data=os.path.join(OUTPUT_DIR, "dl/test_data.h5")
        output:
            model=os.path.join(OUTPUT_DIR, "dl/transformer_model.pt"),
            results=os.path.join(OUTPUT_DIR, "dl/transformer_results.tsv"),
            history=os.path.join(OUTPUT_DIR, "dl/transformer_training_history.png")
        params:
            script_dir=os.path.join(PIPE_DIR, "scripts/dl"),
            epochs=DL_CONFIG.get("epochs", 100),
            batch_size=DL_CONFIG.get("batch_size", 32),
            lr=DL_CONFIG.get("learning_rate", 0.001),
            patience=DL_CONFIG.get("early_stopping_patience", 10)
        threads: config.get("threads", 24)
        resources:
            gpu=1 if DL_CONFIG.get("use_gpu", True) else 0
        log:
            os.path.join(WORK_DIR, "logs/dl/train_transformer.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_dl.yaml")
        shell:
            """
            python {params.script_dir}/train_dl.py \\
                --train-data {input.train_data} \\
                --test-data {input.test_data} \\
                --model-type transformer \\
                --output-model {output.model} \\
                --output-results {output.results} \\
                --output-history {output.history} \\
                --epochs {params.epochs} \\
                --batch-size {params.batch_size} \\
                --learning-rate {params.lr} \\
                --patience {params.patience} \\
                2> {log}
            """

    # -------------------------------------------------------------------------
    # DL Results Summary
    # -------------------------------------------------------------------------
    rule dl_summary:
        """Summarize deep learning results"""
        input:
            cnn_results=os.path.join(OUTPUT_DIR, "dl/cnn_results.tsv"),
            transformer_results=os.path.join(OUTPUT_DIR, "dl/transformer_results.tsv")
        output:
            summary=os.path.join(OUTPUT_DIR, "dl/model_results.tsv"),
            comparison=os.path.join(OUTPUT_DIR, "dl/training_history.png")
        params:
            script_dir=os.path.join(PIPE_DIR, "scripts/dl")
        log:
            os.path.join(WORK_DIR, "logs/dl/summary.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_dl.yaml")
        shell:
            """
            python {params.script_dir}/summarize_dl.py \\
                --cnn-results {input.cnn_results} \\
                --transformer-results {input.transformer_results} \\
                --output-summary {output.summary} \\
                --output-plot {output.comparison} \\
                2> {log}
            """
