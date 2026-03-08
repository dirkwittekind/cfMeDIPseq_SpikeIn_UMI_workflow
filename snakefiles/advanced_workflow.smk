# =============================================================================
# MEDIPIPE Advanced Workflow - Stages 6-7
# =============================================================================
# This Snakefile runs the advanced analysis pipeline:
#   Stage 6: DMR Analysis
#   Stage 7: ML/DL Classification
# =============================================================================

import os
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Determine rule paths
PIPE_DIR = config.get("pipe_dir", "/home/dirk/medipipe_warp")
RULES_DIR = os.path.join(PIPE_DIR, "rules")

# Include common functions and variables
include: os.path.join(RULES_DIR, "core/common.smk")

# Include advanced workflow rules
include: os.path.join(RULES_DIR, "advanced/ml_classification.smk")
include: os.path.join(RULES_DIR, "advanced/dl_models.smk")
include: os.path.join(RULES_DIR, "advanced/ml_discrimination.smk")

# Include AlphaGenome regulatory analysis (optional)
include: os.path.join(RULES_DIR, "advanced/alphagenome.smk")


# =============================================================================
# Main Target Rule
# =============================================================================
rule all:
    """Main target: run all enabled advanced workflow stages"""
    input:
        get_advanced_workflow_targets()


# =============================================================================
# Stage-specific Target Rules
# =============================================================================
rule stage6_dmr:
    """Stage 6: DMR analysis"""
    input:
        os.path.join(OUTPUT_DIR, "dmr/dmr_results.tsv") if config.get("advanced", {}).get("dmr_analysis", False) else []


rule stage7_ml:
    """Stage 7: Machine Learning classification"""
    input:
        os.path.join(OUTPUT_DIR, "ml/model_results.tsv") if config.get("advanced", {}).get("ml_classification", False) else [],
        os.path.join(OUTPUT_DIR, "ml/roc_curve.png") if config.get("advanced", {}).get("ml_classification", False) else []


rule stage7_dl:
    """Stage 7: Deep Learning classification"""
    input:
        os.path.join(OUTPUT_DIR, "dl/model_results.tsv") if config.get("advanced", {}).get("dl_classification", False) else [],
        os.path.join(OUTPUT_DIR, "dl/training_history.png") if config.get("advanced", {}).get("dl_classification", False) else []


# =============================================================================
# ML Discrimination Workflow (Post-Analysis)
# =============================================================================
rule ml_discrimination:
    """Run ML discrimination workflow for group comparison"""
    input:
        os.path.join(OUTPUT_DIR, "ml_discrimination/Outputs/{case}_vs_{ctrl}/ml_discrimination_report.txt".format(
            case=config.get("comparison", {}).get("group1", "AEG"),
            ctrl=config.get("comparison", {}).get("group2", "CTRL")
        )) if config.get("advanced", {}).get("ml_classification", False) else []


# =============================================================================
# Combined ML + DL
# =============================================================================
rule classification:
    """Run both ML and DL classification"""
    input:
        rules.stage7_ml.input,
        rules.stage7_dl.input


# =============================================================================
# AlphaGenome Regulatory Analysis
# =============================================================================
rule alphagenome:
    """Run AlphaGenome regulatory element analysis on DMRs/windows.
    
    Requires: ALPHAGENOME_API_KEY environment variable
    """
    input:
        get_alphagenome_targets() if 'get_alphagenome_targets' in dir() else []
