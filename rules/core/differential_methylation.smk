# =============================================================================
# MEDIPIPE Core Rules - Differential Methylation Analysis
# =============================================================================
# Compares methylation between two groups (e.g., AEG vs CTRL)
# =============================================================================

import os
import pandas as pd

# -----------------------------------------------------------------------------
# Helper functions to get BAMs by group
# -----------------------------------------------------------------------------
def get_group_samples(group_name):
    """Get sample IDs for a specific group"""
    if SAMPLES.empty:
        return []
    if "group" not in SAMPLES.columns:
        return []
    return SAMPLES[SAMPLES["group"] == group_name]["sample_id"].tolist()


def get_group_bams(group_name):
    """Get BAM file paths for a specific group"""
    sample_ids = get_group_samples(group_name)
    bams = []
    for sid in sample_ids:
        if PAIRED_END and USE_UMI:
            bam = os.path.join(OUTPUT_DIR, f"dedup_bam_umi/{sid}_dedup.bam")
        elif PAIRED_END:
            bam = os.path.join(OUTPUT_DIR, f"dedup_bam/{sid}_dedup.bam")
        else:
            bam = os.path.join(OUTPUT_DIR, f"dedup_bam_se/{sid}_dedup.bam")
        bams.append(bam)
    return bams


# Get comparison groups from config
COMPARISON = config.get("comparison", {})
GROUP1_NAME = COMPARISON.get("group1", "")
GROUP2_NAME = COMPARISON.get("group2", "")

# Auto-detect groups if not specified
if not GROUP1_NAME or not GROUP2_NAME:
    if "group" in SAMPLES.columns:
        unique_groups = SAMPLES["group"].unique().tolist()
        # Try to identify case vs control
        ctrl_candidates = ["CTRL", "Control", "Kontrolle", "control", "ctrl", "healthy", "normal"]
        case_groups = [g for g in unique_groups if g not in ctrl_candidates]
        ctrl_groups = [g for g in unique_groups if g in ctrl_candidates]
        
        if case_groups and ctrl_groups:
            GROUP1_NAME = case_groups[0]  # First case group (e.g., AEG)
            GROUP2_NAME = ctrl_groups[0]  # First control group (e.g., CTRL)


# -----------------------------------------------------------------------------
# Differential Methylation Analysis Rule
# -----------------------------------------------------------------------------
if GROUP1_NAME and GROUP2_NAME:
    
    rule medips_differential:
        """Run differential MEDIPS analysis comparing two groups"""
        input:
            bams_group1=lambda w: get_group_bams(GROUP1_NAME),
            bams_group2=lambda w: get_group_bams(GROUP2_NAME)
        output:
            differential=os.path.join(OUTPUT_DIR, f"differential/{GROUP1_NAME}_vs_{GROUP2_NAME}.differential.tsv"),
            dmr_005=os.path.join(OUTPUT_DIR, f"differential/{GROUP1_NAME}_vs_{GROUP2_NAME}.DMRs.padj005.tsv"),
            summary=os.path.join(OUTPUT_DIR, f"differential/{GROUP1_NAME}_vs_{GROUP2_NAME}.summary.txt")
        params:
            genome=GENOME,
            window_size=WINDOW_SIZE,
            group1_name=GROUP1_NAME,
            group2_name=GROUP2_NAME,
            script_dir=os.path.join(PIPE_DIR, "scripts/core"),
            output_prefix=os.path.join(OUTPUT_DIR, f"differential/{GROUP1_NAME}_vs_{GROUP2_NAME}")
        resources:
            mem_mb=120000
        threads: 16
        log:
            os.path.join(WORK_DIR, f"logs/differential/{GROUP1_NAME}_vs_{GROUP2_NAME}.log")
        conda:
            os.path.join(ENV_DIR, "medipipe_r.yaml")
        shell:
            """
            mkdir -p $(dirname {output.differential})
            
            # Build comma-separated BAM lists
            BAM_LIST_1=$(echo {input.bams_group1} | tr ' ' ',')
            BAM_LIST_2=$(echo {input.bams_group2} | tr ' ' ',')
            
            Rscript --vanilla {params.script_dir}/medips_differential.R \\
                --bam-list-1 "$BAM_LIST_1" \\
                --bam-list-2 "$BAM_LIST_2" \\
                --group1-name {params.group1_name} \\
                --group2-name {params.group2_name} \\
                --output {params.output_prefix} \\
                --genome {params.genome} \\
                --window {params.window_size} \\
                --paired \\
                2>&1 | tee {log}
            """
    
    
    # Add differential analysis to workflow targets
    def get_differential_targets():
        """Get targets for differential analysis"""
        targets = []
        if GROUP1_NAME and GROUP2_NAME:
            targets.append(os.path.join(OUTPUT_DIR, f"differential/{GROUP1_NAME}_vs_{GROUP2_NAME}.differential.tsv"))
            targets.append(os.path.join(OUTPUT_DIR, f"differential/{GROUP1_NAME}_vs_{GROUP2_NAME}.summary.txt"))
        return targets
