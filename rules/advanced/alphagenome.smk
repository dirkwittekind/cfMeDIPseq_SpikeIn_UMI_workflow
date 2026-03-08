# =============================================================================
# MEDIPIPE Advanced Rules - AlphaGenome Integration
# =============================================================================
# Optional workflow for regulatory element and motif analysis using Google
# DeepMind's AlphaGenome model. Requires API key.
#
# Features:
#   - Region track prediction (DNASE, CAGE, CHIP_TF)
#   - ISM (In-silico Mutagenesis) hotspot detection
#   - Merged regulatory report generation
#
# Prerequisites:
#   - AlphaGenome API key (set via ALPHAGENOME_API_KEY environment variable)
#   - Input BED file with regions of interest (e.g., DMRs, ML-selected windows)
# =============================================================================

import os

# -----------------------------------------------------------------------------
# AlphaGenome Configuration
# -----------------------------------------------------------------------------
ALPHAGENOME_CONFIG = config.get("alphagenome", {})
USE_ALPHAGENOME = ALPHAGENOME_CONFIG.get("enable", False)

# Only define rules if AlphaGenome is enabled
if USE_ALPHAGENOME:
    
    # -------------------------------------------------------------------------
    # AlphaGenome Region Track Prediction
    # -------------------------------------------------------------------------
    rule alphagenome_region_tracks:
        """Predict regulatory tracks (DNASE, CAGE, CHIP_TF) for genomic regions.
        
        Uses AlphaGenome API to predict chromatin accessibility, gene expression,
        and transcription factor binding for input regions.
        
        Requires: ALPHAGENOME_API_KEY environment variable
        """
        input:
            regions=ALPHAGENOME_CONFIG.get("input_regions", 
                os.path.join(OUTPUT_DIR, "dmr/significant_dmrs.bed"))
        output:
            tsv=os.path.join(OUTPUT_DIR, "alphagenome/region_tracks.tsv")
        log:
            os.path.join(WORK_DIR, "logs/alphagenome/region_tracks.log")
        conda:
            os.path.join(ENV_DIR, "alphagenome.yaml")
        params:
            api_key_env=ALPHAGENOME_CONFIG.get("api_key_env", "ALPHAGENOME_API_KEY"),
            ontology_terms=ALPHAGENOME_CONFIG.get("ontology_terms", [
                "UBERON:0007650",  # Brain
                "UBERON:0001043",  # Esophagus  
                "UBERON:0000945",  # Stomach
                "UBERON:0000178"   # Blood
            ]),
            requested_outputs=ALPHAGENOME_CONFIG.get("requested_outputs", [
                "DNASE", "CAGE", "CHIP_TF"
            ]),
            seq_len_key=ALPHAGENOME_CONFIG.get("sequence_length_key", "SEQUENCE_LENGTH_16KB"),
            top_n=int(ALPHAGENOME_CONFIG.get("top_n_regions", 100)),
            script=os.path.join(PIPE_DIR, "scripts/advanced/alphagenome/alphagenome_region_tracks.py")
        resources:
            mem_mb=16000
        threads: 4
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv})
            mkdir -p $(dirname {log})
            
            # Check for API key
            if [ -z "${{{params.api_key_env}:-}}" ]; then
                echo "[ERROR] AlphaGenome API key not set. Export {params.api_key_env} environment variable." >&2
                exit 1
            fi
            
            python {params.script} \
                --regions {input.regions} \
                --output {output.tsv} \
                --api-key-env {params.api_key_env} \
                --ontology-terms "{params.ontology_terms}" \
                --requested-outputs "{params.requested_outputs}" \
                --sequence-length-key {params.seq_len_key} \
                --top-n {params.top_n} \
                2> {log}
            """


    # -------------------------------------------------------------------------
    # AlphaGenome ISM Hotspot Detection
    # -------------------------------------------------------------------------
    rule alphagenome_ism_hotspots:
        """Identify regulatory hotspots using In-Silico Mutagenesis (ISM).
        
        Performs base-level importance scoring to identify regions where
        sequence changes most affect predicted regulatory activity.
        
        Outputs:
          - hotspots: Contiguous regions of high importance
          - bedgraph: Per-base importance scores
        """
        input:
            regions=ALPHAGENOME_CONFIG.get("input_regions",
                os.path.join(OUTPUT_DIR, "dmr/significant_dmrs.bed"))
        output:
            hotspots=os.path.join(OUTPUT_DIR, "alphagenome/ism_hotspots.tsv"),
            bedgraph=os.path.join(OUTPUT_DIR, "alphagenome/ism_bedgraph.tsv")
        log:
            os.path.join(WORK_DIR, "logs/alphagenome/ism_hotspots.log")
        conda:
            os.path.join(ENV_DIR, "alphagenome.yaml")
        params:
            api_key_env=ALPHAGENOME_CONFIG.get("api_key_env", "ALPHAGENOME_API_KEY"),
            ontology_terms=ALPHAGENOME_CONFIG.get("ontology_terms", [
                "UBERON:0007650", "UBERON:0001043", "UBERON:0000945", "UBERON:0000178"
            ]),
            seq_len_key=ALPHAGENOME_CONFIG.get("sequence_length_key", "SEQUENCE_LENGTH_16KB"),
            ism_output=ALPHAGENOME_CONFIG.get("ism_output", "DNASE"),
            top_n=int(ALPHAGENOME_CONFIG.get("top_n_regions", 100)),
            ism_max_len=int(ALPHAGENOME_CONFIG.get("ism_max_len_bp", 512)),
            q=float(ALPHAGENOME_CONFIG.get("ism_quantile", 0.995)),
            min_hotspot=int(ALPHAGENOME_CONFIG.get("min_hotspot_bp", 8)),
            script=os.path.join(PIPE_DIR, "scripts/advanced/alphagenome/alphagenome_ism_hotspots.py")
        resources:
            mem_mb=32000
        threads: 4
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.hotspots})
            mkdir -p $(dirname {log})
            
            # Check for API key
            if [ -z "${{{params.api_key_env}:-}}" ]; then
                echo "[ERROR] AlphaGenome API key not set. Export {params.api_key_env} environment variable." >&2
                exit 1
            fi
            
            python {params.script} \
                --regions {input.regions} \
                --output-hotspots {output.hotspots} \
                --output-bedgraph {output.bedgraph} \
                --api-key-env {params.api_key_env} \
                --ontology-terms "{params.ontology_terms}" \
                --sequence-length-key {params.seq_len_key} \
                --ism-output {params.ism_output} \
                --top-n {params.top_n} \
                --ism-max-len {params.ism_max_len} \
                --quantile {params.q} \
                --min-hotspot {params.min_hotspot} \
                2> {log}
            """


    # -------------------------------------------------------------------------
    # AlphaGenome Merge Report
    # -------------------------------------------------------------------------
    rule alphagenome_merge_report:
        """Generate merged regulatory analysis report.
        
        Combines region track predictions, ISM hotspots, and input region
        annotations into a comprehensive regulatory analysis report.
        """
        input:
            tracks=os.path.join(OUTPUT_DIR, "alphagenome/region_tracks.tsv"),
            hotspots=os.path.join(OUTPUT_DIR, "alphagenome/ism_hotspots.tsv"),
            regions=ALPHAGENOME_CONFIG.get("input_regions",
                os.path.join(OUTPUT_DIR, "dmr/significant_dmrs.bed"))
        output:
            report=os.path.join(OUTPUT_DIR, "alphagenome/final_report.tsv")
        log:
            os.path.join(WORK_DIR, "logs/alphagenome/merge_report.log")
        conda:
            os.path.join(ENV_DIR, "alphagenome.yaml")
        params:
            script=os.path.join(PIPE_DIR, "scripts/advanced/alphagenome/alphagenome_merge_report.py")
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.report})
            
            python {params.script} \
                --tracks {input.tracks} \
                --hotspots {input.hotspots} \
                --regions {input.regions} \
                --output {output.report} \
                2> {log}
            """


    # -------------------------------------------------------------------------
    # AlphaGenome Analysis Target
    # -------------------------------------------------------------------------
    rule alphagenome_analysis:
        """Run complete AlphaGenome regulatory analysis pipeline."""
        input:
            os.path.join(OUTPUT_DIR, "alphagenome/region_tracks.tsv"),
            os.path.join(OUTPUT_DIR, "alphagenome/ism_hotspots.tsv"),
            os.path.join(OUTPUT_DIR, "alphagenome/ism_bedgraph.tsv"),
            os.path.join(OUTPUT_DIR, "alphagenome/final_report.tsv")


# -----------------------------------------------------------------------------
# Helper function for workflow targets
# -----------------------------------------------------------------------------
def get_alphagenome_targets():
    """Get AlphaGenome analysis targets if enabled."""
    if not USE_ALPHAGENOME:
        return []
    return [
        os.path.join(OUTPUT_DIR, "alphagenome/region_tracks.tsv"),
        os.path.join(OUTPUT_DIR, "alphagenome/ism_hotspots.tsv"),
        os.path.join(OUTPUT_DIR, "alphagenome/final_report.tsv")
    ]
