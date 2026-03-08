#!/bin/bash
# =============================================================================
# Nature Cancer Reference Panel Scoring for cfMeDIP-seq
# Scores samples against cancer-specific DMR panels
# =============================================================================

set -e

BW_MAP="/data/medipipe_data/output/reference_scoring/bw_map.tsv"
PANEL_DIR="/home/dirk/medipipe_warp/resources/Code/15191455/pan_cancer_cfMeDIP_Figures_source_data/pan_cancer_cfMDIP_Figures_source_data/DMRs_Source_Data/BED_files"
OUT_DIR="/data/medipipe_data/output/reference_scoring"
THREADS=12

# Get BigWig list from map
BIGWIGS=$(cut -f3 "$BW_MAP" | tr '\n' ' ')
LABELS=$(cut -f1 "$BW_MAP" | tr '\n' ' ')

echo "=== Nature Cancer Reference Panel Scoring ==="
echo "Samples: $(echo $LABELS | wc -w)"
echo "Output: $OUT_DIR"
echo ""

mkdir -p "$OUT_DIR/panels"

# Define panels to score
declare -A PANELS
PANELS["panCancer_Hyper"]="panCancer_Hyper_DMRs.bed"
PANELS["panCancer_Hypo"]="panCancer_Hypo_DMRs.bed"
PANELS["AML_Hyper"]="AML specific Hyper DMRs.bed"
PANELS["AML_Hypo"]="AML specific Hypo DMRs.bed"
PANELS["Bladder_Hyper"]="Bladder Cancer specific Hyper DMRs.bed"
PANELS["Brain_Hyper"]="Brain Cancer specific Hyper DMRs.bed"
PANELS["Brain_Hypo"]="Brain Cancer specific Hypo DMRs.bed"
PANELS["Breast_Hyper"]="Breast Cancer specific Hyper DMRs.bed"
PANELS["Colorectal_Hyper"]="Colorectal Cancer specific Hyper DMRs.bed"
PANELS["HeadNeck_Hyper"]="Head and Neck Cancer specific Hyper DMRs.bed"
PANELS["HeadNeck_Hypo"]="Head and Neck Cancer specific Hypo DMRs.bed"
PANELS["Lung_Hyper"]="Lung Cancer specific Hyper DMRs.bed"
PANELS["Pancreatic_Hyper"]="Pancreatic Cancer specific Hyper DMRs.bed"
PANELS["Pancreatic_Hypo"]="Pancreatic Cancer specific Hypo DMRs.bed"
PANELS["Prostate_Hyper"]="Prostate Cancer specific Hyper DMRs.bed"
PANELS["Prostate_Hypo"]="Prostate Cancer specific Hypo DMRs.bed"
PANELS["Renal_Hyper"]="Renal Cancer specific Hyper DMRs.bed"
PANELS["Renal_Hypo"]="Renal Cancer specific Hypo DMRs.bed"
PANELS["UvealMelanoma_Hyper"]="Uveal Melanoma specific Hyper DMRs.bed"
PANELS["UvealMelanoma_Hypo"]="Uveal Melanoma specific Hypo DMRs.bed"

echo "Processing ${#PANELS[@]} reference panels..."
echo ""

for panel_id in "${!PANELS[@]}"; do
    bed_file="$PANEL_DIR/${PANELS[$panel_id]}"
    out_raw="$OUT_DIR/panels/${panel_id}.raw.tsv"
    out_npz="$OUT_DIR/panels/${panel_id}.npz"
    
    if [[ ! -f "$bed_file" ]]; then
        echo "WARNING: Panel BED not found: $bed_file"
        continue
    fi
    
    n_regions=$(wc -l < "$bed_file")
    if [[ $n_regions -lt 5 ]]; then
        echo "SKIP: $panel_id ($n_regions regions - too few)"
        continue
    fi
    
    echo "  Processing $panel_id ($n_regions regions)..."
    
    multiBigwigSummary BED-file \
        --bwfiles $BIGWIGS \
        --BED "$bed_file" \
        --labels $LABELS \
        --outFileName "$out_npz" \
        --outRawCounts "$out_raw" \
        --numberOfProcessors $THREADS \
        2>/dev/null || echo "    WARNING: multiBigwigSummary failed for $panel_id"
done

echo ""
echo "Running Python scoring script..."

# Run the Python scoring script
python /home/dirk/medipipe_warp/scripts/advanced/score_reference_panels.py \
    --bw-map "$BW_MAP" \
    --panel-raw-glob "$OUT_DIR/panels/*.raw.tsv" \
    --control-label CTRL \
    --out-prefix "$OUT_DIR/AEG_vs_CTRL"

echo ""
echo "=== Reference Panel Scoring Complete ==="
echo "Results: $OUT_DIR/AEG_vs_CTRL.scores.long.tsv"
echo "Summary: $OUT_DIR/AEG_vs_CTRL.scores.summary.tsv"
