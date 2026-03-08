# =============================================================================
# MEDIPIPE Core Rules - Common Functions and Variables
# =============================================================================

import os
import pandas as pd
from pathlib import Path

# -----------------------------------------------------------------------------
# Directory Configuration
# -----------------------------------------------------------------------------
PIPE_DIR = config.get("pipe_dir", "/home/dirk/medipipe_warp")
WORK_DIR = config.get("work_dir", os.path.join(PIPE_DIR, "work"))
OUTPUT_DIR = config.get("output_dir", "/data/medipipe_data/output")
INPUT_DIR = config.get("input_dir", "/data/raw_data")
RESOURCES_DIR = config.get("resources_dir", os.path.join(PIPE_DIR, "resources"))

# Create work directory if needed
os.makedirs(WORK_DIR, exist_ok=True)

# Environment directory
ENV_DIR = os.path.join(PIPE_DIR, "environments")

# -----------------------------------------------------------------------------
# Reference Files
# -----------------------------------------------------------------------------
REF = config.get("ref_files", {})
GENOME = config.get("genome", "hg38")
BSGENOME = config.get("bsgenome", "BSgenome.Hsapiens.UCSC.hg38")

# -----------------------------------------------------------------------------
# Analysis Options
# -----------------------------------------------------------------------------
PAIRED_END = config.get("paired_end", config.get("paired-end", True))
# Robust UMI toggle (backward compatible with legacy config keys)
USE_UMI = bool(config.get("use_umi", False) or config.get("umi", False) or config.get("add_umi", False))
# Robust spike-in toggle (backward compatible)
USE_SPIKEIN = bool(config.get("use_spikein", False) or config.get("spike_in", False))
WINDOW_SIZE = config.get("window_size", 300)
THREADS = config.get("threads", 24)

# -----------------------------------------------------------------------------
# UMI Configuration
# For Twist UMI Adapter system: 5bp UMI + 2bp skip per read
# Pattern explanation: (?P<umi_N>.{X}) captures X bases as UMI
#                     (?P<discard_N>.{Y}) removes Y bases after UMI
#                     .* matches remaining read
# -----------------------------------------------------------------------------
UMI_ENV = os.path.join(ENV_DIR, "umi_tools.yaml")
UMI_REGEX_R1 = config.get("umi_regex_r1", r"(?P<umi_1>.{5})(?P<discard_1>.{2}).*")
UMI_REGEX_R2 = config.get("umi_regex_r2", r"(?P<umi_2>.{5})(?P<discard_2>.{2}).*")

# -----------------------------------------------------------------------------
# Sample Loading
# -----------------------------------------------------------------------------
def load_samples(samples_file):
    """Load samples from TSV file"""
    if not os.path.exists(samples_file):
        raise ValueError(f"Samples file not found: {samples_file}")
    
    df = pd.read_csv(samples_file, sep="\t")
    
    # Validate required columns
    required_cols = ["sample_id", "R1"]
    if PAIRED_END:
        required_cols.append("R2")
    
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found in samples file")
    
    return df.set_index("sample_id", drop=False).sort_index()


SAMPLES_FILE = config.get("samples", os.path.join(PIPE_DIR, "samplefiles", "samples.tsv"))
if os.path.exists(SAMPLES_FILE):
    SAMPLES = load_samples(SAMPLES_FILE)
    SAMPLE_IDS = SAMPLES["sample_id"].tolist()
else:
    SAMPLES = pd.DataFrame()
    SAMPLE_IDS = []

# -----------------------------------------------------------------------------
# Input Functions
# -----------------------------------------------------------------------------
def get_raw_fastq_r1(wildcards):
    """Get R1 FASTQ file for a sample"""
    return SAMPLES.loc[wildcards.sample, "R1"]


def get_raw_fastq_r2(wildcards):
    """Get R2 FASTQ file for a sample (paired-end only)"""
    if PAIRED_END:
        return SAMPLES.loc[wildcards.sample, "R2"]
    return ""


def get_fastq_for_trim(wildcards):
    """Get FASTQ files for trimming (after optional UMI extraction)"""
    if PAIRED_END and USE_UMI:
        return {
            "r1": os.path.join(WORK_DIR, f"umi_extracted/{wildcards.sample}_R1.fastq.gz"),
            "r2": os.path.join(WORK_DIR, f"umi_extracted/{wildcards.sample}_R2.fastq.gz")
        }
    elif PAIRED_END:
        return {
            "r1": SAMPLES.loc[wildcards.sample, "R1"],
            "r2": SAMPLES.loc[wildcards.sample, "R2"]
        }
    else:
        return {"r1": SAMPLES.loc[wildcards.sample, "R1"]}


def get_trimmed_fastq(wildcards):
    """Get trimmed FASTQ files for alignment"""
    if PAIRED_END:
        return {
            "r1": os.path.join(WORK_DIR, f"trimmed/{wildcards.sample}_R1_val_1.fq.gz"),
            "r2": os.path.join(WORK_DIR, f"trimmed/{wildcards.sample}_R2_val_2.fq.gz")
        }
    else:
        return {"r1": os.path.join(WORK_DIR, f"trimmed/{wildcards.sample}_trimmed.fq.gz")}


def get_dedup_bam(wildcards):
    """Get deduplicated BAM file path"""
    if PAIRED_END and USE_UMI:
        return os.path.join(OUTPUT_DIR, f"dedup_bam_umi/{wildcards.sample}_dedup.bam")
    elif PAIRED_END:
        return os.path.join(OUTPUT_DIR, f"dedup_bam/{wildcards.sample}_dedup.bam")
    else:
        return os.path.join(OUTPUT_DIR, f"dedup_bam_se/{wildcards.sample}_dedup.bam")


def get_bwa_index():
    """Get BWA index path - returns combined index if spike-in enabled"""
    # Spike-in aware: combined index if configured, else genome-only index
    if USE_SPIKEIN and config.get("spike_idx"):
        return config.get("spike_idx")
    return REF.get("bwa_index", os.path.join(RESOURCES_DIR, f"genomes/{GENOME}/bwa_index/{GENOME}.fa"))


def get_spikein_bam(wildcards):
    """Get spike-in BAM file path"""
    if USE_SPIKEIN:
        return os.path.join(OUTPUT_DIR, f"spikein_bam/{wildcards.sample}_spikein.bam")
    return None


# -----------------------------------------------------------------------------
# Output Target Functions
# -----------------------------------------------------------------------------
def get_core_workflow_targets():
    """Get all targets for core workflow (stages 1-5)"""
    targets = []
    
    workflow_config = config.get("workflow", {})
    
    # Stage 1: Preprocessing - always include BAM outputs
    if workflow_config.get("preprocessing", True):
        if USE_UMI:
            targets += expand(os.path.join(OUTPUT_DIR, "dedup_bam_umi/{sample}_dedup.bam"), sample=SAMPLE_IDS)
            targets += expand(os.path.join(OUTPUT_DIR, "dedup_bam_umi/{sample}_dedup.bam.bai"), sample=SAMPLE_IDS)
        else:
            targets += expand(os.path.join(OUTPUT_DIR, "dedup_bam/{sample}_dedup.bam"), sample=SAMPLE_IDS)
            targets += expand(os.path.join(OUTPUT_DIR, "dedup_bam/{sample}_dedup.bam.bai"), sample=SAMPLE_IDS)
    
    # Stage 2: QC Analysis
    if workflow_config.get("qc_analysis", True):
        targets += expand(os.path.join(OUTPUT_DIR, "fastqc/{sample}_R1_fastqc.zip"), sample=SAMPLE_IDS)
        targets += expand(os.path.join(OUTPUT_DIR, "fastqc/{sample}_R2_fastqc.zip"), sample=SAMPLE_IDS)
        targets += expand(os.path.join(OUTPUT_DIR, "stats/{sample}_dedup.stats.txt"), sample=SAMPLE_IDS)
        targets += expand(os.path.join(OUTPUT_DIR, "insert_size/{sample}_insert_size_metrics.txt"), sample=SAMPLE_IDS)
    
    # Stage 3: Methylation Quantification
    if workflow_config.get("methylation_quant", True):
        targets += expand(os.path.join(OUTPUT_DIR, "meth_quant/{sample}_meth_qc.txt"), sample=SAMPLE_IDS)
        targets += expand(os.path.join(OUTPUT_DIR, "meth_quant/{sample}_count.txt"), sample=SAMPLE_IDS)
    
    # Stage 4: Fragment Profile
    if workflow_config.get("fragment_profile", True) and config.get("frag_profile", True):
        targets += expand(os.path.join(OUTPUT_DIR, "fragment_profile/{sample}_fragment_profile.txt"), sample=SAMPLE_IDS)
    
    # Stage 5: Tissue-of-Origin
    if workflow_config.get("tissue_of_origin", True):
        too_config = config.get("tissue_of_origin", {})
        if too_config.get("enable", True):
            targets.append(os.path.join(OUTPUT_DIR, "tissue_of_origin/too_fractions.tsv"))
    
    return targets


def get_advanced_workflow_targets():
    """Get all targets for advanced workflow (stages 6-7)"""
    targets = []
    
    advanced_config = config.get("advanced", {})
    
    # Stage 6: DMR Analysis
    if advanced_config.get("dmr_analysis", False):
        targets.append(os.path.join(OUTPUT_DIR, "dmr/dmr_results.tsv"))
        targets.append(os.path.join(OUTPUT_DIR, "dmr/gsva_scores.tsv"))
    
    # Stage 7: ML Classification
    if advanced_config.get("ml_classification", False):
        targets.append(os.path.join(OUTPUT_DIR, "ml/model_results.tsv"))
        targets.append(os.path.join(OUTPUT_DIR, "ml/roc_curve.png"))
    
    # Stage 7: DL Classification
    if advanced_config.get("dl_classification", False):
        targets.append(os.path.join(OUTPUT_DIR, "dl/model_results.tsv"))
        targets.append(os.path.join(OUTPUT_DIR, "dl/training_history.png"))
    
    return targets
