# MEDIPIPE cfMeDIP-seq Analysis Pipeline

## Project Overview

MEDIPIPE is a comprehensive bioinformatics pipeline for cell-free methylated DNA immunoprecipitation sequencing (cfMeDIP-seq) analysis. It includes preprocessing, quality control, methylation quantification, tissue of origin deconvolution, and machine learning discrimination.

**Project Location:** `/home/dirk/medipipe_warp`
**Data Location:** `/data/medipipe_data`
**Raw Data:** `/data/raw_data/raw_data_temp`

---

## Quick Start Commands

### Launch GUI
```shell
/home/dirk/medipipe_warp/launch_medipipe.sh
```

### Run Core Workflow
```shell
cd /home/dirk/medipipe_warp && snakemake --snakefile snakefiles/core_workflow.smk --configfile work/config.yaml --cores {{cores}} --use-conda --rerun-incomplete --keep-going --conda-frontend conda --latency-wait 60
```

### Dry Run (Preview Jobs)
```shell
cd /home/dirk/medipipe_warp && snakemake --snakefile snakefiles/core_workflow.smk --configfile work/config.yaml --cores 24 --use-conda --dry-run
```

### Monitor Workflow Progress
```shell
watch -n 30 'echo "=== Snakemake Status ===" && ps aux | grep snakemake | grep -v grep | head -1 && echo "" && echo "=== Recent Progress ===" && grep -E "Finished|steps.*done|Error" $(ls -t /home/dirk/medipipe_warp/work/logs/workflow_*.log 2>/dev/null | head -1) 2>/dev/null | tail -5'
```

### Stop Running Workflow
```shell
pkill -f "snakemake.*core_workflow"
```

---

## Directory Structure

```
/home/dirk/medipipe_warp/           # Pipeline code (local SSD)
├── gui/                            # Tkinter GUI application
│   ├── main.py                     # Entry point
│   ├── app.py                      # Main application class
│   └── components/                 # UI components
├── scripts/                        # Analysis scripts
│   ├── core/                       # Core analysis (MEDIPS, NNLS)
│   ├── ml/                         # Machine learning
│   ├── qc/                         # Quality control
│   └── reporting/                  # Report generation
├── rules/                          # Snakemake rules
│   ├── core/                       # Core workflow rules
│   └── qc/                         # QC rules
├── snakefiles/                     # Main Snakemake files
├── environments/                   # Conda environment YAMLs
├── resources/                      # Reference files
├── work/                           # Working directory
│   ├── config.yaml                 # Main configuration
│   └── logs/                       # Workflow logs
└── .warp/workflows/                # Warp workflows

/data/medipipe_data/                # Data storage (Elements drive)
├── output/                         # Analysis outputs
│   ├── fastqc/                     # FastQC reports
│   ├── dedup_bam_umi/              # Deduplicated BAMs
│   ├── meth_quant/                 # Methylation quantification
│   ├── tissue_deconv/              # Tissue of origin
│   ├── bigwig/                     # Coverage tracks
│   └── qc_reports/                 # QC reports
└── temp/                           # Temporary files
```

---

## Workflow Stages

### 1. Preprocessing
- **FastQC:** Raw read quality assessment
- **Trim Galore:** Adapter trimming and quality filtering
- **UMI Extraction:** Extract UMI barcodes (optional)

### 2. Alignment
- **BWA-MEM:** Alignment to reference genome (hg38/hg19)
- **Samtools:** BAM processing and statistics

### 3. Deduplication
- **UMI-tools:** UMI-based deduplication (recommended)
- **Picard:** PCR duplicate marking (alternative)

### 4. Methylation Quantification
- **MEDIPS:** CpG enrichment and methylation calling
- **Window size:** 1000bp (to avoid integer overflow)
- **QSEA:** Optional (disabled by default - very slow)

### 5. Quality Control
- **Fragment profiles:** Insert size distribution
- **Saturation analysis:** Library complexity
- **Spike-in analysis:** IP efficiency (optional)

### 6. Tissue of Origin
- **NNLS deconvolution:** Cell type proportions
- **Reference matrix:** Tissue methylation signatures

---

## QC Metrics and Thresholds

| Metric | Threshold | Description |
|--------|-----------|-------------|
| Dedup reads | >5,000,000 | Minimum deduplicated reads |
| Duplication rate | <90% | Maximum acceptable duplication |
| relH | >1.5, <10.0 | CpG enrichment score |
| GoGe | >1.0 | Observed/expected CpG ratio |
| Spike-in relH | >2.5 | IP efficiency from spike-in |
| Spike-in GoGe | >1.5 | Spike-in specificity |
| Meth ratio | >0.7 | Methylated/total spike-in |

### Generate QC Report
```shell
python3 /home/dirk/medipipe_warp/scripts/qc/generate_qc_report.py --output-dir /data/medipipe_data/output
```

### Collect All QC Data
```shell
python3 /home/dirk/medipipe_warp/scripts/qc/collect_all_qc.py
```

---

## Configuration

### Main Config File
**Location:** `/home/dirk/medipipe_warp/work/config.yaml`

Key parameters:
```yaml
genome: hg38
paired-end: true
use_umi: true
window_size: 300
meth_qc_window_size: 1000  # For MEDIPS (avoids overflow)
threads: 24
output_dir: /data/medipipe_data/output
```

### Samples File
**Location:** `/home/dirk/medipipe_warp/work/samples.tsv`

Format:
```
sample_id	group	R1	R2
2001CSF	CSF	/data/raw_data/raw_data_temp/2001CSF_S1_L004_R1_001.fastq.gz	/data/raw_data/raw_data_temp/2001CSF_S1_L004_R2_001.fastq.gz
```

---

## Machine Learning Analysis

### ML Discrimination
```shell
python3 /home/dirk/medipipe_warp/scripts/ml/ml_exact_permutation.py \
  --matrix {{matrix}} \
  --annotation {{annotation}} \
  --output {{output}} \
  --feature-space {{feature_space}} \
  --case-group {{case_group}} \
  --ctrl-group {{ctrl_group}} \
  --model-type both
```

### ML Visualization
```shell
python3 /home/dirk/medipipe_warp/scripts/ml/ml_visualization.py \
  --matrix "/data/medipipe_data/output/ml_discrimination/matrices/windows/hg38_w2000/all/matrix_autosomal.tsv" \
  --annotation "/data/medipipe_data/output/_upload_package_meta/sample_annotation.active.with_covariates.clean.tsv" \
  --output "/data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/visualizations" \
  --feature-space "autosomal_windows" \
  --case-group "AEG" --ctrl-group "CTRL" --model-type both
```

### NatureCancer Panel Scoring
```shell
python3 /home/dirk/medipipe_warp/scripts/ml/naturecancer_panel_scoring.py
```

---

## Troubleshooting

### Common Issues

#### 1. Integer Overflow in MEDIPS
**Error:** `NAs produced by integer overflow`
**Solution:** Use `meth_qc_window_size: 1000` in config.yaml

#### 2. Filesystem Latency (HDD)
**Error:** `MissingOutputException` after job completes
**Solution:** Add `--latency-wait 60` to snakemake command

#### 3. QSEA Taking Too Long
**Issue:** QSEA whole-genome analysis takes 24+ hours per sample
**Solution:** QSEA is now disabled by default (`--skip_qsea` flag)

#### 4. Conda Environment Issues
**Error:** `libmamba Non-conda folder exists at prefix`
**Solution:** 
```shell
rm -rf /home/dirk/medipipe_warp/.snakemake/conda/*
# Use --conda-frontend conda instead of mamba
```

### Check Workflow Status
```shell
# Is Snakemake running?
ps aux | grep snakemake | grep -v grep

# Latest workflow log
tail -50 $(ls -t /home/dirk/medipipe_warp/work/logs/workflow_*.log | head -1)

# Check specific job logs
tail -20 /home/dirk/medipipe_warp/work/logs/meth_quant/*.log
```

### Restart After Failure
```shell
cd /home/dirk/medipipe_warp
snakemake --unlock --snakefile snakefiles/core_workflow.smk --configfile work/config.yaml
snakemake --snakefile snakefiles/core_workflow.smk --configfile work/config.yaml \
  --cores 24 --use-conda --rerun-incomplete --keep-going --conda-frontend conda --latency-wait 60
```

---

## GUI Features

### Tabs
1. **Workflow:** Select pipeline stages, UMI/spike-in options
2. **Parameters:** Genome, window size, threads, directories
3. **QC Reports:** View sample QC status, metrics, reports
4. **ML Analysis:** Feature selection, model training, visualization

### Launch GUI
```shell
# Option 1: Shell script
/home/dirk/medipipe_warp/launch_medipipe.sh

# Option 2: Direct Python
python3 /home/dirk/medipipe_warp/gui/main.py
```

---

## Scripts Reference

### Core Scripts
| Script | Purpose |
|--------|---------|
| `scripts/core/medips_analysis.R` | MEDIPS methylation quantification |
| `scripts/core/nnls_deconv.py` | Tissue of origin deconvolution |

### QC Scripts
| Script | Purpose |
|--------|---------|
| `scripts/qc/generate_qc_report.py` | Comprehensive QC report |
| `scripts/qc/collect_all_qc.py` | Collect all QC metrics |

### ML Scripts
| Script | Purpose |
|--------|---------|
| `scripts/ml/ml_exact_permutation.py` | Permutation testing |
| `scripts/ml/ml_visualization.py` | Plots and figures |
| `scripts/ml/naturecancer_panel_scoring.py` | Panel z-scores |

---

## Warp Workflows

Available workflows (access via `Ctrl+Shift+R` → search "medipipe"):

1. **MEDIPIPE - Run Core Workflow**
2. **MEDIPIPE - Dry Run**
3. **MEDIPIPE - Monitor Progress**
4. **MEDIPIPE - Stop Workflow**
5. **MEDIPIPE - Generate QC Report**
6. **MEDIPIPE - Collect All QC Data**
7. **MEDIPIPE - ML Visualization**
8. **MEDIPIPE - Launch GUI**

---

## Sample Processing Status

### Current Samples (March 2026)
- 2001CSF - 2010CSF (10 samples)
- Processing: MEDIPS methylation quantification

### Completed Samples
- 148, 149, 150, 152, 154, 155 (older batch)
- C1, C3, C4, C5, C6 (controls)
- Various CSF samples (014-CSF, 015-CSF, etc.)

---

## Version Information

- **Pipeline Version:** MEDIPIPE v2.0
- **Snakemake:** 7.x
- **R:** 4.x (with MEDIPS, BSgenome)
- **Python:** 3.10+
- **Genome:** hg38 (UCSC)

---

## Contact & Support

- **MCP Server:** Available for AI-assisted analysis
- **Documentation:** `/home/dirk/medipipe_warp/mcp/project/scripts.yaml`
- **Warp Workflows:** `/home/dirk/medipipe_warp/.warp/workflows/`
