# cfMeDIP-seq Spike-In UMI Workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive Snakemake workflow for cell-free methylated DNA immunoprecipitation sequencing (cfMeDIP-seq) analysis with **UMI-based deduplication** and **spike-in normalization** support.

## Features

- **UMI Extraction & Deduplication**: Supports Twist UMI adapter system (5bp UMI + 2bp skip)
- **Spike-in Normalization**: Wilson spike-in controls for IP efficiency assessment
- **Methylation Quantification**: MEDIPS-based CpG enrichment analysis
- **Tissue-of-Origin Deconvolution**: NNLS-based cell type composition estimation
- **Quality Control**: Comprehensive QC metrics (relH, GoGe, saturation, fragment profiles)
- **Machine Learning**: Optional ML-based sample discrimination

## Workflow Overview

```
Stage 1: Preprocessing
  └── UMI Extraction → Trimming → BWA Alignment → Deduplication

Stage 2: Quality Control
  └── FastQC → Insert Size → GC Bias → Library Complexity

Stage 3: Methylation Quantification
  └── MEDIPS QC/Quant (genome + spike-in)

Stage 4: Fragment Profile
  └── cfDNA fragment size analysis

Stage 5: Tissue-of-Origin
  └── NNLS deconvolution with reference methylation atlas
```

## Requirements

### Software
- Snakemake ≥7.0
- Conda/Mamba
- R ≥4.0 (with MEDIPS, BSgenome)
- Python ≥3.10

### Reference Data
- Human reference genome (hg38)
- Combined hg38 + spike-in reference (for spike-in analysis)
- Custom BSgenome package for spike-in sequences

## Installation

```bash
# Clone the repository
git clone https://github.com/dirkwittekind2025/cfMeDIPseq_SpikeIn_UMI_workflow.git
cd cfMeDIPseq_SpikeIn_UMI_workflow

# Create conda environments (handled automatically by Snakemake)
snakemake --use-conda --conda-create-envs-only -c1
```

## Quick Start

### 1. Prepare Sample Sheet

Create `samplefiles/samples.tsv`:

```tsv
sample_id	group	R1	R2
Sample1	Case	/path/to/Sample1_R1.fastq.gz	/path/to/Sample1_R2.fastq.gz
Sample2	Control	/path/to/Sample2_R1.fastq.gz	/path/to/Sample2_R2.fastq.gz
```

### 2. Configure Workflow

Edit `work/config.yaml`:

```yaml
# Key configuration options
genome: hg38
paired_end: true
use_umi: true           # Enable UMI deduplication
use_spikein: true       # Enable spike-in analysis
threads: 24
output_dir: /path/to/output

# UMI patterns (Twist adapter)
umi_regex_r1: "(?P<umi_1>.{5})(?P<discard_1>.{2}).*"
umi_regex_r2: "(?P<umi_2>.{5})(?P<discard_2>.{2}).*"
```

### 3. Run Workflow

```bash
# Dry run (preview jobs)
snakemake --snakefile snakefiles/core_workflow.smk \
  --configfile work/config.yaml \
  --cores 24 --use-conda --dry-run

# Execute workflow
snakemake --snakefile snakefiles/core_workflow.smk \
  --configfile work/config.yaml \
  --cores 24 --use-conda --conda-frontend conda \
  --rerun-incomplete --keep-going --latency-wait 60
```

## Configuration Options

### UMI Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_umi` | `true` | Enable UMI-based deduplication |
| `umi_regex_r1` | `(?P<umi_1>.{5})(?P<discard_1>.{2}).*` | UMI pattern for R1 |
| `umi_regex_r2` | `(?P<umi_2>.{5})(?P<discard_2>.{2}).*` | UMI pattern for R2 |

### Spike-in Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_spikein` | `true` | Enable spike-in analysis |
| `spike_idx` | - | Path to combined BWA index |
| `spike_in_chr` | - | Space-separated spike-in contig names |
| `spike_in_bsgenome` | - | Custom BSgenome package name |

### Methylation QC

| Parameter | Default | Description |
|-----------|---------|-------------|
| `window_size` | `300` | Window size for methylation quantification |
| `meth_qc_window_size` | `1000` | Window size for MEDIPS QC (larger to avoid overflow) |

## Output Structure

```
output/
├── fastqc/                     # FastQC reports
├── trim_reports/               # Trimming statistics
├── dedup_bam_umi/             # UMI-deduplicated BAM files
├── stats/                      # Alignment and dedup statistics
├── insert_size/               # Insert size metrics
├── meth_quant/                # Methylation quantification
├── meth_quant_spikein/        # Spike-in methylation QC
├── spikein_bam/               # Spike-in extracted BAMs
├── fragment_profile/          # Fragment size profiles
├── tissue_of_origin/          # ToO deconvolution results
└── qc_reports/                # Summary QC reports
```

## QC Metrics & Thresholds

| Metric | Threshold | Description |
|--------|-----------|-------------|
| Dedup reads | >5,000,000 | Minimum deduplicated reads |
| Duplication rate | <90% | Maximum acceptable duplication |
| relH | >1.5, <10.0 | CpG enrichment score |
| GoGe | >1.0 | Observed/expected CpG ratio |
| Spike-in relH | >2.5 | IP efficiency from spike-in |
| Spike-in GoGe | >1.5 | Spike-in specificity |

## Spike-in Controls

This workflow supports Wilson spike-in controls with:
- **52 synthetic DNA fragments** (80bp, 160bp, 320bp)
- **Varying CpG densities** (35%, 50%, 65% GC content)
- **Methylated and unmethylated variants** for each size/density combination

The spike-in QC provides:
- relH: Relative enrichment of methylated CpGs
- GoGe: Genome-over-gene enrichment ratio
- Meth ratio: Ratio of methylated vs unmethylated spike-ins

## Tissue-of-Origin Deconvolution

Uses NNLS (Non-Negative Least Squares) deconvolution with:
- Reference methylation atlas from Loyfer et al. 2023
- Per-marker min-max scaling adapted for cfMeDIP-seq data
- Output: Cell type fractions per sample

## Directory Structure

```
cfMeDIPseq_SpikeIn_UMI_workflow/
├── snakefiles/              # Main Snakemake workflows
│   ├── core_workflow.smk    # Core analysis pipeline
│   └── advanced_workflow.smk # ML/DL extensions
├── rules/                   # Snakemake rule modules
│   ├── core/               # Core workflow rules
│   │   ├── common.smk      # Shared functions
│   │   ├── reads_qc.smk    # UMI extraction, trimming, FastQC
│   │   ├── mapping.smk     # BWA alignment
│   │   ├── deduplication.smk # UMI/standard dedup
│   │   ├── meth_qc_quant.smk # MEDIPS analysis
│   │   ├── spikein.smk     # Spike-in analysis
│   │   ├── fragment_profile.smk # Fragment analysis
│   │   └── tissue_of_origin.smk # ToO deconvolution
│   ├── qc/                 # QC rules
│   └── advanced/           # ML/DL rules
├── scripts/                # Analysis scripts
│   ├── core/              # R and Python analysis scripts
│   ├── qc/                # QC reporting scripts
│   └── ml/                # Machine learning scripts
├── environments/           # Conda environment definitions
├── configfiles/           # Example configurations
├── samplefiles/           # Sample sheet templates
└── resources/             # Reference files (not in repo)
```

## Citation

If you use this workflow, please cite:

- **MEDIPS**: Lienhard et al. (2014) MEDIPS: genome-wide differential coverage analysis of sequencing data derived from DNA enrichment experiments. Bioinformatics.
- **UMI-tools**: Smith et al. (2017) UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research.
- **Tissue-of-Origin**: Loyfer et al. (2023) A DNA methylation atlas of normal human cell types. Nature.

## License

MIT License - see [LICENSE](LICENSE) file.

## Contributing

Contributions welcome! Please submit issues and pull requests.

## Authors

- Dirk Wittekind

---

*This workflow was developed for cfMeDIP-seq analysis of cell-free DNA methylation.*
