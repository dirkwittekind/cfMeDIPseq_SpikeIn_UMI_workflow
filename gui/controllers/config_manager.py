"""
MEDIPIPE Configuration Manager
Handles loading, saving, and validating configuration files.
"""

import os
import yaml
import shutil
from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field


@dataclass
class MedipipeConfig:
    """Data class for MEDIPIPE configuration"""
    # Directories
    pipe_dir: str = "/home/dirk/medipipe_warp"
    resources_dir: str = "/home/dirk/medipipe_warp/resources"
    work_dir: str = "/home/dirk/medipipe_warp/work"
    input_dir: str = "/data/raw_data"
    output_dir: str = "/data/medipipe_data/output"
    local_output_dir: str = "/home/dirk/medipipe_warp/outputs"
    
    # Sample files
    samples_file: str = ""
    
    # Analysis parameters
    paired_end: bool = True
    use_umi: bool = False
    use_spikein: bool = False
    genome: str = "hg38"
    window_size: int = 300
    threads: int = 24
    
    # UMI patterns (Twist adapter defaults: 5bp UMI + 2bp skip)
    umi_regex_r1: str = r"(?P<umi_1>.{5})(?P<discard_1>.{2}).*"
    umi_regex_r2: str = r"(?P<umi_2>.{5})(?P<discard_2>.{2}).*"
    
    # Spike-in settings
    spike_idx: str = ""  # BWA index for combined genome + spike-in
    spike_in_chr: str = ""  # Space-separated spike-in contig names
    spike_in_bsgenome: str = ""  # BSgenome name for spike-in (e.g., BSgenome.Hsapiens.UCSC.hg38plusspikeins)
    spike_in_bsgenome_pkg: str = ""  # Path to BSgenome tarball
    
    # Workflow stages
    run_preprocessing: bool = True
    run_qc: bool = True
    run_methylation: bool = True
    run_fragment_profile: bool = True
    run_tissue_of_origin: bool = True
    run_dmr: bool = False
    run_ml: bool = False
    run_dl: bool = False
    
    # Group comparison
    group1: str = ""  # Case group
    group2: str = ""  # Control group
    
    # ML Discrimination settings
    ml_input_source: str = "bigwig"  # "bigwig" or "dedup_bam"
    ml_window_size: int = 2000
    ml_kbest: int = 5000
    ml_c: float = 1.0
    ml_mode: str = "methylOnly"  # "methylOnly" or "methylPlusCov"
    ml_covariates: List[str] = field(default_factory=list)
    ml_rskfold_splits: int = 5
    ml_rskfold_repeats: int = 50
    ml_seed: int = 1
    ml_run_dmr: bool = True
    ml_run_univariate: bool = True
    
    # Reference files
    ref_files: Dict[str, str] = field(default_factory=dict)


class ConfigManager:
    """Manages MEDIPIPE configuration"""
    
    def __init__(self, pipe_dir: str = "/home/dirk/medipipe_warp"):
        self.pipe_dir = pipe_dir
        self.config_dir = os.path.join(pipe_dir, "configfiles")
        self.template_path = os.path.join(self.config_dir, "config_template.yaml")
        self.current_config: Optional[MedipipeConfig] = None
        self._config_path: Optional[str] = None
    
    def load_template(self) -> Dict[str, Any]:
        """Load the configuration template"""
        if os.path.exists(self.template_path):
            with open(self.template_path, 'r') as f:
                return yaml.safe_load(f)
        return {}
    
    def load_config(self, config_path: str) -> MedipipeConfig:
        """Load configuration from a YAML file"""
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        self._config_path = config_path
        self.current_config = self._dict_to_config(config_dict)
        return self.current_config
    
    def save_config(self, config: MedipipeConfig, config_path: Optional[str] = None) -> str:
        """Save configuration to a YAML file"""
        if config_path is None:
            config_path = self._config_path or os.path.join(self.config_dir, "config.yaml")
        
        config_dict = self._config_to_dict(config)
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        
        with open(config_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, indent=2, sort_keys=False)
        
        self._config_path = config_path
        self.current_config = config
        return config_path
    
    def create_working_config(self, config: MedipipeConfig) -> str:
        """Create a working config.yaml in the work directory for Snakemake"""
        work_config_path = os.path.join(config.work_dir, "config.yaml")
        
        # Ensure work directory exists
        os.makedirs(config.work_dir, exist_ok=True)
        
        return self.save_config(config, work_config_path)
    
    def _dict_to_config(self, d: Dict[str, Any]) -> MedipipeConfig:
        """Convert dictionary to MedipipeConfig"""
        config = MedipipeConfig()
        
        # Map dictionary keys to config attributes
        config.pipe_dir = d.get("pipe_dir", config.pipe_dir)
        config.resources_dir = d.get("resources_dir", config.resources_dir)
        config.work_dir = d.get("work_dir", config.work_dir)
        config.input_dir = d.get("input_dir", config.input_dir)
        config.output_dir = d.get("output_dir", config.output_dir)
        config.local_output_dir = d.get("local_output_dir", config.local_output_dir)
        
        config.samples_file = d.get("samples", "")
        
        config.paired_end = d.get("paired_end", True)
        config.use_umi = d.get("use_umi", False)
        config.use_spikein = d.get("use_spikein", False)
        config.genome = d.get("genome", "hg38")
        config.window_size = d.get("window_size", 300)
        config.threads = d.get("threads", 24)
        
        # UMI patterns
        config.umi_regex_r1 = d.get("umi_regex_r1", r"(?P<umi_1>.{5})(?P<discard_1>.{2}).*")
        config.umi_regex_r2 = d.get("umi_regex_r2", r"(?P<umi_2>.{5})(?P<discard_2>.{2}).*")
        
        # Spike-in settings
        config.spike_idx = d.get("spike_idx", "")
        config.spike_in_chr = d.get("spike_in_chr", "")
        config.spike_in_bsgenome = d.get("spike_in_bsgenome", "")
        config.spike_in_bsgenome_pkg = d.get("spike_in_bsgenome_pkg", "")
        
        # Workflow stages
        workflow = d.get("workflow", {})
        config.run_preprocessing = workflow.get("preprocessing", True)
        config.run_qc = workflow.get("qc_analysis", True)
        config.run_methylation = workflow.get("methylation_quant", True)
        config.run_fragment_profile = workflow.get("fragment_profile", True)
        config.run_tissue_of_origin = workflow.get("tissue_of_origin", True)
        
        advanced = d.get("advanced", {})
        config.run_dmr = advanced.get("dmr_analysis", False)
        config.run_ml = advanced.get("ml_classification", False)
        config.run_dl = advanced.get("dl_classification", False)
        
        config.ref_files = d.get("ref_files", {})
        
        # Group comparison
        comparison = d.get("comparison", {})
        config.group1 = comparison.get("group1", "")
        config.group2 = comparison.get("group2", "")
        
        # ML Discrimination settings
        ml_discrim = d.get("ml_discrimination", {})
        config.ml_input_source = ml_discrim.get("input_source", "bigwig")
        config.ml_window_size = ml_discrim.get("window_size", 2000)
        config.ml_kbest = ml_discrim.get("kbest", 5000)
        config.ml_c = ml_discrim.get("C", 1.0)
        config.ml_mode = ml_discrim.get("mode", "methylOnly")
        config.ml_covariates = ml_discrim.get("covariates", [])
        config.ml_rskfold_splits = ml_discrim.get("rskfold_splits", 5)
        config.ml_rskfold_repeats = ml_discrim.get("rskfold_repeats", 50)
        config.ml_seed = ml_discrim.get("seed", 1)
        config.ml_run_dmr = ml_discrim.get("run_dmr", True)
        config.ml_run_univariate = ml_discrim.get("run_univariate", True)
        
        return config
    
    def _config_to_dict(self, config: MedipipeConfig) -> Dict[str, Any]:
        """Convert MedipipeConfig to dictionary for YAML"""
        # Get reference files based on genome
        ref_files = config.ref_files or self._get_default_ref_files(config.genome, config.resources_dir)
        
        return {
            "pipe_dir": config.pipe_dir,
            "resources_dir": config.resources_dir,
            "work_dir": config.work_dir,
            "input_dir": config.input_dir,
            "output_dir": config.output_dir,
            "local_output_dir": config.local_output_dir,
            
            "samples": config.samples_file,
            "samples_aggr": config.samples_file.replace(".tsv", "_aggr.tsv") if config.samples_file else "",
            
            "paired_end": config.paired_end,
            "use_umi": config.use_umi,
            "add_umi": config.use_umi,  # Legacy backward compatibility
            "umi": config.use_umi,  # Legacy backward compatibility
            "use_spikein": config.use_spikein,
            "spike_in": config.use_spikein,  # Legacy backward compatibility
            "spike_idx": config.spike_idx,
            "spike_in_chr": config.spike_in_chr,
            "spike_in_bsgenome": config.spike_in_bsgenome,
            "spike_in_bsgenome_pkg": config.spike_in_bsgenome_pkg,
            "umi_regex_r1": config.umi_regex_r1,
            "umi_regex_r2": config.umi_regex_r2,
            "umi_pattern": config.umi_regex_r1,  # Legacy fallback
            "genome": config.genome,
            "bsgenome": f"BSgenome.Hsapiens.UCSC.{config.genome}",
            "window_size": config.window_size,
            "threads": config.threads,
            
            "ref_files": ref_files,
            
            "frag_profile": config.run_fragment_profile,
            "aggregate": False,
            
            "workflow": {
                "preprocessing": config.run_preprocessing,
                "qc_analysis": config.run_qc,
                "methylation_quant": config.run_methylation,
                "fragment_profile": config.run_fragment_profile,
                "tissue_of_origin": config.run_tissue_of_origin,
            },
            
            "advanced": {
                "dmr_analysis": config.run_dmr,
                "ml_classification": config.run_ml,
                "dl_classification": config.run_dl,
            },
            
            "tissue_of_origin": {
                "enable": config.run_tissue_of_origin,
                "signature_tsv": os.path.join(config.resources_dir, "tissue_of_origin/reference_methylation_matrix.tsv"),
                "markers_bed": os.path.join(config.resources_dir, "tissue_of_origin/marker_regions.bed"),
                "scale": "minmax_by_marker",
            },
            
            "qc": {
                "run_fastqc": True,
                "run_multiqc": True,
                "run_coverage": True,
                "run_gc_bias": True,
                "run_overamp": True,
                "generate_bigwig": True,
                "bigwig_bin_size": 25,
                "bigwig_normalization": "CPM",
            },
            
            "ml": {
                "algorithms": ["random_forest", "gradient_boosting", "svm", "logistic_regression"],
                "cv_folds": 5,
                "feature_selection": "variance_threshold",
                "n_features": 1000,
            },
            
            "dl": {
                "models": ["cnn_1d", "transformer"],
                "epochs": 100,
                "batch_size": 32,
                "learning_rate": 0.001,
                "early_stopping_patience": 10,
                "use_gpu": True,
            },
            
            "comparison": {
                "group1": config.group1,
                "group2": config.group2,
            },
            
            "ml_discrimination": {
                "input_source": config.ml_input_source,
                "window_size": config.ml_window_size,
                "kbest": config.ml_kbest,
                "C": config.ml_c,
                "mode": config.ml_mode,
                "covariates": config.ml_covariates,
                "rskfold_splits": config.ml_rskfold_splits,
                "rskfold_repeats": config.ml_rskfold_repeats,
                "seed": config.ml_seed,
                "run_dmr": config.ml_run_dmr,
                "run_univariate": config.ml_run_univariate,
                "effective_genome_size": 2913022398,  # hg38 default
                "promoter_bed": os.path.join(config.resources_dir, "annotations/hg38_promoters_2kb.bed"),
                "gene_bed": os.path.join(config.resources_dir, "annotations/gencode_v46_genes.bed"),
            },
        }
    
    def _get_default_ref_files(self, genome: str, resources_dir: str) -> Dict[str, str]:
        """Get default reference file paths for a genome"""
        genome_dir = os.path.join(resources_dir, "genomes", genome)
        return {
            "fasta": os.path.join(genome_dir, f"{genome}.fa"),
            "bwa_index": os.path.join(genome_dir, f"bwa_index/{genome}.fa"),
            "bowtie2_index": os.path.join(genome_dir, f"bowtie2_index/{genome}"),
            "chrom_sizes": os.path.join(genome_dir, f"{genome}.chrom.sizes"),
            "blacklist": os.path.join(genome_dir, f"{genome}.blacklist.bed"),
            "gc_bias": os.path.join(genome_dir, "gc_bias.tsv"),
            "rmsk": os.path.join(genome_dir, "rmsk.bed"),
        }
    
    def validate_config(self, config: MedipipeConfig) -> List[str]:
        """Validate configuration and return list of errors"""
        errors = []
        
        # Check directories
        if not os.path.exists(config.input_dir):
            errors.append(f"Input directory does not exist: {config.input_dir}")
        
        if not os.path.exists(config.resources_dir):
            errors.append(f"Resources directory does not exist: {config.resources_dir}")
        
        # Check samples file
        if config.samples_file and not os.path.exists(config.samples_file):
            errors.append(f"Samples file does not exist: {config.samples_file}")
        
        # Check reference files
        for name, path in config.ref_files.items():
            if path and not os.path.exists(path):
                errors.append(f"Reference file '{name}' not found: {path}")
        
        return errors
