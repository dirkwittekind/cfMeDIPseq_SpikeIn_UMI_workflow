"""
MEDIPIPE Workflow Panel Component
Panel for selecting and configuring workflow stages.
"""

import tkinter as tk
from tkinter import ttk
from typing import Dict, Callable, Optional


class WorkflowPanel(ttk.LabelFrame):
    """Panel for configuring workflow stages and options"""
    
    def __init__(
        self,
        parent,
        on_change: Optional[Callable[[Dict], None]] = None,
        **kwargs
    ):
        super().__init__(parent, text="Workflow Configuration", padding=10, **kwargs)
        
        self.on_change = on_change
        
        # Variables for workflow options
        self.vars = {
            # Core workflow stages
            "preprocessing": tk.BooleanVar(value=True),
            "qc_analysis": tk.BooleanVar(value=True),
            "methylation_quant": tk.BooleanVar(value=True),
            "fragment_profile": tk.BooleanVar(value=True),
            "tissue_of_origin": tk.BooleanVar(value=True),
            # Advanced workflow stages
            "dmr_analysis": tk.BooleanVar(value=False),
            "ml_classification": tk.BooleanVar(value=False),
            "dl_classification": tk.BooleanVar(value=False),
            # Data options
            "paired_end": tk.BooleanVar(value=True),
            "use_umi": tk.BooleanVar(value=False),
            "use_spikein": tk.BooleanVar(value=False),
        }
        
        self._create_widgets()
        self._setup_bindings()
    
    def _create_widgets(self):
        """Create panel widgets"""
        # Main container
        notebook = ttk.Notebook(self)
        notebook.pack(fill=tk.BOTH, expand=True)
        
        # Tab 1: Core Workflow
        core_frame = ttk.Frame(notebook, padding=10)
        notebook.add(core_frame, text="Core Workflow")
        
        # Workflow stages section
        stages_frame = ttk.LabelFrame(core_frame, text="Workflow Stages", padding=10)
        stages_frame.pack(fill=tk.X, pady=(0, 10))
        
        core_stages = [
            ("preprocessing", "1. Preprocessing", "QC, trimming, UMI extraction"),
            ("qc_analysis", "2. QC Analysis", "Coverage, GC bias, overamplification"),
            ("methylation_quant", "3. Methylation Quantification", "MEDIPS/MeDEStrand/QSEA analysis"),
            ("fragment_profile", "4. Fragment Profile", "cfDNA fragment size analysis"),
            ("tissue_of_origin", "5. Tissue of Origin", "NNLS deconvolution"),
        ]
        
        for i, (key, label, desc) in enumerate(core_stages):
            cb = ttk.Checkbutton(stages_frame, text=label, variable=self.vars[key])
            cb.grid(row=i, column=0, sticky="w", padx=5, pady=3)
            ttk.Label(stages_frame, text=desc, style="Desc.TLabel").grid(row=i, column=1, sticky="w", padx=15)
        
        # Data Options section (moved from separate tab)
        data_frame = ttk.LabelFrame(core_frame, text="Data Options", padding=10)
        data_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Sequencing type
        ttk.Checkbutton(data_frame, text="Paired-end sequencing", 
                        variable=self.vars["paired_end"]).grid(row=0, column=0, sticky="w", padx=5, pady=3)
        ttk.Label(data_frame, text="Uncheck for single-end data", 
                  style="Desc.TLabel").grid(row=0, column=1, sticky="w", padx=15)
        
        # UMI options
        ttk.Checkbutton(data_frame, text="Use UMI-based deduplication", 
                        variable=self.vars["use_umi"]).grid(row=1, column=0, sticky="w", padx=5, pady=3)
        ttk.Label(data_frame, text="Extract UMIs from read sequences (umi_tools)", 
                  style="Desc.TLabel").grid(row=1, column=1, sticky="w", padx=15)
        
        # Spike-in options
        ttk.Checkbutton(data_frame, text="Use spike-in controls", 
                        variable=self.vars["use_spikein"]).grid(row=2, column=0, sticky="w", padx=5, pady=3)
        ttk.Label(data_frame, text="Enable spike-in normalization", 
                  style="Desc.TLabel").grid(row=2, column=1, sticky="w", padx=15)
        
        # UMI Pattern section (for Twist UMI Adapter: 5bp UMI + 2bp skip per read)
        umi_pattern_frame = ttk.LabelFrame(data_frame, text="UMI Patterns (Twist Adapter: 5bp UMI + 2bp skip)", padding=5)
        umi_pattern_frame.grid(row=3, column=0, columnspan=2, sticky="ew", padx=5, pady=(10, 5))
        
        # R1 pattern
        r1_frame = ttk.Frame(umi_pattern_frame)
        r1_frame.pack(fill=tk.X, pady=2)
        ttk.Label(r1_frame, text="R1 Pattern:", width=12).pack(side=tk.LEFT)
        self.umi_pattern_r1_var = tk.StringVar(value=r"(?P<umi_1>.{5})(?P<discard_1>.{2}).*")
        self.umi_pattern_r1_entry = ttk.Entry(r1_frame, textvariable=self.umi_pattern_r1_var, width=40)
        self.umi_pattern_r1_entry.pack(side=tk.LEFT, padx=5)
        
        # R2 pattern
        r2_frame = ttk.Frame(umi_pattern_frame)
        r2_frame.pack(fill=tk.X, pady=2)
        ttk.Label(r2_frame, text="R2 Pattern:", width=12).pack(side=tk.LEFT)
        self.umi_pattern_r2_var = tk.StringVar(value=r"(?P<umi_2>.{5})(?P<discard_2>.{2}).*")
        self.umi_pattern_r2_entry = ttk.Entry(r2_frame, textvariable=self.umi_pattern_r2_var, width=40)
        self.umi_pattern_r2_entry.pack(side=tk.LEFT, padx=5)
        
        # Help text
        ttk.Label(umi_pattern_frame, 
                  text="Pattern: (?P<umi_N>.{X}) = X bp UMI, (?P<discard_N>.{Y}) = Y bp skip",
                  style="Desc.TLabel").pack(anchor="w", pady=(5, 0))
        
        # Tab 2: Advanced Workflow
        adv_frame = ttk.Frame(notebook, padding=10)
        notebook.add(adv_frame, text="Advanced Analysis")
        
        # Analysis stages
        stages_adv_frame = ttk.LabelFrame(adv_frame, text="Analysis Stages", padding=10)
        stages_adv_frame.pack(fill=tk.X, pady=(0, 10))
        
        adv_stages = [
            ("dmr_analysis", "DMR Analysis", "Differential methylation regions"),
            ("ml_classification", "ML Classification", "Random Forest, SVM, Gradient Boosting"),
            ("dl_classification", "DL Classification", "CNN, Transformer models"),
        ]
        
        for i, (key, label, desc) in enumerate(adv_stages):
            cb = ttk.Checkbutton(stages_adv_frame, text=label, variable=self.vars[key])
            cb.grid(row=i, column=0, sticky="w", padx=5, pady=3)
            ttk.Label(stages_adv_frame, text=desc, style="Desc.TLabel").grid(row=i, column=1, sticky="w", padx=15)
        
        # Group Comparison section
        comparison_frame = ttk.LabelFrame(adv_frame, text="Group Comparison", padding=10)
        comparison_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(comparison_frame, text="Select groups to compare for DMR/ML/DL analysis:",
                  font=("Helvetica", 10)).pack(anchor="w", pady=(0, 5))
        
        # Group selection frame
        group_sel_frame = ttk.Frame(comparison_frame)
        group_sel_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(group_sel_frame, text="Group 1 (Case):").grid(row=0, column=0, sticky="w", padx=5)
        self.group1_var = tk.StringVar()
        self.group1_combo = ttk.Combobox(group_sel_frame, textvariable=self.group1_var, 
                                         state="readonly", width=20)
        self.group1_combo.grid(row=0, column=1, sticky="w", padx=10, pady=3)
        
        ttk.Label(group_sel_frame, text="Group 2 (Control):").grid(row=1, column=0, sticky="w", padx=5)
        self.group2_var = tk.StringVar()
        self.group2_combo = ttk.Combobox(group_sel_frame, textvariable=self.group2_var, 
                                         state="readonly", width=20)
        self.group2_combo.grid(row=1, column=1, sticky="w", padx=10, pady=3)
        
        ttk.Label(comparison_frame, text="Groups are defined by the 'group' column in samples file",
                  style="Desc.TLabel").pack(anchor="w", pady=(5, 0))
        
        # ML algorithm selection
        ml_frame = ttk.LabelFrame(adv_frame, text="ML Algorithms", padding=10)
        ml_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.ml_algorithms = {
            "random_forest": tk.BooleanVar(value=True),
            "gradient_boosting": tk.BooleanVar(value=True),
            "svm": tk.BooleanVar(value=True),
            "logistic_regression": tk.BooleanVar(value=True),
        }
        
        ml_alg_frame = ttk.Frame(ml_frame)
        ml_alg_frame.pack(fill=tk.X)
        
        for i, (key, var) in enumerate(self.ml_algorithms.items()):
            label = key.replace("_", " ").title()
            ttk.Checkbutton(ml_alg_frame, text=label, variable=var).grid(row=i//2, column=i%2, sticky="w", padx=10, pady=2)
        
        # DL model selection
        dl_frame = ttk.LabelFrame(adv_frame, text="DL Models", padding=10)
        dl_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.dl_models = {
            "cnn_1d": tk.BooleanVar(value=True),
            "transformer": tk.BooleanVar(value=True),
        }
        
        dl_model_frame = ttk.Frame(dl_frame)
        dl_model_frame.pack(fill=tk.X)
        
        for i, (key, var) in enumerate(self.dl_models.items()):
            label = key.replace("_", " ").upper()
            ttk.Checkbutton(dl_model_frame, text=label, variable=var).grid(row=0, column=i, sticky="w", padx=10, pady=2)
        
        # =====================================================================
        # ML Discrimination Settings (for post-analysis workflow)
        # =====================================================================
        ml_discrim_frame = ttk.LabelFrame(adv_frame, text="ML Discrimination Settings", padding=10)
        ml_discrim_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Input source selection
        input_frame = ttk.Frame(ml_discrim_frame)
        input_frame.pack(fill=tk.X, pady=2)
        ttk.Label(input_frame, text="Input Source:", width=15).pack(side=tk.LEFT)
        self.ml_input_source_var = tk.StringVar(value="bigwig")
        input_combo = ttk.Combobox(input_frame, textvariable=self.ml_input_source_var,
                                   values=["bigwig", "dedup_bam"], state="readonly", width=15)
        input_combo.pack(side=tk.LEFT, padx=5)
        ttk.Label(input_frame, text="Use BigWig or deduplicated BAMs", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
        
        # Window size
        window_frame = ttk.Frame(ml_discrim_frame)
        window_frame.pack(fill=tk.X, pady=2)
        ttk.Label(window_frame, text="Window Size (bp):", width=15).pack(side=tk.LEFT)
        self.ml_window_size_var = tk.StringVar(value="2000")
        window_combo = ttk.Combobox(window_frame, textvariable=self.ml_window_size_var,
                                    values=["500", "1000", "2000", "5000", "10000"], width=10)
        window_combo.pack(side=tk.LEFT, padx=5)
        ttk.Label(window_frame, text="Genome binning window size", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
        
        # K-best features
        kbest_frame = ttk.Frame(ml_discrim_frame)
        kbest_frame.pack(fill=tk.X, pady=2)
        ttk.Label(kbest_frame, text="K-best Features:", width=15).pack(side=tk.LEFT)
        self.ml_kbest_var = tk.StringVar(value="5000")
        kbest_combo = ttk.Combobox(kbest_frame, textvariable=self.ml_kbest_var,
                                   values=["1000", "2000", "5000", "10000", "20000"], width=10)
        kbest_combo.pack(side=tk.LEFT, padx=5)
        ttk.Label(kbest_frame, text="Number of top features for ML", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
        
        # Regularization C parameter
        c_frame = ttk.Frame(ml_discrim_frame)
        c_frame.pack(fill=tk.X, pady=2)
        ttk.Label(c_frame, text="Regularization (C):", width=15).pack(side=tk.LEFT)
        self.ml_c_var = tk.StringVar(value="1.0")
        c_entry = ttk.Entry(c_frame, textvariable=self.ml_c_var, width=10)
        c_entry.pack(side=tk.LEFT, padx=5)
        ttk.Label(c_frame, text="Logistic regression regularization", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
        
        # ML Mode
        mode_frame = ttk.Frame(ml_discrim_frame)
        mode_frame.pack(fill=tk.X, pady=2)
        ttk.Label(mode_frame, text="ML Mode:", width=15).pack(side=tk.LEFT)
        self.ml_mode_var = tk.StringVar(value="methylOnly")
        mode_combo = ttk.Combobox(mode_frame, textvariable=self.ml_mode_var,
                                  values=["methylOnly", "methylPlusCov"], state="readonly", width=15)
        mode_combo.pack(side=tk.LEFT, padx=5)
        ttk.Label(mode_frame, text="methylOnly or include covariates", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
        
        # Covariates (for methylPlusCov mode)
        cov_frame = ttk.Frame(ml_discrim_frame)
        cov_frame.pack(fill=tk.X, pady=2)
        ttk.Label(cov_frame, text="Covariates:", width=15).pack(side=tk.LEFT)
        self.ml_covariates_var = tk.StringVar(value="")
        cov_entry = ttk.Entry(cov_frame, textvariable=self.ml_covariates_var, width=30)
        cov_entry.pack(side=tk.LEFT, padx=5)
        ttk.Label(cov_frame, text="e.g., age sex batch (space-separated)", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
        
        # CV settings
        cv_frame = ttk.Frame(ml_discrim_frame)
        cv_frame.pack(fill=tk.X, pady=2)
        ttk.Label(cv_frame, text="RSKFold:", width=15).pack(side=tk.LEFT)
        self.ml_rskfold_splits_var = tk.StringVar(value="5")
        ttk.Entry(cv_frame, textvariable=self.ml_rskfold_splits_var, width=5).pack(side=tk.LEFT)
        ttk.Label(cv_frame, text="splits ×").pack(side=tk.LEFT, padx=2)
        self.ml_rskfold_repeats_var = tk.StringVar(value="50")
        ttk.Entry(cv_frame, textvariable=self.ml_rskfold_repeats_var, width=5).pack(side=tk.LEFT)
        ttk.Label(cv_frame, text="repeats", style="Desc.TLabel").pack(side=tk.LEFT, padx=5)
        
        # Run DMR analysis checkbox
        dmr_frame = ttk.Frame(ml_discrim_frame)
        dmr_frame.pack(fill=tk.X, pady=2)
        self.ml_run_dmr_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(dmr_frame, text="Run DMR Analysis (limma)", 
                        variable=self.ml_run_dmr_var).pack(side=tk.LEFT)
        ttk.Label(dmr_frame, text="Differential methylation regions", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
        
        # Run univariate drivers checkbox
        univ_frame = ttk.Frame(ml_discrim_frame)
        univ_frame.pack(fill=tk.X, pady=2)
        self.ml_run_univariate_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(univ_frame, text="Run Univariate Feature Drivers", 
                        variable=self.ml_run_univariate_var).pack(side=tk.LEFT)
        ttk.Label(univ_frame, text="Rank features by effect size", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=10)
    
    def _setup_bindings(self):
        """Set up variable trace bindings"""
        for var in self.vars.values():
            var.trace_add("write", self._on_var_change)
    
    def _on_var_change(self, *args):
        """Handle variable changes"""
        if self.on_change:
            self.on_change(self.get_config())
    
    def get_config(self) -> Dict:
        """Get current workflow configuration"""
        # Parse covariates from space-separated string
        covariates_str = self.ml_covariates_var.get().strip()
        covariates = covariates_str.split() if covariates_str else []
        
        return {
            "workflow": {
                "preprocessing": self.vars["preprocessing"].get(),
                "qc_analysis": self.vars["qc_analysis"].get(),
                "methylation_quant": self.vars["methylation_quant"].get(),
                "fragment_profile": self.vars["fragment_profile"].get(),
                "tissue_of_origin": self.vars["tissue_of_origin"].get(),
            },
            "advanced": {
                "dmr_analysis": self.vars["dmr_analysis"].get(),
                "ml_classification": self.vars["ml_classification"].get(),
                "dl_classification": self.vars["dl_classification"].get(),
            },
            "data_options": {
                "paired_end": self.vars["paired_end"].get(),
                "use_umi": self.vars["use_umi"].get(),
                "use_spikein": self.vars["use_spikein"].get(),
                "umi_regex_r1": self.umi_pattern_r1_var.get(),
                "umi_regex_r2": self.umi_pattern_r2_var.get(),
            },
            "comparison": {
                "group1": self.group1_var.get(),
                "group2": self.group2_var.get(),
            },
            "ml": {
                "algorithms": [k for k, v in self.ml_algorithms.items() if v.get()],
            },
            "dl": {
                "models": [k for k, v in self.dl_models.items() if v.get()],
            },
            "ml_discrimination": {
                "input_source": self.ml_input_source_var.get(),
                "window_size": int(self.ml_window_size_var.get()),
                "kbest": int(self.ml_kbest_var.get()),
                "C": float(self.ml_c_var.get()),
                "mode": self.ml_mode_var.get(),
                "covariates": covariates,
                "rskfold_splits": int(self.ml_rskfold_splits_var.get()),
                "rskfold_repeats": int(self.ml_rskfold_repeats_var.get()),
                "run_dmr": self.ml_run_dmr_var.get(),
                "run_univariate": self.ml_run_univariate_var.get(),
            },
        }
    
    def set_config(self, config: Dict):
        """Set workflow configuration from dict"""
        # Workflow stages
        workflow = config.get("workflow", {})
        for key in ["preprocessing", "qc_analysis", "methylation_quant", "fragment_profile", "tissue_of_origin"]:
            if key in workflow:
                self.vars[key].set(workflow[key])
        
        # Advanced stages
        advanced = config.get("advanced", {})
        for key in ["dmr_analysis", "ml_classification", "dl_classification"]:
            if key in advanced:
                self.vars[key].set(advanced[key])
        
        # Data options
        data_opts = config.get("data_options", config)  # Try data_options or root
        for key in ["paired_end", "use_umi", "use_spikein"]:
            if key in data_opts:
                self.vars[key].set(data_opts[key])
        
        # UMI patterns (R1/R2)
        if "umi_regex_r1" in data_opts:
            self.umi_pattern_r1_var.set(data_opts["umi_regex_r1"])
        elif "umi_regex_r1" in config:
            self.umi_pattern_r1_var.set(config["umi_regex_r1"])
        
        if "umi_regex_r2" in data_opts:
            self.umi_pattern_r2_var.set(data_opts["umi_regex_r2"])
        elif "umi_regex_r2" in config:
            self.umi_pattern_r2_var.set(config["umi_regex_r2"])
        
        # Group comparison
        comparison = config.get("comparison", {})
        if "group1" in comparison:
            self.group1_var.set(comparison["group1"])
        if "group2" in comparison:
            self.group2_var.set(comparison["group2"])
        
        # ML algorithms
        ml_config = config.get("ml", {})
        if "algorithms" in ml_config:
            for alg in self.ml_algorithms:
                self.ml_algorithms[alg].set(alg in ml_config["algorithms"])
        
        # DL models
        dl_config = config.get("dl", {})
        if "models" in dl_config:
            for model in self.dl_models:
                self.dl_models[model].set(model in dl_config["models"])
        
        # ML Discrimination settings
        ml_discrim = config.get("ml_discrimination", {})
        if "input_source" in ml_discrim:
            self.ml_input_source_var.set(ml_discrim["input_source"])
        if "window_size" in ml_discrim:
            self.ml_window_size_var.set(str(ml_discrim["window_size"]))
        if "kbest" in ml_discrim:
            self.ml_kbest_var.set(str(ml_discrim["kbest"]))
        if "C" in ml_discrim:
            self.ml_c_var.set(str(ml_discrim["C"]))
        if "mode" in ml_discrim:
            self.ml_mode_var.set(ml_discrim["mode"])
        if "covariates" in ml_discrim:
            self.ml_covariates_var.set(" ".join(ml_discrim["covariates"]))
        if "rskfold_splits" in ml_discrim:
            self.ml_rskfold_splits_var.set(str(ml_discrim["rskfold_splits"]))
        if "rskfold_repeats" in ml_discrim:
            self.ml_rskfold_repeats_var.set(str(ml_discrim["rskfold_repeats"]))
        if "run_dmr" in ml_discrim:
            self.ml_run_dmr_var.set(ml_discrim["run_dmr"])
        if "run_univariate" in ml_discrim:
            self.ml_run_univariate_var.set(ml_discrim["run_univariate"])
    
    def update_available_groups(self, groups: list):
        """Update the available groups in the comparison dropdowns"""
        self.group1_combo['values'] = groups
        self.group2_combo['values'] = groups
        
        # Auto-select first two groups if available
        if len(groups) >= 2:
            if not self.group1_var.get():
                self.group1_var.set(groups[0])
            if not self.group2_var.get():
                self.group2_var.set(groups[1])
        elif len(groups) == 1:
            if not self.group1_var.get():
                self.group1_var.set(groups[0])
    
    def get_comparison_groups(self) -> tuple:
        """Get selected comparison groups"""
        return (self.group1_var.get(), self.group2_var.get())
    
    def is_core_workflow_selected(self) -> bool:
        """Check if any core workflow stage is selected"""
        return any([
            self.vars["preprocessing"].get(),
            self.vars["qc_analysis"].get(),
            self.vars["methylation_quant"].get(),
            self.vars["fragment_profile"].get(),
            self.vars["tissue_of_origin"].get(),
        ])
    
    def is_advanced_workflow_selected(self) -> bool:
        """Check if any advanced workflow stage is selected"""
        return any([
            self.vars["dmr_analysis"].get(),
            self.vars["ml_classification"].get(),
            self.vars["dl_classification"].get(),
        ])
