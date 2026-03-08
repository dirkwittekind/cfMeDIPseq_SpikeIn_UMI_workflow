"""
MEDIPIPE Configuration Panel Component
Panel for setting analysis parameters.
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog
from typing import Dict, Callable, Optional


class ConfigPanel(ttk.LabelFrame):
    """Panel for configuring analysis parameters"""
    
    def __init__(
        self,
        parent,
        on_change: Optional[Callable[[Dict], None]] = None,
        **kwargs
    ):
        super().__init__(parent, text="Analysis Parameters", padding=10, **kwargs)
        
        self.on_change = on_change
        
        # Variables
        self.vars = {
            "genome": tk.StringVar(value="hg38"),
            "window_size": tk.IntVar(value=300),
            "threads": tk.IntVar(value=24),
            "input_dir": tk.StringVar(value="/data/raw_data"),
            "output_dir": tk.StringVar(value="/data/medipipe_data/output"),
            "samples_file": tk.StringVar(value=""),
        }
        
        self._create_widgets()
    
    def _create_widgets(self):
        """Create panel widgets"""
        # Notebook for organization
        notebook = ttk.Notebook(self)
        notebook.pack(fill=tk.BOTH, expand=True)
        
        # Tab 1: General settings
        general_frame = ttk.Frame(notebook, padding=10)
        notebook.add(general_frame, text="General")
        
        row = 0
        
        # Genome selection
        ttk.Label(general_frame, text="Reference Genome:").grid(row=row, column=0, sticky="w", pady=5)
        genome_combo = ttk.Combobox(general_frame, textvariable=self.vars["genome"],
                                    values=["hg38", "hg19", "mm10", "mm39"], state="readonly", width=15)
        genome_combo.grid(row=row, column=1, sticky="w", pady=5, padx=5)
        genome_combo.bind("<<ComboboxSelected>>", self._on_change)
        row += 1
        
        # Window size
        ttk.Label(general_frame, text="Window Size (bp):").grid(row=row, column=0, sticky="w", pady=5)
        window_spin = ttk.Spinbox(general_frame, textvariable=self.vars["window_size"],
                                  from_=100, to=1000, increment=50, width=10)
        window_spin.grid(row=row, column=1, sticky="w", pady=5, padx=5)
        window_spin.bind("<FocusOut>", self._on_change)
        row += 1
        
        # Threads
        ttk.Label(general_frame, text="CPU Threads:").grid(row=row, column=0, sticky="w", pady=5)
        thread_spin = ttk.Spinbox(general_frame, textvariable=self.vars["threads"],
                                  from_=1, to=64, increment=1, width=10)
        thread_spin.grid(row=row, column=1, sticky="w", pady=5, padx=5)
        thread_spin.bind("<FocusOut>", self._on_change)
        row += 1
        
        # Tab 2: Directories
        dir_frame = ttk.Frame(notebook, padding=10)
        notebook.add(dir_frame, text="Directories")
        
        row = 0
        
        # Input directory
        ttk.Label(dir_frame, text="Input Directory:").grid(row=row, column=0, sticky="w", pady=5)
        input_entry = ttk.Entry(dir_frame, textvariable=self.vars["input_dir"], width=40)
        input_entry.grid(row=row, column=1, sticky="ew", pady=5, padx=5)
        ttk.Button(dir_frame, text="Browse", width=8,
                  command=lambda: self._browse_directory("input_dir")).grid(row=row, column=2, pady=5, padx=2)
        row += 1
        
        # Output directory
        ttk.Label(dir_frame, text="Output Directory:").grid(row=row, column=0, sticky="w", pady=5)
        output_entry = ttk.Entry(dir_frame, textvariable=self.vars["output_dir"], width=40)
        output_entry.grid(row=row, column=1, sticky="ew", pady=5, padx=5)
        ttk.Button(dir_frame, text="Browse", width=8,
                  command=lambda: self._browse_directory("output_dir")).grid(row=row, column=2, pady=5, padx=2)
        row += 1
        
        # Samples file
        ttk.Label(dir_frame, text="Samples File:").grid(row=row, column=0, sticky="w", pady=5)
        samples_entry = ttk.Entry(dir_frame, textvariable=self.vars["samples_file"], width=40)
        samples_entry.grid(row=row, column=1, sticky="ew", pady=5, padx=5)
        ttk.Button(dir_frame, text="Browse", width=8,
                  command=self._browse_samples_file).grid(row=row, column=2, pady=5, padx=2)
        row += 1
        
        # Info labels
        ttk.Label(dir_frame, text="Note: Large files (FASTQ, BAM) should be on HDD (/data/)",
                  foreground="gray").grid(row=row, column=0, columnspan=3, sticky="w", pady=10)
        
        dir_frame.columnconfigure(1, weight=1)
        
        # Tab 3: QC Settings
        qc_frame = ttk.Frame(notebook, padding=10)
        notebook.add(qc_frame, text="QC Options")
        
        self.qc_vars = {
            "run_fastqc": tk.BooleanVar(value=True),
            "run_multiqc": tk.BooleanVar(value=True),
            "run_coverage": tk.BooleanVar(value=True),
            "run_gc_bias": tk.BooleanVar(value=True),
            "run_overamp": tk.BooleanVar(value=True),
            "generate_bigwig": tk.BooleanVar(value=True),
        }
        
        row = 0
        for key, var in self.qc_vars.items():
            label = key.replace("_", " ").replace("run ", "").title()
            ttk.Checkbutton(qc_frame, text=label, variable=var).grid(row=row, column=0, sticky="w", pady=2)
            row += 1
        
        # BigWig settings
        ttk.Separator(qc_frame, orient=tk.HORIZONTAL).grid(row=row, column=0, columnspan=2, sticky="ew", pady=10)
        row += 1
        
        ttk.Label(qc_frame, text="BigWig Settings:", font=("", 9, "bold")).grid(row=row, column=0, sticky="w")
        row += 1
        
        self.bigwig_bin_var = tk.IntVar(value=25)
        ttk.Label(qc_frame, text="Bin Size:").grid(row=row, column=0, sticky="w", padx=20)
        ttk.Spinbox(qc_frame, textvariable=self.bigwig_bin_var, from_=10, to=100, width=8).grid(row=row, column=1, sticky="w")
        row += 1
        
        self.bigwig_norm_var = tk.StringVar(value="CPM")
        ttk.Label(qc_frame, text="Normalization:").grid(row=row, column=0, sticky="w", padx=20)
        ttk.Combobox(qc_frame, textvariable=self.bigwig_norm_var,
                     values=["CPM", "RPKM", "BPM", "None"], state="readonly", width=10).grid(row=row, column=1, sticky="w")
        
        # Tab 4: ML/DL Settings
        ml_frame = ttk.Frame(notebook, padding=10)
        notebook.add(ml_frame, text="ML/DL Settings")
        
        # ML settings
        ml_settings_frame = ttk.LabelFrame(ml_frame, text="Machine Learning", padding=10)
        ml_settings_frame.pack(fill=tk.X, pady=5)
        
        row = 0
        
        self.cv_folds_var = tk.IntVar(value=5)
        ttk.Label(ml_settings_frame, text="CV Folds:").grid(row=row, column=0, sticky="w")
        ttk.Spinbox(ml_settings_frame, textvariable=self.cv_folds_var, from_=3, to=10, width=8).grid(row=row, column=1, sticky="w", padx=5)
        row += 1
        
        self.n_features_var = tk.IntVar(value=1000)
        ttk.Label(ml_settings_frame, text="Num Features:").grid(row=row, column=0, sticky="w")
        ttk.Spinbox(ml_settings_frame, textvariable=self.n_features_var, from_=100, to=10000, increment=100, width=8).grid(row=row, column=1, sticky="w", padx=5)
        row += 1
        
        # DL settings
        dl_settings_frame = ttk.LabelFrame(ml_frame, text="Deep Learning", padding=10)
        dl_settings_frame.pack(fill=tk.X, pady=5)
        
        row = 0
        
        self.epochs_var = tk.IntVar(value=100)
        ttk.Label(dl_settings_frame, text="Epochs:").grid(row=row, column=0, sticky="w")
        ttk.Spinbox(dl_settings_frame, textvariable=self.epochs_var, from_=10, to=500, width=8).grid(row=row, column=1, sticky="w", padx=5)
        row += 1
        
        self.batch_size_var = tk.IntVar(value=32)
        ttk.Label(dl_settings_frame, text="Batch Size:").grid(row=row, column=0, sticky="w")
        ttk.Spinbox(dl_settings_frame, textvariable=self.batch_size_var, from_=8, to=128, width=8).grid(row=row, column=1, sticky="w", padx=5)
        row += 1
        
        self.learning_rate_var = tk.DoubleVar(value=0.001)
        ttk.Label(dl_settings_frame, text="Learning Rate:").grid(row=row, column=0, sticky="w")
        ttk.Entry(dl_settings_frame, textvariable=self.learning_rate_var, width=10).grid(row=row, column=1, sticky="w", padx=5)
        row += 1
        
        self.use_gpu_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(dl_settings_frame, text="Use GPU (if available)", variable=self.use_gpu_var).grid(row=row, column=0, columnspan=2, sticky="w")
    
    def _browse_directory(self, var_name: str):
        """Browse for a directory"""
        current = self.vars[var_name].get()
        initial = current if os.path.exists(current) else "/"
        
        path = filedialog.askdirectory(initialdir=initial, title=f"Select {var_name.replace('_', ' ').title()}")
        if path:
            self.vars[var_name].set(path)
            self._on_change()
    
    def _browse_samples_file(self):
        """Browse for samples file"""
        current = self.vars["samples_file"].get()
        initial_dir = os.path.dirname(current) if current else "/home/dirk/medipipe_warp/samplefiles"
        
        path = filedialog.askopenfilename(
            initialdir=initial_dir,
            title="Select Samples File",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if path:
            self.vars["samples_file"].set(path)
            self._on_change()
    
    def _on_change(self, event=None):
        """Handle configuration changes"""
        if self.on_change:
            self.on_change(self.get_config())
    
    def get_config(self) -> Dict:
        """Get current configuration"""
        return {
            "genome": self.vars["genome"].get(),
            "window_size": self.vars["window_size"].get(),
            "threads": self.vars["threads"].get(),
            "input_dir": self.vars["input_dir"].get(),
            "output_dir": self.vars["output_dir"].get(),
            "samples": self.vars["samples_file"].get(),
            
            "qc": {
                "run_fastqc": self.qc_vars["run_fastqc"].get(),
                "run_multiqc": self.qc_vars["run_multiqc"].get(),
                "run_coverage": self.qc_vars["run_coverage"].get(),
                "run_gc_bias": self.qc_vars["run_gc_bias"].get(),
                "run_overamp": self.qc_vars["run_overamp"].get(),
                "generate_bigwig": self.qc_vars["generate_bigwig"].get(),
                "bigwig_bin_size": self.bigwig_bin_var.get(),
                "bigwig_normalization": self.bigwig_norm_var.get(),
            },
            
            "ml": {
                "cv_folds": self.cv_folds_var.get(),
                "n_features": self.n_features_var.get(),
            },
            
            "dl": {
                "epochs": self.epochs_var.get(),
                "batch_size": self.batch_size_var.get(),
                "learning_rate": self.learning_rate_var.get(),
                "use_gpu": self.use_gpu_var.get(),
            },
        }
    
    def set_config(self, config: Dict):
        """Set configuration from dict"""
        # Basic settings
        for key in ["genome", "window_size", "threads", "input_dir", "output_dir"]:
            if key in config and key in self.vars:
                self.vars[key].set(config[key])
        
        if "samples" in config:
            self.vars["samples_file"].set(config["samples"])
        
        # QC settings
        qc_config = config.get("qc", {})
        for key in self.qc_vars:
            if key in qc_config:
                self.qc_vars[key].set(qc_config[key])
        
        if "bigwig_bin_size" in qc_config:
            self.bigwig_bin_var.set(qc_config["bigwig_bin_size"])
        if "bigwig_normalization" in qc_config:
            self.bigwig_norm_var.set(qc_config["bigwig_normalization"])
        
        # ML settings
        ml_config = config.get("ml", {})
        if "cv_folds" in ml_config:
            self.cv_folds_var.set(ml_config["cv_folds"])
        if "n_features" in ml_config:
            self.n_features_var.set(ml_config["n_features"])
        
        # DL settings
        dl_config = config.get("dl", {})
        if "epochs" in dl_config:
            self.epochs_var.set(dl_config["epochs"])
        if "batch_size" in dl_config:
            self.batch_size_var.set(dl_config["batch_size"])
        if "learning_rate" in dl_config:
            self.learning_rate_var.set(dl_config["learning_rate"])
        if "use_gpu" in dl_config:
            self.use_gpu_var.set(dl_config["use_gpu"])
