"""
MEDIPIPE GUI - ML Analysis Panel
Machine learning discrimination analysis with permutation testing and visualization.
"""

import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import subprocess
import threading
from typing import Optional, Callable


class MLAnalysisPanel(ttk.Frame):
    """Panel for ML discrimination analysis"""
    
    def __init__(self, parent, on_log: Optional[Callable] = None, **kwargs):
        super().__init__(parent, **kwargs)
        self.on_log = on_log
        
        # Paths
        self.pipe_dir = "/home/dirk/medipipe_warp"
        self.ml_scripts_dir = os.path.join(self.pipe_dir, "scripts", "ml")
        self.conda_env = "/home/dirk/.conda_envs/medipipe"
        
        # Variables
        self.matrix_path = tk.StringVar()
        self.output_dir = tk.StringVar(value="/data/medipipe_data/output/ml_discrimination")
        self.group1 = tk.StringVar(value="CTRL")
        self.group2 = tk.StringVar(value="AEG")
        self.n_permutations = tk.IntVar(value=1000)
        self.model_type = tk.StringVar(value="L2")
        self.use_exact_perm = tk.BooleanVar(value=True)
        self.run_alphagenome = tk.BooleanVar(value=False)
        
        # Running process
        self._process = None
        self._running = False
        
        self._create_widgets()
    
    def _create_widgets(self):
        """Create the panel widgets"""
        # Main container with padding
        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)
        
        # Title
        ttk.Label(main, text="ML Discrimination Analysis", 
                  font=("Helvetica", 16, "bold")).pack(anchor=tk.W, pady=(0, 15))
        
        # === Input Section ===
        input_frame = ttk.LabelFrame(main, text="Input Data", padding=10)
        input_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Matrix file
        matrix_row = ttk.Frame(input_frame)
        matrix_row.pack(fill=tk.X, pady=5)
        
        ttk.Label(matrix_row, text="Feature Matrix:", width=15, anchor=tk.W).pack(side=tk.LEFT)
        ttk.Entry(matrix_row, textvariable=self.matrix_path, width=50).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(matrix_row, text="Browse...", command=self._browse_matrix, width=10).pack(side=tk.LEFT, padx=5)
        
        # Auto-detect button
        ttk.Button(input_frame, text="🔍 Auto-detect Latest Matrix", 
                   command=self._auto_detect_matrix).pack(anchor=tk.W, pady=(5, 0))
        
        # Output directory
        output_row = ttk.Frame(input_frame)
        output_row.pack(fill=tk.X, pady=5)
        
        ttk.Label(output_row, text="Output Directory:", width=15, anchor=tk.W).pack(side=tk.LEFT)
        ttk.Entry(output_row, textvariable=self.output_dir, width=50).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(output_row, text="Browse...", command=self._browse_output, width=10).pack(side=tk.LEFT, padx=5)
        
        # === Groups Section ===
        groups_frame = ttk.LabelFrame(main, text="Sample Groups", padding=10)
        groups_frame.pack(fill=tk.X, pady=(0, 10))
        
        groups_row = ttk.Frame(groups_frame)
        groups_row.pack(fill=tk.X)
        
        ttk.Label(groups_row, text="Group 1 (Control):", width=15, anchor=tk.W).pack(side=tk.LEFT)
        ttk.Entry(groups_row, textvariable=self.group1, width=15).pack(side=tk.LEFT, padx=5)
        
        ttk.Label(groups_row, text="  vs  ", font=("Helvetica", 12, "bold")).pack(side=tk.LEFT, padx=10)
        
        ttk.Label(groups_row, text="Group 2 (Case):", width=15, anchor=tk.W).pack(side=tk.LEFT)
        ttk.Entry(groups_row, textvariable=self.group2, width=15).pack(side=tk.LEFT, padx=5)
        
        # === Model Settings ===
        model_frame = ttk.LabelFrame(main, text="Model Settings", padding=10)
        model_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Model type
        model_row = ttk.Frame(model_frame)
        model_row.pack(fill=tk.X, pady=5)
        
        ttk.Label(model_row, text="Model Type:", width=15, anchor=tk.W).pack(side=tk.LEFT)
        
        models = [("L2 (Ridge)", "L2"), ("L1 (Lasso)", "L1"), ("ElasticNet", "ElasticNet")]
        for text, value in models:
            ttk.Radiobutton(model_row, text=text, variable=self.model_type, 
                           value=value).pack(side=tk.LEFT, padx=10)
        
        # Permutation settings
        perm_row = ttk.Frame(model_frame)
        perm_row.pack(fill=tk.X, pady=5)
        
        ttk.Label(perm_row, text="Permutations:", width=15, anchor=tk.W).pack(side=tk.LEFT)
        ttk.Spinbox(perm_row, textvariable=self.n_permutations, from_=100, to=10000, 
                    increment=100, width=10).pack(side=tk.LEFT, padx=5)
        
        ttk.Checkbutton(perm_row, text="Use exact permutation (recommended for n<15)", 
                        variable=self.use_exact_perm).pack(side=tk.LEFT, padx=20)
        
        # === Advanced Options ===
        adv_frame = ttk.LabelFrame(main, text="Advanced Options", padding=10)
        adv_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(adv_frame, text="Run AlphaGenome analysis on top features (requires API key)", 
                        variable=self.run_alphagenome).pack(anchor=tk.W)
        
        ttk.Label(adv_frame, text="Note: Set ALPHAGENOME_API_KEY environment variable before running",
                  font=("Helvetica", 10), foreground="gray").pack(anchor=tk.W, pady=(5, 0))
        
        # === Action Buttons ===
        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill=tk.X, pady=15)
        
        self.run_btn = ttk.Button(btn_frame, text="▶  Run ML Analysis", 
                                  command=self._run_analysis, width=20,
                                  style="Run.TButton")
        self.run_btn.pack(side=tk.LEFT, padx=5)
        
        self.stop_btn = ttk.Button(btn_frame, text="⬛ Stop", 
                                   command=self._stop_analysis, width=12,
                                   style="Stop.TButton", state=tk.DISABLED)
        self.stop_btn.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(btn_frame, text="📊 View Results", 
                   command=self._view_results, width=15).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(btn_frame, text="📁 Open Output Folder", 
                   command=self._open_output_folder, width=18).pack(side=tk.RIGHT, padx=5)
        
        # === Log Output ===
        log_frame = ttk.LabelFrame(main, text="Analysis Log", padding=10)
        log_frame.pack(fill=tk.BOTH, expand=True)
        
        self.log_text = tk.Text(log_frame, height=10, wrap=tk.WORD, 
                                font=("Monospace", 11), state=tk.DISABLED)
        log_scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=log_scroll.set)
        
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        log_scroll.pack(side=tk.RIGHT, fill=tk.Y)
    
    def _log(self, message: str):
        """Add message to log"""
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)
        
        if self.on_log:
            self.on_log(message)
    
    def _browse_matrix(self):
        """Browse for matrix file"""
        filepath = filedialog.askopenfilename(
            initialdir="/data/medipipe_data/output/ml_discrimination/matrices",
            title="Select Feature Matrix",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if filepath:
            self.matrix_path.set(filepath)
    
    def _browse_output(self):
        """Browse for output directory"""
        dirpath = filedialog.askdirectory(
            initialdir="/data/medipipe_data/output",
            title="Select Output Directory"
        )
        if dirpath:
            self.output_dir.set(dirpath)
    
    def _auto_detect_matrix(self):
        """Auto-detect the latest feature matrix"""
        base_dir = "/data/medipipe_data/output/ml_discrimination/matrices"
        
        # Look for matrix files
        candidates = []
        for root, dirs, files in os.walk(base_dir):
            for f in files:
                if f.endswith(".tsv") and "matrix" in f.lower():
                    full_path = os.path.join(root, f)
                    candidates.append((os.path.getmtime(full_path), full_path))
        
        if candidates:
            # Sort by modification time (newest first)
            candidates.sort(reverse=True)
            latest = candidates[0][1]
            self.matrix_path.set(latest)
            self._log(f"Auto-detected matrix: {latest}")
        else:
            self._log("No matrix files found in standard location")
            messagebox.showwarning("Not Found", 
                "No matrix files found.\nPlease run the core workflow first to generate feature matrices.")
    
    def _run_analysis(self):
        """Run the ML analysis"""
        # Validate inputs
        if not self.matrix_path.get():
            messagebox.showerror("Error", "Please select a feature matrix file")
            return
        
        if not os.path.exists(self.matrix_path.get()):
            messagebox.showerror("Error", f"Matrix file not found:\n{self.matrix_path.get()}")
            return
        
        # Create output directory
        output_dir = self.output_dir.get()
        os.makedirs(output_dir, exist_ok=True)
        
        # Build command
        script_path = os.path.join(self.ml_scripts_dir, "ml_exact_permutation.py")
        
        if not os.path.exists(script_path):
            messagebox.showerror("Error", f"ML script not found:\n{script_path}")
            return
        
        cmd = [
            "python", script_path,
            "--matrix", self.matrix_path.get(),
            "--output", output_dir,
            "--group1", self.group1.get(),
            "--group2", self.group2.get(),
            "--model", self.model_type.get(),
            "--n_permutations", str(self.n_permutations.get())
        ]
        
        if self.use_exact_perm.get():
            cmd.append("--exact")
        
        self._log(f"Starting ML analysis...")
        self._log(f"Model: {self.model_type.get()}")
        self._log(f"Groups: {self.group1.get()} vs {self.group2.get()}")
        self._log(f"Permutations: {self.n_permutations.get()}")
        self._log("-" * 50)
        
        # Run in background thread
        self._running = True
        self.run_btn.configure(state=tk.DISABLED)
        self.stop_btn.configure(state=tk.NORMAL)
        
        thread = threading.Thread(target=self._run_command, args=(cmd,))
        thread.daemon = True
        thread.start()
    
    def _run_command(self, cmd):
        """Run command in background"""
        try:
            # Activate conda environment
            env = os.environ.copy()
            env["PATH"] = f"{self.conda_env}/bin:" + env.get("PATH", "")
            
            self._process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
                cwd=self.pipe_dir
            )
            
            # Stream output
            for line in self._process.stdout:
                if not self._running:
                    break
                self.after(0, lambda l=line: self._log(l.rstrip()))
            
            self._process.wait()
            
            if self._process.returncode == 0:
                self.after(0, lambda: self._log("\n✓ ML analysis completed successfully!"))
                
                # Run visualization if successful
                if self._running:
                    self.after(0, self._run_visualization)
                
                # Run AlphaGenome if requested
                if self._running and self.run_alphagenome.get():
                    self.after(100, self._run_alphagenome_analysis)
            else:
                self.after(0, lambda: self._log(f"\n✗ Analysis failed (exit code {self._process.returncode})"))
            
        except Exception as e:
            self.after(0, lambda: self._log(f"\nError: {e}"))
        
        finally:
            self._running = False
            self.after(0, self._analysis_complete)
    
    def _run_visualization(self):
        """Run visualization script"""
        viz_script = os.path.join(self.ml_scripts_dir, "ml_visualization.py")
        
        if not os.path.exists(viz_script):
            self._log("Visualization script not found, skipping...")
            return
        
        self._log("\nGenerating visualizations...")
        
        try:
            env = os.environ.copy()
            env["PATH"] = f"{self.conda_env}/bin:" + env.get("PATH", "")
            
            cmd = [
                "python", viz_script,
                "--results_dir", self.output_dir.get()
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=self.pipe_dir)
            
            if result.returncode == 0:
                self._log("✓ Visualizations generated")
            else:
                self._log(f"Visualization warning: {result.stderr[:200]}")
                
        except Exception as e:
            self._log(f"Visualization error: {e}")
    
    def _run_alphagenome_analysis(self):
        """Run AlphaGenome analysis on top features"""
        ag_script = "/home/dirk/MEDIPIPE_PROJECT/work/post_analysis/DMR_ML_Discrimination/scripts/alphagenome/alphagenome_feature_analysis.py"
        
        if not os.path.exists(ag_script):
            self._log("AlphaGenome script not found, skipping...")
            return
        
        if not os.environ.get("ALPHAGENOME_API_KEY"):
            self._log("⚠ ALPHAGENOME_API_KEY not set, skipping AlphaGenome analysis")
            return
        
        self._log("\nRunning AlphaGenome analysis on top features...")
        
        try:
            env = os.environ.copy()
            env["PATH"] = f"{self.conda_env}/bin:" + env.get("PATH", "")
            
            cmd = [
                "python", ag_script,
                "--results_dir", self.output_dir.get(),
                "--top_n", "50"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
            
            if result.returncode == 0:
                self._log("✓ AlphaGenome analysis completed")
            else:
                self._log(f"AlphaGenome warning: {result.stderr[:200]}")
                
        except subprocess.TimeoutExpired:
            self._log("AlphaGenome analysis timed out (10 min limit)")
        except Exception as e:
            self._log(f"AlphaGenome error: {e}")
    
    def _stop_analysis(self):
        """Stop running analysis"""
        self._running = False
        if self._process:
            self._process.terminate()
            self._log("\n⬛ Analysis stopped by user")
    
    def _analysis_complete(self):
        """Reset UI after analysis completes"""
        self.run_btn.configure(state=tk.NORMAL)
        self.stop_btn.configure(state=tk.DISABLED)
    
    def _view_results(self):
        """Open results viewer"""
        output_dir = self.output_dir.get()
        
        # Look for result files
        results_file = os.path.join(output_dir, "ml_results_summary.json")
        
        if os.path.exists(results_file):
            try:
                import json
                with open(results_file) as f:
                    results = json.load(f)
                
                # Show results dialog
                dialog = tk.Toplevel(self)
                dialog.title("ML Analysis Results")
                dialog.geometry("600x500")
                dialog.transient(self.winfo_toplevel())
                
                text = tk.Text(dialog, wrap=tk.WORD, font=("Monospace", 12))
                scroll = ttk.Scrollbar(dialog, orient=tk.VERTICAL, command=text.yview)
                text.configure(yscrollcommand=scroll.set)
                
                text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
                scroll.pack(side=tk.RIGHT, fill=tk.Y)
                
                # Format results
                text.insert(tk.END, "=" * 50 + "\n")
                text.insert(tk.END, "ML DISCRIMINATION RESULTS\n")
                text.insert(tk.END, "=" * 50 + "\n\n")
                
                for model, data in results.items():
                    text.insert(tk.END, f"Model: {model}\n")
                    text.insert(tk.END, f"  AUROC: {data.get('auroc', 'N/A'):.3f}\n")
                    text.insert(tk.END, f"  P-value: {data.get('p_value', 'N/A'):.4f}\n")
                    
                    robustness = data.get('robustness', 'Unknown')
                    if robustness == 'ROBUST':
                        text.insert(tk.END, f"  Status: ✓ ROBUST\n")
                    elif robustness == 'FRAGILE':
                        text.insert(tk.END, f"  Status: ⚠ FRAGILE (may be overfit)\n")
                    else:
                        text.insert(tk.END, f"  Status: {robustness}\n")
                    text.insert(tk.END, "\n")
                
                text.configure(state=tk.DISABLED)
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load results: {e}")
        else:
            messagebox.showinfo("No Results", 
                "No results file found.\nPlease run the analysis first.")
    
    def _open_output_folder(self):
        """Open output folder in file manager"""
        output_dir = self.output_dir.get()
        
        if os.path.exists(output_dir):
            subprocess.Popen(["xdg-open", output_dir])
        else:
            os.makedirs(output_dir, exist_ok=True)
            subprocess.Popen(["xdg-open", output_dir])
    
    def set_matrix_path(self, path: str):
        """Set the matrix path programmatically"""
        self.matrix_path.set(path)
    
    def set_groups(self, group1: str, group2: str):
        """Set the groups programmatically"""
        self.group1.set(group1)
        self.group2.set(group2)
