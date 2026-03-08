"""
MEDIPIPE GUI - QC Panel Component
Displays QC status and reports for samples.
"""

import os
import json
import tkinter as tk
from tkinter import ttk, messagebox
from pathlib import Path
from typing import Dict, List, Optional, Callable


class QCStatus:
    """Container for sample QC status"""
    def __init__(self, sample_id: str):
        self.sample_id = sample_id
        self.has_qc_report = False
        self.qc_pass = None
        self.dedup_reads = None
        self.relH = None
        self.GoGe = None
        self.insert_size = None
        self.avg_quality = None
        self.spikein_status = None
        self.spikein_relH = None
        self.spikein_meth_ratio = None
        # Tissue of origin
        self.too_brain_total = None
        self.too_blood_total = None
        self.too_liver = None
        self.too_status = None
        self.too_top_tissues = None
        self.warnings = []
        self.report_path = None


class QCReportViewer:
    """Dialog to view full QC report for a sample"""
    
    def __init__(self, parent, sample_id: str, report_path: str):
        self.window = tk.Toplevel(parent)
        self.window.title(f"QC Report - {sample_id}")
        self.window.geometry("800x600")
        self.window.transient(parent)
        
        # Text widget with scrollbar
        frame = ttk.Frame(self.window, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(frame, text=f"QC Report: {sample_id}", 
                  font=("Helvetica", 14, "bold")).pack(anchor="w", pady=(0, 10))
        
        text_frame = ttk.Frame(frame)
        text_frame.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(text_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.text = tk.Text(text_frame, wrap=tk.WORD, font=("Courier", 11),
                           yscrollcommand=scrollbar.set)
        self.text.pack(fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.text.yview)
        
        # Load report
        try:
            with open(report_path) as f:
                content = f.read()
            self.text.insert("1.0", content)
        except Exception as e:
            self.text.insert("1.0", f"Error loading report: {e}")
        
        self.text.config(state=tk.DISABLED)
        
        # Close button
        ttk.Button(frame, text="Close", command=self.window.destroy,
                   width=15).pack(pady=10)


class QCPanel(ttk.Frame):
    """Panel displaying QC status for all samples"""
    
    def __init__(self, parent, output_dir: str = "/data/medipipe_data/output",
                 on_refresh: Optional[Callable] = None):
        super().__init__(parent)
        self.output_dir = Path(output_dir)
        self.on_refresh = on_refresh
        self.qc_statuses: Dict[str, QCStatus] = {}
        
        self._create_widgets()
    
    def _create_widgets(self):
        """Create the QC panel widgets"""
        # Header with refresh button
        header = ttk.Frame(self)
        header.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(header, text="Sample QC Status", 
                  font=("Helvetica", 14, "bold")).pack(side=tk.LEFT)
        
        ttk.Button(header, text="🔄 Refresh", command=self.refresh_qc_status,
                   width=12).pack(side=tk.RIGHT)
        ttk.Button(header, text="📊 View Summary", command=self._view_summary,
                   width=14).pack(side=tk.RIGHT, padx=5)
        
        # QC Treeview
        columns = ("sample", "status", "reads", "relH", "GoGe", "quality", "spikein", "too")
        self.tree = ttk.Treeview(self, columns=columns, show="headings", 
                                 selectmode="browse")
        
        # Configure columns
        self.tree.heading("sample", text="Sample")
        self.tree.heading("status", text="QC Status")
        self.tree.heading("reads", text="Dedup Reads")
        self.tree.heading("relH", text="relH")
        self.tree.heading("GoGe", text="GoGe")
        self.tree.heading("quality", text="Avg Quality")
        self.tree.heading("spikein", text="Spike-in")
        self.tree.heading("too", text="Tissue Origin")
        
        self.tree.column("sample", width=120, minwidth=80)
        self.tree.column("status", width=70, minwidth=50)
        self.tree.column("reads", width=90, minwidth=70)
        self.tree.column("relH", width=50, minwidth=40)
        self.tree.column("GoGe", width=50, minwidth=40)
        self.tree.column("quality", width=60, minwidth=50)
        self.tree.column("spikein", width=80, minwidth=60)
        self.tree.column("too", width=120, minwidth=100)
        
        # Scrollbar
        scrollbar = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Tags for coloring
        self.tree.tag_configure("pass", foreground="#198754")  # Green
        self.tree.tag_configure("fail", foreground="#dc3545")  # Red
        self.tree.tag_configure("warning", foreground="#fd7e14")  # Orange
        self.tree.tag_configure("pending", foreground="#6c757d")  # Gray
        
        # Bind double-click to view report
        self.tree.bind("<Double-1>", self._on_double_click)
    
    def set_output_dir(self, output_dir: str):
        """Set the output directory and refresh"""
        self.output_dir = Path(output_dir)
        self.refresh_qc_status()
    
    def refresh_qc_status(self, samples: Optional[List[str]] = None):
        """Refresh QC status for all samples"""
        self.tree.delete(*self.tree.get_children())
        self.qc_statuses.clear()
        
        # Find samples from QC reports or dedup stats
        if samples:
            sample_list = samples
        else:
            sample_list = self._detect_samples()
        
        # Load QC data
        self._load_summary_qc()
        
        # Load per-sample data for any missing
        for sample_id in sample_list:
            if sample_id not in self.qc_statuses:
                self._load_sample_qc(sample_id)
        
        # Populate tree
        for sample_id in sorted(sample_list):
            self._add_sample_row(sample_id)
        
        if self.on_refresh:
            self.on_refresh(self.qc_statuses)
    
    def _detect_samples(self) -> List[str]:
        """Detect samples from output directory"""
        samples = set()
        
        # Check dedup_bam_pe folder
        dedup_dir = self.output_dir / "dedup_bam_pe"
        if dedup_dir.exists():
            for f in dedup_dir.glob("*_dedup.bam.stats.txt"):
                sample = f.stem.replace("_dedup.bam.stats", "")
                samples.add(sample)
        
        # Also check QC reports folder
        qc_dir = self.output_dir / "qc_reports" / "per_sample"
        if qc_dir.exists():
            for f in qc_dir.glob("*_qc.json"):
                sample = f.stem.replace("_qc", "")
                samples.add(sample)
        
        return sorted(samples)
    
    def _load_summary_qc(self):
        """Load data from summary QC report"""
        json_path = self.output_dir / "qc_reports" / "qc_report.json"
        
        # Also check comprehensive report
        if not json_path.exists():
            json_path = self.output_dir / "qc_report_comprehensive.json"
        
        if not json_path.exists():
            return
        
        try:
            with open(json_path) as f:
                data = json.load(f)
            
            for sample_data in data.get("samples", []):
                sample_id = sample_data.get("sample_id")
                if not sample_id:
                    continue
                
                qc = QCStatus(sample_id)
                qc.has_qc_report = True
                qc.qc_pass = sample_data.get("qc_pass", True)
                qc.dedup_reads = sample_data.get("dedup_reads")
                qc.relH = sample_data.get("relH")
                qc.GoGe = sample_data.get("GoGe")
                qc.insert_size = sample_data.get("insert_size_avg")
                qc.avg_quality = sample_data.get("avg_quality")
                qc.spikein_status = sample_data.get("spikein_status")
                qc.spikein_relH = sample_data.get("spikein_relH")
                qc.spikein_meth_ratio = sample_data.get("spikein_meth_ratio")
                # Tissue of origin
                qc.too_brain_total = sample_data.get("too_brain_total")
                qc.too_blood_total = sample_data.get("too_blood_total")
                qc.too_liver = sample_data.get("too_liver")
                qc.too_status = sample_data.get("too_status")
                qc.too_top_tissues = sample_data.get("too_top_tissues")
                
                warnings = sample_data.get("qc_warnings", "")
                if warnings:
                    qc.warnings = warnings.split("; ")
                
                # Find report path
                per_sample = self.output_dir / "qc_reports" / "per_sample" / f"{sample_id}_qc.txt"
                if per_sample.exists():
                    qc.report_path = str(per_sample)
                
                self.qc_statuses[sample_id] = qc
                
        except Exception as e:
            print(f"Error loading QC summary: {e}")
    
    def _load_sample_qc(self, sample_id: str):
        """Load QC data for a single sample"""
        qc = QCStatus(sample_id)
        
        # Check for per-sample JSON
        json_path = self.output_dir / "qc_reports" / "per_sample" / f"{sample_id}_qc.json"
        if json_path.exists():
            try:
                with open(json_path) as f:
                    data = json.load(f)
                
                if data.get("samples"):
                    sample_data = data["samples"][0]
                    qc.has_qc_report = True
                    qc.qc_pass = sample_data.get("qc_pass", True)
                    qc.dedup_reads = sample_data.get("dedup_reads")
                    qc.relH = sample_data.get("relH")
                    qc.GoGe = sample_data.get("GoGe")
                    qc.avg_quality = sample_data.get("avg_quality")
                    qc.spikein_status = sample_data.get("spikein_status")
                    qc.report_path = str(json_path.with_suffix(".txt"))
            except Exception as e:
                print(f"Error loading QC for {sample_id}: {e}")
        
        # Check if stats file exists but no QC report yet
        stats_path = self.output_dir / "dedup_bam_pe" / f"{sample_id}_dedup.bam.stats.txt"
        if stats_path.exists() and not qc.has_qc_report:
            qc.has_qc_report = False  # Stats exist but no QC report generated
        
        self.qc_statuses[sample_id] = qc
    
    def _add_sample_row(self, sample_id: str):
        """Add a sample row to the tree"""
        qc = self.qc_statuses.get(sample_id, QCStatus(sample_id))
        
        # Determine status display
        if not qc.has_qc_report:
            status = "Pending"
            tag = "pending"
        elif qc.qc_pass:
            status = "✅ PASS"
            tag = "pass"
        else:
            status = "❌ FAIL"
            tag = "fail"
        
        # Format values
        reads = f"{qc.dedup_reads:,}" if qc.dedup_reads else "—"
        relH = f"{qc.relH:.2f}" if qc.relH else "—"
        GoGe = f"{qc.GoGe:.2f}" if qc.GoGe else "—"
        quality = f"{qc.avg_quality:.1f}" if qc.avg_quality else "—"
        
        # Spike-in status
        if qc.spikein_status:
            if qc.spikein_status == "PASS":
                spikein = "✅ PASS"
            elif qc.spikein_status.startswith("WARNING"):
                spikein = "⚠️ WARN"
                if tag == "pass":
                    tag = "warning"
            else:
                spikein = qc.spikein_status[:15]
        elif qc.spikein_relH:
            spikein = f"relH: {qc.spikein_relH:.1f}"
        else:
            spikein = "—"
        
        # Tissue of origin
        if qc.too_brain_total is not None:
            brain_pct = f"{qc.too_brain_total*100:.0f}%"
            blood_pct = f"{qc.too_blood_total*100:.0f}%" if qc.too_blood_total else "0%"
            too = f"🧠{brain_pct} 🩸{blood_pct}"
            if qc.too_status and qc.too_status != "NORMAL":
                too = f"⚠️ {too}"
        else:
            too = "—"
        
        self.tree.insert("", "end", sample_id, 
                        values=(sample_id, status, reads, relH, GoGe, quality, spikein, too),
                        tags=(tag,))
    
    def _on_double_click(self, event):
        """Handle double-click to view full report"""
        item = self.tree.selection()
        if not item:
            return
        
        sample_id = item[0]
        qc = self.qc_statuses.get(sample_id)
        
        if qc and qc.report_path and os.path.exists(qc.report_path):
            QCReportViewer(self.winfo_toplevel(), sample_id, qc.report_path)
        else:
            # Try to find the report
            paths = [
                self.output_dir / "qc_reports" / "per_sample" / f"{sample_id}_qc.txt",
                self.output_dir / "qc_reports" / "qc_report.txt",
                self.output_dir / "qc_report_comprehensive.txt"
            ]
            
            for path in paths:
                if path.exists():
                    QCReportViewer(self.winfo_toplevel(), sample_id, str(path))
                    return
            
            messagebox.showinfo("No Report", 
                f"No QC report found for sample {sample_id}.\n\n"
                "Run the workflow to generate QC reports.")
    
    def _view_summary(self):
        """View the summary QC report"""
        paths = [
            self.output_dir / "qc_reports" / "qc_report.txt",
            self.output_dir / "qc_report_comprehensive.txt"
        ]
        
        for path in paths:
            if path.exists():
                QCReportViewer(self.winfo_toplevel(), "Summary", str(path))
                return
        
        messagebox.showinfo("No Summary Report", 
            "No summary QC report found.\n\n"
            "Run the workflow with all samples to generate the summary report.")
    
    def get_qc_status(self, sample_id: str) -> Optional[QCStatus]:
        """Get QC status for a specific sample"""
        return self.qc_statuses.get(sample_id)
    
    def get_all_passing_samples(self) -> List[str]:
        """Get list of all samples that pass QC"""
        return [sid for sid, qc in self.qc_statuses.items() 
                if qc.has_qc_report and qc.qc_pass]
    
    def get_all_failing_samples(self) -> List[str]:
        """Get list of all samples that fail QC"""
        return [sid for sid, qc in self.qc_statuses.items() 
                if qc.has_qc_report and not qc.qc_pass]
