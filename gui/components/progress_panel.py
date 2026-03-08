"""
MEDIPIPE Progress Panel Component
Panel for monitoring workflow execution progress.
"""

import tkinter as tk
from tkinter import ttk
from typing import Optional
import threading
import queue


class ProgressPanel(ttk.LabelFrame):
    """Panel for displaying workflow progress and logs"""
    
    def __init__(self, parent, **kwargs):
        super().__init__(parent, text="Execution Progress", padding=10, **kwargs)
        
        self._log_queue = queue.Queue()
        self._create_widgets()
        self._start_log_polling()
    
    def _create_widgets(self):
        """Create panel widgets"""
        # Top status bar
        status_frame = ttk.Frame(self)
        status_frame.pack(fill=tk.X, pady=(5, 15))
        
        # Status label
        self.status_var = tk.StringVar(value="Ready")
        status_label = ttk.Label(status_frame, textvariable=self.status_var, 
                                 font=("Helvetica", 13, "bold"))
        status_label.pack(side=tk.LEFT)
        
        # Current task label
        self.task_var = tk.StringVar(value="")
        task_label = ttk.Label(status_frame, textvariable=self.task_var, 
                               foreground="gray", font=("Helvetica", 11))
        task_label.pack(side=tk.LEFT, padx=25)
        
        # Progress bar
        progress_frame = ttk.Frame(self)
        progress_frame.pack(fill=tk.X, pady=(0, 15))
        
        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(progress_frame, variable=self.progress_var,
                                            maximum=100, mode="determinate")
        self.progress_bar.pack(side=tk.LEFT, fill=tk.X, expand=True, ipady=3)
        
        self.progress_label = ttk.Label(progress_frame, text="0%", width=8, 
                                        font=("Helvetica", 11, "bold"))
        self.progress_label.pack(side=tk.LEFT, padx=10)
        
        # Job counter
        counter_frame = ttk.Frame(self)
        counter_frame.pack(fill=tk.X, pady=(0, 15))
        
        ttk.Label(counter_frame, text="Jobs:", font=("Helvetica", 11)).pack(side=tk.LEFT)
        self.jobs_var = tk.StringVar(value="0 / 0")
        ttk.Label(counter_frame, textvariable=self.jobs_var, 
                  font=("Helvetica", 11, "bold")).pack(side=tk.LEFT, padx=8)
        
        ttk.Label(counter_frame, text="Current Rule:", 
                  font=("Helvetica", 11)).pack(side=tk.LEFT, padx=(25, 0))
        self.rule_var = tk.StringVar(value="-")
        ttk.Label(counter_frame, textvariable=self.rule_var, 
                  font=("Helvetica", 11)).pack(side=tk.LEFT, padx=8)
        
        ttk.Label(counter_frame, text="Sample:", 
                  font=("Helvetica", 11)).pack(side=tk.LEFT, padx=(25, 0))
        self.sample_var = tk.StringVar(value="-")
        ttk.Label(counter_frame, textvariable=self.sample_var, 
                  font=("Helvetica", 11)).pack(side=tk.LEFT, padx=8)
        
        # Log area
        log_frame = ttk.LabelFrame(self, text="Execution Log", padding=8)
        log_frame.pack(fill=tk.BOTH, expand=True)
        
        self.log_text = tk.Text(log_frame, height=12, width=80, state=tk.DISABLED,
                               wrap=tk.WORD, font=("Consolas", 10))
        log_scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=log_scroll.set)
        
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        log_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Configure text tags for coloring
        self.log_text.tag_configure("error", foreground="red", font=("Consolas", 10, "bold"))
        self.log_text.tag_configure("warning", foreground="orange")
        self.log_text.tag_configure("success", foreground="green", font=("Consolas", 10, "bold"))
        self.log_text.tag_configure("info", foreground="blue")
        self.log_text.tag_configure("rule", foreground="purple", font=("Consolas", 10, "bold"))
        
        # Control buttons
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill=tk.X, pady=(12, 5))
        
        self.clear_btn = ttk.Button(btn_frame, text="Clear Log", command=self.clear_log, width=12)
        self.clear_btn.pack(side=tk.LEFT, padx=4)
        
        self.autoscroll_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(btn_frame, text="Auto-scroll", variable=self.autoscroll_var).pack(side=tk.LEFT, padx=20)
        
        ttk.Button(btn_frame, text="Export Log...", command=self._export_log, width=12).pack(side=tk.RIGHT, padx=4)
    
    def _start_log_polling(self):
        """Start polling for log messages"""
        self._poll_logs()
    
    def _poll_logs(self):
        """Poll log queue and update display"""
        try:
            while True:
                msg = self._log_queue.get_nowait()
                self._append_log(msg)
        except queue.Empty:
            pass
        
        # Schedule next poll
        self.after(100, self._poll_logs)
    
    def _append_log(self, message: str):
        """Append a message to the log"""
        self.log_text.configure(state=tk.NORMAL)
        
        # Determine tag based on content
        tag = None
        msg_lower = message.lower()
        if "error" in msg_lower or "failed" in msg_lower:
            tag = "error"
        elif "warning" in msg_lower:
            tag = "warning"
        elif "complete" in msg_lower or "success" in msg_lower or "done" in msg_lower:
            tag = "success"
        elif message.startswith("rule "):
            tag = "rule"
        
        if tag:
            self.log_text.insert(tk.END, message + "\n", tag)
        else:
            self.log_text.insert(tk.END, message + "\n")
        
        self.log_text.configure(state=tk.DISABLED)
        
        # Auto-scroll if enabled
        if self.autoscroll_var.get():
            self.log_text.see(tk.END)
    
    def add_log(self, message: str):
        """Add a log message (thread-safe)"""
        self._log_queue.put(message)
    
    def clear_log(self):
        """Clear the log display"""
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.delete("1.0", tk.END)
        self.log_text.configure(state=tk.DISABLED)
    
    def set_status(self, status: str):
        """Set the status text"""
        self.status_var.set(status)
    
    def set_task(self, task: str):
        """Set the current task text"""
        self.task_var.set(task)
    
    def set_progress(self, percent: float, completed: int = 0, total: int = 0):
        """Set progress bar value"""
        self.progress_var.set(percent)
        self.progress_label.configure(text=f"{percent:.0f}%")
        self.jobs_var.set(f"{completed} / {total}")
    
    def set_current_rule(self, rule: str):
        """Set the current rule name"""
        self.rule_var.set(rule)
    
    def set_current_sample(self, sample: str):
        """Set the current sample name"""
        self.sample_var.set(sample)
    
    def update_from_progress(self, progress):
        """Update from a WorkflowProgress object"""
        self.set_status(progress.status.value.title())
        self.set_task(progress.message)
        self.set_progress(progress.progress_percent, progress.completed_jobs, progress.total_jobs)
        self.set_current_rule(progress.current_rule or "-")
        self.set_current_sample(progress.current_sample or "-")
    
    def reset(self):
        """Reset all progress indicators"""
        self.set_status("Ready")
        self.set_task("")
        self.set_progress(0, 0, 0)
        self.set_current_rule("-")
        self.set_current_sample("-")
        self.clear_log()
    
    def set_indeterminate(self, indeterminate: bool = True):
        """Switch progress bar to indeterminate mode"""
        if indeterminate:
            self.progress_bar.configure(mode="indeterminate")
            self.progress_bar.start(10)
        else:
            self.progress_bar.stop()
            self.progress_bar.configure(mode="determinate")
    
    def _export_log(self):
        """Export log to a file"""
        from tkinter import filedialog
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("Log files", "*.log"), ("All files", "*.*")],
            title="Export Log"
        )
        
        if filepath:
            content = self.log_text.get("1.0", tk.END)
            with open(filepath, "w") as f:
                f.write(content)
            self.add_log(f"Log exported to: {filepath}")


class MiniProgressBar(ttk.Frame):
    """A compact progress bar widget for use in other panels"""
    
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)
        
        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(self, variable=self.progress_var,
                                            maximum=100, mode="determinate", length=200)
        self.progress_bar.pack(side=tk.LEFT)
        
        self.label = ttk.Label(self, text="0%", width=5)
        self.label.pack(side=tk.LEFT, padx=5)
    
    def set_progress(self, percent: float):
        """Set progress value"""
        self.progress_var.set(percent)
        self.label.configure(text=f"{percent:.0f}%")
    
    def reset(self):
        """Reset to 0%"""
        self.set_progress(0)
