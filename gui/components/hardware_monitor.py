"""
MEDIPIPE Hardware Monitor Component
Real-time monitoring of CPU, RAM, and GPU usage.
"""

import tkinter as tk
from tkinter import ttk
from typing import Optional, Callable
import threading
import time

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# Try to import GPU monitoring
try:
    import subprocess
    def get_nvidia_gpu_info():
        """Get NVIDIA GPU info using nvidia-smi"""
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=utilization.gpu,memory.used,memory.total,name', 
                 '--format=csv,noheader,nounits'],
                capture_output=True, text=True, timeout=2
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                gpus = []
                for line in lines:
                    parts = [p.strip() for p in line.split(',')]
                    if len(parts) >= 4:
                        gpus.append({
                            'utilization': float(parts[0]),
                            'memory_used': float(parts[1]),
                            'memory_total': float(parts[2]),
                            'name': parts[3]
                        })
                return gpus
        except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
            pass
        return []
    
    GPU_AVAILABLE = len(get_nvidia_gpu_info()) > 0
except:
    GPU_AVAILABLE = False
    def get_nvidia_gpu_info():
        return []


class HardwareUsageBar(ttk.Frame):
    """A single hardware usage bar with label and percentage"""
    
    def __init__(self, parent, label: str, color: str = "#0d6efd", show_details: bool = True, **kwargs):
        super().__init__(parent, **kwargs)
        
        self.label_text = label
        self.color = color
        self.show_details = show_details
        
        self._create_widgets()
    
    def _create_widgets(self):
        """Create the usage bar widgets"""
        # Label
        self.label = ttk.Label(self, text=self.label_text, width=8, font=("Helvetica", 10))
        self.label.pack(side=tk.LEFT, padx=(0, 5))
        
        # Progress bar container (for custom coloring)
        self.bar_frame = ttk.Frame(self)
        self.bar_frame.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        # Progress bar
        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(
            self.bar_frame, 
            variable=self.progress_var,
            maximum=100, 
            mode="determinate",
            length=150
        )
        self.progress_bar.pack(fill=tk.X, expand=True, ipady=2)
        
        # Percentage label
        self.percent_label = ttk.Label(self, text="0%", width=5, font=("Helvetica", 10, "bold"))
        self.percent_label.pack(side=tk.LEFT, padx=(5, 0))
        
        # Details label (e.g., "8.5 / 16.0 GB")
        if self.show_details:
            self.detail_label = ttk.Label(self, text="", width=14, font=("Helvetica", 9), foreground="gray")
            self.detail_label.pack(side=tk.LEFT, padx=(5, 0))
    
    def set_value(self, percent: float, detail_text: str = ""):
        """Set the usage value"""
        self.progress_var.set(min(100, max(0, percent)))
        self.percent_label.configure(text=f"{percent:.0f}%")
        
        # Color based on usage level
        if percent >= 90:
            self.percent_label.configure(foreground="red")
        elif percent >= 70:
            self.percent_label.configure(foreground="orange")
        else:
            self.percent_label.configure(foreground="green")
        
        if self.show_details and hasattr(self, 'detail_label'):
            self.detail_label.configure(text=detail_text)


class HardwareMonitor(ttk.LabelFrame):
    """Panel for monitoring system hardware usage"""
    
    def __init__(self, parent, on_cores_change: Optional[Callable[[int], None]] = None, **kwargs):
        super().__init__(parent, text="System Resources", padding=10, **kwargs)
        
        self.on_cores_change = on_cores_change
        self._monitoring = False
        self._monitor_thread = None
        
        self._create_widgets()
        self.start_monitoring()
    
    def _create_widgets(self):
        """Create monitor widgets"""
        # Top row: Cores selection
        cores_frame = ttk.Frame(self)
        cores_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(cores_frame, text="CPU Cores to use:", font=("Helvetica", 10)).pack(side=tk.LEFT)
        
        # Get available cores
        if PSUTIL_AVAILABLE:
            max_cores = psutil.cpu_count(logical=True) or 8
            physical_cores = psutil.cpu_count(logical=False) or 4
        else:
            max_cores = 8
            physical_cores = 4
        
        self.cores_var = tk.IntVar(value=min(24, max_cores))
        
        # Cores spinbox
        self.cores_spinbox = ttk.Spinbox(
            cores_frame, 
            from_=1, 
            to=max_cores, 
            textvariable=self.cores_var,
            width=5,
            command=self._on_cores_changed
        )
        self.cores_spinbox.pack(side=tk.LEFT, padx=10)
        
        # Quick select buttons
        ttk.Button(cores_frame, text=f"Half ({physical_cores})", width=10,
                   command=lambda: self._set_cores(physical_cores)).pack(side=tk.LEFT, padx=2)
        ttk.Button(cores_frame, text=f"All ({max_cores})", width=10,
                   command=lambda: self._set_cores(max_cores)).pack(side=tk.LEFT, padx=2)
        
        # Info label
        self.cores_info = ttk.Label(cores_frame, 
                                     text=f"({physical_cores} physical, {max_cores} logical)",
                                     font=("Helvetica", 9), foreground="gray")
        self.cores_info.pack(side=tk.LEFT, padx=10)
        
        # Separator
        ttk.Separator(self, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=5)
        
        # Usage bars frame
        usage_frame = ttk.Frame(self)
        usage_frame.pack(fill=tk.X, pady=5)
        
        # CPU Usage
        cpu_frame = ttk.Frame(usage_frame)
        cpu_frame.pack(fill=tk.X, pady=3)
        self.cpu_bar = HardwareUsageBar(cpu_frame, "CPU:", color="#0d6efd", show_details=True)
        self.cpu_bar.pack(fill=tk.X)
        
        # RAM Usage
        ram_frame = ttk.Frame(usage_frame)
        ram_frame.pack(fill=tk.X, pady=3)
        self.ram_bar = HardwareUsageBar(ram_frame, "RAM:", color="#198754", show_details=True)
        self.ram_bar.pack(fill=tk.X)
        
        # GPU Usage (if available)
        if GPU_AVAILABLE:
            gpu_frame = ttk.Frame(usage_frame)
            gpu_frame.pack(fill=tk.X, pady=3)
            self.gpu_bar = HardwareUsageBar(gpu_frame, "GPU:", color="#6f42c1", show_details=True)
            self.gpu_bar.pack(fill=tk.X)
            
            # GPU Memory
            gpu_mem_frame = ttk.Frame(usage_frame)
            gpu_mem_frame.pack(fill=tk.X, pady=3)
            self.gpu_mem_bar = HardwareUsageBar(gpu_mem_frame, "VRAM:", color="#d63384", show_details=True)
            self.gpu_mem_bar.pack(fill=tk.X)
        else:
            self.gpu_bar = None
            self.gpu_mem_bar = None
            
            # Show GPU not available message
            ttk.Label(usage_frame, text="GPU: Not detected or nvidia-smi not available",
                     font=("Helvetica", 9), foreground="gray").pack(fill=tk.X, pady=3)
        
        # Status line
        self.status_frame = ttk.Frame(self)
        self.status_frame.pack(fill=tk.X, pady=(5, 0))
        
        self.status_var = tk.StringVar(value="Monitoring system resources...")
        ttk.Label(self.status_frame, textvariable=self.status_var, 
                 font=("Helvetica", 9), foreground="gray").pack(side=tk.LEFT)
    
    def _set_cores(self, cores: int):
        """Set the number of cores"""
        self.cores_var.set(cores)
        self._on_cores_changed()
    
    def _on_cores_changed(self):
        """Handle cores value change"""
        if self.on_cores_change:
            self.on_cores_change(self.cores_var.get())
    
    def get_cores(self) -> int:
        """Get the selected number of cores"""
        return self.cores_var.get()
    
    def set_cores(self, cores: int):
        """Set the number of cores from external code"""
        self.cores_var.set(cores)
    
    def start_monitoring(self):
        """Start the monitoring thread"""
        if self._monitoring:
            return
        
        self._monitoring = True
        self._monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self._monitor_thread.start()
    
    def stop_monitoring(self):
        """Stop the monitoring thread"""
        self._monitoring = False
    
    def _monitor_loop(self):
        """Background monitoring loop"""
        while self._monitoring:
            try:
                self._update_stats()
            except Exception as e:
                pass
            time.sleep(1)  # Update every second
    
    def _update_stats(self):
        """Update hardware statistics (called from background thread)"""
        if not PSUTIL_AVAILABLE:
            self.after(0, lambda: self.status_var.set("psutil not installed - monitoring unavailable"))
            return
        
        # CPU
        cpu_percent = psutil.cpu_percent(interval=None)
        cpu_count = psutil.cpu_count()
        cpu_freq = psutil.cpu_freq()
        freq_str = f"{cpu_freq.current/1000:.1f} GHz" if cpu_freq else ""
        
        # RAM
        mem = psutil.virtual_memory()
        ram_percent = mem.percent
        ram_used = mem.used / (1024**3)
        ram_total = mem.total / (1024**3)
        ram_detail = f"{ram_used:.1f} / {ram_total:.1f} GB"
        
        # Schedule UI update on main thread
        self.after(0, lambda: self._update_ui(cpu_percent, freq_str, ram_percent, ram_detail))
        
        # GPU (if available)
        if self.gpu_bar:
            gpus = get_nvidia_gpu_info()
            if gpus:
                gpu = gpus[0]  # Use first GPU
                gpu_util = gpu['utilization']
                gpu_mem_used = gpu['memory_used'] / 1024  # MB to GB
                gpu_mem_total = gpu['memory_total'] / 1024
                gpu_mem_percent = (gpu['memory_used'] / gpu['memory_total']) * 100 if gpu['memory_total'] > 0 else 0
                gpu_detail = gpu['name'][:20] if len(gpu['name']) > 20 else gpu['name']
                gpu_mem_detail = f"{gpu_mem_used:.1f} / {gpu_mem_total:.1f} GB"
                
                self.after(0, lambda u=gpu_util, d=gpu_detail, mp=gpu_mem_percent, md=gpu_mem_detail: 
                          self._update_gpu_ui(u, d, mp, md))
    
    def _update_ui(self, cpu_percent: float, cpu_detail: str, ram_percent: float, ram_detail: str):
        """Update UI elements (must be called on main thread)"""
        self.cpu_bar.set_value(cpu_percent, cpu_detail)
        self.ram_bar.set_value(ram_percent, ram_detail)
    
    def _update_gpu_ui(self, gpu_util: float, gpu_detail: str, gpu_mem_percent: float, gpu_mem_detail: str):
        """Update GPU UI elements (must be called on main thread)"""
        if self.gpu_bar:
            self.gpu_bar.set_value(gpu_util, gpu_detail)
        if self.gpu_mem_bar:
            self.gpu_mem_bar.set_value(gpu_mem_percent, gpu_mem_detail)
    
    def destroy(self):
        """Clean up when widget is destroyed"""
        self.stop_monitoring()
        super().destroy()


class CompactHardwareMonitor(ttk.Frame):
    """A compact version of the hardware monitor for embedding in status bars"""
    
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)
        
        self._monitoring = False
        self._create_widgets()
        self.start_monitoring()
    
    def _create_widgets(self):
        """Create compact monitor widgets"""
        # CPU
        ttk.Label(self, text="CPU:", font=("Helvetica", 9)).pack(side=tk.LEFT, padx=(0, 2))
        self.cpu_var = tk.StringVar(value="0%")
        self.cpu_label = ttk.Label(self, textvariable=self.cpu_var, width=4, 
                                   font=("Helvetica", 9, "bold"))
        self.cpu_label.pack(side=tk.LEFT, padx=(0, 10))
        
        # RAM
        ttk.Label(self, text="RAM:", font=("Helvetica", 9)).pack(side=tk.LEFT, padx=(0, 2))
        self.ram_var = tk.StringVar(value="0%")
        self.ram_label = ttk.Label(self, textvariable=self.ram_var, width=4,
                                   font=("Helvetica", 9, "bold"))
        self.ram_label.pack(side=tk.LEFT, padx=(0, 10))
        
        # GPU (if available)
        if GPU_AVAILABLE:
            ttk.Label(self, text="GPU:", font=("Helvetica", 9)).pack(side=tk.LEFT, padx=(0, 2))
            self.gpu_var = tk.StringVar(value="0%")
            self.gpu_label = ttk.Label(self, textvariable=self.gpu_var, width=4,
                                       font=("Helvetica", 9, "bold"))
            self.gpu_label.pack(side=tk.LEFT)
        else:
            self.gpu_var = None
            self.gpu_label = None
    
    def start_monitoring(self):
        """Start monitoring"""
        if self._monitoring:
            return
        self._monitoring = True
        self._update()
    
    def stop_monitoring(self):
        """Stop monitoring"""
        self._monitoring = False
    
    def _update(self):
        """Update stats"""
        if not self._monitoring:
            return
        
        if PSUTIL_AVAILABLE:
            cpu = psutil.cpu_percent(interval=None)
            ram = psutil.virtual_memory().percent
            
            self.cpu_var.set(f"{cpu:.0f}%")
            self.ram_var.set(f"{ram:.0f}%")
            
            # Color coding
            self.cpu_label.configure(foreground="red" if cpu >= 90 else "orange" if cpu >= 70 else "green")
            self.ram_label.configure(foreground="red" if ram >= 90 else "orange" if ram >= 70 else "green")
            
            if GPU_AVAILABLE and self.gpu_var:
                gpus = get_nvidia_gpu_info()
                if gpus:
                    gpu_util = gpus[0]['utilization']
                    self.gpu_var.set(f"{gpu_util:.0f}%")
                    self.gpu_label.configure(foreground="red" if gpu_util >= 90 else "orange" if gpu_util >= 70 else "green")
        
        # Schedule next update
        self.after(1000, self._update)
    
    def destroy(self):
        """Clean up"""
        self.stop_monitoring()
        super().destroy()
