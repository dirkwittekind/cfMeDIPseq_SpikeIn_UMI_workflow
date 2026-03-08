"""
MEDIPIPE Workflow Controller
Handles Snakemake workflow execution and monitoring.
"""

import os
import subprocess
import threading
import queue
import signal
import time
import shutil
from pathlib import Path
from typing import Optional, Callable, List, Dict, Any
from dataclasses import dataclass
from enum import Enum


class WorkflowStatus(Enum):
    """Workflow execution status"""
    IDLE = "idle"
    RUNNING = "running"
    PAUSED = "paused"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class WorkflowProgress:
    """Progress information for a running workflow"""
    status: WorkflowStatus = WorkflowStatus.IDLE
    current_rule: str = ""
    current_sample: str = ""
    completed_jobs: int = 0
    total_jobs: int = 0
    progress_percent: float = 0.0
    message: str = ""
    error: str = ""


class WorkflowController:
    """Controls Snakemake workflow execution"""
    
    def __init__(self, pipe_dir: str = "/home/dirk/medipipe_warp"):
        self.pipe_dir = pipe_dir
        self.snakefiles_dir = os.path.join(pipe_dir, "snakefiles")
        self.work_dir = os.path.join(pipe_dir, "work")
        
        # Conda environment for running snakemake
        self.conda_env = "/home/dirk/.conda_envs/medipipe"
        self.conda_prefix = os.path.expanduser("~/miniconda3")  # or ~/anaconda3
        
        self.process: Optional[subprocess.Popen] = None
        self.status = WorkflowStatus.IDLE
        self.progress = WorkflowProgress()
        
        self._log_queue: queue.Queue = queue.Queue()
        self._log_thread: Optional[threading.Thread] = None
        self._stop_event = threading.Event()
        
        # Callbacks
        self._on_progress: Optional[Callable[[WorkflowProgress], None]] = None
        self._on_log: Optional[Callable[[str], None]] = None
        self._on_complete: Optional[Callable[[bool], None]] = None
    
    def set_callbacks(
        self,
        on_progress: Optional[Callable[[WorkflowProgress], None]] = None,
        on_log: Optional[Callable[[str], None]] = None,
        on_complete: Optional[Callable[[bool], None]] = None
    ):
        """Set callback functions for progress updates"""
        self._on_progress = on_progress
        self._on_log = on_log
        self._on_complete = on_complete
    
    def _get_conda_run_prefix(self) -> List[str]:
        """Get the command prefix to run within conda environment"""
        # Use conda run to execute commands in the environment
        return [
            "conda", "run", "-p", self.conda_env, "--no-capture-output"
        ]
    
    def _get_shell_command(self, cmd: List[str]) -> str:
        """Build a shell command that activates conda and runs the command"""
        # Build command string for bash -c execution
        cmd_str = " ".join(cmd)
        return f"source {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.conda_env} && {cmd_str}"
    
    def run_dry_run(
        self,
        config_path: str,
        workflow: str = "core",
        targets: Optional[List[str]] = None
    ) -> bool:
        """
        Run Snakemake in dry-run mode.
        
        Args:
            config_path: Path to config.yaml
            workflow: "core" or "advanced"
            targets: Optional specific targets to build
        
        Returns:
            True if dry-run succeeds
        """
        snakefile = self._get_snakefile(workflow)
        
        snakemake_cmd = [
            "snakemake",
            "--snakefile", snakefile,
            "--configfile", config_path,
            "--dry-run",
            "--printshellcmds",
            "--reason",
        ]
        
        if targets:
            snakemake_cmd.extend(targets)
        
        # Build full shell command with conda activation
        shell_cmd = self._get_shell_command(snakemake_cmd)
        
        try:
            result = subprocess.run(
                ["bash", "-c", shell_cmd],
                cwd=self.work_dir,
                capture_output=True,
                text=True,
                timeout=120
            )
            
            if self._on_log:
                self._on_log(result.stdout)
                if result.stderr:
                    self._on_log(f"STDERR: {result.stderr}")
            
            return result.returncode == 0
        
        except subprocess.TimeoutExpired:
            if self._on_log:
                self._on_log("Dry run timed out")
            return False
        except Exception as e:
            if self._on_log:
                self._on_log(f"Dry run error: {e}")
            return False
    
    def _cleanup_failed_conda_envs(self):
        """Clean up any failed/incomplete conda environments before starting workflow.
        
        This removes empty directories that mamba sees as 'non-conda folders',
        which can happen when previous runs were interrupted.
        """
        # Check both possible conda locations
        conda_dirs = [
            os.path.join(self.work_dir, ".snakemake/conda"),
            os.path.join(self.pipe_dir, ".snakemake/conda")
        ]
        
        for conda_dir in conda_dirs:
            if not os.path.exists(conda_dir):
                continue
            
            for item in os.listdir(conda_dir):
                item_path = os.path.join(conda_dir, item)
                if os.path.isdir(item_path):
                    # Check if this is a valid conda environment (has conda-meta)
                    conda_meta = os.path.join(item_path, "conda-meta")
                    if not os.path.exists(conda_meta):
                        # Not a valid conda env - remove it
                        try:
                            shutil.rmtree(item_path)
                            if self._on_log:
                                self._on_log(f"Cleaned up incomplete conda env: {item}")
                        except Exception as e:
                            if self._on_log:
                                self._on_log(f"Warning: Could not remove {item_path}: {e}")
    
    def start_workflow(
        self,
        config_path: str,
        workflow: str = "core",
        targets: Optional[List[str]] = None,
        cores: int = 24,
        use_conda: bool = True
    ) -> bool:
        """
        Start the Snakemake workflow.
        
        Args:
            config_path: Path to config.yaml
            workflow: "core" or "advanced"
            targets: Optional specific targets
            cores: Number of cores to use
            use_conda: Whether to use conda environments
        
        Returns:
            True if workflow started successfully
        """
        if self.status == WorkflowStatus.RUNNING:
            return False
        
        # Clean up any failed conda environments from previous runs
        if use_conda:
            self._cleanup_failed_conda_envs()
        
        snakefile = self._get_snakefile(workflow)
        
        snakemake_cmd = [
            "snakemake",
            "--snakefile", snakefile,
            "--configfile", config_path,
            "--cores", str(cores),
            "--keep-going",
            "--rerun-incomplete",
            "--printshellcmds",
        ]
        
        if use_conda:
            # Use conda instead of mamba due to mamba bug with "Non-conda folder exists at prefix"
            snakemake_cmd.extend(["--use-conda", "--conda-frontend", "conda"])
        
        if targets:
            snakemake_cmd.extend(targets)
        
        # Build full shell command with conda activation
        shell_cmd = self._get_shell_command(snakemake_cmd)
        
        try:
            # Ensure work directory exists
            os.makedirs(self.work_dir, exist_ok=True)
            
            # Start the process with bash -c to activate conda
            self.process = subprocess.Popen(
                ["bash", "-c", shell_cmd],
                cwd=self.work_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                preexec_fn=os.setsid  # Create new process group
            )
            
            self.status = WorkflowStatus.RUNNING
            self.progress.status = WorkflowStatus.RUNNING
            self.progress.message = "Workflow started"
            self._stop_event.clear()
            
            # Start log monitoring thread
            self._log_thread = threading.Thread(target=self._monitor_output, daemon=True)
            self._log_thread.start()
            
            if self._on_progress:
                self._on_progress(self.progress)
            
            return True
        
        except Exception as e:
            self.status = WorkflowStatus.FAILED
            self.progress.status = WorkflowStatus.FAILED
            self.progress.error = str(e)
            
            if self._on_log:
                self._on_log(f"Failed to start workflow: {e}")
            
            return False
    
    def stop_workflow(self) -> bool:
        """Stop the running workflow"""
        if self.process is None or self.status != WorkflowStatus.RUNNING:
            return False
        
        try:
            # Send SIGTERM to the process group
            os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            
            # Wait a bit for graceful shutdown
            time.sleep(2)
            
            # Force kill if still running
            if self.process.poll() is None:
                os.killpg(os.getpgid(self.process.pid), signal.SIGKILL)
            
            self._stop_event.set()
            self.status = WorkflowStatus.CANCELLED
            self.progress.status = WorkflowStatus.CANCELLED
            self.progress.message = "Workflow cancelled by user"
            
            if self._on_progress:
                self._on_progress(self.progress)
            
            if self._on_complete:
                self._on_complete(False)
            
            return True
        
        except Exception as e:
            if self._on_log:
                self._on_log(f"Error stopping workflow: {e}")
            return False
    
    def _get_snakefile(self, workflow: str) -> str:
        """Get the path to the appropriate Snakefile"""
        if workflow == "core":
            return os.path.join(self.snakefiles_dir, "core_workflow.smk")
        elif workflow == "advanced":
            return os.path.join(self.snakefiles_dir, "advanced_workflow.smk")
        else:
            return os.path.join(self.snakefiles_dir, "core_workflow.smk")
    
    def _monitor_output(self):
        """Monitor Snakemake output in a separate thread"""
        if self.process is None:
            return
        
        try:
            for line in iter(self.process.stdout.readline, ''):
                if self._stop_event.is_set():
                    break
                
                line = line.strip()
                if not line:
                    continue
                
                # Parse Snakemake output for progress
                self._parse_snakemake_output(line)
                
                # Send to log callback
                if self._on_log:
                    self._on_log(line)
            
            # Process completed
            return_code = self.process.wait()
            
            if return_code == 0:
                self.status = WorkflowStatus.COMPLETED
                self.progress.status = WorkflowStatus.COMPLETED
                self.progress.message = "Workflow completed successfully"
                self.progress.progress_percent = 100.0
            else:
                self.status = WorkflowStatus.FAILED
                self.progress.status = WorkflowStatus.FAILED
                self.progress.error = f"Workflow failed with return code {return_code}"
            
            if self._on_progress:
                self._on_progress(self.progress)
            
            if self._on_complete:
                self._on_complete(return_code == 0)
        
        except Exception as e:
            self.status = WorkflowStatus.FAILED
            self.progress.status = WorkflowStatus.FAILED
            self.progress.error = str(e)
            
            if self._on_log:
                self._on_log(f"Monitoring error: {e}")
    
    def _parse_snakemake_output(self, line: str):
        """Parse Snakemake output to extract progress information"""
        # Look for rule execution
        if line.startswith("rule ") or "rule " in line:
            parts = line.split()
            if len(parts) >= 2:
                self.progress.current_rule = parts[1].rstrip(":")
        
        # Look for job counts
        if "of " in line and "steps" in line:
            # Example: "1 of 10 steps (10%) done"
            import re
            match = re.search(r"(\d+) of (\d+) steps \((\d+)%\)", line)
            if match:
                self.progress.completed_jobs = int(match.group(1))
                self.progress.total_jobs = int(match.group(2))
                self.progress.progress_percent = float(match.group(3))
        
        # Look for sample being processed
        if "wildcards:" in line:
            parts = line.split("wildcards:")
            if len(parts) > 1:
                wildcards_str = parts[1].strip()
                if "sample=" in wildcards_str:
                    sample = wildcards_str.split("sample=")[1].split(",")[0].split()[0]
                    self.progress.current_sample = sample
        
        self.progress.message = f"Running: {self.progress.current_rule}"
        
        if self._on_progress:
            self._on_progress(self.progress)
    
    def get_status(self) -> WorkflowStatus:
        """Get current workflow status"""
        return self.status
    
    def get_progress(self) -> WorkflowProgress:
        """Get current progress information"""
        return self.progress
    
    def is_running(self) -> bool:
        """Check if workflow is currently running"""
        return self.status == WorkflowStatus.RUNNING
