"""
MEDIPIPE GUI - Main Application
Full Tkinter desktop GUI for cfMeDIP-seq analysis.
"""

import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from tkinter import font as tkfont
from typing import Optional, List

# Fix GTK file dialog font size on Linux
os.environ['GDK_SCALE'] = '1'
os.environ['GDK_DPI_SCALE'] = '1.2'  # Make dialogs slightly larger

# Import components
from gui.components.file_browser import FileBrowser, FileSelectionDialog
from gui.components.workflow_panel import WorkflowPanel
from gui.components.config_panel import ConfigPanel
from gui.components.progress_panel import ProgressPanel
from gui.components.hardware_monitor import HardwareMonitor, CompactHardwareMonitor
from gui.components.ml_panel import MLAnalysisPanel
from gui.components.qc_panel import QCPanel

# Import controllers
from gui.controllers.config_manager import ConfigManager, MedipipeConfig
from gui.controllers.sample_manager import SampleManager, Sample
from gui.controllers.workflow_controller import WorkflowController, WorkflowStatus


class MedipipeApp:
    """Main MEDIPIPE GUI Application"""
    
    VERSION = "1.0.0"
    PIPE_DIR = "/home/dirk/medipipe_warp"
    
    def __init__(self):
        self.root = tk.Tk()
        self.root.title(f"MEDIPIPE - cfMeDIP-seq Analysis Pipeline v{self.VERSION}")
        self.root.geometry("1200x800")
        self.root.minsize(1000, 700)
        
        # Controllers
        self.config_manager = ConfigManager(self.PIPE_DIR)
        self.sample_manager = SampleManager()
        self.workflow_controller = WorkflowController(self.PIPE_DIR)
        
        # Set up workflow callbacks
        self.workflow_controller.set_callbacks(
            on_progress=self._on_workflow_progress,
            on_log=self._on_workflow_log,
            on_complete=self._on_workflow_complete
        )
        
        # Current configuration
        self.current_config: Optional[MedipipeConfig] = None
        
        # Apply styles first (before creating widgets)
        self._apply_styles()
        
        # Create UI
        self._create_menu()
        self._create_header()
        self._create_toolbar()
        self._create_main_layout()
        self._create_statusbar()
        
        # Load default config
        self._load_default_config()
    
    def _apply_styles(self):
        """Apply custom styles with modern look"""
        style = ttk.Style()
        style.theme_use("clam")
        
        # Color scheme - modern and professional
        self.colors = {
            'bg': '#f8f9fa',           # Light gray background
            'fg': '#212529',           # Dark text
            'accent': '#0d6efd',       # Blue accent
            'accent_dark': '#0a58ca',  # Darker blue for hover
            'success': '#198754',      # Green for success/run
            'danger': '#dc3545',       # Red for stop/error
            'warning': '#ffc107',      # Yellow for warnings
            'secondary': '#6c757d',    # Gray for secondary text
            'border': '#dee2e6',       # Light border
            'card_bg': '#ffffff',      # White card background
            'heading': '#343a40',      # Dark heading
        }
        
        # Set default Tk fonts to be larger (increased for better readability)
        default_font = tkfont.nametofont("TkDefaultFont")
        default_font.configure(family="Helvetica", size=14)
        
        text_font = tkfont.nametofont("TkTextFont")
        text_font.configure(family="Helvetica", size=14)
        
        fixed_font = tkfont.nametofont("TkFixedFont")
        fixed_font.configure(family="Monospace", size=13)
        
        # Font definitions (increased for better readability)
        self.fonts = {
            'default': ("Helvetica", 14),
            'heading': ("Helvetica", 16, "bold"),
            'subheading': ("Helvetica", 15, "bold"),
            'small': ("Helvetica", 12),
            'button': ("Helvetica", 14),
            'title': ("Helvetica", 20, "bold"),
        }
        
        # Configure root window
        self.root.configure(bg=self.colors['bg'])
        self.root.option_add('*Font', 'Helvetica 14')
        self.root.option_add('*Dialog.msg.font', 'Helvetica 14')
        
        # Base styles
        style.configure(".", 
                       font=self.fonts['default'], 
                       foreground=self.colors['fg'],
                       background=self.colors['bg'])
        
        # Frame styles
        style.configure("TFrame", background=self.colors['bg'])
        style.configure("Card.TFrame", background=self.colors['card_bg'], relief="solid", borderwidth=1)
        
        # Label styles
        style.configure("TLabel", 
                       font=self.fonts['default'], 
                       padding=3, 
                       foreground=self.colors['fg'],
                       background=self.colors['bg'])
        style.configure("Heading.TLabel", 
                       font=self.fonts['heading'], 
                       foreground=self.colors['heading'],
                       padding=(0, 5))
        style.configure("Title.TLabel", 
                       font=self.fonts['title'], 
                       foreground=self.colors['accent'],
                       padding=(0, 10))
        style.configure("Desc.TLabel", 
                       font=self.fonts['small'], 
                       foreground=self.colors['secondary'],
                       padding=2)
        
        # Button styles
        style.configure("TButton", 
                       font=self.fonts['button'], 
                       padding=(12, 8),
                       foreground=self.colors['fg'])
        style.map("TButton",
                 foreground=[('active', self.colors['accent'])],
                 background=[('active', self.colors['border'])])
        
        style.configure("Accent.TButton", 
                       font=self.fonts['button'], 
                       foreground="white",
                       background=self.colors['accent'],
                       padding=(15, 10))
        
        style.configure("Run.TButton", 
                       font=("Helvetica", 12, "bold"), 
                       foreground="white",
                       background=self.colors['success'],
                       padding=(20, 10))
        style.map("Run.TButton",
                 foreground=[('disabled', self.colors['secondary'])])
        
        style.configure("Stop.TButton", 
                       font=("Helvetica", 12, "bold"), 
                       foreground="white",
                       background=self.colors['danger'],
                       padding=(20, 10))
        
        style.configure("Toolbar.TButton", 
                       font=self.fonts['small'], 
                       padding=(10, 6),
                       foreground=self.colors['fg'])
        
        # Checkbutton and Radiobutton
        style.configure("TCheckbutton", 
                       font=self.fonts['default'], 
                       padding=5,
                       foreground=self.colors['fg'],
                       background=self.colors['bg'])
        style.configure("TRadiobutton", 
                       font=self.fonts['default'], 
                       padding=5,
                       foreground=self.colors['fg'],
                       background=self.colors['bg'])
        
        # Entry and Combobox
        style.configure("TEntry", 
                       font=self.fonts['default'], 
                       padding=8,
                       foreground=self.colors['fg'])
        style.configure("TCombobox", 
                       font=self.fonts['default'], 
                       padding=8,
                       foreground=self.colors['fg'])
        style.configure("TSpinbox", 
                       font=self.fonts['default'], 
                       padding=8,
                       foreground=self.colors['fg'])
        
        # LabelFrame
        style.configure("TLabelframe", 
                       font=self.fonts['subheading'], 
                       padding=15,
                       background=self.colors['card_bg'],
                       relief="solid",
                       borderwidth=1)
        style.configure("TLabelframe.Label", 
                       font=self.fonts['subheading'], 
                       foreground=self.colors['heading'],
                       background=self.colors['card_bg'])
        
        # Notebook tabs
        style.configure("TNotebook", background=self.colors['bg'])
        style.configure("TNotebook.Tab", 
                       font=self.fonts['default'], 
                       padding=(20, 10),
                       foreground=self.colors['fg'])
        style.map("TNotebook.Tab",
                 foreground=[('selected', self.colors['accent'])],
                 background=[('selected', self.colors['card_bg'])])
        
        # Treeview (sample list) - note: font for content must be set via option_add
        # because style.configure doesn't reliably set treeview item fonts on Linux
        treeview_font = tkfont.Font(family="Helvetica", size=14)
        self.root.option_add('*Treeview*Font', treeview_font)
        
        style.configure("Treeview", 
                       font=self.fonts['default'], 
                       rowheight=42,  # Increase row height for better readability
                       foreground=self.colors['fg'],
                       background=self.colors['card_bg'],
                       fieldbackground=self.colors['card_bg'])
        style.configure("Treeview.Heading", 
                       font=("Helvetica", 15, "bold"),  # Larger heading font
                       padding=12,
                       foreground=self.colors['heading'],
                       background=self.colors['bg'])
        style.map("Treeview",
                 background=[('selected', self.colors['accent'])],
                 foreground=[('selected', 'white')])
        
        # Force treeview font globally (Linux workaround)
        style.layout("Treeview", [("Treeview.treearea", {"sticky": "nswe"})])
        
        # Progressbar
        style.configure("TProgressbar",
                       troughcolor=self.colors['border'],
                       background=self.colors['accent'],
                       thickness=20)
        
        # Separator
        style.configure("TSeparator", background=self.colors['border'])
        
        # Status bar
        style.configure("Status.TLabel", 
                       font=self.fonts['small'], 
                       padding=8,
                       foreground=self.colors['fg'],
                       background=self.colors['border'])
        style.configure("StatusBar.TFrame", background=self.colors['border'])
    
    def _create_menu(self):
        """Create the menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="New Configuration", command=self._new_config)
        file_menu.add_command(label="Open Configuration...", command=self._open_config)
        file_menu.add_command(label="Save Configuration", command=self._save_config)
        file_menu.add_command(label="Save Configuration As...", command=self._save_config_as)
        file_menu.add_separator()
        file_menu.add_command(label="Import Samples...", command=self._import_samples)
        file_menu.add_command(label="Export Samples...", command=self._export_samples)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self._on_close)
        
        # Workflow menu
        workflow_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Workflow", menu=workflow_menu)
        workflow_menu.add_command(label="Run Core Workflow", command=lambda: self._run_workflow("core"))
        workflow_menu.add_command(label="Run Advanced Workflow", command=lambda: self._run_workflow("advanced"))
        workflow_menu.add_separator()
        workflow_menu.add_command(label="Dry Run", command=self._dry_run)
        workflow_menu.add_command(label="Stop Workflow", command=self._stop_workflow)
        
        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="Documentation", command=self._show_docs)
        help_menu.add_command(label="About", command=self._show_about)
    
    def _create_header(self):
        """Create the application header with title and quick actions"""
        header = ttk.Frame(self.root, padding=(20, 15))
        header.pack(fill=tk.X)
        
        # Left side: Logo and title
        title_frame = ttk.Frame(header)
        title_frame.pack(side=tk.LEFT)
        
        ttk.Label(title_frame, text="🧬 MEDIPIPE", style="Title.TLabel").pack(side=tk.LEFT)
        ttk.Label(title_frame, text="cfMeDIP-seq Analysis Pipeline", 
                  style="Desc.TLabel").pack(side=tk.LEFT, padx=(10, 0), pady=(8, 0))
        
        # Right side: Quick actions
        actions_frame = ttk.Frame(header)
        actions_frame.pack(side=tk.RIGHT)
        
        # Run button - prominent
        self.run_btn = ttk.Button(actions_frame, text="▶  Run Analysis", 
                                  command=self._run_workflow_default,
                                  style="Run.TButton", width=16)
        self.run_btn.pack(side=tk.LEFT, padx=(0, 10))
        
        # Stop button
        self.stop_btn = ttk.Button(actions_frame, text="⬛  Stop", 
                                   command=self._stop_workflow,
                                   style="Stop.TButton", width=10, state=tk.DISABLED)
        self.stop_btn.pack(side=tk.LEFT, padx=(0, 10))
        
        # Separator line
        ttk.Separator(self.root, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=10)
    
    def _create_toolbar(self):
        """Create the toolbar with organized button groups"""
        toolbar = ttk.Frame(self.root, padding=(15, 10))
        toolbar.pack(fill=tk.X)
        
        # === Configuration Group ===
        config_group = ttk.LabelFrame(toolbar, text="Configuration", padding=(10, 5))
        config_group.pack(side=tk.LEFT, padx=(0, 15))
        
        ttk.Button(config_group, text="📂 Open", command=self._open_config, 
                   width=10).pack(side=tk.LEFT, padx=2)
        ttk.Button(config_group, text="💾 Save", command=self._save_config, 
                   width=10).pack(side=tk.LEFT, padx=2)
        ttk.Button(config_group, text="📝 New", command=self._new_config, 
                   width=10).pack(side=tk.LEFT, padx=2)
        
        # === Samples Group ===
        samples_group = ttk.LabelFrame(toolbar, text="Samples", padding=(10, 5))
        samples_group.pack(side=tk.LEFT, padx=(0, 15))
        
        ttk.Button(samples_group, text="📁 Add Files", command=self._add_files, 
                   width=12).pack(side=tk.LEFT, padx=2)
        ttk.Button(samples_group, text="🔍 Auto-Discover", command=self._discover_samples, 
                   width=14).pack(side=tk.LEFT, padx=2)
        ttk.Button(samples_group, text="📋 Import TSV", command=self._import_samples, 
                   width=12).pack(side=tk.LEFT, padx=2)
        
        # === Workflow Group ===
        workflow_group = ttk.LabelFrame(toolbar, text="Workflow", padding=(10, 5))
        workflow_group.pack(side=tk.LEFT, padx=(0, 15))
        
        ttk.Button(workflow_group, text="🔎 Dry Run", command=self._dry_run, 
                   width=12).pack(side=tk.LEFT, padx=2)
        ttk.Button(workflow_group, text="📊 View DAG", command=self._show_dag, 
                   width=12).pack(side=tk.LEFT, padx=2)
    
    def _create_main_layout(self):
        """Create the main application layout"""
        # Main paned window
        main_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel: Samples and File Browser
        left_frame = ttk.Frame(main_paned)
        main_paned.add(left_frame, weight=1)
        
        self._create_samples_panel(left_frame)
        
        # Right panel: Configuration and Progress
        right_paned = ttk.PanedWindow(main_paned, orient=tk.VERTICAL)
        main_paned.add(right_paned, weight=2)
        
        # Top right: Configuration
        config_frame = ttk.Frame(right_paned)
        right_paned.add(config_frame, weight=1)
        
        self._create_config_area(config_frame)
        
        # Bottom right: Hardware Monitor + Progress
        progress_frame = ttk.Frame(right_paned)
        right_paned.add(progress_frame, weight=1)
        
        # Hardware monitor at the top of progress section
        self.hardware_monitor = HardwareMonitor(
            progress_frame, 
            on_cores_change=self._on_cores_changed
        )
        self.hardware_monitor.pack(fill=tk.X, pady=(0, 5))
        
        # Progress panel below hardware monitor
        self.progress_panel = ProgressPanel(progress_frame)
        self.progress_panel.pack(fill=tk.BOTH, expand=True)
    
    def _create_samples_panel(self, parent):
        """Create the samples panel"""
        samples_frame = ttk.LabelFrame(parent, text="Samples", padding=10)
        samples_frame.pack(fill=tk.BOTH, expand=True)
        
        # Sample list
        list_frame = ttk.Frame(samples_frame)
        list_frame.pack(fill=tk.BOTH, expand=True)
        
        columns = ("sample_id", "group", "status")
        self.sample_list = ttk.Treeview(list_frame, columns=columns, show="headings", selectmode="extended")
        
        # Configure font tag for readable text (Linux fix)
        self.sample_list.tag_configure('default', font=self.fonts['default'])
        
        self.sample_list.heading("sample_id", text="Sample ID")
        self.sample_list.heading("group", text="Group")
        self.sample_list.heading("status", text="Status")
        
        self.sample_list.column("sample_id", width=200, minwidth=150)
        self.sample_list.column("group", width=120, minwidth=100)
        self.sample_list.column("status", width=80, minwidth=60)
        
        sample_scroll = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self.sample_list.yview)
        self.sample_list.configure(yscrollcommand=sample_scroll.set)
        
        self.sample_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sample_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Group selector frame
        group_frame = ttk.Frame(samples_frame)
        group_frame.pack(fill=tk.X, pady=(10, 5))
        
        ttk.Label(group_frame, text="Assign to Group:", font=("Helvetica", 10)).pack(side=tk.LEFT, padx=(0, 5))
        
        # Group dropdown
        self.group_select_var = tk.StringVar(value="")
        self.group_combo = ttk.Combobox(group_frame, textvariable=self.group_select_var, 
                                        width=15, state="normal")
        self.group_combo['values'] = ["default", "CTRL", "CASE", "AEG", "Treatment", "Control"]
        self.group_combo.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(group_frame, text="Apply to Selected", 
                   command=self._apply_group_to_selected, width=14).pack(side=tk.LEFT, padx=5)
        
        ttk.Separator(group_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        ttk.Button(group_frame, text="+ New Group", 
                   command=self._create_new_group, width=12).pack(side=tk.LEFT, padx=5)
        
        # Sample buttons
        btn_frame = ttk.Frame(samples_frame)
        btn_frame.pack(fill=tk.X, pady=(10, 5))
        
        ttk.Button(btn_frame, text="Add...", command=self._add_files, width=12).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Remove", command=self._remove_samples, width=12).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Validate", command=self._validate_samples, width=12).pack(side=tk.LEFT, padx=4)
        
        # Sample count label
        self.sample_count_var = tk.StringVar(value="0 samples")
        ttk.Label(btn_frame, textvariable=self.sample_count_var, font=("Helvetica", 11, "bold")).pack(side=tk.RIGHT, padx=5)
    
    def _create_config_area(self, parent):
        """Create the configuration area"""
        notebook = ttk.Notebook(parent)
        notebook.pack(fill=tk.BOTH, expand=True)
        
        # Tab 1: Workflow configuration
        workflow_frame = ttk.Frame(notebook, padding=5)
        notebook.add(workflow_frame, text="Workflow")
        
        self.workflow_panel = WorkflowPanel(workflow_frame, on_change=self._on_workflow_config_change)
        self.workflow_panel.pack(fill=tk.BOTH, expand=True)
        
        # Tab 2: Analysis parameters
        config_frame = ttk.Frame(notebook, padding=5)
        notebook.add(config_frame, text="Parameters")
        
        self.config_panel = ConfigPanel(config_frame, on_change=self._on_config_change)
        self.config_panel.pack(fill=tk.BOTH, expand=True)
        
        # Tab 3: QC Reports
        qc_frame = ttk.Frame(notebook, padding=5)
        notebook.add(qc_frame, text="QC Reports")
        
        self.qc_panel = QCPanel(qc_frame, output_dir=self.current_config.output_dir if self.current_config else "/data/medipipe_data/output")
        self.qc_panel.pack(fill=tk.BOTH, expand=True)
        
        # Tab 4: ML Analysis
        ml_frame = ttk.Frame(notebook, padding=5)
        notebook.add(ml_frame, text="ML Analysis")
        
        self.ml_panel = MLAnalysisPanel(ml_frame, on_log=self._on_ml_log)
        self.ml_panel.pack(fill=tk.BOTH, expand=True)
    
    def _create_statusbar(self):
        """Create the status bar"""
        statusbar = ttk.Frame(self.root, padding=(10, 8))
        statusbar.pack(fill=tk.X, side=tk.BOTTOM)
        
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(statusbar, textvariable=self.status_var, style="Status.TLabel").pack(side=tk.LEFT, padx=10)
        
        ttk.Separator(statusbar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        self.workflow_status_var = tk.StringVar(value="Idle")
        ttk.Label(statusbar, textvariable=self.workflow_status_var, 
                  style="Status.TLabel", font=("Helvetica", 10, "bold")).pack(side=tk.LEFT, padx=10)
        
        ttk.Separator(statusbar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        # Compact hardware monitor in status bar
        self.compact_hw_monitor = CompactHardwareMonitor(statusbar)
        self.compact_hw_monitor.pack(side=tk.LEFT, padx=10)
        
        # Right side info
        ttk.Label(statusbar, text=f"v{self.VERSION}", style="Status.TLabel").pack(side=tk.RIGHT, padx=10)
    
    def _load_default_config(self):
        """Load default configuration"""
        self.current_config = MedipipeConfig()
        self._update_ui_from_config()
        self.status_var.set("Loaded default configuration")
    
    def _update_ui_from_config(self):
        """Update UI components from current config"""
        if not self.current_config:
            return
        
        # Update config panel
        config_dict = {
            "genome": self.current_config.genome,
            "window_size": self.current_config.window_size,
            "threads": self.current_config.threads,
            "input_dir": self.current_config.input_dir,
            "output_dir": self.current_config.output_dir,
            "samples": self.current_config.samples_file,
        }
        self.config_panel.set_config(config_dict)
        
        # Update hardware monitor cores
        if hasattr(self, 'hardware_monitor'):
            self.hardware_monitor.set_cores(self.current_config.threads)
        
        # Update workflow panel
        workflow_dict = {
            "workflow": {
                "preprocessing": self.current_config.run_preprocessing,
                "qc_analysis": self.current_config.run_qc,
                "methylation_quant": self.current_config.run_methylation,
                "fragment_profile": self.current_config.run_fragment_profile,
                "tissue_of_origin": self.current_config.run_tissue_of_origin,
            },
            "advanced": {
                "dmr_analysis": self.current_config.run_dmr,
                "ml_classification": self.current_config.run_ml,
                "dl_classification": self.current_config.run_dl,
            },
            "data_options": {
                "paired_end": self.current_config.paired_end,
                "use_umi": self.current_config.use_umi,
                "use_spikein": self.current_config.use_spikein,
                "umi_regex_r1": self.current_config.umi_regex_r1,
                "umi_regex_r2": self.current_config.umi_regex_r2,
            },
        }
        self.workflow_panel.set_config(workflow_dict)
        
        # Update QC panel output directory
        if hasattr(self, 'qc_panel'):
            self.qc_panel.set_output_dir(self.current_config.output_dir)
    
    def _on_cores_changed(self, cores: int):
        """Handle cores change from hardware monitor"""
        if self.current_config:
            self.current_config.threads = cores
            # Also update config panel if it has a threads field
            try:
                self.config_panel.set_config({"threads": cores})
            except:
                pass
            self.status_var.set(f"CPU cores set to {cores}")
    
    def _update_config_from_ui(self):
        """Update current config from UI components"""
        if not self.current_config:
            self.current_config = MedipipeConfig()
        
        # Get from config panel
        config_dict = self.config_panel.get_config()
        self.current_config.genome = config_dict.get("genome", "hg38")
        self.current_config.window_size = config_dict.get("window_size", 300)
        
        # Get threads from hardware monitor (takes precedence)
        if hasattr(self, 'hardware_monitor'):
            self.current_config.threads = self.hardware_monitor.get_cores()
        else:
            self.current_config.threads = config_dict.get("threads", 24)
        self.current_config.input_dir = config_dict.get("input_dir", "/data/raw_data")
        self.current_config.output_dir = config_dict.get("output_dir", "/data/medipipe_data/output")
        self.current_config.samples_file = config_dict.get("samples", "")
        
        # Get from workflow panel
        workflow_dict = self.workflow_panel.get_config()
        workflow = workflow_dict.get("workflow", {})
        self.current_config.run_preprocessing = workflow.get("preprocessing", True)
        self.current_config.run_qc = workflow.get("qc_analysis", True)
        self.current_config.run_methylation = workflow.get("methylation_quant", True)
        self.current_config.run_fragment_profile = workflow.get("fragment_profile", True)
        self.current_config.run_tissue_of_origin = workflow.get("tissue_of_origin", True)
        
        advanced = workflow_dict.get("advanced", {})
        self.current_config.run_dmr = advanced.get("dmr_analysis", False)
        self.current_config.run_ml = advanced.get("ml_classification", False)
        self.current_config.run_dl = advanced.get("dl_classification", False)
        
        data_opts = workflow_dict.get("data_options", {})
        self.current_config.paired_end = data_opts.get("paired_end", True)
        self.current_config.use_umi = data_opts.get("use_umi", False)
        self.current_config.use_spikein = data_opts.get("use_spikein", False)
        self.current_config.umi_regex_r1 = data_opts.get("umi_regex_r1", r"(?P<umi_1>.{5})(?P<discard_1>.{2}).*")
        self.current_config.umi_regex_r2 = data_opts.get("umi_regex_r2", r"(?P<umi_2>.{5})(?P<discard_2>.{2}).*")
    
    def _update_sample_list(self):
        """Update the sample list display"""
        # Clear current items
        for item in self.sample_list.get_children():
            self.sample_list.delete(item)
        
        # Add samples
        for sample in self.sample_manager.samples:
            status = "✓" if sample.is_valid else "✗"
            self.sample_list.insert("", "end", sample.sample_id,
                                   values=(sample.sample_id, sample.group, status),
                                   tags=('default',))
        
        # Update count
        valid = len(self.sample_manager.get_valid_samples())
        total = len(self.sample_manager.samples)
        self.sample_count_var.set(f"{valid}/{total} valid samples")
        
        # Update available groups in workflow panel and sample panel dropdown
        groups = self.sample_manager.get_groups()
        self.workflow_panel.update_available_groups(groups)
        
        # Update group dropdown in sample panel
        if hasattr(self, 'group_combo'):
            self._update_group_combo()
    
    # === File Operations ===
    
    def _new_config(self):
        """Create a new configuration"""
        if messagebox.askyesno("New Configuration", "Create new configuration? Unsaved changes will be lost."):
            self.current_config = MedipipeConfig()
            self._update_ui_from_config()
            self.sample_manager.samples = []
            self._update_sample_list()
            self.status_var.set("New configuration created")
    
    def _open_config(self):
        """Open a configuration file"""
        filepath = filedialog.askopenfilename(
            initialdir=os.path.join(self.PIPE_DIR, "configfiles"),
            title="Open Configuration",
            filetypes=[("YAML files", "*.yaml *.yml"), ("All files", "*.*")]
        )
        
        if filepath:
            try:
                self.current_config = self.config_manager.load_config(filepath)
                self._update_ui_from_config()
                self.status_var.set(f"Loaded: {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load configuration:\n{e}")
    
    def _save_config(self):
        """Save the current configuration"""
        self._update_config_from_ui()
        
        if self.config_manager._config_path:
            try:
                self.config_manager.save_config(self.current_config)
                self.status_var.set(f"Saved: {os.path.basename(self.config_manager._config_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save configuration:\n{e}")
        else:
            self._save_config_as()
    
    def _save_config_as(self):
        """Save configuration to a new file"""
        self._update_config_from_ui()
        
        filepath = filedialog.asksaveasfilename(
            initialdir=os.path.join(self.PIPE_DIR, "configfiles"),
            title="Save Configuration As",
            defaultextension=".yaml",
            filetypes=[("YAML files", "*.yaml"), ("All files", "*.*")]
        )
        
        if filepath:
            try:
                self.config_manager.save_config(self.current_config, filepath)
                self.status_var.set(f"Saved: {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save configuration:\n{e}")
    
    # === Sample Operations ===
    
    def _add_files(self):
        """Add FASTQ files using file browser dialog"""
        dialog = FileSelectionDialog(
            self.root,
            title="Select FASTQ Files",
            initial_dir=self.current_config.input_dir if self.current_config else "/data/raw_data",
            file_filter="fastq",
            multi_select=True
        )
        
        files = dialog.show()
        
        if files:
            # Add ONLY the selected files as samples (not the whole directory)
            new_samples = self.sample_manager.add_files(
                files,
                paired_end=self.current_config.paired_end if self.current_config else True
            )
            self._update_sample_list()
            self.status_var.set(f"Added {len(new_samples)} samples from {len(files)} files")
    
    def _discover_samples(self):
        """Auto-discover samples in input directory"""
        input_dir = self.current_config.input_dir if self.current_config else "/data/raw_data"
        
        if not os.path.exists(input_dir):
            messagebox.showerror("Error", f"Input directory not found:\n{input_dir}")
            return
        
        samples = self.sample_manager.discover_samples(
            input_dir,
            paired_end=self.current_config.paired_end if self.current_config else True
        )
        
        self._update_sample_list()
        self.status_var.set(f"Discovered {len(samples)} samples")
    
    def _remove_samples(self):
        """Remove selected samples"""
        selection = self.sample_list.selection()
        if not selection:
            return
        
        for item in selection:
            self.sample_manager.remove_sample(item)
        
        self._update_sample_list()
    
    def _apply_group_to_selected(self):
        """Apply selected group to selected samples"""
        selection = self.sample_list.selection()
        if not selection:
            messagebox.showwarning("No Selection", "Please select samples to assign to a group")
            return
        
        group = self.group_select_var.get().strip()
        if not group:
            messagebox.showwarning("No Group", "Please select or enter a group name")
            return
        
        for item in selection:
            self.sample_manager.update_sample_group(item, group)
        
        self._update_sample_list()
        self._update_group_combo()
        self.status_var.set(f"Assigned {len(selection)} samples to group '{group}'")
    
    def _create_new_group(self):
        """Create a new group via dialog"""
        dialog = tk.Toplevel(self.root)
        dialog.title("Create New Group")
        dialog.geometry("350x150")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Center dialog
        dialog.update_idletasks()
        x = (dialog.winfo_screenwidth() - dialog.winfo_width()) // 2
        y = (dialog.winfo_screenheight() - dialog.winfo_height()) // 2
        dialog.geometry(f"+{x}+{y}")
        
        ttk.Label(dialog, text="Enter new group name:", font=("Helvetica", 11)).pack(pady=(15, 5))
        
        group_var = tk.StringVar()
        entry = ttk.Entry(dialog, textvariable=group_var, width=30, font=("Helvetica", 11))
        entry.pack(pady=10)
        entry.focus()
        
        def add_group():
            group = group_var.get().strip()
            if group:
                # Add to combo box values
                current_values = list(self.group_combo['values'])
                if group not in current_values:
                    current_values.append(group)
                    self.group_combo['values'] = current_values
                self.group_select_var.set(group)
                self.status_var.set(f"Created new group '{group}'")
            dialog.destroy()
        
        entry.bind("<Return>", lambda e: add_group())
        
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Create", command=add_group, width=12).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=dialog.destroy, width=12).pack(side=tk.LEFT, padx=5)
    
    def _update_group_combo(self):
        """Update group combo box with available groups from samples"""
        existing_groups = self.sample_manager.get_groups()
        default_groups = ["default", "CTRL", "CASE", "AEG", "Treatment", "Control"]
        all_groups = list(set(default_groups + existing_groups))
        all_groups.sort()
        self.group_combo['values'] = all_groups
    
    def _validate_samples(self):
        """Validate all samples and show status"""
        valid, invalid = self.sample_manager.validate_all(
            paired_end=self.current_config.paired_end if self.current_config else True
        )
        self._update_sample_list()
        
        if invalid == 0:
            messagebox.showinfo("Validation", f"All {valid} samples are valid! ✓")
        else:
            messagebox.showwarning("Validation", f"Validation complete:\n• {valid} valid samples\n• {invalid} invalid samples\n\nCheck the Status column for details.")
    
    def _set_sample_group(self):
        """Set group for selected samples (legacy method - now uses dropdown)"""
        self._apply_group_to_selected()
    
    def _import_samples(self):
        """Import samples from TSV file"""
        filepath = filedialog.askopenfilename(
            initialdir=os.path.join(self.PIPE_DIR, "samplefiles"),
            title="Import Samples",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        
        if filepath:
            try:
                self.sample_manager.load_samples_tsv(
                    filepath,
                    paired_end=self.current_config.paired_end if self.current_config else True
                )
                self._update_sample_list()
                self.status_var.set(f"Imported samples from {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to import samples:\n{e}")
    
    def _export_samples(self):
        """Export samples to TSV file"""
        if not self.sample_manager.samples:
            messagebox.showwarning("Warning", "No samples to export")
            return
        
        filepath = filedialog.asksaveasfilename(
            initialdir=os.path.join(self.PIPE_DIR, "samplefiles"),
            title="Export Samples",
            defaultextension=".tsv",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        
        if filepath:
            try:
                self.sample_manager.save_samples_tsv(filepath)
                self.status_var.set(f"Exported samples to {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to export samples:\n{e}")
    
    # === Workflow Operations ===
    
    def _run_workflow_default(self):
        """Run workflow based on what's selected"""
        if self.workflow_panel.is_core_workflow_selected():
            self._run_workflow("core")
        elif self.workflow_panel.is_advanced_workflow_selected():
            self._run_workflow("advanced")
        else:
            messagebox.showwarning("Warning", "No workflow stages selected")
    
    def _run_workflow(self, workflow: str = "core"):
        """Run the specified workflow"""
        # Validate
        if not self.sample_manager.get_valid_samples():
            messagebox.showerror("Error", "No valid samples. Please add samples before running.")
            return
        
        # Update and save config
        self._update_config_from_ui()
        
        # Save samples file
        samples_path = os.path.join(self.PIPE_DIR, "samplefiles", "current_samples.tsv")
        self.sample_manager.save_samples_tsv(samples_path)
        self.current_config.samples_file = samples_path
        
        # Create working config
        config_path = self.config_manager.create_working_config(self.current_config)
        
        # Start workflow
        self.progress_panel.reset()
        self.progress_panel.set_status("Starting...")
        
        success = self.workflow_controller.start_workflow(
            config_path,
            workflow=workflow,
            cores=self.current_config.threads,
            use_conda=True
        )
        
        if success:
            self.run_btn.configure(state=tk.DISABLED)
            self.stop_btn.configure(state=tk.NORMAL)
            self.workflow_status_var.set(f"Running: {workflow}")
            self.status_var.set(f"Started {workflow} workflow")
        else:
            messagebox.showerror("Error", "Failed to start workflow")
    
    def _dry_run(self):
        """Perform a dry run"""
        self._update_config_from_ui()
        
        if not self.sample_manager.samples:
            messagebox.showwarning("Warning", "No samples loaded")
            return
        
        # Save samples and config
        samples_path = os.path.join(self.PIPE_DIR, "samplefiles", "current_samples.tsv")
        self.sample_manager.save_samples_tsv(samples_path)
        self.current_config.samples_file = samples_path
        config_path = self.config_manager.create_working_config(self.current_config)
        
        # Determine workflow
        workflow = "core" if self.workflow_panel.is_core_workflow_selected() else "advanced"
        
        # Run dry-run
        self.progress_panel.clear_log()
        self.progress_panel.add_log(f"=== Dry Run ({workflow} workflow) ===")
        
        success = self.workflow_controller.run_dry_run(config_path, workflow)
        
        if success:
            self.progress_panel.add_log("Dry run completed successfully")
            self.status_var.set("Dry run successful")
        else:
            self.progress_panel.add_log("Dry run failed")
            self.status_var.set("Dry run failed")
    
    def _stop_workflow(self):
        """Stop the running workflow"""
        if messagebox.askyesno("Stop Workflow", "Are you sure you want to stop the workflow?"):
            self.workflow_controller.stop_workflow()
            self.run_btn.configure(state=tk.NORMAL)
            self.stop_btn.configure(state=tk.DISABLED)
    
    def _show_dag(self):
        """Show the workflow DAG (directed acyclic graph)"""
        self._update_config_from_ui()
        
        if not self.sample_manager.samples:
            messagebox.showwarning("Warning", "No samples loaded. Add samples to visualize the workflow.")
            return
        
        # Save samples and config
        samples_path = os.path.join(self.PIPE_DIR, "samplefiles", "current_samples.tsv")
        try:
            self.sample_manager.save_samples_tsv(samples_path)
            self.current_config.samples_file = samples_path
            config_path = self.config_manager.create_working_config(self.current_config)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to prepare config: {e}")
            return
        
        # Determine workflow
        workflow = "core" if self.workflow_panel.is_core_workflow_selected() else "advanced"
        snakefile = f"{workflow}_workflow.smk"
        
        self.progress_panel.add_log(f"Generating DAG for {workflow} workflow...")
        
        # Run Snakemake with --dag
        import subprocess
        dag_path = os.path.join(self.current_config.work_dir, "workflow_dag.svg")
        
        try:
            cmd = [
                "snakemake", "--dag",
                "-s", os.path.join(self.PIPE_DIR, "snakefiles", snakefile),
                "--configfile", config_path
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.PIPE_DIR)
            
            if result.returncode == 0:
                # Convert dot to SVG
                dot_cmd = ["dot", "-Tsvg", "-o", dag_path]
                dot_result = subprocess.run(dot_cmd, input=result.stdout, capture_output=True, text=True)
                
                if dot_result.returncode == 0 and os.path.exists(dag_path):
                    # Try to open with default viewer
                    subprocess.Popen(["xdg-open", dag_path])
                    self.progress_panel.add_log(f"DAG saved to {dag_path}")
                else:
                    # Just show the text output
                    self.progress_panel.add_log("DAG (text):\n" + result.stdout[:1000])
            else:
                self.progress_panel.add_log(f"DAG generation failed: {result.stderr}")
        except Exception as e:
            self.progress_panel.add_log(f"Error generating DAG: {e}")
    
    # === Workflow Callbacks ===
    
    def _on_workflow_progress(self, progress):
        """Handle workflow progress updates"""
        self.root.after(0, lambda: self.progress_panel.update_from_progress(progress))
        self.root.after(0, lambda: self.workflow_status_var.set(progress.status.value.title()))
    
    def _on_workflow_log(self, message: str):
        """Handle workflow log messages"""
        self.root.after(0, lambda: self.progress_panel.add_log(message))
    
    def _on_workflow_complete(self, success: bool):
        """Handle workflow completion"""
        def update():
            self.run_btn.configure(state=tk.NORMAL)
            self.stop_btn.configure(state=tk.DISABLED)
            
            if success:
                self.workflow_status_var.set("Completed")
                self.status_var.set("Workflow completed successfully")
                # Refresh QC reports after successful workflow
                if hasattr(self, 'qc_panel'):
                    self.qc_panel.refresh_qc_status()
                messagebox.showinfo("Success", "Workflow completed successfully!")
            else:
                self.workflow_status_var.set("Failed")
                self.status_var.set("Workflow failed")
                messagebox.showerror("Error", "Workflow failed. Check the log for details.")
        
        self.root.after(0, update)
    
    # === Config Change Handlers ===
    
    def _on_workflow_config_change(self, config: dict):
        """Handle workflow configuration changes"""
        pass  # Config is read when needed
    
    def _on_config_change(self, config: dict):
        """Handle parameter configuration changes"""
        # Update QC panel if output_dir changed
        if 'output_dir' in config and hasattr(self, 'qc_panel'):
            self.qc_panel.set_output_dir(config['output_dir'])
    
    def _on_ml_log(self, message: str):
        """Handle ML panel log messages - forward to progress panel"""
        if hasattr(self, 'progress_panel'):
            self.progress_panel.add_log(f"[ML] {message}")
    
    # === Help Menu ===
    
    def _show_docs(self):
        """Show documentation"""
        messagebox.showinfo("Documentation", 
            "MEDIPIPE Documentation\n\n"
            "For detailed documentation, please visit:\n"
            "https://github.com/medipipe/docs\n\n"
            "Or check the README.md in the project directory.")
    
    def _show_about(self):
        """Show about dialog"""
        messagebox.showinfo("About MEDIPIPE",
            f"MEDIPIPE v{self.VERSION}\n\n"
            "cfMeDIP-seq Analysis Pipeline\n\n"
            "A comprehensive pipeline for cell-free methylated DNA\n"
            "immunoprecipitation sequencing analysis.\n\n"
            "Features:\n"
            "• Quality control and preprocessing\n"
            "• Methylation quantification\n"
            "• Fragment profile analysis\n"
            "• Tissue of origin deconvolution\n"
            "• Machine learning classification\n"
            "• Deep learning models")
    
    def _on_close(self):
        """Handle window close"""
        if self.workflow_controller.is_running():
            if not messagebox.askyesno("Confirm Exit", 
                "A workflow is running. Are you sure you want to exit?"):
                return
            self.workflow_controller.stop_workflow()
        
        self.root.destroy()
    
    def run(self):
        """Run the application"""
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
        self.root.mainloop()


if __name__ == "__main__":
    app = MedipipeApp()
    app.run()
