#!/usr/bin/env python3
"""
MEDIPIPE ML Analysis GUI - Enhanced Version
Modern UI with file browser and preview panel for cfMeDIP-seq analysis

Features:
- Modern styling with larger fonts
- Left panel: Output file tree browser
- Right panel: File preview (images, text, tables)
- ML model configuration (L2, ElasticNet, both)
- Regularization parameter tuning (C, alpha)

Author: MEDIPIPE
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import subprocess
import threading
import os
from pathlib import Path
from datetime import datetime
import json

# Try to import PIL for image preview
try:
    from PIL import Image, ImageTk
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    print("Note: Install Pillow (pip install Pillow) for image preview support")

# Try to import pandas for table preview
try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


class ModernStyle:
    """Modern styling constants"""
    # Colors - Dark theme
    BG_DARK = "#1e1e2e"
    BG_MEDIUM = "#313244"
    BG_LIGHT = "#45475a"
    BG_INPUT = "#11111b"
    FG_PRIMARY = "#cdd6f4"
    FG_SECONDARY = "#a6adc8"
    FG_DIM = "#6c7086"
    ACCENT = "#89b4fa"
    ACCENT_HOVER = "#b4befe"
    SUCCESS = "#a6e3a1"
    WARNING = "#f9e2af"
    ERROR = "#f38ba8"
    
    # Fonts - larger sizes
    FONT_FAMILY = "Ubuntu"
    FONT_SIZE_SMALL = 11
    FONT_SIZE_NORMAL = 12
    FONT_SIZE_LARGE = 14
    FONT_SIZE_TITLE = 18
    FONT_SIZE_HEADER = 16
    
    @classmethod
    def configure_styles(cls, root):
        """Configure ttk styles for modern look"""
        style = ttk.Style(root)
        
        # Try to use clam theme as base
        try:
            style.theme_use('clam')
        except:
            pass
        
        # Configure default styles
        style.configure(".",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            background=cls.BG_MEDIUM,
            foreground=cls.FG_PRIMARY,
            borderwidth=0)
        
        # Frame styles
        style.configure("TFrame", background=cls.BG_MEDIUM)
        style.configure("Dark.TFrame", background=cls.BG_DARK)
        style.configure("Card.TFrame", background=cls.BG_LIGHT, relief="flat")
        
        # Label styles
        style.configure("TLabel",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            background=cls.BG_MEDIUM,
            foreground=cls.FG_PRIMARY,
            padding=2)
        style.configure("Title.TLabel",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_TITLE, "bold"),
            foreground=cls.ACCENT,
            padding=5)
        style.configure("Header.TLabel",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_HEADER, "bold"),
            foreground=cls.ACCENT,
            padding=3)
        style.configure("Dim.TLabel",
            foreground=cls.FG_DIM)
        
        # Button styles
        style.configure("TButton",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            padding=(16, 10),
            background=cls.BG_LIGHT)
        style.map("TButton",
            background=[('active', cls.BG_DARK), ('pressed', cls.BG_DARK)])
        
        style.configure("Accent.TButton",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_LARGE, "bold"),
            padding=(20, 12),
            background=cls.ACCENT)
        style.map("Accent.TButton",
            background=[('active', cls.ACCENT_HOVER)])
        
        style.configure("Small.TButton",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_SMALL),
            padding=(8, 4))
        
        # Entry styles
        style.configure("TEntry",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            padding=8,
            fieldbackground=cls.BG_INPUT,
            foreground=cls.FG_PRIMARY)
        
        # Combobox styles
        style.configure("TCombobox",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            padding=6,
            fieldbackground=cls.BG_INPUT)
        
        # Spinbox styles
        style.configure("TSpinbox",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            padding=6,
            fieldbackground=cls.BG_INPUT)
        
        # Notebook styles
        style.configure("TNotebook",
            background=cls.BG_DARK,
            borderwidth=0)
        style.configure("TNotebook.Tab",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_LARGE),
            padding=(25, 12),
            background=cls.BG_MEDIUM)
        style.map("TNotebook.Tab",
            background=[('selected', cls.BG_LIGHT)],
            foreground=[('selected', cls.ACCENT)])
        
        # Treeview styles
        style.configure("Treeview",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            rowheight=32,
            background=cls.BG_DARK,
            foreground=cls.FG_PRIMARY,
            fieldbackground=cls.BG_DARK,
            borderwidth=0)
        style.configure("Treeview.Heading",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL, "bold"),
            background=cls.BG_MEDIUM,
            foreground=cls.ACCENT)
        style.map("Treeview",
            background=[('selected', cls.BG_LIGHT)],
            foreground=[('selected', cls.ACCENT)])
        
        # LabelFrame styles
        style.configure("TLabelframe",
            background=cls.BG_MEDIUM,
            borderwidth=2,
            relief="flat")
        style.configure("TLabelframe.Label",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_LARGE, "bold"),
            foreground=cls.ACCENT,
            background=cls.BG_MEDIUM,
            padding=5)
        
        # Scale styles
        style.configure("TScale",
            background=cls.BG_MEDIUM,
            troughcolor=cls.BG_DARK)
        
        # Checkbutton styles  
        style.configure("TCheckbutton",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            background=cls.BG_MEDIUM,
            foreground=cls.FG_PRIMARY,
            padding=5)
        
        # Radiobutton styles
        style.configure("TRadiobutton",
            font=(cls.FONT_FAMILY, cls.FONT_SIZE_NORMAL),
            background=cls.BG_MEDIUM,
            foreground=cls.FG_PRIMARY,
            padding=5)
        
        # Progressbar
        style.configure("TProgressbar",
            background=cls.ACCENT,
            troughcolor=cls.BG_DARK,
            borderwidth=0,
            thickness=8)
        
        # Scrollbar
        style.configure("TScrollbar",
            background=cls.BG_LIGHT,
            troughcolor=cls.BG_DARK,
            borderwidth=0)


class FileTreePanel(ttk.Frame):
    """Left panel with output file tree browser"""
    
    def __init__(self, parent, on_file_select=None):
        super().__init__(parent)
        self.on_file_select = on_file_select
        self.current_dir = None
        self.setup_ui()
        
    def setup_ui(self):
        self.configure(style="Dark.TFrame")
        
        # Header
        header = ttk.Frame(self, style="Dark.TFrame")
        header.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Label(header, text="📁 Output Files", style="Header.TLabel").pack(side=tk.LEFT)
        ttk.Button(header, text="🔄", width=4, style="Small.TButton",
                   command=self.refresh).pack(side=tk.RIGHT)
        
        # Directory selector
        dir_frame = ttk.Frame(self, style="Dark.TFrame")
        dir_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.dir_var = tk.StringVar()
        self.dir_entry = ttk.Entry(dir_frame, textvariable=self.dir_var)
        self.dir_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))
        ttk.Button(dir_frame, text="...", width=4, style="Small.TButton",
                   command=self.browse_dir).pack(side=tk.RIGHT)
        
        # File tree
        tree_frame = ttk.Frame(self, style="Dark.TFrame")
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.tree = ttk.Treeview(tree_frame, show="tree", selectmode="browse")
        self.tree.column("#0", width=280, stretch=True)
        
        # Scrollbars
        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        
        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)
        
        # Bind selection
        self.tree.bind("<<TreeviewSelect>>", self._on_select)
        self.tree.bind("<Double-1>", self._on_double_click)
        
        # File type icons (using emoji)
        self.icons = {
            '.png': '🖼️', '.pdf': '📊', '.jpg': '🖼️', '.jpeg': '🖼️',
            '.tsv': '📋', '.csv': '📋', '.txt': '📝', '.json': '⚙️',
            '.log': '📜', '.html': '🌐', 'dir': '📁'
        }
        
    def browse_dir(self):
        directory = filedialog.askdirectory(
            title="Select Output Directory",
            initialdir=self.dir_var.get() or os.path.expanduser("~")
        )
        if directory:
            self.set_directory(directory)
            
    def set_directory(self, directory):
        self.current_dir = directory
        self.dir_var.set(directory)
        self.refresh()
        
    def refresh(self):
        """Refresh the file tree"""
        for item in self.tree.get_children():
            self.tree.delete(item)
            
        if not self.current_dir or not os.path.exists(self.current_dir):
            return
            
        self._populate_tree("", self.current_dir)
        
    def _populate_tree(self, parent, path, depth=0):
        """Recursively populate tree with files"""
        if depth > 3:
            return
            
        try:
            items = sorted(os.listdir(path))
        except PermissionError:
            return
            
        dirs = [i for i in items if os.path.isdir(os.path.join(path, i)) and not i.startswith('.')]
        files = [i for i in items if os.path.isfile(os.path.join(path, i)) and not i.startswith('.')]
        
        for name in dirs:
            full_path = os.path.join(path, name)
            icon = self.icons.get('dir', '📁')
            node = self.tree.insert(parent, "end", text=f"{icon} {name}",
                                   values=(full_path,), open=False)
            self._populate_tree(node, full_path, depth + 1)
            
        for name in files:
            full_path = os.path.join(path, name)
            ext = os.path.splitext(name)[1].lower()
            icon = self.icons.get(ext, '📄')
            self.tree.insert(parent, "end", text=f"{icon} {name}", values=(full_path,))
            
    def _on_select(self, event):
        selection = self.tree.selection()
        if selection and self.on_file_select:
            item = selection[0]
            values = self.tree.item(item, "values")
            if values:
                file_path = values[0]
                if os.path.isfile(file_path):
                    self.on_file_select(file_path)
                    
    def _on_double_click(self, event):
        """Open file with system default application"""
        selection = self.tree.selection()
        if selection:
            item = selection[0]
            values = self.tree.item(item, "values")
            if values:
                file_path = values[0]
                if os.path.isfile(file_path):
                    try:
                        subprocess.Popen(['xdg-open', file_path])
                    except:
                        pass


class PreviewPanel(ttk.Frame):
    """Right panel for file preview"""
    
    def __init__(self, parent):
        super().__init__(parent)
        self.current_image = None
        self.setup_ui()
        
    def setup_ui(self):
        self.configure(style="Dark.TFrame")
        
        # Header
        header = ttk.Frame(self, style="Dark.TFrame")
        header.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Label(header, text="👁️ Preview", style="Header.TLabel").pack(side=tk.LEFT)
        self.filename_var = tk.StringVar(value="No file selected")
        ttk.Label(header, textvariable=self.filename_var,
                 style="Dim.TLabel").pack(side=tk.RIGHT)
        
        # Preview area
        self.preview_frame = ttk.Frame(self, style="Dark.TFrame")
        self.preview_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Image preview label
        self.image_label = ttk.Label(self.preview_frame,
                                     text="Select a file from the tree\nto preview it here",
                                     anchor="center", justify="center")
        self.image_label.pack(fill=tk.BOTH, expand=True)
        
        # Text preview (hidden by default)
        self.text_frame = ttk.Frame(self.preview_frame, style="Dark.TFrame")
        self.text_preview = scrolledtext.ScrolledText(
            self.text_frame,
            wrap=tk.WORD,
            font=("Consolas", 11),
            bg="#11111b",
            fg="#cdd6f4",
            insertbackground="#cdd6f4",
            relief="flat",
            borderwidth=0,
            padx=10,
            pady=10
        )
        self.text_preview.pack(fill=tk.BOTH, expand=True)
        
        # Table preview
        self.table_frame = ttk.Frame(self.preview_frame, style="Dark.TFrame")
        
    def preview_file(self, file_path):
        """Preview the selected file"""
        if not os.path.exists(file_path):
            return
            
        self.filename_var.set(os.path.basename(file_path))
        ext = os.path.splitext(file_path)[1].lower()
        
        # Hide all preview widgets
        self.image_label.pack_forget()
        self.text_frame.pack_forget()
        self.table_frame.pack_forget()
        
        if ext in ['.png', '.jpg', '.jpeg', '.gif', '.bmp']:
            self._preview_image(file_path)
        elif ext in ['.tsv', '.csv']:
            self._preview_table(file_path)
        elif ext in ['.txt', '.log', '.json', '.yaml', '.yml', '.md', '.html']:
            self._preview_text(file_path)
        elif ext == '.pdf':
            self._show_placeholder("📊 PDF File\n\nDouble-click to open")
        else:
            self._show_file_info(file_path)
            
    def _preview_image(self, file_path):
        """Preview an image file"""
        if not HAS_PIL:
            self._show_placeholder("Install Pillow for image preview:\npip install Pillow")
            return
            
        try:
            img = Image.open(file_path)
            
            # Get frame size
            self.preview_frame.update_idletasks()
            max_w = max(400, self.preview_frame.winfo_width() - 40)
            max_h = max(400, self.preview_frame.winfo_height() - 40)
            
            # Resize maintaining aspect ratio
            img.thumbnail((max_w, max_h), Image.Resampling.LANCZOS)
            
            self.current_image = ImageTk.PhotoImage(img)
            self.image_label.configure(image=self.current_image, text="")
            self.image_label.pack(fill=tk.BOTH, expand=True)
            
        except Exception as e:
            self._show_placeholder(f"Could not load image:\n{e}")
            
    def _preview_table(self, file_path):
        """Preview a TSV/CSV file"""
        if not HAS_PANDAS:
            self._preview_text(file_path)
            return
            
        try:
            sep = '\t' if file_path.endswith('.tsv') else ','
            df = pd.read_csv(file_path, sep=sep, nrows=100)
            
            # Clear previous table
            for widget in self.table_frame.winfo_children():
                widget.destroy()
                
            # Create treeview for table
            tree = ttk.Treeview(self.table_frame, show="headings", height=20)
            
            # Configure columns
            cols = list(df.columns)[:10]  # Limit columns
            tree["columns"] = cols
            for col in cols:
                tree.heading(col, text=str(col)[:15])
                tree.column(col, width=100, minwidth=60)
                
            # Add data
            for _, row in df.head(50).iterrows():
                values = [str(v)[:20] for v in row[cols].values]
                tree.insert("", "end", values=values)
                
            # Scrollbars
            vsb = ttk.Scrollbar(self.table_frame, orient="vertical", command=tree.yview)
            hsb = ttk.Scrollbar(self.table_frame, orient="horizontal", command=tree.xview)
            tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
            
            tree.grid(row=0, column=0, sticky="nsew")
            vsb.grid(row=0, column=1, sticky="ns")
            hsb.grid(row=1, column=0, sticky="ew")
            
            self.table_frame.grid_rowconfigure(0, weight=1)
            self.table_frame.grid_columnconfigure(0, weight=1)
            self.table_frame.pack(fill=tk.BOTH, expand=True)
            
        except Exception as e:
            self._preview_text(file_path)
            
    def _preview_text(self, file_path):
        """Preview a text file"""
        try:
            with open(file_path, 'r', errors='replace') as f:
                content = f.read(50000)
                
            self.text_preview.delete(1.0, tk.END)
            self.text_preview.insert(1.0, content)
            self.text_frame.pack(fill=tk.BOTH, expand=True)
            
        except Exception as e:
            self._show_placeholder(f"Could not read file:\n{e}")
            
    def _show_placeholder(self, text):
        """Show placeholder text"""
        self.image_label.configure(text=text, image="")
        self.image_label.pack(fill=tk.BOTH, expand=True)
        
    def _show_file_info(self, file_path):
        """Show file info for unknown types"""
        size = os.path.getsize(file_path)
        if size > 1024*1024:
            size_str = f"{size/1024/1024:.1f} MB"
        elif size > 1024:
            size_str = f"{size/1024:.1f} KB"
        else:
            size_str = f"{size:,} bytes"
            
        self._show_placeholder(f"📄 {os.path.basename(file_path)}\n\nSize: {size_str}\n\nDouble-click to open")


class MedipipeMLGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("MEDIPIPE ML Discrimination Analysis")
        self.root.geometry("1700x1000")
        self.root.minsize(1400, 900)
        
        # Apply modern styling
        ModernStyle.configure_styles(root)
        
        # Configure root background
        self.root.configure(bg=ModernStyle.BG_DARK)
        
        # Configuration variables
        self.setup_variables()
        
        # Build UI
        self.setup_ui()
        
        # Log startup
        self.log_message("🚀 MEDIPIPE ML Analysis GUI initialized")
        
    def setup_variables(self):
        """Initialize all configuration variables"""
        self.matrix_path = tk.StringVar(value="")
        self.annotation_path = tk.StringVar(value="")
        self.output_dir = tk.StringVar(value="")
        
        self.case_group = tk.StringVar(value="AEG")
        self.ctrl_group = tk.StringVar(value="CTRL")
        self.feature_space = tk.StringVar(value="autosomal_windows")
        
        self.model_type = tk.StringVar(value="both")
        self.top_k = tk.IntVar(value=5000)
        self.C_param = tk.DoubleVar(value=1.0)
        self.alpha_param = tk.DoubleVar(value=0.01)
        self.l1_ratio = tk.DoubleVar(value=0.5)
        self.n_jobs = tk.IntVar(value=-1)
        
        self.top_n_heatmap = tk.IntVar(value=20)
        self.top_n_barplot = tk.IntVar(value=30)
        
        self.run_permutation = tk.BooleanVar(value=True)
        self.run_visualization = tk.BooleanVar(value=True)
        self.run_alphagenome = tk.BooleanVar(value=False)  # AlphaGenome analysis
        self.alphagenome_filter_intergenic = tk.BooleanVar(value=True)
        self.alphagenome_top_n = tk.IntVar(value=10)
        
        self.current_process = None
        self.is_running = tk.BooleanVar(value=False)
        
    def setup_ui(self):
        """Build the main user interface with 3-panel layout"""
        # Main container
        main = ttk.Frame(self.root, style="Dark.TFrame")
        main.pack(fill=tk.BOTH, expand=True)
        
        # Create PanedWindow for resizable panels
        self.paned = ttk.PanedWindow(main, orient=tk.HORIZONTAL)
        self.paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel - File tree
        left_frame = ttk.Frame(self.paned, style="Dark.TFrame", width=300)
        self.file_panel = FileTreePanel(left_frame, on_file_select=self.on_file_selected)
        self.file_panel.pack(fill=tk.BOTH, expand=True)
        self.paned.add(left_frame, weight=1)
        
        # Center panel - Configuration
        center_frame = ttk.Frame(self.paned, style="Dark.TFrame")
        self.notebook = ttk.Notebook(center_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.create_data_tab()
        self.create_ml_params_tab()
        self.create_run_tab()
        
        self.paned.add(center_frame, weight=3)
        
        # Right panel - Preview
        right_frame = ttk.Frame(self.paned, style="Dark.TFrame", width=400)
        self.preview_panel = PreviewPanel(right_frame)
        self.preview_panel.pack(fill=tk.BOTH, expand=True)
        self.paned.add(right_frame, weight=2)
        
    def create_data_tab(self):
        """Create Data Selection tab"""
        tab = ttk.Frame(self.notebook, padding=20)
        self.notebook.add(tab, text="  📂 Data  ")
        
        # Paths section
        paths = ttk.LabelFrame(tab, text="Input Files", padding=15)
        paths.pack(fill=tk.X, pady=10)
        
        for label_text, var, browse_cmd in [
            ("Feature Matrix:", self.matrix_path, self.browse_matrix),
            ("Sample Annotation:", self.annotation_path, self.browse_annotation),
            ("Output Directory:", self.output_dir, self.browse_output)
        ]:
            row = ttk.Frame(paths)
            row.pack(fill=tk.X, pady=8)
            ttk.Label(row, text=label_text, width=18).pack(side=tk.LEFT)
            ttk.Entry(row, textvariable=var).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)
            ttk.Button(row, text="Browse", command=browse_cmd).pack(side=tk.RIGHT)
        
        # Groups section
        groups = ttk.LabelFrame(tab, text="Sample Groups", padding=15)
        groups.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(groups)
        row.pack(fill=tk.X, pady=8)
        
        ttk.Label(row, text="Case Group:").pack(side=tk.LEFT)
        ttk.Entry(row, textvariable=self.case_group, width=12).pack(side=tk.LEFT, padx=(5, 30))
        
        ttk.Label(row, text="Control Group:").pack(side=tk.LEFT)
        ttk.Entry(row, textvariable=self.ctrl_group, width=12).pack(side=tk.LEFT, padx=(5, 30))
        
        ttk.Label(row, text="Feature Space:").pack(side=tk.LEFT)
        ttk.Combobox(row, textvariable=self.feature_space, width=20,
                    values=["autosomal_windows", "promoters", "gene_bodies",
                           "cpg_islands", "enhancers"]).pack(side=tk.LEFT, padx=5)
        
        # Quick actions
        actions = ttk.LabelFrame(tab, text="Quick Actions", padding=15)
        actions.pack(fill=tk.X, pady=10)
        
        btn_row = ttk.Frame(actions)
        btn_row.pack(fill=tk.X)
        
        ttk.Button(btn_row, text="🚀 Load AEG vs CTRL",
                   command=self.load_default).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_row, text="🔍 Find Matrices",
                   command=self.detect_matrices).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_row, text="📂 Browse Output",
                   command=self.browse_and_show_output).pack(side=tk.LEFT, padx=5)
        
    def create_ml_params_tab(self):
        """Create ML Parameters tab"""
        tab = ttk.Frame(self.notebook, padding=20)
        self.notebook.add(tab, text="  ⚙️ Parameters  ")
        
        # Model selection
        model = ttk.LabelFrame(tab, text="Model Selection", padding=15)
        model.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(model)
        row.pack(fill=tk.X, pady=5)
        for text, val in [("L2 (Ridge)", "L2"), ("Elastic-Net", "ElasticNet"), ("Both ✓", "both")]:
            ttk.Radiobutton(row, text=text, variable=self.model_type, value=val).pack(side=tk.LEFT, padx=15)
        
        # Regularization
        reg = ttk.LabelFrame(tab, text="Regularization", padding=15)
        reg.pack(fill=tk.X, pady=10)
        
        for label, var, hint in [
            ("L2 Regularization (C):", self.C_param, "↑ = less regularization"),
            ("ElasticNet Alpha:", self.alpha_param, "↓ = less regularization"),
            ("Top-K Features:", self.top_k, "ANOVA F-score selection")
        ]:
            row = ttk.Frame(reg)
            row.pack(fill=tk.X, pady=8)
            ttk.Label(row, text=label, width=20).pack(side=tk.LEFT)
            ttk.Entry(row, textvariable=var, width=10).pack(side=tk.LEFT, padx=10)
            ttk.Label(row, text=hint, style="Dim.TLabel").pack(side=tk.LEFT, padx=10)
        
        # Visualization
        viz = ttk.LabelFrame(tab, text="Visualization Options", padding=15)
        viz.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(viz)
        row.pack(fill=tk.X, pady=5)
        ttk.Label(row, text="Top N Heatmap:").pack(side=tk.LEFT)
        ttk.Spinbox(row, from_=5, to=100, textvariable=self.top_n_heatmap, width=8).pack(side=tk.LEFT, padx=(5, 30))
        ttk.Label(row, text="Top N Bar Plot:").pack(side=tk.LEFT)
        ttk.Spinbox(row, from_=10, to=100, textvariable=self.top_n_barplot, width=8).pack(side=tk.LEFT, padx=5)
        
        # Computation
        comp = ttk.LabelFrame(tab, text="Computation", padding=15)
        comp.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(comp)
        row.pack(fill=tk.X, pady=5)
        ttk.Label(row, text="Parallel Jobs:").pack(side=tk.LEFT)
        ttk.Spinbox(row, from_=-1, to=128, textvariable=self.n_jobs, width=8).pack(side=tk.LEFT, padx=10)
        ttk.Label(row, text="(-1 = all cores)", style="Dim.TLabel").pack(side=tk.LEFT)
        
        # Presets
        presets = ttk.LabelFrame(tab, text="Presets", padding=15)
        presets.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(presets)
        row.pack(fill=tk.X)
        for text, cmd in [("Default", self.preset_default), ("High Reg.", self.preset_high_reg),
                          ("Low Reg.", self.preset_low_reg), ("Sparse", self.preset_sparse)]:
            ttk.Button(row, text=text, command=cmd).pack(side=tk.LEFT, padx=5)
        
    def create_run_tab(self):
        """Create Run tab"""
        tab = ttk.Frame(self.notebook, padding=20)
        self.notebook.add(tab, text="  ▶️ Run  ")
        
        # Analysis selection
        sel = ttk.LabelFrame(tab, text="Analysis Selection", padding=15)
        sel.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(sel)
        row.pack(fill=tk.X)
        ttk.Checkbutton(row, text="Exact Permutation Test",
                        variable=self.run_permutation).pack(side=tk.LEFT, padx=15)
        ttk.Checkbutton(row, text="Visualizations (Volcano, Heatmap)",
                        variable=self.run_visualization).pack(side=tk.LEFT, padx=15)
        
        # AlphaGenome section
        ag = ttk.LabelFrame(tab, text="🧬 AlphaGenome Analysis (Advanced)", padding=15)
        ag.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(ag)
        row.pack(fill=tk.X, pady=5)
        ttk.Checkbutton(row, text="Run AlphaGenome on intergenic driver regions",
                        variable=self.run_alphagenome).pack(side=tk.LEFT)
        
        row2 = ttk.Frame(ag)
        row2.pack(fill=tk.X, pady=5)
        ttk.Checkbutton(row2, text="Filter to intergenic only",
                        variable=self.alphagenome_filter_intergenic).pack(side=tk.LEFT, padx=15)
        ttk.Label(row2, text="Top N regions:").pack(side=tk.LEFT, padx=(20, 5))
        ttk.Spinbox(row2, from_=1, to=50, textvariable=self.alphagenome_top_n, width=8).pack(side=tk.LEFT)
        
        ag_note = ttk.Label(ag, text="⚠️ Requires ALPHAGENOME_API_KEY environment variable", style="Dim.TLabel")
        ag_note.pack(anchor=tk.W, pady=(5, 0))
        
        # Controls
        ctrl = ttk.LabelFrame(tab, text="Controls", padding=15)
        ctrl.pack(fill=tk.X, pady=10)
        
        row = ttk.Frame(ctrl)
        row.pack(fill=tk.X)
        
        self.run_btn = ttk.Button(row, text="▶️ Run Analysis",
                                  command=self.run_analysis, style="Accent.TButton")
        self.run_btn.pack(side=tk.LEFT, padx=5)
        
        self.stop_btn = ttk.Button(row, text="⏹️ Stop",
                                   command=self.stop_analysis, state=tk.DISABLED)
        self.stop_btn.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(row, text="📋 Show Command", command=self.show_command).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="🗑️ Clear", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        
        # Progress
        prog = ttk.LabelFrame(tab, text="Progress", padding=15)
        prog.pack(fill=tk.X, pady=10)
        
        self.progress_var = tk.StringVar(value="Ready")
        ttk.Label(prog, textvariable=self.progress_var).pack(anchor=tk.W)
        self.progress_bar = ttk.Progressbar(prog, mode='indeterminate')
        self.progress_bar.pack(fill=tk.X, pady=5)
        
        # Log
        log = ttk.LabelFrame(tab, text="Log", padding=10)
        log.pack(fill=tk.BOTH, expand=True, pady=10)
        
        self.log_text = scrolledtext.ScrolledText(
            log, height=12, wrap=tk.WORD,
            font=("Consolas", 11),
            bg="#11111b", fg="#cdd6f4",
            insertbackground="#cdd6f4",
            relief="flat", borderwidth=0
        )
        self.log_text.pack(fill=tk.BOTH, expand=True)
        
    # Event handlers
    def on_file_selected(self, file_path):
        self.preview_panel.preview_file(file_path)
        
    def browse_and_show_output(self):
        if self.output_dir.get():
            self.file_panel.set_directory(self.output_dir.get())
        else:
            self.browse_output()
            if self.output_dir.get():
                self.file_panel.set_directory(self.output_dir.get())
        
    # Browse functions
    def _get_initial_dir(self, dirs):
        for d in dirs:
            if os.path.exists(d):
                return d
        return os.path.expanduser("~")
    
    def browse_matrix(self):
        f = filedialog.askopenfilename(
            title="Select Feature Matrix",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")],
            initialdir=self._get_initial_dir([
                "/data/medipipe_data/output/ml_discrimination/matrices",
                os.path.expanduser("~")])
        )
        if f:
            self.matrix_path.set(f)
            self.log_message(f"📊 Matrix: {f}")
            
    def browse_annotation(self):
        f = filedialog.askopenfilename(
            title="Select Sample Annotation",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")],
            initialdir=self._get_initial_dir([
                "/home/dirk/medipipe_warp/samplefiles",
                os.path.expanduser("~")])
        )
        if f:
            self.annotation_path.set(f)
            self.log_message(f"📋 Annotation: {f}")
            
    def browse_output(self):
        d = filedialog.askdirectory(
            title="Select Output Directory",
            initialdir=self._get_initial_dir([
                "/data/medipipe_data/output/ml_discrimination/Outputs",
                os.path.expanduser("~")])
        )
        if d:
            self.output_dir.set(d)
            self.log_message(f"📁 Output: {d}")
            self.file_panel.set_directory(d)
    
    # Quick actions
    def load_default(self):
        self.matrix_path.set("/data/medipipe_data/output/ml_discrimination/matrices/windows/hg38_w2000/all/matrix.tsv")
        self.annotation_path.set("/home/dirk/medipipe_warp/samplefiles/current_samples.tsv")
        self.output_dir.set("/data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL")
        self.case_group.set("AEG")
        self.ctrl_group.set("CTRL")
        self.feature_space.set("autosomal_windows")
        self.log_message("✅ Loaded AEG vs CTRL configuration")
        if os.path.exists(self.output_dir.get()):
            self.file_panel.set_directory(self.output_dir.get())
        
    def detect_matrices(self):
        base = Path("/data/medipipe_data/output/ml_discrimination/matrices")
        if base.exists():
            matrices = list(base.glob("**/*.tsv"))
            self.log_message(f"Found {len(matrices)} matrices")
            for m in matrices[:5]:
                self.log_message(f"  📊 {m}")
        else:
            self.log_message("❌ Matrix directory not found")
    
    # Presets
    def preset_default(self):
        self.C_param.set(1.0); self.alpha_param.set(0.01); self.top_k.set(5000)
        self.log_message("✅ Default preset")
        
    def preset_high_reg(self):
        self.C_param.set(0.1); self.alpha_param.set(0.1); self.top_k.set(1000)
        self.log_message("✅ High regularization preset")
        
    def preset_low_reg(self):
        self.C_param.set(10.0); self.alpha_param.set(0.001); self.top_k.set(10000)
        self.log_message("✅ Low regularization preset")
        
    def preset_sparse(self):
        self.C_param.set(1.0); self.alpha_param.set(0.05); self.top_k.set(2000)
        self.log_message("✅ Sparse preset")
    
    # Commands
    def build_permutation_cmd(self):
        script = "/home/dirk/medipipe_warp/scripts/ml/ml_exact_permutation.py"
        return ["python3", script,
                "--matrix", self.matrix_path.get(),
                "--annotation", self.annotation_path.get(),
                "--output", os.path.join(self.output_dir.get(), "exact_permutation"),
                "--feature-space", self.feature_space.get(),
                "--case-group", self.case_group.get(),
                "--ctrl-group", self.ctrl_group.get(),
                "--top-k", str(self.top_k.get()),
                "--C", str(self.C_param.get()),
                "--alpha", str(self.alpha_param.get()),
                "--model-type", self.model_type.get(),
                "--n-jobs", str(self.n_jobs.get())]
    
    def build_viz_cmd(self):
        script = "/home/dirk/medipipe_warp/scripts/ml/ml_visualization.py"
        return ["python3", script,
                "--matrix", self.matrix_path.get(),
                "--annotation", self.annotation_path.get(),
                "--output", os.path.join(self.output_dir.get(), "visualizations"),
                "--feature-space", self.feature_space.get(),
                "--case-group", self.case_group.get(),
                "--ctrl-group", self.ctrl_group.get(),
                "--top-k", str(self.top_k.get()),
                "--C", str(self.C_param.get()),
                "--alpha", str(self.alpha_param.get()),
                "--model-type", self.model_type.get(),
                "--top-n-heatmap", str(self.top_n_heatmap.get()),
                "--top-n-barplot", str(self.top_n_barplot.get())]
    
    def build_alphagenome_cmd(self):
        """Build command to run AlphaGenome on intergenic driver regions"""
        script = "/home/dirk/medipipe_warp/scripts/ml/run_alphagenome_drivers.py"
        drivers_file = os.path.join(self.output_dir.get(), "visualizations", "top10_L2_drivers_table.tsv")
        ag_output = os.path.join(self.output_dir.get(), "alphagenome")
        
        cmd = ["python3", script,
               "--drivers", drivers_file,
               "--output-dir", ag_output,
               "--top-n", str(self.alphagenome_top_n.get())]
        
        if self.alphagenome_filter_intergenic.get():
            cmd.append("--filter-intergenic")
        
        return cmd
    
    def show_command(self):
        self.log_message("\n" + "="*50)
        if self.run_permutation.get():
            self.log_message("Permutation: " + " ".join(self.build_permutation_cmd()))
        if self.run_visualization.get():
            self.log_message("Visualization: " + " ".join(self.build_viz_cmd()))
        if self.run_alphagenome.get():
            self.log_message("AlphaGenome: " + " ".join(self.build_alphagenome_cmd()))
        self.log_message("="*50)
    
    # Run
    def validate(self):
        errors = []
        if not self.matrix_path.get() or not os.path.exists(self.matrix_path.get()):
            errors.append("Matrix file missing")
        if not self.annotation_path.get() or not os.path.exists(self.annotation_path.get()):
            errors.append("Annotation file missing")
        if not self.output_dir.get():
            errors.append("Output directory not set")
        if not self.run_permutation.get() and not self.run_visualization.get() and not self.run_alphagenome.get():
            errors.append("No analysis selected")
        if self.run_alphagenome.get() and not os.environ.get("ALPHAGENOME_API_KEY"):
            errors.append("ALPHAGENOME_API_KEY environment variable not set")
        return errors
    
    def run_analysis(self):
        errors = self.validate()
        if errors:
            messagebox.showerror("Error", "\n".join(errors))
            return
        
        os.makedirs(self.output_dir.get(), exist_ok=True)
        
        self.is_running.set(True)
        self.run_btn.config(state=tk.DISABLED)
        self.stop_btn.config(state=tk.NORMAL)
        self.progress_bar.start()
        
        threading.Thread(target=self._run_thread, daemon=True).start()
    
    def _run_thread(self):
        try:
            if self.run_permutation.get():
                self.progress_var.set("Running permutation test...")
                self.log_message("\n🚀 STARTING PERMUTATION TEST")
                cmd = self.build_permutation_cmd()
                self.current_process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                                       stderr=subprocess.STDOUT, text=True)
                for line in self.current_process.stdout:
                    self.log_message(line.rstrip())
                    if not self.is_running.get():
                        break
                self.current_process.wait()
                self.log_message("✅ Permutation test done" if self.current_process.returncode == 0 else "❌ Failed")
            
            if self.run_visualization.get() and self.is_running.get():
                self.progress_var.set("Generating visualizations...")
                self.log_message("\n🎨 GENERATING VISUALIZATIONS")
                cmd = self.build_viz_cmd()
                self.current_process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                                       stderr=subprocess.STDOUT, text=True)
                for line in self.current_process.stdout:
                    self.log_message(line.rstrip())
                    if not self.is_running.get():
                        break
                self.current_process.wait()
                self.log_message("✅ Visualizations done" if self.current_process.returncode == 0 else "❌ Failed")
            
            if self.run_alphagenome.get() and self.is_running.get():
                self.progress_var.set("Running AlphaGenome analysis...")
                self.log_message("\n🧬 RUNNING ALPHAGENOME ANALYSIS")
                self.log_message("Analyzing intergenic driver regions for regulatory activity...")
                cmd = self.build_alphagenome_cmd()
                self.current_process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                                       stderr=subprocess.STDOUT, text=True)
                for line in self.current_process.stdout:
                    self.log_message(line.rstrip())
                    if not self.is_running.get():
                        break
                self.current_process.wait()
                if self.current_process.returncode == 0:
                    self.log_message("✅ AlphaGenome analysis done")
                    self.log_message("   Results: " + os.path.join(self.output_dir.get(), "alphagenome", "Outputs", "tables"))
                else:
                    self.log_message("❌ AlphaGenome failed - check ALPHAGENOME_API_KEY")
            
            self.log_message("\n🎉 COMPLETE")
            self.root.after(0, self.file_panel.refresh)
            
        except Exception as e:
            self.log_message(f"❌ Error: {e}")
        finally:
            self.root.after(0, self._reset_ui)
    
    def stop_analysis(self):
        self.is_running.set(False)
        if self.current_process:
            self.current_process.terminate()
            self.log_message("⏹️ Stopped")
        self._reset_ui()
    
    def _reset_ui(self):
        self.is_running.set(False)
        self.run_btn.config(state=tk.NORMAL)
        self.stop_btn.config(state=tk.DISABLED)
        self.progress_bar.stop()
        self.progress_var.set("Ready")
        self.current_process = None
    
    # Logging
    def log_message(self, msg):
        ts = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(tk.END, f"[{ts}] {msg}\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()
    
    def clear_log(self):
        self.log_text.delete(1.0, tk.END)


def main():
    root = tk.Tk()
    app = MedipipeMLGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
