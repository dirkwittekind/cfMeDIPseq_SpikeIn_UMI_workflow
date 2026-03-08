"""
MEDIPIPE File Browser Component
Word-like file browser for selecting FASTQ files.
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog
from tkinter import font as tkfont
from typing import List, Callable, Optional


class FileBrowser(ttk.Frame):
    """
    A file browser component with folder navigation and file selection.
    Styled like a Word/Office file browser dialog.
    """
    
    FASTQ_EXTENSIONS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    
    def __init__(
        self,
        parent,
        initial_dir: str = "/data/raw_data",
        on_selection_change: Optional[Callable[[List[str]], None]] = None,
        file_filter: str = "fastq",
        multi_select: bool = True
    ):
        super().__init__(parent)
        
        self.current_dir = initial_dir
        self.on_selection_change = on_selection_change
        self.file_filter = file_filter
        self.multi_select = multi_select
        self.selected_files: List[str] = []
        
        self._create_widgets()
        self._refresh_view()
    
    def _create_widgets(self):
        """Create the file browser widgets"""
        # Configure font for Treeviews (ensures readable font size on Linux)
        self.tree_font = tkfont.Font(family="Helvetica", size=11)
        self.tree_heading_font = tkfont.Font(family="Helvetica", size=11, weight="bold")
        
        # Top toolbar with navigation
        toolbar = ttk.Frame(self)
        toolbar.pack(fill=tk.X, padx=5, pady=5)
        
        # Back button
        self.back_btn = ttk.Button(toolbar, text="⬆ Up", width=8, command=self._go_up)
        self.back_btn.pack(side=tk.LEFT, padx=2)
        
        # Path entry
        self.path_var = tk.StringVar(value=self.current_dir)
        self.path_entry = ttk.Entry(toolbar, textvariable=self.path_var, width=50)
        self.path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        self.path_entry.bind("<Return>", lambda e: self._navigate_to_path())
        
        # Go button
        ttk.Button(toolbar, text="Go", width=6, command=self._navigate_to_path).pack(side=tk.LEFT, padx=2)
        
        # Browse button
        ttk.Button(toolbar, text="Browse...", width=10, command=self._browse_dialog).pack(side=tk.LEFT, padx=2)
        
        # Main content frame with folder tree and file list
        content = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        content.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel: Folder tree
        tree_frame = ttk.LabelFrame(content, text="Folders", padding=5)
        
        self.folder_tree = ttk.Treeview(tree_frame, selectmode="browse", show="tree")
        self.folder_tree.tag_configure('default', font=self.tree_font)
        folder_scroll = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=self.folder_tree.yview)
        self.folder_tree.configure(yscrollcommand=folder_scroll.set)
        
        self.folder_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        folder_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.folder_tree.bind("<<TreeviewSelect>>", self._on_folder_select)
        
        content.add(tree_frame, weight=1)
        
        # Right panel: File list
        file_frame = ttk.LabelFrame(content, text="Files", padding=5)
        
        # Column headers
        columns = ("name", "size", "modified")
        self.file_list = ttk.Treeview(file_frame, columns=columns, show="headings",
                                      selectmode="extended" if self.multi_select else "browse")
        self.file_list.tag_configure('default', font=self.tree_font)
        
        self.file_list.heading("name", text="Name", anchor=tk.W)
        self.file_list.heading("size", text="Size", anchor=tk.E)
        self.file_list.heading("modified", text="Modified", anchor=tk.W)
        
        self.file_list.column("name", width=350, minwidth=200)
        self.file_list.column("size", width=100, minwidth=80)
        self.file_list.column("modified", width=180, minwidth=120)
        
        file_scroll_y = ttk.Scrollbar(file_frame, orient=tk.VERTICAL, command=self.file_list.yview)
        file_scroll_x = ttk.Scrollbar(file_frame, orient=tk.HORIZONTAL, command=self.file_list.xview)
        self.file_list.configure(yscrollcommand=file_scroll_y.set, xscrollcommand=file_scroll_x.set)
        
        self.file_list.grid(row=0, column=0, sticky="nsew")
        file_scroll_y.grid(row=0, column=1, sticky="ns")
        file_scroll_x.grid(row=1, column=0, sticky="ew")
        
        file_frame.columnconfigure(0, weight=1)
        file_frame.rowconfigure(0, weight=1)
        
        self.file_list.bind("<<TreeviewSelect>>", self._on_file_select)
        self.file_list.bind("<Double-1>", self._on_double_click)
        
        content.add(file_frame, weight=3)
        
        # Bottom status bar
        status_frame = ttk.Frame(self)
        status_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(status_frame, textvariable=self.status_var).pack(side=tk.LEFT)
        
        self.selection_var = tk.StringVar(value="0 files selected")
        ttk.Label(status_frame, textvariable=self.selection_var).pack(side=tk.RIGHT)
    
    def _refresh_view(self):
        """Refresh the folder tree and file list"""
        self._populate_folder_tree()
        self._populate_file_list()
    
    def _populate_folder_tree(self):
        """Populate the folder tree view"""
        # Clear existing items
        for item in self.folder_tree.get_children():
            self.folder_tree.delete(item)
        
        if not os.path.exists(self.current_dir):
            return
        
        # Add common locations
        self.folder_tree.insert("", "end", "raw_data", text="📁 Raw Data (/data/raw_data)", tags=('default',))
        self.folder_tree.insert("", "end", "home", text="📁 Home", tags=('default',))
        self.folder_tree.insert("", "end", "current", text="📁 Current Directory", tags=('default',))
        
        # Add subfolders of current directory
        try:
            for item in sorted(os.listdir(self.current_dir)):
                item_path = os.path.join(self.current_dir, item)
                if os.path.isdir(item_path):
                    self.folder_tree.insert("", "end", item_path, text=f"  📂 {item}", tags=('default',))
        except PermissionError:
            pass
    
    def _populate_file_list(self):
        """Populate the file list view"""
        # Clear existing items
        for item in self.file_list.get_children():
            self.file_list.delete(item)
        
        if not os.path.exists(self.current_dir):
            self.status_var.set(f"Directory not found: {self.current_dir}")
            return
        
        file_count = 0
        dir_count = 0
        
        try:
            items = sorted(os.listdir(self.current_dir))
            
            # Add directories first
            for item in items:
                item_path = os.path.join(self.current_dir, item)
                if os.path.isdir(item_path):
                    dir_count += 1
                    self.file_list.insert("", "end", item_path, values=(f"📂 {item}", "<DIR>", ""), tags=('default',))
            
            # Then add files
            for item in items:
                item_path = os.path.join(self.current_dir, item)
                if os.path.isfile(item_path):
                    # Apply filter
                    if self.file_filter == "fastq":
                        if not item.lower().endswith(self.FASTQ_EXTENSIONS):
                            continue
                    
                    file_count += 1
                    size = self._format_size(os.path.getsize(item_path))
                    mtime = self._format_time(os.path.getmtime(item_path))
                    self.file_list.insert("", "end", item_path, values=(f"📄 {item}", size, mtime), tags=('default',))
            
            self.status_var.set(f"{dir_count} folders, {file_count} files")
        
        except PermissionError:
            self.status_var.set("Permission denied")
    
    def _format_size(self, size: int) -> str:
        """Format file size in human-readable format"""
        for unit in ["B", "KB", "MB", "GB"]:
            if size < 1024:
                return f"{size:.1f} {unit}"
            size /= 1024
        return f"{size:.1f} TB"
    
    def _format_time(self, timestamp: float) -> str:
        """Format timestamp"""
        from datetime import datetime
        return datetime.fromtimestamp(timestamp).strftime("%Y-%m-%d %H:%M")
    
    def _go_up(self):
        """Navigate to parent directory"""
        parent = os.path.dirname(self.current_dir)
        if parent and parent != self.current_dir:
            self.current_dir = parent
            self.path_var.set(self.current_dir)
            self._refresh_view()
    
    def _navigate_to_path(self):
        """Navigate to the path in the entry field"""
        path = self.path_var.get().strip()
        if os.path.isdir(path):
            self.current_dir = path
            self._refresh_view()
        else:
            self.status_var.set(f"Invalid directory: {path}")
    
    def _browse_dialog(self):
        """Open a folder selection dialog"""
        path = filedialog.askdirectory(initialdir=self.current_dir, title="Select Directory")
        if path:
            self.current_dir = path
            self.path_var.set(path)
            self._refresh_view()
    
    def _on_folder_select(self, event):
        """Handle folder tree selection"""
        selection = self.folder_tree.selection()
        if not selection:
            return
        
        item_id = selection[0]
        
        # Handle special locations
        if item_id == "raw_data":
            self.current_dir = "/data/raw_data"
        elif item_id == "home":
            self.current_dir = os.path.expanduser("~")
        elif item_id == "current":
            self.current_dir = os.getcwd()
        elif os.path.isdir(item_id):
            self.current_dir = item_id
        
        self.path_var.set(self.current_dir)
        self._refresh_view()
    
    def _on_file_select(self, event):
        """Handle file list selection"""
        selection = self.file_list.selection()
        self.selected_files = []
        
        for item in selection:
            if os.path.isfile(item):
                self.selected_files.append(item)
        
        self.selection_var.set(f"{len(self.selected_files)} files selected")
        
        if self.on_selection_change:
            self.on_selection_change(self.selected_files)
    
    def _on_double_click(self, event):
        """Handle double-click on item"""
        selection = self.file_list.selection()
        if not selection:
            return
        
        item_path = selection[0]
        if os.path.isdir(item_path):
            self.current_dir = item_path
            self.path_var.set(item_path)
            self._refresh_view()
    
    def get_selected_files(self) -> List[str]:
        """Get list of selected file paths"""
        return self.selected_files
    
    def set_directory(self, path: str):
        """Set the current directory"""
        if os.path.isdir(path):
            self.current_dir = path
            self.path_var.set(path)
            self._refresh_view()
    
    def clear_selection(self):
        """Clear the current selection"""
        self.file_list.selection_remove(*self.file_list.selection())
        self.selected_files = []
        self.selection_var.set("0 files selected")


class FileSelectionDialog(tk.Toplevel):
    """
    A popup dialog for file selection.
    Similar to File > Open in Word.
    """
    
    def __init__(
        self,
        parent,
        title: str = "Select Files",
        initial_dir: str = "/data/raw_data",
        file_filter: str = "fastq",
        multi_select: bool = True
    ):
        super().__init__(parent)
        self.title(title)
        self.geometry("900x600")
        self.transient(parent)
        self.grab_set()
        
        self.selected_files: List[str] = []
        self.result: Optional[List[str]] = None
        
        # File browser
        self.browser = FileBrowser(
            self,
            initial_dir=initial_dir,
            file_filter=file_filter,
            multi_select=multi_select
        )
        self.browser.pack(fill=tk.BOTH, expand=True)
        
        # Button frame
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Button(btn_frame, text="Cancel", command=self._cancel, width=12).pack(side=tk.RIGHT, padx=5)
        ttk.Button(btn_frame, text="Select", command=self._select, width=12).pack(side=tk.RIGHT, padx=5)
        
        # Center the dialog
        self.update_idletasks()
        x = (self.winfo_screenwidth() - self.winfo_width()) // 2
        y = (self.winfo_screenheight() - self.winfo_height()) // 2
        self.geometry(f"+{x}+{y}")
    
    def _select(self):
        """Handle Select button click"""
        self.result = self.browser.get_selected_files()
        self.destroy()
    
    def _cancel(self):
        """Handle Cancel button click"""
        self.result = None
        self.destroy()
    
    def show(self) -> Optional[List[str]]:
        """Show the dialog and return the result"""
        self.wait_window()
        return self.result
