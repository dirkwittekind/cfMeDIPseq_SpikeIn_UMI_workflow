#!/usr/bin/env python3
"""
MEDIPIPE GUI - Main Entry Point
cfMeDIP-seq Analysis Pipeline with GUI
"""

import os
import sys

# Add the parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gui.app import MedipipeApp


def main():
    """Main entry point for MEDIPIPE GUI"""
    app = MedipipeApp()
    app.run()


if __name__ == "__main__":
    main()
