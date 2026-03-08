#!/bin/bash
# MEDIPIPE Launcher Script
# Launches the MEDIPIPE GUI application

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    zenity --error --title="MEDIPIPE Error" --text="Python 3 is not installed or not in PATH."
    exit 1
fi

# Check for required Python packages
python3 -c "import tkinter" 2>/dev/null || {
    zenity --error --title="MEDIPIPE Error" --text="Python tkinter is not installed.\n\nInstall with: sudo apt install python3-tk"
    exit 1
}

python3 -c "import yaml" 2>/dev/null || {
    zenity --error --title="MEDIPIPE Error" --text="PyYAML is not installed.\n\nInstall with: pip install pyyaml"
    exit 1
}

# GTK settings for better font rendering in file dialogs
export GDK_SCALE=1
export GDK_DPI_SCALE=1.2

# Set GTK theme settings for larger fonts in dialogs
export GTK_THEME=Adwaita
if [ -f ~/.config/gtk-3.0/settings.ini ]; then
    : # Use existing settings
else
    mkdir -p ~/.config/gtk-3.0
    cat > ~/.config/gtk-3.0/settings.ini << 'EOF'
[Settings]
gtk-font-name=Sans 12
gtk-application-prefer-dark-theme=0
EOF
fi

# Launch the application
exec python3 "$SCRIPT_DIR/gui/main.py" "$@"
