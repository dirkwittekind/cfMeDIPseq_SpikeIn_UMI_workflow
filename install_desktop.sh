#!/bin/bash
# MEDIPIPE Desktop Installer
# Installs the MEDIPIPE desktop launcher for the current user

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DESKTOP_FILE="$SCRIPT_DIR/medipipe.desktop"
ICON_SVG="$SCRIPT_DIR/resources/icons/medipipe_icon.svg"
ICON_PNG="$SCRIPT_DIR/resources/icons/medipipe_icon.png"
LAUNCH_SCRIPT="$SCRIPT_DIR/launch_medipipe.sh"

echo "=== MEDIPIPE Desktop Installer ==="
echo ""

# Make launch script executable
chmod +x "$LAUNCH_SCRIPT"
echo "✓ Made launch script executable"

# Convert SVG to PNG if needed
if [ ! -f "$ICON_PNG" ]; then
    echo "Converting icon to PNG..."
    if command -v rsvg-convert &> /dev/null; then
        rsvg-convert -w 256 -h 256 "$ICON_SVG" -o "$ICON_PNG"
        echo "✓ Icon converted using rsvg-convert"
    elif command -v convert &> /dev/null; then
        convert -background none -density 256 "$ICON_SVG" -resize 256x256 "$ICON_PNG"
        echo "✓ Icon converted using ImageMagick"
    elif command -v inkscape &> /dev/null; then
        inkscape --export-type=png --export-filename="$ICON_PNG" -w 256 -h 256 "$ICON_SVG"
        echo "✓ Icon converted using Inkscape"
    else
        echo "⚠ No SVG converter found. Install one of: librsvg2-bin, imagemagick, or inkscape"
        echo "  Using SVG icon directly (may not display in all file managers)"
        # Update desktop file to use SVG
        sed -i "s|Icon=.*|Icon=$ICON_SVG|" "$DESKTOP_FILE"
    fi
fi

# Update paths in desktop file
sed -i "s|Exec=.*|Exec=$LAUNCH_SCRIPT|" "$DESKTOP_FILE"
if [ -f "$ICON_PNG" ]; then
    sed -i "s|Icon=.*|Icon=$ICON_PNG|" "$DESKTOP_FILE"
fi

# Install to user's applications directory
APPLICATIONS_DIR="$HOME/.local/share/applications"
mkdir -p "$APPLICATIONS_DIR"
cp "$DESKTOP_FILE" "$APPLICATIONS_DIR/medipipe.desktop"
echo "✓ Installed to $APPLICATIONS_DIR/medipipe.desktop"

# Also copy to Desktop if it exists
if [ -d "$HOME/Desktop" ]; then
    cp "$DESKTOP_FILE" "$HOME/Desktop/medipipe.desktop"
    chmod +x "$HOME/Desktop/medipipe.desktop"
    # Mark as trusted (for GNOME)
    gio set "$HOME/Desktop/medipipe.desktop" metadata::trusted true 2>/dev/null || true
    echo "✓ Created desktop shortcut at $HOME/Desktop/medipipe.desktop"
fi

# Update desktop database
if command -v update-desktop-database &> /dev/null; then
    update-desktop-database "$APPLICATIONS_DIR" 2>/dev/null || true
    echo "✓ Updated desktop database"
fi

echo ""
echo "=== Installation Complete ==="
echo ""
echo "You can now launch MEDIPIPE by:"
echo "  1. Clicking the desktop icon (if created)"
echo "  2. Searching for 'MEDIPIPE' in your application menu"
echo "  3. Running: $LAUNCH_SCRIPT"
echo ""
