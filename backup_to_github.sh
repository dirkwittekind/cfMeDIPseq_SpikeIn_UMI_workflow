#!/bin/bash
# =============================================================================
# cfMeDIP-seq Workflow - Automatic Private Backup Script
# =============================================================================
# This script automatically backs up the entire workflow to a private GitHub
# repository. Runs via cron for regular backups.
# 
# Usage: ./backup_to_github.sh [commit message]
# Cron:  0 */6 * * * /home/dirk/medipipe_warp/backup_to_github.sh "Automatic backup"
# =============================================================================

set -e

# Configuration
BACKUP_REPO_URL="https://github.com/dirkwittekind/cfMeDIPseq_backup_private.git"
BRANCH="main"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BACKUP_DIR="/tmp/medipipe_backup_$$"
LOG_FILE="/home/dirk/medipipe_warp/work/logs/backup.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Logging function
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "[$timestamp] $1" | tee -a "$LOG_FILE"
}

# Create log directory if needed
mkdir -p "$(dirname "$LOG_FILE")"

log "${GREEN}=== cfMeDIP-seq Private Backup ===${NC}"

# Get commit message
if [ -n "$1" ]; then
    COMMIT_MSG="$1"
else
    COMMIT_MSG="Automatic backup - $(date '+%Y-%m-%d %H:%M:%S')"
fi

# Create temporary backup directory
log "${YELLOW}Creating backup directory...${NC}"
rm -rf "$BACKUP_DIR"
mkdir -p "$BACKUP_DIR"

# Clone existing repo or init new one
cd "$BACKUP_DIR"
if git clone "$BACKUP_REPO_URL" repo 2>/dev/null; then
    cd repo
    log "Cloned existing repository"
else
    mkdir repo && cd repo
    git init
    git remote add origin "$BACKUP_REPO_URL"
    git branch -M main
    log "Initialized new repository"
fi

# Configure git
git config user.email "dirkwittekind@github.com"
git config user.name "Dirk Wittekind"

# Copy files to backup (excluding large data files)
log "${YELLOW}Copying workflow files...${NC}"

# Core workflow files
rsync -av --exclude='.git' \
    --exclude='.snakemake' \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='*.bam' \
    --exclude='*.bam.bai' \
    --exclude='*.fastq*' \
    --exclude='*.fq*' \
    --exclude='*.fa' \
    --exclude='*.fasta' \
    --exclude='*.bed' \
    --exclude='*.bigwig' \
    --exclude='*.bw' \
    --exclude='*.RData' \
    --exclude='*.rds' \
    --exclude='work/trimmed' \
    --exclude='work/aligned' \
    --exclude='work/sorted' \
    --exclude='work/umi_extracted' \
    --exclude='work/spikein_bsgenome' \
    --exclude='resources/genomes' \
    --exclude='resources/bsgenome' \
    --exclude='outputs' \
    --exclude='*.log' \
    --exclude='*.tmp' \
    --exclude='*.zip' \
    "$SCRIPT_DIR/" ./ 2>/dev/null || true

# Copy specific config files (with sensitive paths)
cp -f "$SCRIPT_DIR/work/config.yaml" ./work/ 2>/dev/null || true

# Copy current samples file
cp -f "$SCRIPT_DIR/samplefiles/current_samples.tsv" ./samplefiles/ 2>/dev/null || true

# Create backup metadata
cat > BACKUP_INFO.md << EOF
# cfMeDIP-seq Workflow Backup

**Backup Date:** $(date '+%Y-%m-%d %H:%M:%S')
**Source:** $SCRIPT_DIR
**Branch:** $BRANCH

## Contents
- Snakemake workflows and rules
- Analysis scripts (Python, R, shell)
- Conda environment definitions
- Configuration files
- GUI application
- MCP server documentation

## Excluded
- Raw data files (FASTQ, BAM)
- Reference genomes
- Large output files
- Temporary files
- Log files

## Restore
\`\`\`bash
git clone https://github.com/dirkwittekind/cfMeDIPseq_backup_private.git
# Then restore resources/genomes from original location
\`\`\`
EOF

# Stage all files
log "${YELLOW}Staging files...${NC}"
git add -A

# Check if there are changes
if git diff --cached --quiet; then
    log "${GREEN}No changes to backup. Repository is up to date.${NC}"
    rm -rf "$BACKUP_DIR"
    exit 0
fi

# Commit and push
log "${YELLOW}Committing changes...${NC}"
git commit -m "$COMMIT_MSG

Backup includes:
- Snakefiles and rules
- Scripts (Python, R, shell)
- Environments
- Configurations
- GUI and MCP

Co-Authored-By: Oz <oz-agent@warp.dev>"

log "${YELLOW}Pushing to private repository...${NC}"
git push -u origin "$BRANCH" --force 2>&1 | tee -a "$LOG_FILE"

# Cleanup
rm -rf "$BACKUP_DIR"

log "${GREEN}✓ Backup completed successfully!${NC}"
log "  Repository: $BACKUP_REPO_URL"
log "  Branch: $BRANCH"
