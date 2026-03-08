#!/bin/bash
# =============================================================================
# cfMeDIP-seq Workflow - GitHub Sync Script
# =============================================================================
# This script automatically syncs workflow changes to GitHub
# Usage: ./sync_to_github.sh [commit message]
# =============================================================================

set -e

# Configuration
REPO_URL="https://github.com/dirkwittekind/cfMeDIPseq_SpikeIn_UMI_workflow.git"
BRANCH="main"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo -e "${GREEN}=== cfMeDIP-seq Workflow GitHub Sync ===${NC}"
echo ""

# Check if git is initialized
if [ ! -d ".git" ]; then
    echo -e "${YELLOW}Initializing git repository...${NC}"
    git init
    git remote add origin "$REPO_URL" 2>/dev/null || git remote set-url origin "$REPO_URL"
    git branch -M main
fi

# Ensure remote is set correctly
git remote set-url origin "$REPO_URL" 2>/dev/null || git remote add origin "$REPO_URL"

# Get commit message
if [ -n "$1" ]; then
    COMMIT_MSG="$1"
else
    COMMIT_MSG="Update workflow - $(date '+%Y-%m-%d %H:%M:%S')"
fi

# Files to include in the repository
echo -e "${YELLOW}Adding workflow files...${NC}"

# Add core workflow files
git add -f snakefiles/*.smk 2>/dev/null || true
git add -f rules/**/*.smk 2>/dev/null || true
git add -f scripts/**/*.py 2>/dev/null || true
git add -f scripts/**/*.R 2>/dev/null || true
git add -f scripts/**/*.sh 2>/dev/null || true
git add -f environments/*.yaml 2>/dev/null || true
git add -f configfiles/*.yaml 2>/dev/null || true
git add -f samplefiles/*_template.tsv 2>/dev/null || true

# Add documentation
git add -f README.md 2>/dev/null || true
git add -f LICENSE 2>/dev/null || true
git add -f .gitignore 2>/dev/null || true
git add -f docs/*.md 2>/dev/null || true

# Add GitHub workflows
git add -f .github/workflows/*.yml 2>/dev/null || true

# Add this sync script
git add -f sync_to_github.sh 2>/dev/null || true

# Show what's staged
echo ""
echo -e "${YELLOW}Staged files:${NC}"
git diff --cached --name-status

# Check if there are changes to commit
if git diff --cached --quiet; then
    echo ""
    echo -e "${GREEN}No changes to commit. Repository is up to date.${NC}"
    exit 0
fi

# Commit changes
echo ""
echo -e "${YELLOW}Committing changes...${NC}"
git commit -m "$COMMIT_MSG

Co-Authored-By: Oz <oz-agent@warp.dev>"

# Push to GitHub
echo ""
echo -e "${YELLOW}Pushing to GitHub...${NC}"
git push -u origin "$BRANCH" --force-with-lease 2>/dev/null || git push -u origin "$BRANCH"

echo ""
echo -e "${GREEN}✓ Successfully synced to GitHub!${NC}"
echo -e "  Repository: $REPO_URL"
echo -e "  Branch: $BRANCH"
echo ""
