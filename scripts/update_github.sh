#!/bin/bash
# Update GitHub repository with new documentation
# Run from: ~/Documents/Github/noise_decorrelation_HIV/

set -e

echo "📚 Updating GitHub repository with documentation..."
echo ""

# Check we're in the right directory
if [ ! -d ".git" ]; then
    echo "❌ Error: Not in a git repository!"
    echo "Run this from: ~/Documents/Github/noise_decorrelation_HIV/"
    exit 1
fi

# Show what we're about to add
echo "Files to add:"
ls -1 README.md PROJECT_STRUCTURE.md QUICKSTART.md .gitignore LICENSE data/README.md 2>/dev/null || echo "  (Some files may be missing - download them first)"

echo ""
read -p "Continue with git add/commit/push? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

# Add all documentation files
echo ""
echo "Adding files to git..."
git add README.md
git add PROJECT_STRUCTURE.md
git add QUICKSTART.md
git add .gitignore
git add LICENSE
git add data/README.md

# Show status
echo ""
echo "Git status:"
git status

# Commit
echo ""
echo "Committing..."
git commit -m "docs: Add comprehensive project documentation

- Updated README with manuscript-ready overview
- Added PROJECT_STRUCTURE.md (complete repository guide)
- Added QUICKSTART.md (reproducibility instructions)
- Added data/README.md (data directory reference)
- Added MIT LICENSE (Copyright 2025 A.C. Demidont, DO)
- Configured .gitignore for Python/results files

Ready for Nature Communications submission."

# Push
echo ""
echo "Pushing to GitHub..."
git push origin main

echo ""
echo "✅ Documentation updated on GitHub!"
echo ""
echo "View your repository:"
echo "https://github.com/TheonlyqueenAC/noise_decorrelation_HIV"
