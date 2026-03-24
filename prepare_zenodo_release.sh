#!/bin/bash
# Prepare repository for Zenodo release
# Run from: ~/Documents/Github/noise_decorrelation_hiv/

set -e

echo "🎯 Preparing repository for Zenodo DOI..."
echo ""

# Check we're in the right directory
if [ ! -d ".git" ]; then
    echo "❌ Error: Not in a git repository!"
    echo "Run this from: ~/Documents/Github/noise_decorrelation_hiv/"
    exit 1
fi

# Add CITATION.cff file
echo "📝 Adding CITATION.cff file..."
if [ -f "CITATION.cff" ]; then
    echo "✓ CITATION.cff already exists"
else
    echo "❌ Please download CITATION.cff from Claude first"
    exit 1
fi

# Add LICENSE if not already there
echo ""
echo "📄 Checking LICENSE file..."
if [ -f "LICENSE" ]; then
    echo "✓ LICENSE exists"
else
    echo "⚠️  LICENSE file missing - download LICENSE.txt from Claude"
    echo "   and rename it to LICENSE (no extension)"
fi

# Commit these files
echo ""
echo "💾 Committing citation files..."
git add CITATION.cff
git add LICENSE 2>/dev/null || echo "  (LICENSE not found, skipping)"

if git diff --cached --quiet; then
    echo "  No new files to commit"
else
    git commit -m "chore: Add CITATION.cff for Zenodo and GitHub citation"
    echo "✓ Committed"
fi

# Push to GitHub
echo ""
echo "📤 Pushing to GitHub..."
git push origin main

echo ""
echo "✅ Repository prepared for Zenodo!"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "NEXT STEPS:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "1. CREATE GITHUB RELEASE:"
echo "   • Go to: https://github.com/Nyx-Dynamics/noise_decorrelation_hiv/releases"
echo "   • Click 'Create a new release'"
echo "   • Tag: v1.0.0"
echo "   • Title: v1.0.0 - Initial Manuscript Submission Release"
echo "   • Copy release notes from RELEASE_NOTES_v1.0.0.md"
echo "   • Publish release"
echo ""
echo "2. CONNECT TO ZENODO:"
echo "   • Go to: https://zenodo.org"
echo "   • Log in with GitHub"
echo "   • Go to GitHub settings (your profile → GitHub)"
echo "   • Find 'noise_decorrelation_hiv' and toggle ON"
echo "   • Wait 5-10 minutes for DOI to be generated"
echo ""
echo "3. GET YOUR DOI:"
echo "   • Go to Zenodo dashboard"
echo "   • Find your repository archive"
echo "   • Copy the DOI (10.5281/zenodo.XXXXXXX)"
echo ""
echo "4. UPDATE REPOSITORY WITH DOI:"
echo "   • Edit README.md and add DOI badge"
echo "   • Update CITATION.cff with real DOI"
echo "   • Update Data Availability section in manuscript"
echo ""
echo "See ZENODO_SETUP_GUIDE.md for detailed instructions!"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
