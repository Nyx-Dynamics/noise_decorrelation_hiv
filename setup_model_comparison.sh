#!/bin/bash
# Setup Model Comparison Analysis
# Run from: ~/Documents/Github/noise_decorrelation_HIV/

set -e

echo "🔬 Setting up Model Comparison Analysis"
echo ""

# Check we're in the right place
if [ ! -d "quantum" ]; then
    echo "❌ Error: quantum/ directory not found"
    echo "Run this from: ~/Documents/Github/noise_decorrelation_HIV/"
    exit 1
fi

# Create results directories
echo "📁 Creating results directories..."
mkdir -p quantum/results/model_comparison_group
mkdir -p quantum/results/model_comparison_individual
echo "✓ Created"

# Copy scripts to quantum directory
echo ""
echo "📋 Copying comparison scripts..."

if [ -f "model_comparison_GROUP_LEVEL.py" ]; then
    cp model_comparison_GROUP_LEVEL.py quantum/
    echo "  ✓ Copied group-level script"
else
    echo "  ⚠️  model_comparison_GROUP_LEVEL.py not found - download from Claude"
fi

if [ -f "model_comparison_INDIVIDUAL_LEVEL.py" ]; then
    cp model_comparison_INDIVIDUAL_LEVEL.py quantum/
    echo "  ✓ Copied individual-level script"
else
    echo "  ⚠️  model_comparison_INDIVIDUAL_LEVEL.py not found - download from Claude"
fi

# Check for data files
echo ""
echo "📊 Checking data files..."

if [ -f "CRITICAL_STUDIES_COMPLETE_DATA.csv" ] || [ -f "data/extracted/CRITICAL_STUDIES_COMPLETE_DATA.csv" ]; then
    echo "  ✓ Group-level data found"
else
    echo "  ❌ CRITICAL_STUDIES_COMPLETE_DATA.csv not found"
    echo "     Need this for group-level analysis!"
fi

if [ -f "VALCOUR_2015_INDIVIDUAL_PATIENTS.csv" ] || [ -f "data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv" ]; then
    echo "  ✓ Individual patient data found"
else
    echo "  ❌ VALCOUR_2015_INDIVIDUAL_PATIENTS.csv not found"
    echo "     Need this for individual-level analysis!"
fi

# Check Python environment
echo ""
echo "🐍 Checking Python environment..."

if command -v python &> /dev/null; then
    python_version=$(python --version 2>&1)
    echo "  ✓ Python found: $python_version"
    
    # Check key packages
    python -c "import pymc" 2>/dev/null && echo "  ✓ PyMC installed" || echo "  ❌ PyMC not installed"
    python -c "import arviz" 2>/dev/null && echo "  ✓ ArviZ installed" || echo "  ❌ ArviZ not installed"
else
    echo "  ❌ Python not found in PATH"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "READY TO RUN MODEL COMPARISON"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Step 1: Run group-level analysis (main manuscript)"
echo "  cd quantum/"
echo "  python model_comparison_GROUP_LEVEL.py"
echo ""
echo "Step 2: Run individual-level analysis (validation)"
echo "  cd quantum/"
echo "  python model_comparison_INDIVIDUAL_LEVEL.py"
echo ""
echo "See MODEL_COMPARISON_GUIDE.md for detailed instructions!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
