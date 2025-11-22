#!/bin/bash

# Check if data directory exists
if [ -d "data/analysis" ]; then
    echo "Cleaning up previous analysis results..."
    rm -rf data/analysis
fi

# Create a results directory where we'll copy the analysis files
mkdir -p results

# Run the improved data generation script
echo "Generating fresh sample data..."
./generate_data.sh

# Run the analysis
echo "Running the quantum coherence analysis..."
python hiv_quantum_stats.py

# Copy important results to the results directory
echo "Copying results to the results directory..."
cp -r data/analysis/figures results/
cp data/analysis/*.json results/
cp data/analysis/*.md results/

echo "Analysis complete! Results are available in the 'results' directory."
echo "Key files:"
echo "  - results/quantum_coherence_analysis_report.md - Main analysis report"
echo "  - results/figures/ - Visualizations of the analysis"
echo "  - results/*.json - Detailed statistical results"