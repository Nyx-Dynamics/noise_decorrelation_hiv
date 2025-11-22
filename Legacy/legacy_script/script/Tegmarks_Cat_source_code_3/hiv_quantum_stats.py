if 'temperature_profile' in temp and 'statistics' in temp['temperature_profile']:
    stats = temp['temperature_profile']['statistics']
    f.write(f"- Temperature range: **{stats['min_temperature']:.2f}°C - {stats['max_temperature']:.2f}°C**\n")
    f.write(f"- Average temperature: {stats['average_temperature']:.2f}°C\n")
    f.write(f"- Temperature standard deviation: {stats['temperature_std']:.2f}°C\n")
    f.write(
        f"- Maximum temperature of {stats['max_temperature']:.2f}°C occurs at t = {stats['time_of_max_temperature']:.2f}\n")

    # Fever metrics
if 'temperature_profile' in temp and 'fever_metrics' in temp['temperature_profile']:
    fever = temp['temperature_profile']['fever_metrics']
    f.write(f"\n- Temperature exceeds the fever threshold of {fever['threshold']:.1f}°C for:")
    f.write(f"\n  - Duration: {fever['duration']:.2f} time units\n")
    f.write(f"  - Percentage: {fever['percentage_of_time']:.1f}% of the simulation\n")

    # Cyclical analysis
if 'temperature_profile' in temp and 'cyclical_analysis' in temp['temperature_profile']:
    cycle = temp['temperature_profile']['cyclical_analysis']
    f.write(f"\n- Dominant temperature cycle: **{cycle['dominant_period']:.2f} time units**\n")

    # Integrity analysis
f.write("\n### Temperature Effects on Structural Integrity\n\n")
if 'integrity_analysis' in temp and 'regression_analysis' in temp['integrity_analysis']:
    reg = temp['integrity_analysis']['regression_analysis']

    # Regular grid
    if 'Regular' in reg:
        r_reg = reg['Regular']
        f.write(f"- **Regular Grid**:\n")
        f.write(f"  - Temperature sensitivity: {r_reg['slope']:.4f} integrity units per °C\n")
        f.write(f"  - R²: {r_reg['r_squared']:.4f}, p-value: {r_reg['p_value']:.4e}\n")
        f.write(f"  - {r_reg['interpretation']}\n")

    # Fibonacci grid
    if 'Fibonacci' in reg:
        r_fib = reg['Fibonacci']
        f.write(f"\n- **Fibonacci Grid**:\n")
        f.write(f"  - Temperature sensitivity: {r_fib['slope']:.4f} integrity units per °C\n")
        f.write(f"  - R²: {r_fib['r_squared']:.4f}, p-value: {r_fib['p_value']:.4e}\n")
        f.write(f"  - {r_fib['interpretation']}\n")

    # Comparison
    if 'comparison' in reg:
        comp = reg['comparison']
        f.write(f"\n- **Comparison**:\n")
        f.write(f"  - {comp['interpretation']}\n")

# Resonance analysis
if 'resonance_analysis' in temp and 'correlation_analysis' in temp['resonance_analysis']:
    f.write("\n### Temperature Effects on Phi-Resonance\n\n")
    corr = temp['resonance_analysis']['correlation_analysis']

    # Regular grid
    if 'Regular' in corr:
        c_reg = corr['Regular']
        f.write(f"- **Regular Grid**:\n")
        f.write(f"  - Correlation with temperature: r = {c_reg['correlation']:.4f}, p-value: {c_reg['p_value']:.4e}\n")
        f.write(f"  - {c_reg['interpretation']}\n")

    # Fibonacci grid
    if 'Fibonacci' in corr:
        c_fib = corr['Fibonacci']
        f.write(f"\n- **Fibonacci Grid**:\n")
        f.write(f"  - Correlation with temperature: r = {c_fib['correlation']:.4f}, p-value: {c_fib['p_value']:.4e}\n")
        f.write(f"  - {c_fib['interpretation']}\n")

    # Comparison
    if 'comparison' in corr:
        comp = corr['comparison']
        f.write(f"\n- **Comparison**:\n")
        f.write(f"  - {comp['interpretation']}\n")
else:
    f.write("Temperature effects analysis was not performed due to missing data.\n")

f.write("\n")

# Methodology section
f.write("## 4. Methodology\n\n")

f.write("### Statistical Methods\n\n")
f.write("This analysis employed the following statistical techniques:\n\n")
f.write("- **Power Law Fitting**: Nonlinear least squares regression to determine decay exponents\n")
f.write("- **Statistical Significance Testing**: Paired t-tests to assess differences between grid types\n")
f.write("- **Correlation Analysis**: Pearson correlation to quantify relationships between variables\n")
f.write("- **Regression Analysis**: Linear regression to model temperature dependencies\n")
f.write("- **Gaussian Curve Fitting**: To characterize the sensitivity to scaling factor variations\n")

f.write("\n### Software\n\n")
f.write("Analysis was performed using Python with the following libraries:\n\n")
f.write("- **pandas**: Data manipulation and analysis\n")
f.write("- **numpy**: Numerical computing\n")
f.write("- **scipy**: Statistical tests and curve fitting\n")
f.write("- **matplotlib**: Data visualization\n")
if HAVE_STATSMODELS:
    f.write("- **statsmodels**: Statistical modeling\n")

f.write("\n### Significance Criteria\n\n")
f.write("- Statistical significance determined using a threshold of p < 0.05\n")
f.write("- R-squared values used to assess goodness-of-fit\n")
f.write("- Correlation strengths interpreted according to standard conventions\n")

f.write("\n")

# Data Sources section
f.write("## 5. Data Sources\n\n")

f.write("The following data files were used in this analysis:\n\n")

for name in self.data.keys():
    if isinstance(self.data[name], pd.DataFrame):
        rows = len(self.data[name])
        cols = len(self.data[name].columns)
        f.write(f"- **{name}**: {rows} rows × {cols} columns\n")
    else:
        f.write(f"- **{name}**\n")

f.write("\n")

# Footer
f.write("---\n\n")
f.write("Report generated by HIVQuantumStatAnalyzer\n")

print(f"Report generated at {report_path}")
return report_path
except Exception as e:
print(f"Error generating report: {e}")
return None


def analyze_all(self):
    """Run all available analyses."""
    print("\nRunning all analyses...")

    # First, make sure we have data
    if not self.data:
        print("No data loaded. Loading data files...")
        self.load_data("*.csv")
        self.load_data("*.json")

    # Run basic statistics
    self.compute_basic_statistics()

    # Run specific analyses
    print("\nRunning specialized analyses...")
    self.analyze_core_dynamics()
    self.analyze_golden_ratio()
    self.analyze_temperature_effects()

    # Generate report
    report_path = self.generate_report()

    print("\nAll analyses complete!")
    if report_path:
        print(f"Comprehensive report saved to {report_path}")

    return self.stats


def main():
    """Parse arguments and run analysis."""
    parser = argparse.ArgumentParser(description='Statistical analysis of HIV quantum coherence simulation data')

    parser.add_argument('--data-dir', type=str, required=True,
                        help='Directory containing simulation data files')

    parser.add_argument('--output-dir', type=str, default=None,
                        help='Directory where analysis results will be saved')

    parser.add_argument('--report', action='store_true',
                        help='Generate a comprehensive analysis report')

    parser.add_argument('--plot', action='store_true',
                        help='Show plots during analysis')

    args = parser.parse_args()

    # Print banner
    print("=" * 80)
    print("HIV Quantum Coherence Statistical Analyzer")
    print("=" * 80)

    # Create analyzer
    analyzer = HIVQuantumStatAnalyzer(args.data_dir, args.output_dir)

    # Load all data
    analyzer.load_data("*.csv")
    analyzer.load_data("*.json")

    if args.report:
        # Run all analyses and generate report
        analyzer.analyze_all()
    else:
        # Run basic statistics
        analyzer.compute_basic_statistics()

        # Individual analyses
        analyzer.analyze_core_dynamics()
        analyzer.analyze_golden_ratio()
        analyzer.analyze_temperature_effects()

    # Show plots if requested
    if args.plot:
        plt.show()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)