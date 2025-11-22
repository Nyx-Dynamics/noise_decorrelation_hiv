import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class TegmarkDecoherenceModel:
    """
    Implementation of Tegmark's quantum decoherence model with specific focus on
    differential decay between regular and Fibonacci grid patterns.

    Reference: Tegmark's "Importance of quantum decoherence in brain processes"
    with extensions for HIV protein quantum resonance.
    """

    def __init__(self):
        """Initialize the decoherence model with default parameters."""
        # Physical constants
        self.hbar = 1.0545718e-34  # Reduced Planck constant (J·s)
        self.k_B = 1.380649e-23  # Boltzmann constant (J/K)

        # Default parameters
        self.mass = 1.67e-27  # Typical protein mass (kg)
        self.length_scale = 1e-9  # Typical length scale (m)
        self.temp = 310.0  # Body temperature (K)

        # Grid-specific parameters
        self.regular_decay_exponent = 1.0  # Exponential decay for regular grid
        self.fibonacci_decay_exponent = -10.1  # Power law decay for Fibonacci grid

        # Environment coupling parameters
        self.env_coupling_regular = 0.01  # Environmental coupling for regular grid
        self.env_coupling_fibonacci = 0.001  # Environmental coupling for Fibonacci grid

        # Decoherence parameters
        self.lambda_thermal = None  # Thermal decoherence rate
        self.lambda_regular = None  # Regular grid additional decoherence
        self.lambda_fibonacci = None  # Fibonacci grid adjusted decoherence

        self._calculate_decoherence_rates()

    def _calculate_decoherence_rates(self):
        """Calculate decoherence rates based on Tegmark's formalism."""
        # Thermal decoherence rate - from Tegmark's equation
        # λ_thermal ≈ (k_B·T·Δx²) / (ħ²/m)
        self.lambda_thermal = (self.k_B * self.temp * self.length_scale ** 2) / \
                              (self.hbar ** 2 / self.mass)

        # For computational convenience, we scale to arbitrary units
        # where the thermal rate becomes 1.0
        self.lambda_thermal_scaled = 1.0

        # Regular grid follows exponential decay as per thermal models
        self.lambda_regular = self.lambda_thermal_scaled * (1.0 + self.env_coupling_regular)

        # Fibonacci grid follows power-law decay due to quantum protection
        # This is implemented in the differential equations
        self.lambda_fibonacci = self.lambda_thermal_scaled * self.env_coupling_fibonacci

    def set_temperature(self, temp):
        """
        Set the temperature of the system.

        Parameters:
        -----------
        temp : float
            Temperature in Kelvin
        """
        self.temp = temp
        self._calculate_decoherence_rates()

    def decoherence_function_regular(self, t, coherence):
        """
        Differential equation for regular grid decoherence.

        Parameters:
        -----------
        t : float
            Time
        coherence : float
            Current coherence value

        Returns:
        --------
        float
            Rate of change of coherence
        """
        # Exponential decay: dC/dt = -λ·C
        return -self.lambda_regular * coherence

    def decoherence_function_fibonacci(self, t, coherence):
        """
        Differential equation for Fibonacci grid decoherence.

        Parameters:
        -----------
        t : float
            Time
        coherence : float
            Current coherence value

        Returns:
        --------
        float
            Rate of change of coherence
        """
        # Power law behavior: dC/dt = -λ·C / (1 + t)^α
        # This creates power-law decay rather than exponential
        protection_factor = (1 + t) ** self.fibonacci_decay_exponent

        # Ensure protection factor doesn't go below minimum threshold
        min_protection = 1e-10
        if protection_factor < min_protection:
            protection_factor = min_protection

        return -self.lambda_fibonacci * coherence * protection_factor

    def simulate(self, t_span=(0, 10), t_points=1000, initial_coherence=1.0):
        """
        Simulate decoherence over time for both regular and Fibonacci grids.

        Parameters:
        -----------
        t_span : tuple
            Start and end time
        t_points : int
            Number of time points to simulate
        initial_coherence : float
            Initial coherence value

        Returns:
        --------
        tuple
            (time_points, regular_coherence, fibonacci_coherence)
        """
        # Time points
        t_eval = np.linspace(t_span[0], t_span[1], t_points)

        # Simulate regular grid decoherence
        solution_regular = solve_ivp(
            self.decoherence_function_regular,
            t_span,
            [initial_coherence],
            t_eval=t_eval,
            method='RK45'
        )

        # Simulate Fibonacci grid decoherence
        solution_fibonacci = solve_ivp(
            self.decoherence_function_fibonacci,
            t_span,
            [initial_coherence],
            t_eval=t_eval,
            method='RK45'
        )

        # Extract results
        time_points = solution_regular.t
        regular_coherence = solution_regular.y[0]
        fibonacci_coherence = solution_fibonacci.y[0]

        return time_points, regular_coherence, fibonacci_coherence

    def analytic_regular(self, t):
        """
        Analytical solution for regular grid coherence.

        Parameters:
        -----------
        t : ndarray
            Time points

        Returns:
        --------
        ndarray
            Coherence values
        """
        # Regular grid follows simple exponential decay
        return np.exp(-self.lambda_regular * t)

    def analytic_fibonacci(self, t):
        """
        Analytical solution for Fibonacci grid coherence.

        Parameters:
        -----------
        t : ndarray
            Time points

        Returns:
        --------
        ndarray
            Coherence values
        """
        # For Fibonacci grid with power-law protection
        if self.fibonacci_decay_exponent != -1:
            # For exponent != -1, we get this form
            exponent = self.fibonacci_decay_exponent + 1
            return np.exp(-self.lambda_fibonacci * ((1 + t) ** exponent - 1) / exponent)
        else:
            # For exponent = -1, we get logarithmic form
            return np.exp(-self.lambda_fibonacci * np.log(1 + t))

    def plot_decoherence(self, t_span=(0, 10), t_points=1000, use_analytic=True,
                         ax=None, log_scale=True, show_ratio=True):
        """
        Plot decoherence over time for both grids.

        Parameters:
        -----------
        t_span : tuple
            Start and end time
        t_points : int
            Number of time points
        use_analytic : bool
            Use analytical solutions instead of numerical
        ax : matplotlib.axes.Axes, optional
            Axes to plot on
        log_scale : bool
            Use logarithmic scale for y-axis
        show_ratio : bool
            Show coherence ratio plot

        Returns:
        --------
        matplotlib.figure.Figure
            Figure object
        """
        # Create time points
        t = np.linspace(t_span[0], t_span[1], t_points)

        # Calculate coherence values
        if use_analytic:
            regular_coherence = self.analytic_regular(t)
            fibonacci_coherence = self.analytic_fibonacci(t)
        else:
            t, regular_coherence, fibonacci_coherence = self.simulate(
                t_span, t_points
            )

        # Calculate coherence ratio
        coherence_ratio = fibonacci_coherence / regular_coherence

        # Create figure if not provided
        if ax is None:
            if show_ratio:
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
            else:
                fig, ax1 = plt.subplots(figsize=(10, 8))
                ax2 = None
        else:
            ax1 = ax
            ax2 = None
            fig = ax1.figure

        # Plot coherence values
        ax1.plot(t, regular_coherence, 'b-', linewidth=2, label='Regular Grid')
        ax1.plot(t, fibonacci_coherence, 'r-', linewidth=2, label='Fibonacci Grid')

        # Set logarithmic scale for y-axis if requested
        if log_scale:
            ax1.set_yscale('log')

        # Add labels and legend
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Quantum Coherence')
        ax1.set_title('Tegmark Decoherence Model: Regular vs. Fibonacci Grid')
        ax1.legend()
        ax1.grid(True)

        # Add coherence ratio plot if requested
        if show_ratio and ax2 is not None:
            ax2.plot(t, coherence_ratio, 'g-', linewidth=2)
            ax2.set_xlabel('Time')
            ax2.set_ylabel('Fibonacci/Regular Ratio')
            ax2.set_title('Coherence Advantage Ratio')
            ax2.axhline(y=1, color='k', linestyle='--')
            ax2.grid(True)

            # Annotate maximum advantage
            max_idx = np.argmax(coherence_ratio)
            max_time = t[max_idx]
            max_ratio = coherence_ratio[max_idx]

            ax2.plot(max_time, max_ratio, 'ro', markersize=8)
            ax2.annotate(
                f'Max: {max_ratio:.2f}x at t={max_time:.2f}',
                xy=(max_time, max_ratio),
                xytext=(max_time * 1.1, max_ratio * 0.9),
                arrowprops=dict(facecolor='red', shrink=0.05),
                ha='center'
            )

        plt.tight_layout()
        return fig

    def get_tegmark_timescales(self):
        """
        Calculate key timescales in Tegmark's model for both grid types.

        Returns:
        --------
        dict
            Dictionary of timescales in seconds
        """
        # Thermal decoherence timescale (1/λ)
        thermal_timescale = 1.0 / self.lambda_thermal

        # Regular grid decoherence timescale
        regular_timescale = 1.0 / self.lambda_regular

        # Fibonacci grid has a time-dependent timescale due to power-law decay
        # We'll calculate effective timescale at t=1 for reference
        t_ref = 1.0
        protection_factor = (1 + t_ref) ** self.fibonacci_decay_exponent
        fibonacci_timescale = 1.0 / (self.lambda_fibonacci * protection_factor)

        # Calculate long-term advantage factor
        # For t → ∞, the advantage scales with t^|fibonacci_exponent|
        long_term_advantage = 10 ** abs(self.fibonacci_decay_exponent)

        return {
            'thermal_decoherence': thermal_timescale,
            'regular_grid': regular_timescale,
            'fibonacci_grid': fibonacci_timescale,
            'theoretical_advantage': long_term_advantage
        }


def integrate_with_extended_analyzer(extended_analyzer):
    """
    Integrate Tegmark's decoherence model with the extended analyzer.

    Parameters:
    -----------
    extended_analyzer : HIVQuantumExtendedAnalyzer
        The extended analyzer object

    Returns:
    --------
    function
        Function to add to extended analyzer
    """

    def analyze_tegmark_decoherence(
            self, temperature=310.0, t_span=(0, 3), t_points=301,
            save_to_file=True, output_file="data/tegmark_decoherence.csv"):
        """
        Analyze quantum decoherence using Tegmark's model.

        Parameters:
        -----------
        temperature : float
            Temperature in Kelvin
        t_span : tuple
            Start and end time for simulation
        t_points : int
            Number of time points
        save_to_file : bool
            Save results to CSV file
        output_file : str
            Output file path

        Returns:
        --------
        dict
            Decoherence analysis results
        """
        print("\nAnalyzing quantum decoherence using Tegmark's model...")

        # Create Tegmark model
        model = TegmarkDecoherenceModel()
        model.set_temperature(temperature)

        # Get timescales
        timescales = model.get_tegmark_timescales()
        print(f"Thermal decoherence timescale: {timescales['thermal_decoherence']:.2e} s")
        print(f"Regular grid decoherence timescale: {timescales['regular_grid']:.2e} s")
        print(f"Fibonacci grid decoherence timescale: {timescales['fibonacci_grid']:.2e} s")
        print(f"Theoretical long-term advantage: {timescales['theoretical_advantage']:.2e}x")

        # Simulate decoherence
        time_points, regular_coherence, fibonacci_coherence = model.simulate(
            t_span, t_points
        )

        # Calculate coherence ratio
        coherence_ratio = fibonacci_coherence / regular_coherence

        # Create dataframe
        import pandas as pd
        df = pd.DataFrame({
            'Time': time_points,
            'Regular_Coherence': regular_coherence,
            'Fibonacci_Coherence': fibonacci_coherence,
            'Coherence_Ratio': coherence_ratio
        })

        # Save to file if requested
        if save_to_file:
            import os
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            df.to_csv(output_file, index=False)
            print(f"Saved Tegmark decoherence data to {output_file}")

        # Create visualization
        try:
            fig = model.plot_decoherence(t_span, t_points)
            if fig and hasattr(self, 'fig_dir') and self.fig_dir:
                fig_path = os.path.join(self.fig_dir, 'tegmark_decoherence_analysis.png')
                fig.savefig(fig_path, dpi=300, bbox_inches='tight')
                print(f"Saved figure to {fig_path}")
        except Exception as e:
            print(f"Error creating Tegmark decoherence plot: {e}")

        # Calculate key metrics
        # Find maximum advantage
        max_idx = np.argmax(coherence_ratio)
        max_time = float(time_points[max_idx])
        max_ratio = float(coherence_ratio[max_idx])

        # Calculate half-lives
        half_idx_regular = np.where(regular_coherence <= 0.5)[0]
        half_idx_fibonacci = np.where(fibonacci_coherence <= 0.5)[0]

        if len(half_idx_regular) > 0:
            half_time_regular = float(time_points[half_idx_regular[0]])
        else:
            half_time_regular = float('inf')

        if len(half_idx_fibonacci) > 0:
            half_time_fibonacci = float(time_points[half_idx_fibonacci[0]])
        else:
            half_time_fibonacci = float('inf')

        half_life_ratio = half_time_fibonacci / half_time_regular if half_time_regular > 0 else float('inf')

        # Calculate decoherence rates (negative log slope)
        # Using first and last 10% of points to estimate rates
        n_points = len(time_points)
        early_idx = int(n_points * 0.1)
        late_idx = int(n_points * 0.9)

        if early_idx < late_idx:
            # Regular grid - should be approximately constant (exponential decay)
            log_reg_early = np.log(regular_coherence[:early_idx + 1])
            time_early = time_points[:early_idx + 1]
            from scipy.stats import linregress
            reg_early_slope, _, _, _, _ = linregress(time_early, log_reg_early)

            log_reg_late = np.log(regular_coherence[late_idx:])
            time_late = time_points[late_idx:]
            reg_late_slope, _, _, _, _ = linregress(time_late, log_reg_late)

            # Fibonacci grid - should decrease over time (power law)
            log_fib_early = np.log(fibonacci_coherence[:early_idx + 1])
            fib_early_slope, _, _, _, _ = linregress(time_early, log_fib_early)

            log_fib_late = np.log(fibonacci_coherence[late_idx:])
            fib_late_slope, _, _, _, _ = linregress(time_late, log_fib_late)
        else:
            # Not enough points for regression
            reg_early_slope = reg_late_slope = fib_early_slope = fib_late_slope = 0

        # Prepare results
        results = {
            'timescales': timescales,
            'max_advantage': {
                'value': max_ratio,
                'time': max_time,
                'interpretation': f"Maximum coherence advantage of {max_ratio:.4f}x occurs at t={max_time:.2f}"
            },
            'half_life': {
                'regular': half_time_regular,
                'fibonacci': half_time_fibonacci,
                'ratio': half_life_ratio,
                'interpretation': f"Fibonacci grid half-life is {half_life_ratio:.4f}x longer than regular grid"
            },
            'decoherence_rates': {
                'regular_early': -reg_early_slope,
                'regular_late': -reg_late_slope,
                'fibonacci_early': -fib_early_slope,
                'fibonacci_late': -fib_late_slope,
                'interpretation': f"Regular grid shows constant decay rate (~{-reg_early_slope:.4f}), while Fibonacci grid shows decreasing rate (from {-fib_early_slope:.4f} to {-fib_late_slope:.4f})"
            },
            'model_parameters': {
                'regular_exponent': model.regular_decay_exponent,
                'fibonacci_exponent': model.fibonacci_decay_exponent,
                'temperature': temperature
            }
        }

        # Save results
        if hasattr(self, 'base') and self.base is not None:
            self.base.save_stats('tegmark_decoherence_analysis', results)
        elif hasattr(self, 'save_stats'):
            self.save_stats('tegmark_decoherence_analysis', results)

        return results, df

    # Add the method to the extended analyzer
    extended_analyzer.analyze_tegmark_decoherence = analyze_tegmark_decoherence.__get__(extended_analyzer)

    return analyze_tegmark_decoherence


# Direct usage example
if __name__ == "__main__":
    # Simple demonstration of the decoherence model
    model = TegmarkDecoherenceModel()

    # Set different exponents
    model.regular_decay_exponent = 1.0  # Exponential decay
    model.fibonacci_decay_exponent = -10.1  # Power law protection

    # Calculate decoherence over time
    t, reg_coherence, fib_coherence = model.simulate(t_span=(0, 10), t_points=1000)

    # Calculate coherence ratio
    ratio = fib_coherence / reg_coherence

    # Print key results
    print("Tegmark Decoherence Model Results")
    print("================================")
    print(f"Regular grid decay exponent: {model.regular_decay_exponent}")
    print(f"Fibonacci grid decay exponent: {model.fibonacci_decay_exponent}")

    # Calculate decoherence at specific times
    times = [0.1, 1.0, 10.0]
    for time in times:
        idx = np.abs(t - time).argmin()
        print(f"\nAt t = {t[idx]:.2f}:")
        print(f"  Regular coherence: {reg_coherence[idx]:.6e}")
        print(f"  Fibonacci coherence: {fib_coherence[idx]:.6e}")
        print(f"  Coherence ratio: {ratio[idx]:.4f}x")

    # Plot results
    fig = model.plot_decoherence(t_span=(0, 10), t_points=1000)
    plt.show()