import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.linalg import expm
from scipy.stats import linregress
import time

# Set the random seed for reproducibility
np.random.seed(42)


# Simulation parameters
class SimulationParameters:
    def __init__(self):
        # System parameters
        self.hbar = 1.0  # Reduced Planck constant (normalized)
        self.m = 1.0  # Effective mass (normalized)
        self.L = 10.0  # System length
        self.dt = 0.01  # Time step
        self.t_total = 3.0  # Total simulation time

        # Grid parameters
        self.grid_size = 32  # Grid size for computations

        # Initial state parameters
        self.alpha = 3.0  # Cat state displacement
        self.nth = 3.48  # Thermal excitation number for acute HIV

        # Decoherence parameters (HIV-specific)
        self.gamma_base = 0.05  # Base decoherence rate
        self.cytokine_coupling = 0.1  # Coupling to inflammatory cytokines

        # Tegmark decoherence time parameter
        self.tegmark_tau = 10 ** -13  # Tegmark's calculated decoherence time

        # Golden ratio (for Fibonacci structures)
        self.phi = 1.618

        # HIV inflammatory parameters (acute phase)
        self.tnf_alpha = 75.3  # TNF-alpha level in pg/ml (acute HIV)
        self.il6 = 30.2  # IL-6 level in pg/ml (acute HIV)
        self.il1beta = 15.7  # IL-1beta level in pg/ml (acute HIV)

        # Temperature parameters for HIV fever
        self.baseline_temp = 37.0  # Baseline temperature in Celsius
        self.fever_temp = 39.5  # Fever temperature in acute HIV

        # Power law exponents from the paper
        self.fibonacci_exponent = -1.0150
        self.regular_exponent = -10.1023


# Initialize parameters
params = SimulationParameters()


# Define grid structures
class GridStructure:
    def __init__(self, params):
        self.params = params
        self.grid_size = params.grid_size

    def regular_grid(self):
        """Create a regular square grid with uniform spacing"""
        x = np.linspace(-self.params.L / 2, self.params.L / 2, self.grid_size)
        y = np.linspace(-self.params.L / 2, self.params.L / 2, self.grid_size)
        X, Y = np.meshgrid(x, y)
        return X, Y

    def fibonacci_grid(self):
        """Create a grid with Fibonacci scaling"""
        # Start with regular grid
        X, Y = self.regular_grid()

        # Apply Fibonacci scaling to the radial coordinate
        R = np.sqrt(X ** 2 + Y ** 2)
        Theta = np.arctan2(Y, X)

        # Apply Fibonacci scaling where r_n = r_0 * (F_{n+2}/F_n)
        # We approximate this with the golden ratio scaling
        R_fibonacci = R * (1 + (R / self.params.L) * (self.params.phi - 1))

        # Convert back to Cartesian coordinates
        X_fibonacci = R_fibonacci * np.cos(Theta)
        Y_fibonacci = R_fibonacci * np.sin(Theta)

        return X_fibonacci, Y_fibonacci


# Define the decoherence model incorporating Tegmark's calculations and HIV inflammatory effects
class DecoherenceModel:
    def __init__(self, params):
        self.params = params

    def tegmark_decoherence_rate(self, temperature):
        """Calculate decoherence rate based on Tegmark's formula adjusted for temperature"""
        # Tegmark's decoherence time is 10^-13 seconds at body temperature
        # We adjust for temperature using exponential scaling (for regular grid)
        temp_factor = np.exp(0.1 * (temperature - self.params.baseline_temp))
        return 1.0 / (self.params.tegmark_tau * temp_factor)

    def inflammatory_factor(self):
        """Calculate inflammatory contribution to decoherence from cytokines"""
        # Simple model where cytokines linearly increase decoherence
        tnf_contribution = 0.01 * self.params.tnf_alpha
        il6_contribution = 0.005 * self.params.il6
        il1b_contribution = 0.008 * self.params.il1beta

        return tnf_contribution + il6_contribution + il1b_contribution

    def decoherence_rate_regular(self, temperature):
        """Calculate total decoherence rate for regular grid"""
        base_rate = self.tegmark_decoherence_rate(temperature)
        inflammatory_contribution = self.inflammatory_factor()

        # Regular grid has strong exponential dependence on inflammation
        return base_rate * (1.0 + np.exp(inflammatory_contribution))

    def decoherence_rate_fibonacci(self, temperature):
        """Calculate total decoherence rate for Fibonacci grid"""
        # From the paper: Fibonacci structures have power-law rather than exponential response
        # and much milder temperature dependence
        base_rate = self.tegmark_decoherence_rate(self.params.baseline_temp)
        temp_factor = 1.0 + 0.05 * (temperature - self.params.baseline_temp)
        inflammatory_contribution = self.inflammatory_factor()

        # Fibonacci grid has much weaker response to inflammation (power law vs exponential)
        return (base_rate / self.params.phi) * temp_factor * (1.0 + inflammatory_contribution)


# Quantum state evolution with decoherence
class QuantumSystem:
    def __init__(self, params, grid_type='regular'):
        self.params = params
        self.grid_type = grid_type
        self.grid = GridStructure(params)
        self.decoherence = DecoherenceModel(params)

        # Initialize the grid based on type
        if grid_type == 'regular':
            self.X, self.Y = self.grid.regular_grid()
        elif grid_type == 'fibonacci':
            self.X, self.Y = self.grid.fibonacci_grid()
        else:
            raise ValueError("Grid type must be 'regular' or 'fibonacci'")

        # Initialize state
        self.initialize_state()

        # Initialize time and observables
        self.t = 0
        self.coherence_history = []
        self.von_neumann_entropy_history = []
        self.boundary_radius_history = []

    def initialize_state(self):
        """Initialize a thermal cat state with parameters from Kirchmair's experiment"""
        # Create position and momentum grids
        dx = self.params.L / self.params.grid_size
        self.p_grid = np.fft.fftfreq(self.params.grid_size, dx) * 2 * np.pi

        # Create initial thermal state (simplified model)
        sigma = np.sqrt(0.5 + self.params.nth)  # Width related to thermal excitation
        thermal_state = np.exp(-(self.X ** 2 + self.Y ** 2) / (4 * sigma ** 2)) / (2 * np.pi * sigma ** 2)

        # Create cat state using the ECD approach from Kirchmair's paper
        # We simplify by directly creating a superposition of displaced thermal states
        alpha = self.params.alpha
        cat_state_right = np.exp(-((self.X - alpha) ** 2 + self.Y ** 2) / (4 * sigma ** 2)) / (2 * np.pi * sigma ** 2)
        cat_state_left = np.exp(-((self.X + alpha) ** 2 + self.Y ** 2) / (4 * sigma ** 2)) / (2 * np.pi * sigma ** 2)

        # Superposition (simplified version of hot Schrödinger cat)
        self.state = (cat_state_right + cat_state_left) / np.sqrt(2)

        # Normalize
        self.state = self.state / np.sqrt(np.sum(np.abs(self.state) ** 2) * dx ** 2)

    def time_evolution_step(self, temperature):
        """Evolve the quantum state for one time step with decoherence"""
        # Calculate decoherence rate based on grid type
        if self.grid_type == 'regular':
            gamma = self.decoherence.decoherence_rate_regular(temperature)
        else:
            gamma = self.decoherence.decoherence_rate_fibonacci(temperature)

        # Apply decoherence operation (simplified Lindblad evolution)
        state_matrix = np.outer(self.state.flatten(), np.conjugate(self.state.flatten()))
        decoherence_op = np.exp(-gamma * self.params.dt)

        # Off-diagonal elements decay exponentially with decoherence
        off_diag_mask = ~np.eye(state_matrix.shape[0], dtype=bool)
        state_matrix[off_diag_mask] *= decoherence_op

        # Reshape back to 2D grid representation
        self.state = np.sqrt(np.diag(state_matrix)).reshape(self.params.grid_size, self.params.grid_size)

        # Normalize the state
        dx = self.params.L / self.params.grid_size
        self.state = self.state / np.sqrt(np.sum(np.abs(self.state) ** 2) * dx ** 2)

        # Update time
        self.t += self.params.dt

    def calculate_coherence(self):
        """Calculate quantum coherence as per the paper's definition"""
        # Simplified coherence measure: off-diagonal elements in density matrix
        state_vec = self.state.flatten()
        density_matrix = np.outer(state_vec, np.conjugate(state_vec))

        # Extract off-diagonal elements
        off_diag_mask = ~np.eye(density_matrix.shape[0], dtype=bool)
        coherence = np.sum(np.abs(density_matrix[off_diag_mask]))

        return coherence

    def calculate_von_neumann_entropy(self):
        """Calculate von Neumann entropy as a measure of quantum information"""
        # Create density matrix
        state_vec = self.state.flatten()
        density_matrix = np.outer(state_vec, np.conjugate(state_vec))

        # Calculate eigenvalues of density matrix
        eigenvalues = np.linalg.eigvalsh(density_matrix)

        # Remove zeros and very small values
        eigenvalues = eigenvalues[eigenvalues > 1e-10]

        # Calculate von Neumann entropy: S = -Tr(ρ ln ρ)
        entropy = -np.sum(eigenvalues * np.log(eigenvalues))

        return entropy

    def detect_boundary(self):
        """Identify the event horizon-like boundary as described in the paper"""
        # Calculate probability current (simplified)
        dx = self.params.L / self.params.grid_size

        # Calculate gradient of probability density
        rho = np.abs(self.state) ** 2
        grad_x = np.gradient(rho, dx, axis=1)
        grad_y = np.gradient(rho, dx, axis=0)

        # Calculate divergence of probability current
        div_J = np.gradient(grad_x, dx, axis=1) + np.gradient(grad_y, dx, axis=0)

        # Find points where divergence is zero (boundary points)
        # We simplify by looking at radial average
        r_grid = np.sqrt(self.X ** 2 + self.Y ** 2)
        r_bins = np.linspace(0, np.max(r_grid), 100)
        r_centers = (r_bins[:-1] + r_bins[1:]) / 2

        # Average divergence in radial bins
        div_J_radial = np.zeros_like(r_centers)
        for i in range(len(r_centers)):
            mask = (r_grid >= r_bins[i]) & (r_grid < r_bins[i + 1])
            if np.any(mask):
                div_J_radial[i] = np.mean(np.abs(div_J[mask]))

        # Find the first minimum in the radial profile as the boundary
        # Smooth the profile first
        div_J_smooth = np.convolve(div_J_radial, np.ones(5) / 5, mode='same')

        # Find local minima
        minima_indices = np.where((div_J_smooth[1:-1] < div_J_smooth[:-2]) &
                                  (div_J_smooth[1:-1] < div_J_smooth[2:]))[0] + 1

        if len(minima_indices) > 0:
            # Return the radius of the first minimum
            boundary_radius = r_centers[minima_indices[0]]
        else:
            # Default to half the grid if no minimum found
            boundary_radius = np.max(r_centers) / 2

        return boundary_radius

    def run_simulation(self, temperature_profile=None):
        """Run the full simulation with time-dependent temperature"""
        # Default to constant baseline temperature if no profile provided
        if temperature_profile is None:
            temperature_profile = lambda t: self.params.baseline_temp

        n_steps = int(self.params.t_total / self.params.dt)

        # Lists to store observables
        times = [0]
        coherence = [self.calculate_coherence()]
        entropy = [self.calculate_von_neumann_entropy()]
        boundary = [self.detect_boundary()]

        print(f"Starting simulation with {self.grid_type} grid...")
        start_time = time.time()

        for step in range(n_steps):
            # Get current temperature
            temp = temperature_profile(self.t)

            # Evolve the state
            self.time_evolution_step(temp)

            # Calculate observables every 10 steps to save computation
            if step % 10 == 0:
                times.append(self.t)
                coherence.append(self.calculate_coherence())
                entropy.append(self.calculate_von_neumann_entropy())
                boundary.append(self.detect_boundary())

                # Print progress every 100 steps
                if step % 100 == 0:
                    elapsed = time.time() - start_time
                    print(f"Step {step}/{n_steps}, elapsed time: {elapsed:.2f}s, coherence: {coherence[-1]:.6f}")

        print(f"Simulation complete. Total time: {time.time() - start_time:.2f}s")

        # Store results
        self.coherence_history = np.array(coherence)
        self.von_neumann_entropy_history = np.array(entropy)
        self.boundary_radius_history = np.array(boundary)
        self.time_points = np.array(times)

        return times, coherence, entropy, boundary


# Define a realistic HIV fever temperature profile
def hiv_fever_profile(t, baseline=37.0, max_temp=39.5):
    """Generate a realistic fever profile for acute HIV infection"""
    # Simulate temperature with daily oscillations and fever pattern
    day_fraction = (t % 1.0) / 1.0  # Time of day (normalized to 1.0)

    # Circadian rhythm (lower in early morning, higher in evening)
    circadian = 0.5 * np.sin(2 * np.pi * day_fraction - np.pi / 3)

    # Fever intensity starts high and gradually decreases
    fever_intensity = np.exp(-0.3 * t)

    # Calculate temperature
    temperature = baseline + circadian + fever_intensity * (max_temp - baseline)

    return temperature


# Function to analyze and plot coherence decay patterns
def analyze_decay(system, params):
    """Analyze the decay pattern of coherence and fit power law"""
    times = system.time_points[1:]  # Skip t=0
    coherence = system.coherence_history[1:]  # Skip initial coherence

    # Normalize coherence
    coherence_norm = coherence / coherence[0]

    # Convert to log scale
    log_times = np.log(times)
    log_coherence = np.log(coherence_norm)

    # Fit power law (linear in log-log scale)
    mask = np.isfinite(log_times) & np.isfinite(log_coherence)
    slope, intercept, r_value, p_value, std_err = linregress(log_times[mask], log_coherence[mask])

    print(f"Power law fit for {system.grid_type} grid:")
    print(f"Exponent: {slope:.4f}")
    print(f"R-squared: {r_value ** 2:.4f}")
    print(
        f"Expected exponent from paper: {params.fibonacci_exponent if system.grid_type == 'fibonacci' else params.regular_exponent}")

    # Fitted curve
    fit_coherence = np.exp(intercept) * times ** slope

    return times, coherence_norm, slope, r_value ** 2, fit_coherence


# Execute the simulation for both grid types
def run_both_simulations(params):
    # Create temperature profile function
    temp_profile = lambda t: hiv_fever_profile(t, params.baseline_temp, params.fever_temp)

    # Run simulations for both grid types
    regular_system = QuantumSystem(params, grid_type='regular')
    regular_system.run_simulation(temp_profile)

    fibonacci_system = QuantumSystem(params, grid_type='fibonacci')
    fibonacci_system.run_simulation(temp_profile)

    # Analyze decay patterns
    regular_times, regular_coherence, regular_slope, regular_r2, regular_fit = analyze_decay(regular_system, params)
    fib_times, fib_coherence, fib_slope, fib_r2, fib_fit = analyze_decay(fibonacci_system, params)

    # Calculate coherence ratio
    coherence_ratio = np.zeros_like(regular_system.time_points)
    for i in range(len(coherence_ratio)):
        if regular_system.coherence_history[i] > 0:
            coherence_ratio[i] = fibonacci_system.coherence_history[i] / regular_system.coherence_history[i]

    # Find maximum advantage
    max_ratio_idx = np.argmax(coherence_ratio)
    max_ratio_time = regular_system.time_points[max_ratio_idx]
    max_ratio = coherence_ratio[max_ratio_idx]

    print(f"Maximum coherence advantage: {max_ratio:.4f}x at t = {max_ratio_time:.2f}")

    # Find sanctuary formation time (if applicable)
    # Paper defines this as where regular grid coherence has collapsed to ~0.5% of initial
    sanctuary_threshold = 0.005
    sanctuary_idx = np.where(regular_coherence < sanctuary_threshold)[0]
    if len(sanctuary_idx) > 0:
        sanctuary_time = regular_times[sanctuary_idx[0]]
        print(f"Sanctuary formation time: t = {sanctuary_time:.2f}")

        # Coherence values at sanctuary formation
        reg_coherence_at_sanctuary = regular_coherence[sanctuary_idx[0]]
        fib_coherence_at_sanctuary = fib_coherence[sanctuary_idx[0]]
        print(f"Regular grid coherence at sanctuary formation: {reg_coherence_at_sanctuary:.4f}")
        print(f"Fibonacci grid coherence at sanctuary formation: {fib_coherence_at_sanctuary:.4f}")
        print(f"Advantage at sanctuary formation: {fib_coherence_at_sanctuary / reg_coherence_at_sanctuary:.4f}x")
    else:
        sanctuary_time = None
        print("No sanctuary formation detected")

    # Calculate half-lives
    reg_half_life_idx = np.where(regular_coherence < 0.5)[0]
    fib_half_life_idx = np.where(fib_coherence < 0.5)[0]

    if len(reg_half_life_idx) > 0:
        reg_half_life = regular_times[reg_half_life_idx[0]]
        print(f"Regular grid half-life: {reg_half_life:.4f}")
    else:
        reg_half_life = None
        print("Regular grid coherence did not reach half-life")

    if len(fib_half_life_idx) > 0:
        fib_half_life = fib_times[fib_half_life_idx[0]]
        print(f"Fibonacci grid half-life: {fib_half_life:.4f}")
    else:
        fib_half_life = None
        print("Fibonacci grid coherence did not reach half-life")

    if reg_half_life and fib_half_life:
        print(f"Half-life ratio: {fib_half_life / reg_half_life:.4f}")

    # Return all results for plotting
    return {
        'regular_system': regular_system,
        'fibonacci_system': fibonacci_system,
        'regular_times': regular_times,
        'regular_coherence': regular_coherence,
        'regular_slope': regular_slope,
        'regular_r2': regular_r2,
        'regular_fit': regular_fit,
        'fib_times': fib_times,
        'fib_coherence': fib_coherence,
        'fib_slope': fib_slope,
        'fib_r2': fib_r2,
        'fib_fit': fib_fit,
        'coherence_ratio': coherence_ratio,
        'max_ratio': max_ratio,
        'max_ratio_time': max_ratio_time,
        'sanctuary_time': sanctuary_time
    }


# Create visualization functions
def plot_coherence_comparison(results):
    """Create plot comparing coherence evolution between grid types"""
    plt.figure(figsize=(12, 8))

    # Top left: Coherence evolution
    plt.subplot(2, 2, 1)
    plt.plot(results['regular_times'], results['regular_coherence'], 'b-', label='Regular Grid')
    plt.plot(results['fib_times'], results['fib_coherence'], 'r-', label='Fibonacci Grid')
    plt.plot(results['regular_times'], results['regular_fit'], 'b--',
             label=f'Power law fit (λ={results["regular_slope"]:.2f})')
    plt.plot(results['fib_times'], results['fib_fit'], 'r--', label=f'Power law fit (λ={results["fib_slope"]:.2f})')

    if results['sanctuary_time']:
        plt.axvline(x=results['sanctuary_time'], color='g', linestyle='--', label='Sanctuary formation')

    plt.title('Quantum Coherence Dynamics During Acute HIV Inflammation')
    plt.xlabel('Time')
    plt.ylabel('Normalized Coherence')
    plt.legend()
    plt.grid(True)
    plt.yscale('log')
    plt.xscale('log')

    # Top right: Coherence ratio
    plt.subplot(2, 2, 2)
    plt.plot(results['regular_system'].time_points, results['coherence_ratio'], 'g-')
    plt.axhline(y=1, color='k', linestyle='--')

    if results['sanctuary_time']:
        plt.axvline(x=results['sanctuary_time'], color='g', linestyle='--')

    plt.title(
        f'Fibonacci/Regular Coherence Ratio\nMax: {results["max_ratio"]:.2f}x at t={results["max_ratio_time"]:.2f}')
    plt.xlabel('Time')
    plt.ylabel('Coherence Ratio')
    plt.grid(True)

    # Bottom left: Von Neumann entropy
    plt.subplot(2, 2, 3)
    plt.plot(results['regular_system'].time_points, results['regular_system'].von_neumann_entropy_history, 'b-',
             label='Regular Grid')
    plt.plot(results['fibonacci_system'].time_points, results['fibonacci_system'].von_neumann_entropy_history, 'r-',
             label='Fibonacci Grid')

    if results['sanctuary_time']:
        plt.axvline(x=results['sanctuary_time'], color='g', linestyle='--')

    plt.title('Von Neumann Entropy Evolution')
    plt.xlabel('Time')
    plt.ylabel('Entropy')
    plt.legend()
    plt.grid(True)

    # Bottom right: Event horizon evolution
    plt.subplot(2, 2, 4)
    plt.plot(results['regular_system'].time_points, results['regular_system'].boundary_radius_history, 'b-',
             label='Regular Grid')
    plt.plot(results['fibonacci_system'].time_points, results['fibonacci_system'].boundary_radius_history, 'r-',
             label='Fibonacci Grid')

    if results['sanctuary_time']:
        plt.axvline(x=results['sanctuary_time'], color='g', linestyle='--')

    plt.title('Event Horizon-like Boundary Evolution')
    plt.xlabel('Time')
    plt.ylabel('Boundary Radius')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('quantum_coherence_comparison.png', dpi=300)
    plt.show()


def plot_final_state_comparison(results):
    """Plot final states for both grid types"""
    plt.figure(figsize=(12, 6))

    # Create custom colormap with zero at white
    cmap = plt.cm.RdBu_r

    # Plot regular grid final state
    plt.subplot(1, 2, 1)
    regular_state = np.abs(results['regular_system'].state) ** 2
    plt.pcolormesh(results['regular_system'].X, results['regular_system'].Y, regular_state, cmap=cmap)
    plt.colorbar(label='Probability Density')
    plt.title(f'Regular Grid Final State\nCoherence: {results["regular_coherence"][-1]:.6f}')
    plt.xlabel('x')
    plt.ylabel('y')

    # Plot fibonacci grid final state
    plt.subplot(1, 2, 2)
    fibonacci_state = np.abs(results['fibonacci_system'].state) ** 2
    plt.pcolormesh(results['fibonacci_system'].X, results['fibonacci_system'].Y, fibonacci_state, cmap=cmap)
    plt.colorbar(label='Probability Density')
    plt.title(f'Fibonacci Grid Final State\nCoherence: {results["fib_coherence"][-1]:.6f}')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.tight_layout()
    plt.savefig('final_state_comparison.png', dpi=300)
    plt.show()


def plot_system_evolution_animation(results, grid_type='regular', num_frames=8):
    """Create a panel of system states at different time points"""
    system = results['regular_system'] if grid_type == 'regular' else results['fibonacci_system']

    plt.figure(figsize=(15, 10))

    # Create custom colormap with zero at white
    cmap = plt.cm.RdBu_r

    # Select time points to show
    time_indices = np.linspace(0, len(system.time_points) - 1, num_frames, dtype=int)

    for i, idx in enumerate(time_indices):
        plt.subplot(2, 4, i + 1)

        state = np.abs(system.state) ** 2 if i == num_frames - 1 else np.abs(
            QuantumSystem(params, grid_type=grid_type).state) ** 2

        plt.pcolormesh(system.X, system.Y, state, cmap=cmap)
        plt.colorbar(label='Probability Density')
        plt.title(f't = {system.time_points[idx]:.2f}\nCoherence: {system.coherence_history[idx]:.4f}')
        plt.xlabel('x')
        plt.ylabel('y')

    plt.tight_layout()
    plt.savefig(f'{grid_type}_evolution.png', dpi=300)
    plt.show()


def plot_golden_ratio_sensitivity(params):
    """Analyze sensitivity to deviations from the golden ratio"""
    # Range of scaling factors to test
    scaling_factors = np.linspace(1.0, 2.0, 21)
    coherence_ratios = []

    # Create temperature profile function
    temp_profile = lambda t: hiv_fever_profile(t, params.baseline_temp, params.fever_temp)

    # Run simulation with regular grid once
    regular_system = QuantumSystem(params, grid_type='regular')
    regular_system.run_simulation(temp_profile)

    # Run simulations for each scaling factor
    for phi in scaling_factors:
        print(f"Testing scaling factor φ = {phi:.3f}")

        # Modify the parameters
        test_params = SimulationParameters()
        test_params.phi = phi

        # Run simulation
        test_system = QuantumSystem(test_params, grid_type='fibonacci')
        test_system.run_simulation(temp_profile)

        # Calculate final coherence ratio
        final_coherence_ratio = test_system.coherence_history[-1] / regular_system.coherence_history[-1]
        coherence_ratios.append(final_coherence_ratio)

        print(f"Final coherence ratio: {final_coherence_ratio:.4f}")

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(scaling_factors, coherence_ratios, 'o-')
    plt.axvline(x=1.618, color='r', linestyle='--', label='Golden Ratio (φ = 1.618)')

    # Find and mark the maximum
    max_idx = np.argmax(coherence_ratios)
    max_phi = scaling_factors[max_idx]
    max_ratio = coherence_ratios[max_idx]

    plt.plot(max_phi, max_ratio, 'ro', markersize=10)
    plt.annotate(f'Max: {max_ratio:.2f} at φ = {max_phi:.3f}',
                 xy=(max_phi, max_ratio), xytext=(max_phi - 0.3, max_ratio + 20),
                 arrowprops=dict(facecolor='black', shrink=0.05))

    plt.title('Final Coherence Ratio vs. φ-Deviation')
    plt.xlabel('Scaling Factor (φ)')
    plt.ylabel('Fibonacci/Regular Coherence Ratio')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig('golden_ratio_sensitivity.png', dpi=300)
    plt.show()

    return scaling_factors, coherence_ratios


def run_extended_analysis(params, t_total=120.0):
    """Run extended analysis with multiplicative decay model as described in paper"""
    # Save original t_total
    original_t_total = params.t_total
    params.t_total = t_total

    # Define the extended decay rates based on paper
    reg_degradation_factor = 0.1  # Much faster decay for regular structures
    fib_degradation_factor = 0.001  # Much slower decay for Fibonacci structures

    # Time points
    n_steps = 1000
    times = np.linspace(0, t_total, n_steps)

    # Calculate coherence values
    # Regular grid: true exponential decay C(t) = C₀ · e^(-decay_factor·t)
    regular_coherence = np.exp(-reg_degradation_factor * times)

    # Fibonacci grid: multiplicative decay C(t) = C₀ · (1 - fib_degradation_factor)^t
    fibonacci_coherence = (1 - fib_degradation_factor) ** times

    # Calculate coherence ratio
    coherence_ratio = fibonacci_coherence / regular_coherence

    # Find the final ratio value
    final_ratio = coherence_ratio[-1]

    print(f"Extended analysis complete")
    print(f"Final coherence values:")
    print(f"Regular grid: {regular_coherence[-1]:.6e}")
    print(f"Fibonacci grid: {fibonacci_coherence[-1]:.6e}")
    print(f"Final advantage ratio: {final_ratio:.6e} ({np.log10(final_ratio):.2f} orders of magnitude)")

    # Restore original t_total
    params.t_total = original_t_total

    # Plot the results
    plt.figure(figsize=(12, 8))

    # Coherence comparison
    plt.subplot(2, 1, 1)
    plt.plot(times, regular_coherence, 'b-', label='Regular Grid')
    plt.plot(times, fibonacci_coherence, 'r-', label='Fibonacci Grid')
    plt.title('Extended Coherence Analysis (0-120 time units)')
    plt.xlabel('Time')
    plt.ylabel('Coherence')
    plt.legend()
    plt.grid(True)
    plt.yscale('log')

    # Coherence ratio
    plt.subplot(2, 1, 2)
    plt.plot(times, coherence_ratio, 'g-')
    plt.title(f'Coherence Ratio\nFinal advantage: {final_ratio:.2e}')
    plt.xlabel('Time')
    plt.ylabel('Fibonacci/Regular Ratio')
    plt.grid(True)
    plt.yscale('log')

    plt.tight_layout()
    plt.savefig('extended_analysis.png', dpi=300)
    plt.show()

    return times, regular_coherence, fibonacci_coherence, coherence_ratio


# Main execution
if __name__ == "__main__":
    print("Starting Fibonacci vs. Regular Grid Coherence Simulation Under HIV Inflammatory Conditions")

    # Run simulations and generate plots
    results = run_both_simulations(params)
    plot_coherence_comparison(results)
    plot_final_state_comparison(results)

    # Run Golden Ratio sensitivity analysis (optional - computationally intensive)
    # phi_values, coherence_ratios = plot_golden_ratio_sensitivity(params)

    # Run extended analysis
    # ext_times, ext_reg_coherence, ext_fib_coherence, ext_ratio = run_extended_analysis(params, t_total=120.0)

    print("Simulation complete!")