import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import sys
import datetime
import json
from scipy import integrate
from scipy import signal

# Force matplotlib to use the Agg backend (non-interactive)
plt.switch_backend('Agg')

# Define the absolute output directory path
output_dir = "/Users/acdmbp4/Desktop/Microtubule_Simulation/script/Core Scripts/comparison_results"
if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        sys.exit(1)

# Constants (same as in acute_hiv.py)
hbar = 1.0  # Reduced Planck's constant
m = 1.0  # Effective mass
L = 10.0  # Axial length of microtubule
R_inner = 7.0  # Inner radius of microtubule (arbitrary units)
R_outer = 12.5  # Outer radius of microtubule (arbitrary units)
N_r = 100  # Number of radial grid points
N_z = 100  # Number of axial grid points
dr = (R_outer - R_inner) / N_r  # Radial step size
dz = L / N_z  # Axial step size
dt = 0.01  # Time step size
time_steps = 300  # Total time steps

# Simulation parameters - MODIFIED for chronic HIV
V_0 = 4.5  # Peak cytokine potential (slightly lower than acute)
Gamma_0 = 0.04  # Baseline decoherence rate (lower than acute)
alpha_c = 0.08  # Scaling factor for cytokine-induced decoherence (lower than acute)

# Target time points for visualization
target_times = [0, 1.5, 2.5, 3.0]
target_frames = [int(t / dt) for t in target_times]
print(f"Target frames for visualization: {target_frames}")

# Create spatial grids
r = np.linspace(R_inner, R_outer, N_r)
z = np.linspace(0, L, N_z)
R, Z = np.meshgrid(r, z)  # 2D grid for visualization


# Fibonacci sequence generation
def generate_fibonacci_sequence(n):
    """Generate Fibonacci sequence up to n terms"""
    fib = [0, 1]
    for i in range(2, n):
        fib.append(fib[i - 1] + fib[i - 2])
    return np.array(fib)


# Generate and normalize Fibonacci sequence for scaling
fib_seq = generate_fibonacci_sequence(max(N_r, N_z) + 10)
fib_r = fib_seq[-(N_r):]  # Use last N_r elements
fib_z = fib_seq[-(N_z):]  # Use last N_z elements

# Normalize sequences to physical dimensions
fib_r_scaled = R_inner + (fib_r / np.max(fib_r)) * (R_outer - R_inner)
fib_z_scaled = (fib_z / np.max(fib_z)) * L

# Create Fibonacci-scaled grid
R_fib = np.zeros((N_z, N_r))
Z_fib = np.zeros((N_z, N_r))
for i in range(N_z):
    for j in range(N_r):
        R_fib[i, j] = fib_r_scaled[j]
        Z_fib[i, j] = fib_z_scaled[i]

# Initialize wavefunctions
sigma_z = L / 10
z0 = L / 2

# Regular wavefunction
Psi_reg = np.exp(-0.5 * ((Z - z0) / sigma_z) ** 2) * np.sin(np.pi * (R - R_inner) / (R_outer - R_inner))
Psi_reg /= np.sqrt(np.sum(np.abs(Psi_reg) ** 2) * dr * dz)

# Fibonacci wavefunction
Psi_fib = np.exp(-0.5 * ((Z_fib - z0) / sigma_z) ** 2) * np.sin(np.pi * (R_fib - R_inner) / (R_outer - R_inner))
Psi_fib /= np.sqrt(np.sum(np.abs(Psi_fib) ** 2) * dr * dz)

# Verify initial normalization
initial_norm_reg = np.sum(np.abs(Psi_reg) ** 2) * dr * dz
initial_norm_fib = np.sum(np.abs(Psi_fib) ** 2) * dr * dz
print(f"Initial normalization - Regular: {initial_norm_reg:.6f}, Fibonacci: {initial_norm_fib:.6f}")

# Store initial coherence values for adaptive threshold calculation
initial_max_prob_reg = np.max(np.abs(Psi_reg) ** 2)
initial_max_prob_fib = np.max(np.abs(Psi_fib) ** 2)

# Coherence history tracking for adaptive threshold
coherence_history_reg = []  # For tracking coherence decay in regular grid
coherence_history_fib = []  # For tracking coherence decay in Fibonacci grid
coherence_history_window = 20  # Window size for chronic phase (longer than acute)


# IMPROVEMENT: Custom initial cytokine concentration for chronic HIV
# Characterized by widespread, more diffuse inflammation
def initialize_chronic_cytokines():
    """
    Initialize cytokine field for chronic HIV with:
    1. More widespread, diffuse inflammation
    2. Lower intensity peaks but higher baseline
    3. Pattern consistent with chronic neuroinflammation
    """
    cytokines = np.zeros((N_z, N_r))

    # More anatomically accurate initialization
    # Diffuse distribution with deeper penetration characteristic of chronic state
    for i in range(N_z):
        for j in range(N_r):
            # Distance from outer boundary (normalized)
            outer_dist = (R_outer - r[j]) / (R_outer - R_inner)

            # Distance from center axially (normalized)
            center_dist = abs(z[i] - L / 2) / (L / 2)

            # CHRONIC PHASE: Has more widespread, diffuse cytokine distribution
            # Deeper penetration than acute phase
            if outer_dist < 0.5:  # Deeper penetration (50% of radius)
                # More uniform profile with less steep dropoff
                cytokines[i, j] = 0.6 * np.exp(-(center_dist ** 2) / 0.8)
            else:
                # Gentler decay toward inner regions
                # Characteristic of chronic phase with greater penetration
                cytokines[i, j] = 0.3 * np.exp(-((outer_dist - 0.5) / 0.5) ** 2) * np.exp(-(center_dist ** 2) / 1.0)

    # Add diffuse regions representing sustained inflammation
    # Chronic HIV has more widespread, less intense inflammatory regions
    hotspot_count = 4  # More hotspots than acute
    hotspot_locations = [
        (int(N_z * 0.3), int(N_r * 0.75)),  # Lower region
        (int(N_z * 0.5), int(N_r * 0.8)),  # Middle region
        (int(N_z * 0.7), int(N_r * 0.7)),  # Upper region
        (int(N_z * 0.4), int(N_r * 0.6))  # Inner region (deeper penetration)
    ]

    for h_i, h_j in hotspot_locations:
        hotspot_size = 10 + np.random.randint(5)  # Larger size characteristic of chronic phase
        hotspot_intensity = 0.6  # Lower intensity but more widespread
        for i in range(N_z):
            for j in range(N_r):
                dist_sq = ((i - h_i) / hotspot_size) ** 2 + ((j - h_j) / hotspot_size) ** 2
                if dist_sq < 1.0:
                    # Gentler intensity profile characteristic of chronic inflammation
                    intensity_factor = hotspot_intensity * np.exp(-1.5 * dist_sq)
                    cytokines[i, j] = min(1.0, cytokines[i, j] + intensity_factor)

    # Add a higher baseline level throughout the domain (chronic background inflammation)
    baseline = 0.15  # Higher baseline characteristic of chronic state
    cytokines = np.clip(cytokines + baseline, 0, 1)

    return cytokines


# Generate initial cytokine field for chronic HIV
C = initialize_chronic_cytokines()

# Potentials with cytokine influence
V_base = 5.0 * np.cos(2 * np.pi * Z / L)  # Base potential (tubulin periodicity)
V_walls = np.zeros_like(R)
V_walls[R < R_inner + dr] = 1e6  # Confinement at inner wall
V_walls[R > R_outer - dr] = 1e6  # Confinement at outer wall
V_reg = V_base + V_walls + V_0 * C  # Add cytokine potential
V_fib = V_base + V_walls + V_0 * C  # Add cytokine potential


# Modified decoherence term for chronic HIV
def calculate_decoherence(cytokine_field, baseline_rate, alpha, state="chronic"):
    """
    Calculate decoherence rate with exponential scaling and state-specific parameters
    Different HIV states have different decoherence characteristics
    """
    if state == "chronic":
        # Chronic HIV has less intense but more persistent decoherence effect
        alpha_modified = alpha * 0.9  # Lower impact compared to acute
        baseline_modified = baseline_rate * 1.1  # Higher baseline for chronic
        # More uniform decoherence with higher baseline - characteristic of chronic phase
        return baseline_modified * np.exp(alpha_modified * cytokine_field)
    else:
        # Default exponential relationship
        return baseline_rate * np.exp(alpha * cytokine_field)


# Adaptive event horizon threshold specific to chronic HIV
def calculate_event_horizon(Psi, initial_max_prob, coherence_level, grid_type="regular", min_threshold=0.05,
                            max_threshold=0.2):
    """
    Enhanced event horizon detection with adaptive threshold based on:
    1. Coherence level
    2. Grid type (regular vs Fibonacci)
    3. Specific adjustments for chronic HIV state
    """
    prob_density = np.abs(Psi) ** 2
    current_max_prob = np.max(prob_density)

    # Different threshold parameters for regular vs Fibonacci grids
    # Chronic HIV has more diffuse boundaries and slower but sustained decoherence
    if grid_type == "regular":
        # Regular grid needs lower threshold in chronic state due to more spread out patterns
        threshold_ratio = max(min_threshold, min(max_threshold, coherence_level * 0.18))
    else:
        # Fibonacci grid still performs better for boundary detection
        threshold_ratio = max(min_threshold, min(max_threshold, coherence_level * 0.15))

    threshold = threshold_ratio * current_max_prob

    # Calculate horizon for each z position
    horizon_r_indices = np.zeros(N_z, dtype=int)
    for i in range(N_z):
        # Find where probability drops below threshold
        # Moving from outer boundary inward
        for j in range(N_r - 1, 0, -1):
            if prob_density[i, j] > threshold:
                horizon_r_indices[i] = j
                break

    # Convert indices to radial positions
    return r[horizon_r_indices]


# Enhanced dispersion measures
def calculate_dispersion_metrics(wavefunction, R_grid, Z_grid, dr, dz):
    """
    Calculate multiple dispersion metrics beyond just variance:
    1. Spatial variance (standard)
    2. Entropy
    3. Kurtosis (measure of "tailedness")
    """
    prob = np.abs(wavefunction) ** 2

    # Calculate variance (standard dispersion measure)
    r_mean = np.sum(R_grid * prob) * dr * dz
    z_mean = np.sum(Z_grid * prob) * dr * dz
    var_r = np.sum((R_grid - r_mean) ** 2 * prob) * dr * dz
    var_z = np.sum((Z_grid - z_mean) ** 2 * prob) * dr * dz
    total_variance = var_r + var_z

    # Calculate entropy (information-theoretic measure of dispersion)
    # Add small epsilon to avoid log(0)
    epsilon = 1e-10
    entropy = -np.sum(prob * np.log(prob + epsilon)) * dr * dz

    # Calculate kurtosis (measure of the "tailedness" of the distribution)
    # Higher values indicate more extreme deviations
    r_kurtosis = np.sum((R_grid - r_mean) ** 4 * prob) * dr * dz / (var_r ** 2)
    z_kurtosis = np.sum((Z_grid - z_mean) ** 4 * prob) * dr * dz / (var_z ** 2)
    total_kurtosis = r_kurtosis + z_kurtosis

    return {
        'variance': total_variance,
        'entropy': entropy,
        'kurtosis': total_kurtosis
    }


# Time evolution function with decoherence and improved normalization
def evolve_wavefunction(Psi, V, Gamma, dr, dz, dt, R_grid):
    # Interior points for finite difference
    Psi_new = np.zeros_like(Psi, dtype=complex)

    for i in range(1, N_z - 1):
        for j in range(1, N_r - 1):
            # Second derivative in r
            d2r = (Psi[i, j + 1] - 2 * Psi[i, j] + Psi[i, j - 1]) / dr ** 2

            # First derivative in r (for cylindrical term)
            dr_term = (Psi[i, j + 1] - Psi[i, j - 1]) / (2 * dr)

            # Second derivative in z
            d2z = (Psi[i + 1, j] - 2 * Psi[i, j] + Psi[i - 1, j]) / dz ** 2

            # Cylindrical coordinate correction term
            cyl_term = dr_term / R_grid[i, j]

            # Laplacian in cylindrical coordinates
            laplacian = d2r + cyl_term + d2z

            # Time evolution step
            Psi_new[i, j] = Psi[i, j] + 1j * hbar * dt / (2 * m) * laplacian - 1j * dt / hbar * V[i, j] * Psi[i, j] - \
                            Gamma[i, j] * dt * Psi[i, j]

    # Apply boundary conditions
    Psi_new[0, :] = 0  # z = 0
    Psi_new[-1, :] = 0  # z = L
    Psi_new[:, 0] = 0  # r = R_inner
    Psi_new[:, -1] = 0  # r = R_outer

    # Enforce normalization at every step to prevent numerical drift
    norm = np.sqrt(np.sum(np.abs(Psi_new) ** 2) * dr * dz)
    if norm > 0:
        Psi_new /= norm
    else:
        print("WARNING: Zero norm encountered during evolution!")

    return Psi_new


# Initial decoherence rates with specific chronic phase parameters
Gamma_cytokine = calculate_decoherence(C, Gamma_0, alpha_c, state="chronic")

# Store probability densities and calculated metrics
Psi_reg_list = []
Psi_fib_list = []
horizon_reg_list = []
horizon_fib_list = []
variance_reg_list = []
variance_fib_list = []
entropy_reg_list = []
entropy_fib_list = []
kurtosis_reg_list = []
kurtosis_fib_list = []
coherence_reg_list = []
coherence_fib_list = []
cytokine_list = []
time_list = []
norm_check_list = []  # To track normalization drift

# Initial state
Psi_reg_prob = np.abs(Psi_reg) ** 2
Psi_fib_prob = np.abs(Psi_fib) ** 2
Psi_reg_list.append(Psi_reg_prob)
Psi_fib_list.append(Psi_fib_prob)
horizon_reg_list.append(calculate_event_horizon(Psi_reg, initial_max_prob_reg, 1.0, grid_type="regular"))
horizon_fib_list.append(calculate_event_horizon(Psi_fib, initial_max_prob_fib, 1.0, grid_type="fibonacci"))

# Store all dispersion metrics
disp_reg = calculate_dispersion_metrics(Psi_reg, R, Z, dr, dz)
disp_fib = calculate_dispersion_metrics(Psi_fib, R_fib, Z_fib, dr, dz)
variance_reg_list.append(disp_reg['variance'])
variance_fib_list.append(disp_fib['variance'])
entropy_reg_list.append(disp_reg['entropy'])
entropy_fib_list.append(disp_fib['entropy'])
kurtosis_reg_list.append(disp_reg['kurtosis'])
kurtosis_fib_list.append(disp_fib['kurtosis'])

coherence_reg_list.append(1.0)  # Initial state has perfect coherence
coherence_fib_list.append(1.0)  # Initial state has perfect coherence
coherence_history_reg.append(1.0)
coherence_history_fib.append(1.0)
cytokine_list.append(C.copy())
time_list.append(0)
norm_check_list.append((1.0, 1.0))  # Initial normalization check

# Save the initial state visualization immediately
try:
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig)

    # Regular grid probability density
    ax1 = fig.add_subplot(gs[0, 0])
    cont1 = ax1.contourf(Z, R, Psi_reg_prob, levels=50, cmap='viridis')
    ax1.plot(z, horizon_reg_list[-1], 'r--', linewidth=2, label='Event Horizon')
    ax1.set_title(f'Regular Grid (t=0.0)')
    ax1.set_xlabel('Axial Position (z)')
    ax1.set_ylabel('Radial Position (r)')
    ax1.legend()
    plt.colorbar(cont1, ax=ax1, label='Probability Density')

    # Fibonacci grid probability density
    ax2 = fig.add_subplot(gs[0, 1])
    cont2 = ax2.contourf(Z, R, Psi_fib_prob, levels=50, cmap='viridis')
    ax2.plot(z, horizon_fib_list[-1], 'r--', linewidth=2, label='Event Horizon')
    ax2.set_title(f'Fibonacci Grid (t=0.0)')
    ax2.set_xlabel('Axial Position (z)')
    ax2.set_ylabel('Radial Position (r)')
    ax2.legend()
    plt.colorbar(cont2, ax=ax2, label='Probability Density')

    # Cytokine concentration
    ax3 = fig.add_subplot(gs[0, 2])
    cont3 = ax3.contourf(Z, R, C, levels=20, cmap='plasma')
    ax3.set_title(f'Cytokine Field (t=0.0) - Chronic HIV')
    ax3.set_xlabel('Axial Position (z)')
    ax3.set_ylabel('Radial Position (r)')
    plt.colorbar(cont3, ax=ax3, label='Concentration')

    plt.tight_layout()
    save_path = os.path.join(output_dir, "chronic_hiv_t0.0.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Initial state visualization saved to {save_path}")
except Exception as e:
    print(f"Error creating initial visualization: {e}")

# Time evolution loop
print("Starting time evolution simulation for chronic HIV state...")
for step in range(1, time_steps + 1):
    current_time = step * dt

    # Enhanced cytokine dynamics for chronic HIV
    # Characterized by slower diffusion and decay than acute
    laplacian_C = np.zeros_like(C)
    for i in range(1, N_z - 1):
        for j in range(1, N_r - 1):
            laplacian_C[i, j] = (C[i + 1, j] + C[i - 1, j] + C[i, j + 1] + C[i, j - 1] - 4 * C[i, j]) / (dr * dz)

    # More realistic cytokine dynamics for chronic phase
    # Lower diffusion rate but slower decay rate (more persistent)
    diffusion_rate = 0.03  # Lower diffusion rate for chronic vs. acute

    # Almost constant decay rate characteristic of chronic inflammation
    # Chronic phase has more stable, persistent inflammation
    base_decay_rate = 0.01  # Lower decay rate than acute phase
    time_factor = 0.02 * current_time  # Slower increase over time
    decay_rate = base_decay_rate * (1 + time_factor)

    # Small random fluctuations to model chronic dynamics
    fluctuation_amplitude = 0.005  # Smaller fluctuations in chronic phase
    fluctuations = fluctuation_amplitude * (np.random.random(C.shape) - 0.5)

    # In chronic phase, cytokines are more stable with slower changes
    C = C + dt * (diffusion_rate * laplacian_C - decay_rate * C + fluctuations)
    C = np.clip(C, 0, 1)  # Keep within bounds

    # Update potentials and decoherence rates based on new cytokine field
    V_reg = V_base + V_walls + V_0 * C
    V_fib = V_base + V_walls + V_0 * C

    # State-specific decoherence calculation
    Gamma_cytokine = calculate_decoherence(C, Gamma_0, alpha_c, state="chronic")

    # Periodic normalization check
    pre_norm_reg = np.sum(np.abs(Psi_reg) ** 2) * dr * dz
    pre_norm_fib = np.sum(np.abs(Psi_fib) ** 2) * dr * dz

    # Evolve wavefunctions
    Psi_reg = evolve_wavefunction(Psi_reg, V_reg, Gamma_cytokine, dr, dz, dt, R)
    Psi_fib = evolve_wavefunction(Psi_fib, V_fib, Gamma_cytokine, dr, dz, dt, R_fib)

    # Check normalization after evolution
    post_norm_reg = np.sum(np.abs(Psi_reg) ** 2) * dr * dz
    post_norm_fib = np.sum(np.abs(Psi_fib) ** 2) * dr * dz

    # Calculate probability densities
    Psi_reg_prob = np.abs(Psi_reg) ** 2
    Psi_fib_prob = np.abs(Psi_fib) ** 2

    # Calculate coherence (correlation with initial state)
    initial_reg = Psi_reg_list[0]
    initial_fib = Psi_fib_list[0]
    coh_reg = np.sum(np.sqrt(initial_reg * Psi_reg_prob)) * dr * dz
    coh_fib = np.sum(np.sqrt(initial_fib * Psi_fib_prob)) * dr * dz

    # Update rolling average of coherence for adaptive threshold
    coherence_history_reg.append(coh_reg)
    coherence_history_fib.append(coh_fib)

    if len(coherence_history_reg) > coherence_history_window:
        coherence_history_reg.pop(0)
        coherence_history_fib.pop(0)

    rolling_coherence_reg = sum(coherence_history_reg) / len(coherence_history_reg)
    rolling_coherence_fib = sum(coherence_history_fib) / len(coherence_history_fib)

    # Store results periodically
    if step % 10 == 0 or step in target_frames:
        Psi_reg_list.append(Psi_reg_prob)
        Psi_fib_list.append(Psi_fib_prob)

        # Grid-specific threshold calculations
        horizon_reg_list.append(calculate_event_horizon(Psi_reg, initial_max_prob_reg,
                                                        rolling_coherence_reg, grid_type="regular"))
        horizon_fib_list.append(calculate_event_horizon(Psi_fib, initial_max_prob_fib,
                                                        rolling_coherence_fib, grid_type="fibonacci"))

        # Store comprehensive dispersion metrics
        disp_reg = calculate_dispersion_metrics(Psi_reg, R, Z, dr, dz)
        disp_fib = calculate_dispersion_metrics(Psi_fib, R_fib, Z_fib, dr, dz)
        variance_reg_list.append(disp_reg['variance'])
        variance_fib_list.append(disp_fib['variance'])
        entropy_reg_list.append(disp_reg['entropy'])
        entropy_fib_list.append(disp_fib['entropy'])
        kurtosis_reg_list.append(disp_reg['kurtosis'])
        kurtosis_fib_list.append(disp_fib['kurtosis'])

        coherence_reg_list.append(coh_reg)
        coherence_fib_list.append(coh_fib)
        cytokine_list.append(C.copy())
        time_list.append(current_time)
        norm_check_list.append((post_norm_reg, post_norm_fib))

        # Print progress update
        if step % 50 == 0:
            print(f"Completed step {step}/{time_steps} ({step / time_steps * 100:.1f}%)")
            print(f"  Time: {current_time:.2f}")
            print(f"  Standard variance: {disp_reg['variance']:.4f}")
            print(f"  Fibonacci variance: {disp_fib['variance']:.4f}")
            print(f"  Coherence standard: {coh_reg:.4f}")
            print(f"  Coherence fibonacci: {coh_fib:.4f}")
            print(f"  Rolling coherence reg: {rolling_coherence_reg:.4f}")
            print(f"  Rolling coherence fib: {rolling_coherence_fib:.4f}")
            print(f"  Normalization (reg, fib): ({post_norm_reg:.6f}, {post_norm_fib:.6f})")

    # Generate visualizations at target time points
    if step in target_frames:
        t_index = target_frames.index(step)
        target_time = target_times[t_index]

        try:
            print(f"Creating visualization for t = {target_time}")

            fig = plt.figure(figsize=(18, 12))
            gs = GridSpec(2, 3, figure=fig)

            # Regular grid probability density
            ax1 = fig.add_subplot(gs[0, 0])
            cont1 = ax1.contourf(Z, R, Psi_reg_prob, levels=50, cmap='viridis')
            ax1.plot(z, horizon_reg_list[-1], 'r--', linewidth=2, label='Event Horizon')
            ax1.set_title(f'Regular Grid (t={target_time:.1f})')
            ax1.set_xlabel('Axial Position (z)')
            ax1.set_ylabel('Radial Position (r)')
            ax1.legend()
            plt.colorbar(cont1, ax=ax1, label='Probability Density')

            # Fibonacci grid probability density
            ax2 = fig.add_subplot(gs[0, 1])
            cont2 = ax2.contourf(Z, R, Psi_fib_prob, levels=50, cmap='viridis')
            ax2.plot(z, horizon_fib_list[-1], 'r--', linewidth=2, label='Event Horizon')
            ax2.set_title(f'Fibonacci Grid (t={target_time:.1f})')
            ax2.set_xlabel('Axial Position (z)')
            ax2.set_ylabel('Radial Position (r)')
            ax2.legend()
            plt.colorbar(cont2, ax=ax2, label='Probability Density')

            # Cytokine concentration
            ax3 = fig.add_subplot(gs[0, 2])
            cont3 = ax3.contourf(Z, R, C, levels=20, cmap='plasma')
            ax3.set_title(f'Cytokine Field (t={target_time:.1f}) - Chronic HIV')
            ax3.set_xlabel('Axial Position (z)')
            ax3.set_ylabel('Radial Position (r)')
            plt.colorbar(cont3, ax=ax3, label='Concentration')

            # Coherence comparison
            ax4 = fig.add_subplot(gs[1, 0:2])
            times = np.array(time_list)
            coh_reg_array = np.array(coherence_reg_list)
            coh_fib_array = np.array(coherence_fib_list)
            ax4.plot(times, coh_reg_array, 'b-', label='Regular Grid')
            ax4.plot(times, coh_fib_array, 'r-', label='Fibonacci Grid')
            ax4.set_xlabel('Time')
            ax4.set_ylabel('Coherence Measure')
            ax4.set_title('Coherence')
            ax4.axvline(x=target_time, color='k', linestyle='--')
            ax4.legend()
            ax4.grid(True)

            # Dispersion comparison (variance)
            ax5 = fig.add_subplot(gs[1, 2])
            var_reg_array = np.array(variance_reg_list)
            var_fib_array = np.array(variance_fib_list)
            ax5.plot(times, var_reg_array, 'b-', label='Regular Grid')
            ax5.plot(times, var_fib_array, 'r-', label='Fibonacci Grid')
            ax5.set_xlabel('Time')
            ax5.set_ylabel('Spatial Variance')
            ax5.set_title('Wavefunction Dispersion')
            ax5.axvline(x=target_time, color='k', linestyle='--')
            ax5.legend()
            ax5.grid(True)

            plt.tight_layout()
            save_path = os.path.join(output_dir, f"chronic_hiv_t{target_time:.1f}.png")
            plt.savefig(save_path, dpi=300)
            plt.close()
            print(f"Visualization saved to {save_path}")

        except Exception as e:
            print(f"Error creating visualization for t={target_time}: {e}")

# Final analysis and metrics calculation
print("\nSimulation completed. Calculating final metrics...")

# Calculate integrated coherence preservation (area under coherence curve)
times_array = np.array(time_list)
coherence_reg_array = np.array(coherence_reg_list)
coherence_fib_array = np.array(coherence_fib_list)

# Trapezoidal integration of coherence curves
integrated_coherence_reg = integrate.trapz(coherence_reg_array, times_array)
integrated_coherence_fib = integrate.trapz(coherence_fib_array, times_array)
coherence_advantage = (integrated_coherence_fib / integrated_coherence_reg - 1) * 100

# Calculate mean event horizon radii
horizon_reg_array = np.array(horizon_reg_list)
horizon_fib_array = np.array(horizon_fib_list)
mean_horizon_reg = np.mean(horizon_reg_array, axis=1)
mean_horizon_fib = np.mean(horizon_fib_array, axis=1)

# Calculate final dispersion statistics
final_var_reg = variance_reg_list[-1]
final_var_fib = variance_fib_list[-1]
final_entropy_reg = entropy_reg_list[-1]
final_entropy_fib = entropy_fib_list[-1]
final_kurtosis_reg = kurtosis_reg_list[-1]
final_kurtosis_fib = kurtosis_fib_list[-1]

# Calculate spectral analysis of coherence fluctuations
# This can reveal characteristic frequencies in the decoherence process
freq_reg, psd_reg = signal.welch(coherence_reg_array, 1 / dt, nperseg=min(len(coherence_reg_array), 128))
freq_fib, psd_fib = signal.welch(coherence_fib_array, 1 / dt, nperseg=min(len(coherence_fib_array), 128))

# Find dominant frequencies
dom_freq_reg = freq_reg[np.argmax(psd_reg)]
dom_freq_fib = freq_fib[np.argmax(psd_fib)]

# Phase space trajectory analysis
# This can reveal attractor-like behavior in the quantum dynamics
phase_space_reg = []
phase_space_fib = []

for i in range(1, len(variance_reg_list)):
    phase_space_reg.append((variance_reg_list[i - 1], variance_reg_list[i]))
    phase_space_fib.append((variance_fib_list[i - 1], variance_fib_list[i]))

# Generate summary figures for the entire simulation
try:
    # Figure 1: Main summary plot
    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(3, 3, figure=fig)

    # Coherence comparison over time
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax1.plot(time_list, coherence_reg_list, 'b-', label='Regular Grid')
    ax1.plot(time_list, coherence_fib_list, 'r-', label='Fibonacci Grid')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Coherence Measure')
    ax1.set_title('Quantum Coherence Preservation - Chronic HIV')
    ax1.legend()
    ax1.grid(True)

    # Variance/dispersion over time
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.plot(time_list, variance_reg_list, 'b-', label='Regular Grid')
    ax2.plot(time_list, variance_fib_list, 'r-', label='Fibonacci Grid')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Spatial Variance')
    ax2.set_title('Wavefunction Dispersion')
    ax2.legend()
    ax2.grid(True)

    # Entropy over time
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(time_list, entropy_reg_list, 'b-', label='Regular Grid')
    ax3.plot(time_list, entropy_fib_list, 'r-', label='Fibonacci Grid')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Entropy')
    ax3.set_title('Information Entropy')
    ax3.legend()
    ax3.grid(True)

    # Kurtosis over time
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(time_list, kurtosis_reg_list, 'b-', label='Regular Grid')
    ax4.plot(time_list, kurtosis_fib_list, 'r-', label='Fibonacci Grid')
    ax4.set_xlabel('Time')
    ax4.set_ylabel('Kurtosis')
    ax4.set_title('Distribution Kurtosis')
    ax4.legend()
    ax4.grid(True)

    # Plot mean event horizon radius over time
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.plot(time_list, [np.mean(h) for h in horizon_reg_list], 'b-', label='Regular Grid')
    ax5.plot(time_list, [np.mean(h) for h in horizon_fib_list], 'r-', label='Fibonacci Grid')
    ax5.set_xlabel('Time')
    ax5.set_ylabel('Mean Radius')
    ax5.set_title('Event Horizon Evolution')
    ax5.legend()
    ax5.grid(True)

    # Power spectral density of coherence fluctuations
    ax6 = fig.add_subplot(gs[2, 0])
    ax6.semilogy(freq_reg, psd_reg, 'b-', label='Regular Grid')
    ax6.semilogy(freq_fib, psd_fib, 'r-', label='Fibonacci Grid')
    ax6.set_xlabel('Frequency (Hz)')
    ax6.set_ylabel('PSD')
    ax6.set_title('Spectral Analysis of Coherence')
    ax6.legend()
    ax6.grid(True)

    # Phase space plot
    ax7 = fig.add_subplot(gs[2, 1])
    ax7.plot([p[0] for p in phase_space_reg], [p[1] for p in phase_space_reg], 'b.-', label='Regular Grid')
    ax7.plot([p[0] for p in phase_space_fib], [p[1] for p in phase_space_fib], 'r.-', label='Fibonacci Grid')
    ax7.set_xlabel('Variance(t)')
    ax7.set_ylabel('Variance(t+1)')
    ax7.set_title('Phase Space Trajectory')
    ax7.legend()
    ax7.grid(True)

    # Normalization drift check
    ax8 = fig.add_subplot(gs[2, 2])
    norm_reg = [n[0] for n in norm_check_list]
    norm_fib = [n[1] for n in norm_check_list]
    ax8.plot(time_list, norm_reg, 'b-', label='Regular Grid')
    ax8.plot(time_list, norm_fib, 'r-', label='Fibonacci Grid')
    ax8.set_xlabel('Time')
    ax8.set_ylabel('Normalization')
    ax8.set_title('Wavefunction Normalization')
    ax8.axhline(y=1.0, color='k', linestyle='--')
    ax8.legend()
    ax8.grid(True)

    plt.tight_layout()
    summary_path = os.path.join(output_dir, "chronic_hiv_summary.png")
    plt.savefig(summary_path, dpi=300)
    plt.close()
    print(f"Summary visualization saved to {summary_path}")

    # Figure 2: Final state comparison
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig)

    # Final regular grid probability density
    ax1 = fig.add_subplot(gs[0, 0])
    cont1 = ax1.contourf(Z, R, Psi_reg_list[-1], levels=50, cmap='viridis')
    ax1.plot(z, horizon_reg_list[-1], 'r--', linewidth=2, label='Event Horizon')
    ax1.set_title(f'Final Regular Grid (t={time_list[-1]:.1f})')
    ax1.set_xlabel('Axial Position (z)')
    ax1.set_ylabel('Radial Position (r)')
    ax1.legend()
    plt.colorbar(cont1, ax=ax1, label='Probability Density')

    # Final Fibonacci grid probability density
    ax2 = fig.add_subplot(gs[0, 1])
    cont2 = ax2.contourf(Z, R, Psi_fib_list[-1], levels=50, cmap='viridis')
    ax2.plot(z, horizon_fib_list[-1], 'r--', linewidth=2, label='Event Horizon')
    ax2.set_title(f'Final Fibonacci Grid (t={time_list[-1]:.1f})')
    ax2.set_xlabel('Axial Position (z)')
    ax2.set_ylabel('Radial Position (r)')
    ax2.legend()
    plt.colorbar(cont2, ax=ax2, label='Probability Density')

    # Final cytokine concentration
    ax3 = fig.add_subplot(gs[0, 2])
    cont3 = ax3.contourf(Z, R, cytokine_list[-1], levels=20, cmap='plasma')
    ax3.set_title(f'Final Cytokine Field (t={time_list[-1]:.1f}) - Chronic HIV')
    ax3.set_xlabel('Axial Position (z)')
    ax3.set_ylabel('Radial Position (r)')
    plt.colorbar(cont3, ax=ax3, label='Concentration')

    # Initial vs Final probability profile along axial midpoint
    ax4 = fig.add_subplot(gs[1, 0:])
    mid_r = N_r // 2
    ax4.plot(z, Psi_reg_list[0][:, mid_r], 'b--', label='Initial Regular')
    ax4.plot(z, Psi_fib_list[0][:, mid_r], 'r--', label='Initial Fibonacci')
    ax4.plot(z, Psi_reg_list[-1][:, mid_r], 'b-', label='Final Regular')
    ax4.plot(z, Psi_fib_list[-1][:, mid_r], 'r-', label='Final Fibonacci')
    ax4.set_xlabel('Axial Position (z)')
    ax4.set_ylabel('Probability Density')
    ax4.set_title('Initial vs Final State Comparison (Radial Midpoint)')
    ax4.legend()
    ax4.grid(True)

    plt.tight_layout()
    final_path = os.path.join(output_dir, "chronic_hiv_final_state.png")
    plt.savefig(final_path, dpi=300)
    plt.close()
    print(f"Final state visualization saved to {final_path}")

except Exception as e:
    print(f"Error creating summary visualizations: {e}")

# Save numerical results as JSON for later analysis and comparison
results = {
    'simulation_parameters': {
        'hbar': hbar,
        'm': m,
        'L': L,
        'R_inner': R_inner,
        'R_outer': R_outer,
        'V_0': V_0,
        'Gamma_0': Gamma_0,
    'alpha_c': alpha_c,
    'IL6_coupling_strength': 0.28,
    'time_steps': time_steps,
        'dt': dt,
        'total_time': time_steps * dt
    },
    'coherence_metrics': {
        'final_coherence_regular': coherence_reg_list[-1],
        'final_coherence_fibonacci': coherence_fib_list[-1],
        'integrated_coherence_regular': float(integrated_coherence_reg),
        'integrated_coherence_fibonacci': float(integrated_coherence_fib),
        'coherence_advantage_percent': float(coherence_advantage),
        'dominant_frequency_regular': float(dom_freq_reg),
        'dominant_frequency_fibonacci': float(dom_freq_fib)
    },
    'dispersion_metrics': {
        'final_variance_regular': float(final_var_reg),
        'final_variance_fibonacci': float(final_var_fib),
        'final_entropy_regular': float(final_entropy_reg),
        'final_entropy_fibonacci': float(final_entropy_fib),
        'final_kurtosis_regular': float(final_kurtosis_reg),
        'final_kurtosis_fibonacci': float(final_kurtosis_fib)
    },
    'event_horizon': {
        'final_mean_radius_regular': float(np.mean(horizon_reg_list[-1])),
        'final_mean_radius_fibonacci': float(np.mean(horizon_fib_list[-1]))
    },
    'state': "chronic_hiv",
    'timestamp': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "final_coherence_fibonacci": coherence_fib_list[-1],
    "final_coherence_regular": coherence_reg_list[-1],
    "half_life_fibonacci": float(coherence_fib_array.sum() / len(coherence_fib_array)),
    "half_life_regular": float(coherence_reg_array.sum() / len(coherence_reg_array)),
    "auc_fibonacci": float(integrated_coherence_fib),
    "auc_regular": float(integrated_coherence_reg)
}

# Save results to JSON file
json_path = os.path.join(output_dir, "chronic_hiv_results.json")
try:
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=4)
    print(f"Results saved to {json_path}")
except Exception as e:
    print(f"Error saving results to JSON: {e}")

# Save coherence & variance arrays to .npz for plotting
timestamp_str = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
npz_path = os.path.join(output_dir, "chronic_hiv_data.npz")
np.savez(npz_path,
         time=np.array(time_list),
         coherence_reg=np.array(coherence_reg_list),
         coherence_fib=np.array(coherence_fib_list),
         variance_reg=np.array(variance_reg_list),
         variance_fib=np.array(variance_fib_list))
print(f"✅ Saved .npz data for visualization: {npz_path}")
print("\nChronic HIV simulation completed successfully!")
print(f"Final coherence - Regular: {coherence_reg_list[-1]:.4f}, Fibonacci: {coherence_fib_list[-1]:.4f}")
print(f"Coherence advantage: {coherence_advantage:.2f}%")
print(f"Final variance - Regular: {final_var_reg:.4f}, Fibonacci: {final_var_fib:.4f}")
print(f"Dominant frequency - Regular: {dom_freq_reg:.4f} Hz, Fibonacci: {dom_freq_fib:.4f} Hz")
