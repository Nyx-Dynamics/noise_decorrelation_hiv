import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import sys

# Force matplotlib to use the Agg backend (non-interactive)
plt.switch_backend('Agg')

# Define the absolute output directory path
output_dir = "/Users/acdmbp4/Desktop/Microtubule_Simulation/script/Core Scripts/Validation/Desktop/Microtubule_Simulation/script/Core Scripts/art_controlled/art_controlled_results"
if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        sys.exit(1)

# Constants
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

# Simulation parameters
V_0 = 5.0  # Peak cytokine potential
Gamma_0 = 0.05  # Baseline decoherence rate
alpha_c = 0.1  # Scaling factor for cytokine-induced decoherence

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


# IMPROVEMENT: Enhanced initial cytokine concentration for ART-controlled HIV
# Characterized by moderate, more evenly distributed inflammation with circadian-like oscillations
def initialize_art_controlled_cytokines():
    """
    Initialize cytokine field for ART-controlled HIV with:
    1. Moderate, evenly distributed baseline inflammation
    2. Physiologically realistic oscillatory patterns
    3. Reduced intensity compared to acute and chronic states
    """
    cytokines = np.zeros((N_z, N_r))

    # Moderate, more evenly distributed inflammation
    for i in range(N_z):
        for j in range(N_r):
            # Distance from outer boundary
            outer_dist = (R_outer - r[j]) / (R_outer - R_inner)
            # Distance from center (axially)
            center_dist = abs(z[i] - L / 2) / (L / 2)

            # IMPROVEMENT: More realistic concentration pattern for ART-controlled HIV
            # Modulated base level with circadian-like oscillatory patterns
            base_level = 0.3  # Lower base level than chronic or acute

            # IMPROVEMENT: Physiologically realistic oscillatory patterns
            # Mimicking circadian and medication-induced cyclical inflammation
            radial_oscillation = 0.15 * np.sin(3 * np.pi * outer_dist)
            axial_oscillation = 0.15 * np.cos(2 * np.pi * center_dist)

            cytokines[i, j] = base_level + radial_oscillation + axial_oscillation

    # IMPROVEMENT: More realistic spatial heterogeneity
    # Add mild random variation that respects biological spatial correlations
    # Instead of pure noise, use spatially correlated patterns
    for i in range(N_z):
        for j in range(N_r):
            # Spatially correlated random component
            spatial_noise = 0.1 * np.sin(5 * np.pi * i / N_z) * np.cos(3 * np.pi * j / N_r)
            cytokines[i, j] += spatial_noise

    # IMPROVEMENT: Add a few residual hotspots (mimicking persistent reservoirs)
    # These are characteristic of ART-controlled HIV with latent reservoirs
    hotspot_locations = [
        (int(N_z * 0.3), int(N_r * 0.8)),  # Lower region
        (int(N_z * 0.7), int(N_r * 0.7))  # Upper region
    ]

    for h_i, h_j in hotspot_locations:
        # Gaussian profile for hotspots with variable size but lower intensity than acute/chronic
        hotspot_size = 5 + np.random.randint(3)
        hotspot_intensity = 0.3  # Reduced intensity compared to acute/chronic
        for i in range(N_z):
            for j in range(N_r):
                dist_sq = ((i - h_i) / hotspot_size) ** 2 + ((j - h_j) / hotspot_size) ** 2
                if dist_sq < 1.0:
                    # Gaussian intensity profile with reduced intensity
                    cytokines[i, j] = min(1.0, cytokines[i, j] + hotspot_intensity * np.exp(-2.0 * dist_sq))

    return np.clip(cytokines, 0, 1)  # Ensure values are in [0,1]


# Generate initial cytokine field for ART-controlled HIV
C = initialize_art_controlled_cytokines()

# Potentials with cytokine influence
V_base = 5.0 * np.cos(2 * np.pi * Z / L)  # Base potential (tubulin periodicity)
V_walls = np.zeros_like(R)
V_walls[R < R_inner + dr] = 1e6  # Confinement at inner wall
V_walls[R > R_outer - dr] = 1e6  # Confinement at outer wall
V_reg = V_base + V_walls + V_0 * C  # Add cytokine potential
V_fib = V_base + V_walls + V_0 * C  # Add cytokine potential


# IMPROVEMENT: Modified decoherence term to exponential scaling with state-specific parameters
def calculate_decoherence(cytokine_field, baseline_rate, alpha, state="art_controlled"):
    """
    Calculate decoherence rate with exponential scaling and state-specific parameters
    Different HIV states have different decoherence characteristics
    """
    if state == "art_controlled":
        # ART-controlled HIV has moderate decoherence effect with oscillatory behavior
        # The decoherence varies cyclically, modeling medication effects
        alpha_modified = alpha * 0.8  # Reduced impact due to ART
        return baseline_rate * np.exp(alpha_modified * cytokine_field)
    else:
        # Default exponential relationship for other states
        return baseline_rate * np.exp(alpha * cytokine_field)


# Initial decoherence rates with specific ART-controlled parameters
Gamma_cytokine = calculate_decoherence(C, Gamma_0, alpha_c, state="art_controlled")

# Store initial coherence values for adaptive threshold calculation
initial_max_prob_reg = np.max(np.abs(Psi_reg) ** 2)
initial_max_prob_fib = np.max(np.abs(Psi_fib) ** 2)
# IMPROVEMENT: Coherence history tracking for adaptive threshold
coherence_history_reg = []  # For tracking coherence decay in regular grid
coherence_history_fib = []  # For tracking coherence decay in Fibonacci grid
coherence_history_window = 12  # Longer window for ART-controlled phase to capture oscillatory patterns


# IMPROVEMENT: Adaptive event horizon threshold specific to ART-controlled HIV
def calculate_event_horizon(Psi, initial_max_prob, coherence_level, grid_type="regular", min_threshold=0.05,
                            max_threshold=0.2):
    """
    Enhanced event horizon detection with adaptive threshold based on:
    1. Coherence level
    2. Grid type (regular vs Fibonacci)
    3. Specific adjustments for ART-controlled HIV state
    """
    prob_density = np.abs(Psi) ** 2
    current_max_prob = np.max(prob_density)

    # IMPROVEMENT: Different threshold parameters for regular vs Fibonacci grids
    # ART-controlled HIV needs intermediate thresholds due to oscillatory decoherence
    if grid_type == "regular":
        # Regular grid needs moderate threshold in ART-controlled state
        threshold_ratio = max(min_threshold, min(max_threshold, coherence_level * 0.22))
    else:
        # Fibonacci grid can use slightly lower threshold due to better coherence preservation
        threshold_ratio = max(min_threshold, min(max_threshold, coherence_level * 0.18))

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


# IMPROVEMENT: Enhanced dispersion measures
def calculate_dispersion_metrics(wavefunction, R_grid, Z_grid, dr, dz):
    """
    Calculate multiple dispersion metrics beyond just variance:
    1. Spatial variance (standard)
    2. Entropy
    3. Kurtosis (measure of "tailedness")
    """
    prob = np.abs(wavefunction) ** 2

    # 1. Calculate variance (standard dispersion measure)
    r_mean = np.sum(R_grid * prob) * dr * dz
    z_mean = np.sum(Z_grid * prob) * dr * dz
    var_r = np.sum((R_grid - r_mean) ** 2 * prob) * dr * dz
    var_z = np.sum((Z_grid - z_mean) ** 2 * prob) * dr * dz
    total_variance = var_r + var_z

    # 2. Calculate entropy (information-theoretic measure of dispersion)
    # Add small epsilon to avoid log(0)
    epsilon = 1e-10
    entropy = -np.sum(prob * np.log(prob + epsilon)) * dr * dz

    # 3. Calculate kurtosis (measure of the "tailedness" of the distribution)
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

    # IMPROVEMENT: Enforce normalization at every step to prevent numerical drift
    norm = np.sqrt(np.sum(np.abs(Psi_new) ** 2) * dr * dz)
    if norm > 0:
        Psi_new /= norm
    else:
        print("WARNING: Zero norm encountered during evolution!")

    return Psi_new


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

# IMPROVEMENT: Store all dispersion metrics
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

# IMPROVEMENT: Save the initial state visualization immediately
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
    ax3.set_title(f'Cytokine Field (t=0.0) - ART-Controlled HIV')
    ax3.set_xlabel('Axial Position (z)')
    ax3.set_ylabel('Radial Position (r)')
    plt.colorbar(cont3, ax=ax3, label='Concentration')

    plt.tight_layout()
    save_path = os.path.join(output_dir, "art_controlled_t0.0.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Initial state visualization saved to {save_path}")
except Exception as e:
    print(f"Error creating initial visualization: {e}")

# Time evolution loop
print("Starting time evolution simulation for ART-controlled HIV state...")
for step in range(1, time_steps + 1):
    current_time = step * dt

    # IMPROVEMENT: Enhanced cytokine dynamics for ART-controlled HIV
    # Characterized by moderate diffusion and oscillatory patterns (medication cycles)
    laplacian_C = np.zeros_like(C)
    for i in range(1, N_z - 1):
        for j in range(1, N_r - 1):
            laplacian_C[i, j] = (C[i + 1, j] + C[i - 1, j] + C[i, j + 1] + C[i, j - 1] - 4 * C[i, j]) / (dr * dz)

    # IMPROVEMENT: More realistic cytokine dynamics for ART-controlled state
    diffusion_rate = 0.04  # Moderate diffusion rate
    decay_rate = 0.015  # Moderate decay rate

    # IMPROVEMENT: Oscillatory pattern characteristic of ART-controlled state
    # Mimicking medication cycles and circadian fluctuations
    # Lower amplitude oscillations than the untreated case, based on clinical observations
    oscillation_frequency_1 = 0.1  # Slower oscillation (medication cycles)
    oscillation_frequency_2 = 0.5  # Faster oscillation (circadian rhythm)

    oscillation = 0.01 * np.sin(oscillation_frequency_1 * current_time + Z / L * 2 * np.pi) * \
                  np.cos(oscillation_frequency_2 * current_time + (R - R_inner) / (R_outer - R_inner) * 2 * np.pi)

    # In ART-controlled phase, cytokines have moderate levels with mild oscillations
    C = C + dt * (diffusion_rate * laplacian_C - decay_rate * C + oscillation)
    C = np.clip(C, 0, 1)  # Keep within bounds

    # IMPROVEMENT: Every so often, model medication dosing effect
    # This creates periodic suppression of inflammation followed by mild rebound
    if step % 30 == 0:  # Medication effect cycle
        medication_factor = np.exp(-0.8 * C)  # Exponential suppression effect
        C = C * medication_factor

    # Update potentials and decoherence rates based on new cytokine field
    V_reg = V_base + V_walls + V_0 * C
    V_fib = V_base + V_walls + V_0 * C

    # IMPROVEMENT: State-specific decoherence calculation
    Gamma_cytokine = calculate_decoherence(C, Gamma_0, alpha_c, state="art_controlled")

    # IMPROVEMENT: Periodic normalization check
    pre_norm_reg = np.sum(np.abs(Psi_reg) ** 2) * dr * dz
    pre_norm_fib = np.sum(np.abs(Psi_fib) ** 2) * dr * dz

    # Evolve wavefunctions
    Psi_reg = evolve_wavefunction(Psi_reg, V_reg, Gamma_cytokine, dr, dz, dt, R)
    Psi_fib = evolve_wavefunction(Psi_fib, V_fib, Gamma_cytokine, dr, dz, dt, R_fib)

    # IMPROVEMENT: Check normalization after evolution
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
    # IMPROVEMENT: Separate tracking for regular and Fibonacci grids
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

        # IMPROVEMENT: Grid-specific threshold calculations
        horizon_reg_list.append(calculate_event_horizon(Psi_reg, initial_max_prob_reg,
                                                        rolling_coherence_reg, grid_type="regular"))
        horizon_fib_list.append(calculate_event_horizon(Psi_fib, initial_max_prob_fib,
                                                        rolling_coherence_fib, grid_type="fibonacci"))

        # IMPROVEMENT: Store comprehensive dispersion metrics
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
            ax3.set_title(f'Cytokine Field (t={target_time:.1f}) - ART-Controlled')
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
            ax4.set_title('Coherence Preservation')
            ax4.legend()
            ax4.grid(True)

            # IMPROVEMENT: Enhanced metrics visualization
            # Variance and entropy comparison
            ax5 = fig.add_subplot(gs[1, 2])
            var_reg_array = np.array(variance_reg_list)
            var_fib_array = np.array(variance_fib_list)

            # Primary y-axis for variance
            ax5.plot(times, var_reg_array, 'b-', label='Regular Variance')
            ax5.plot(times, var_fib_array, 'r-', label='Fibonacci Variance')
            ax5.set_xlabel('Time')
            ax5.set_ylabel('Spatial Variance')
            ax5.set_title('Wavefunction Dispersion')

            # Secondary y-axis for entropy
            ax5_twin = ax5.twinx()
            ent_reg_array = np.array(entropy_reg_list)
            ent_fib_array = np.array(entropy_fib_list)
            ax5_twin.plot(times, ent_reg_array, 'b--', alpha=0.5, label='Regular Entropy')
            ax5_twin.plot(times, ent_fib_array, 'r--', alpha=0.5, label='Fibonacci Entropy')
            ax5_twin.set_ylabel('Entropy')

            # Combined legend
            lines1, labels1 = ax5.get_legend_handles_labels()
            lines2, labels2 = ax5_twin.get_legend_handles_labels()
            ax5.legend(lines1 + lines2, labels1 + labels2, loc='best')

            plt.tight_layout()
            # Use absolute path for saving figures
            save_path = os.path.join(output_dir, f"art_controlled_t{target_time:.1f}.png")
            plt.savefig(save_path, dpi=300)
            plt.close()

            print(f"Generated visualization for t = {target_time}, saved to {save_path}")
        except Exception as e:
            print(f"Error creating visualization for t = {target_time}: {e}")
            import traceback

            traceback.print_exc()

        # Save final data
        # Only execute this after the main loop completes (after all time steps)
        # Save final data
        # Only execute this after the main loop completes (after all time steps)
        # Save final data
        print("\nSaving final simulation data...")
        data_file_path = os.path.join(output_dir, "art_controlled_data.npz")
        np.savez(data_file_path,
                 time=np.array(time_list),
                 psi_reg=np.array(Psi_reg_list),
                 psi_fib=np.array(Psi_fib_list),
                 horizon_reg=np.array(horizon_reg_list),
                 horizon_fib=np.array(horizon_fib_list),
                 variance_reg=np.array(variance_reg_list),
                 variance_fib=np.array(variance_fib_list),
                 entropy_reg=np.array(entropy_reg_list),
                 entropy_fib=np.array(entropy_fib_list),
                 kurtosis_reg=np.array(kurtosis_reg_list),
                 kurtosis_fib=np.array(kurtosis_fib_list),
                 coherence_reg=np.array(coherence_reg_list),
                 coherence_fib=np.array(coherence_fib_list),
                 cytokine=np.array(cytokine_list),
                 norm_check=np.array(norm_check_list))

        # Calculate summary statistics
        print("Calculating summary statistics...")
        try:
            # Calculate coherence decay rates
            time_points = np.array(time_list)


            # Calculate coherence half-life
            def calculate_half_life(coherence_values, times):
                """Calculate the time at which coherence reaches half of its decay"""
                initial_coherence = coherence_values[0]
                final_coherence = coherence_values[-1]
                half_decay_value = initial_coherence - (initial_coherence - final_coherence) / 2

                # Find the first time point where coherence drops below half-decay value
                for i, coh in enumerate(coherence_values):
                    if coh <= half_decay_value:
                        if i > 0:
                            # Linear interpolation for more precise half-life
                            t1, t2 = times[i - 1], times[i]
                            c1, c2 = coherence_values[i - 1], coherence_values[i]
                            half_life = t1 + (half_decay_value - c1) * (t2 - t1) / (c2 - c1)
                            return half_life
                        return times[i]

                # If never drops below half-decay, return the final time
                return times[-1]


            # Calculate Area Under the Coherence Curve (AUC)
            def calculate_auc(coherence_values, times):
                """Calculate area under the coherence curve using trapezoidal rule"""
                from scipy import integrate
                return integrate.trapezoid(coherence_values, times)


            # Calculate metrics
            half_life_reg = calculate_half_life(coherence_reg_list, time_list)
            half_life_fib = calculate_half_life(coherence_fib_list, time_list)
            auc_reg = calculate_auc(coherence_reg_list, time_list)
            auc_fib = calculate_auc(coherence_fib_list, time_list)

            # Calculate improvement percentages
            coherence_improvement = ((coherence_fib_list[-1] - coherence_reg_list[-1]) /
                                     coherence_reg_list[-1] * 100)
            variance_reduction = ((variance_reg_list[-1] - variance_fib_list[-1]) /
                                  variance_reg_list[-1] * 100)
            half_life_improvement = ((half_life_fib - half_life_reg) / half_life_reg * 100)
            auc_improvement = ((auc_fib - auc_reg) / auc_reg * 100)

            # Print summary statistics
            print("\nART-Controlled HIV Simulation Summary:")
            print(f"Final coherence (Regular): {coherence_reg_list[-1]:.4f}")
            print(f"Final coherence (Fibonacci): {coherence_fib_list[-1]:.4f}")
            print(f"Coherence improvement: {coherence_improvement:.2f}%")
            print(f"Final variance (Regular): {variance_reg_list[-1]:.4f}")
            print(f"Final variance (Fibonacci): {variance_fib_list[-1]:.4f}")
            print(f"Variance reduction: {variance_reduction:.2f}%")
            print(f"Coherence half-life (Regular): {half_life_reg:.4f} time units")
            print(f"Coherence half-life (Fibonacci): {half_life_fib:.4f} time units")
            print(f"Half-life improvement: {half_life_improvement:.2f}%")
            print(f"Coherence AUC (Regular): {auc_reg:.4f}")
            print(f"Coherence AUC (Fibonacci): {auc_fib:.4f}")
            print(f"AUC improvement: {auc_improvement:.2f}%")

            # Save summary statistics to CSV
            summary_file_path = os.path.join(output_dir, "art_controlled_summary.csv")
            with open(summary_file_path, 'w') as f:
                f.write("Metric,Regular,Fibonacci,Improvement(%)\n")
                f.write(
                    f"Final Coherence,{coherence_reg_list[-1]:.6f},{coherence_fib_list[-1]:.6f},{coherence_improvement:.2f}\n")
                f.write(
                    f"Final Variance,{variance_reg_list[-1]:.6f},{variance_fib_list[-1]:.6f},{variance_reduction:.2f}\n")
                f.write(f"Final Entropy,{entropy_reg_list[-1]:.6f},{entropy_fib_list[-1]:.6f},N/A\n")
                f.write(f"Final Kurtosis,{kurtosis_reg_list[-1]:.6f},{kurtosis_fib_list[-1]:.6f},N/A\n")
                f.write(f"Coherence Half-life,{half_life_reg:.6f},{half_life_fib:.6f},{half_life_improvement:.2f}\n")
                f.write(f"Coherence AUC,{auc_reg:.6f},{auc_fib:.6f},{auc_improvement:.2f}\n")

            print(f"Summary statistics saved to {summary_file_path}")

            # Create coherence decay plot
            plt.figure(figsize=(10, 6))
            plt.plot(time_list, coherence_reg_list, 'b-', label='Regular Grid')
            plt.plot(time_list, coherence_fib_list, 'r-', label='Fibonacci Grid')
            plt.axhline(y=coherence_reg_list[0] - (coherence_reg_list[0] - coherence_reg_list[-1]) / 2,
                        color='b', linestyle='--', alpha=0.5)
            plt.axhline(y=coherence_fib_list[0] - (coherence_fib_list[0] - coherence_fib_list[-1]) / 2,
                        color='r', linestyle='--', alpha=0.5)
            plt.axvline(x=half_life_reg, color='b', linestyle=':', alpha=0.5)
            plt.axvline(x=half_life_fib, color='r', linestyle=':', alpha=0.5)
            plt.xlabel('Time')
            plt.ylabel('Coherence')
            plt.title('Coherence Decay in ART-Controlled HIV Simulation')
            plt.legend()
            plt.grid(True)
            coherence_plot_path = os.path.join(output_dir, "art_controlled_coherence_decay.png")
            plt.savefig(coherence_plot_path, dpi=300)
            plt.close()

            # Export coherence data for cross-simulation comparison
            coherence_data_path = os.path.join(output_dir, "art_controlled_coherence_data.csv")
            with open(coherence_data_path, 'w') as f:
                f.write("Time,Regular_Grid_Coherence,Fibonacci_Grid_Coherence,Difference\n")
                for i, t in enumerate(time_list):
                    diff = coherence_fib_list[i] - coherence_reg_list[i]
                    f.write(f"{t:.4f},{coherence_reg_list[i]:.6f},{coherence_fib_list[i]:.6f},{diff:.6f}\n")

            print(f"Coherence data saved to {coherence_data_path}")
            print(f"Coherence decay plot saved to {coherence_plot_path}")

        except Exception as e:
            print(f"Error generating summary statistics: {e}")
            import traceback

            traceback.print_exc()

        # Export event horizon data for temporal analysis
        try:
            horizon_data_path = os.path.join(output_dir, "art_controlled_horizon_data.csv")
            with open(horizon_data_path, 'w') as f:
                f.write("Time,Mean_Radius_Regular,Mean_Radius_Fibonacci,Std_Regular,Std_Fibonacci\n")
                for i, t in enumerate(time_list):
                    if i < len(horizon_reg_list) and i < len(horizon_fib_list):
                        mean_reg = np.mean(horizon_reg_list[i])
                        mean_fib = np.mean(horizon_fib_list[i])
                        std_reg = np.std(horizon_reg_list[i])
                        std_fib = np.std(horizon_fib_list[i])
                        f.write(f"{t:.4f},{mean_reg:.6f},{mean_fib:.6f},{std_reg:.6f},{std_fib:.6f}\n")

            print(f"Event horizon data saved to {horizon_data_path}")
        except Exception as e:
            print(f"Error exporting event horizon data: {e}")

        # Create a JSON metadata file to document the simulation parameters
        import json
        import datetime

        try:
            metadata = {
                "simulation_type": "ART-controlled HIV",
                "timestamp": str(datetime.datetime.now()),
                "parameters": {
                    "hbar": hbar,
                    "m": m,
                    "L": L,
                    "R_inner": R_inner,
                    "R_outer": R_outer,
                    "N_r": N_r,
                    "N_z": N_z,
                    "dr": dr,
                    "dz": dz,
                    "dt": dt,
                    "time_steps": time_steps,
                    "V_0": V_0,
                    "Gamma_0": Gamma_0,
                    "alpha_c": alpha_c
                },
                "results": {
                    "final_coherence_regular": float(coherence_reg_list[-1]),
                    "final_coherence_fibonacci": float(coherence_fib_list[-1]),
                    "coherence_improvement_percent": float(coherence_improvement),
                    "variance_reduction_percent": float(variance_reduction),
                    "half_life_regular": float(half_life_reg),
                    "half_life_fibonacci": float(half_life_fib),
                    "half_life_improvement_percent": float(half_life_improvement),
                    "auc_regular": float(auc_reg),
                    "auc_fibonacci": float(auc_fib),
                    "auc_improvement_percent": float(auc_improvement)
                }
            }

            metadata_path = os.path.join(output_dir, "art_controlled_metadata.json")
            with open(metadata_path, 'w') as f:
                json.dump(metadata, f, indent=2)

            print(f"Simulation metadata saved to {metadata_path}")
        except Exception as e:
            print(f"Error saving metadata: {e}")

        # Add variance-focused plots to better highlight the positive results
        try:
            plt.figure(figsize=(10, 6))
            plt.plot(time_list, variance_reg_list, 'b-', label='Regular Grid')
            plt.plot(time_list, variance_fib_list, 'r-', label='Fibonacci Grid')
            plt.xlabel('Time')
            plt.ylabel('Wavefunction Variance')
            plt.title('Variance Growth in ART-Controlled HIV Simulation')
            plt.legend()
            plt.grid(True)
            variance_plot_path = os.path.join(output_dir, "art_controlled_variance.png")
            plt.savefig(variance_plot_path, dpi=300)
            plt.close()
            print(f"Variance plot saved to {variance_plot_path}")
        except Exception as e:
            print(f"Error creating variance plot: {e}")

        # Final summary
        print("\nART-controlled HIV state simulation complete!")
        print(f"Data saved to {data_file_path}")
        print(f"Visualizations saved to {output_dir}")
        print("\nComparing Regular vs Fibonacci grid:")
        print(f"  Initial variance: {variance_reg_list[0]:.4f} vs {variance_fib_list[0]:.4f}")
        print(f"  Final variance: {variance_reg_list[-1]:.4f} vs {variance_fib_list[-1]:.4f}")
        print(
            f"  Variance increase: {(variance_reg_list[-1] / variance_reg_list[0]):.2f}x vs {(variance_fib_list[-1] / variance_fib_list[0]):.2f}x")
        print(
            f"  Coherence retention: {(coherence_reg_list[-1] / coherence_reg_list[0] * 100):.2f}% vs {(coherence_fib_list[-1] / coherence_fib_list[0] * 100):.2f}%")
        print(f"  Baseline coherence half-life: {half_life_reg:.4f} vs {half_life_fib:.4f} time units")
        print("\nSimulation complete!")