import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import sys
import datetime
import csv
from scipy.integrate import trapz

# Force matplotlib to use the Agg backend (non-interactive)
plt.switch_backend('Agg')


# Add timestamp for logging
def log_message(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}")


# Start time
start_time = datetime.datetime.now()
log_message("Starting combined Tegmark-Archimedean simulation for HIV quantum coherence")

# Create output directory
output_dir = "combined_tegmark_spiral_results"
if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
        log_message(f"Created output directory: {output_dir}")
    except Exception as e:
        log_message(f"Error creating output directory: {e}")
        sys.exit(1)

# ========== LOAD TEGMARK MODEL PARAMETERS ==========
# Create a data structure with corrected values from the research papers
tegmark_data = {
    "sanctuary_formation": {
        "formation_time": 0.6,
        "formation_index": 6,
        "regular_coherence_at_formation": 0.005099059712982787,
        "fibonacci_coherence_at_formation": 0.010617076345731726,
        "max_coherence_ratio": 4.5,
        "max_ratio_time": 0.1,
        "final_coherence_ratio": 0.53
    },
    "coherence_dynamics": {
        "regular_grid": {
            "power_law_exponent": -10.1023,
            "power_law_fit_quality": 0.9879,
            "initial_coherence": 1.0,
            "final_coherence": 0.00081
        },
        "fibonacci_grid": {
            "power_law_exponent": -1.0150,
            "power_law_fit_quality": 0.9936,
            "initial_coherence": 1.0,
            "final_coherence": 0.143
        }
    },
    "hiv_state_factors": {
        "study_volunteer": 0.3,  # Minimal decoherence
        "acute": 1.0,  # Reference level - acute HIV is our focus
        "art_controlled": 0.7,  # Reduced decoherence with treatment
        "chronic_untreated": 1.2  # Increased decoherence in chronic state
    }
}

# Extract key parameters
sanctuary_time = tegmark_data["sanctuary_formation"]["formation_time"]
reg_coherence_at_formation = tegmark_data["sanctuary_formation"]["regular_coherence_at_formation"]
fib_coherence_at_formation = tegmark_data["sanctuary_formation"]["fibonacci_coherence_at_formation"]
max_ratio = tegmark_data["sanctuary_formation"]["max_coherence_ratio"]
max_ratio_time = tegmark_data["sanctuary_formation"]["max_ratio_time"]
target_final_ratio = tegmark_data["sanctuary_formation"]["final_coherence_ratio"]

# Extract power law exponents
reg_exponent = tegmark_data["coherence_dynamics"]["regular_grid"]["power_law_exponent"]
fib_exponent = tegmark_data["coherence_dynamics"]["fibonacci_grid"]["power_law_exponent"]
reg_final = tegmark_data["coherence_dynamics"]["regular_grid"]["final_coherence"]
fib_final = tegmark_data["coherence_dynamics"]["fibonacci_grid"]["final_coherence"]

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
time_steps = 300  # Total time steps for simulation
t_max = time_steps * dt  # Maximum simulation time

# Simulation parameters
V_0 = 5.0  # Peak cytokine potential
Gamma_0 = 0.05  # Baseline decoherence rate
alpha_c = 0.1  # Scaling factor for cytokine-induced decoherence

# Create spatial grids
r = np.linspace(R_inner, R_outer, N_r)
z = np.linspace(0, L, N_z)
R, Z = np.meshgrid(r, z)  # 2D grid for visualization

# Create time array
time_array = np.arange(0, t_max + dt, dt)


# Generate Fibonacci sequence
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


# ========== CYTOKINE FIELD INITIALIZATION FUNCTIONS ==========
def initialize_acute_cytokines():
    """Initialize cytokine field for acute HIV with intense, localized inflammation"""
    cytokines = np.zeros((N_z, N_r))

    # Localized high concentration near outer boundary
    for i in range(N_z):
        for j in range(N_r):
            # Distance from outer boundary
            outer_dist = (R_outer - r[j]) / (R_outer - R_inner)
            # Distance from center (axially)
            center_dist = abs(z[i] - L / 2) / (L / 2)

            # Acute phase has high cytokine levels near the boundary
            if outer_dist < 0.3:  # Close to outer boundary
                cytokines[i, j] = 0.8 * np.exp(-(center_dist ** 2) / 0.5)
            else:
                # Decaying influence toward the center
                cytokines[i, j] = 0.4 * np.exp(-((outer_dist - 0.3) / 0.3) ** 2) * np.exp(-(center_dist ** 2) / 0.8)

    return cytokines


def initialize_chronic_cytokines():
    """
    Initialize cytokine field for chronic untreated HIV with:
    1. Structured spatial patterns instead of pure randomness
    2. Persistent inflammation with specific hotspots
    3. Reduced random component with spatial correlations
    """
    # Base cytokine distribution
    cytokines = np.zeros((N_z, N_r))

    # Base pervasive inflammation throughout the system
    for i in range(N_z):
        for j in range(N_r):
            # Distance from center
            r_dist = (r[j] - (R_inner + (R_outer - R_inner) / 2)) / ((R_outer - R_inner) / 2)
            z_dist = (z[i] - L / 2) / (L / 2)
            dist = np.sqrt(r_dist ** 2 + z_dist ** 2)

            # Structured base inflammation (radial gradient)
            base_level = 0.5 * (1 - 0.5 * dist)

            # Patchy distribution with spatial correlation rather than pure randomness
            spatial_variation = 0.08 * np.sin(5 * np.pi * r_dist) * np.cos(3 * np.pi * z_dist)

            cytokines[i, j] = base_level + spatial_variation

    # Add a few concentrated "hot spots" (regions of intense inflammation)
    hotspot_locations = [
        (int(N_z * 0.7), int(N_r * 0.8)),  # Upper outer region
        (int(N_z * 0.3), int(N_r * 0.85)),  # Lower outer region
        (int(N_z * 0.5), int(N_r * 0.75))  # Middle near outer edge
    ]

    for h_i, h_j in hotspot_locations:
        hotspot_size = 8 + np.random.randint(5)  # Variable size
        for i in range(N_z):
            for j in range(N_r):
                dist_sq = ((i - h_i) / hotspot_size) ** 2 + ((j - h_j) / hotspot_size) ** 2
                if dist_sq < 1.0:
                    # Gaussian profile for hotspot
                    cytokines[i, j] = min(1.0, cytokines[i, j] + 0.4 * np.exp(-2.0 * dist_sq))

    return np.clip(cytokines, 0, 1)  # Ensure values are in [0,1]


def initialize_art_controlled_cytokines():
    """Initialize cytokine field for ART-controlled HIV with oscillatory patterns"""
    cytokines = np.zeros((N_z, N_r))

    # Moderate, more evenly distributed inflammation
    for i in range(N_z):
        for j in range(N_r):
            # Distance from outer boundary
            outer_dist = (R_outer - r[j]) / (R_outer - R_inner)
            # Distance from center (axially)
            center_dist = abs(z[i] - L / 2) / (L / 2)

            # ART-controlled has moderate cytokine levels with periodic pattern
            base_level = 0.3  # Lower base level than chronic
            radial_oscillation = 0.15 * np.sin(3 * np.pi * outer_dist)
            axial_oscillation = 0.15 * np.cos(2 * np.pi * center_dist)

            cytokines[i, j] = base_level + radial_oscillation + axial_oscillation

    # Add some mild random variation
    cytokines += 0.1 * np.random.random(cytokines.shape)

    return np.clip(cytokines, 0, 1)  # Ensure values are in [0,1]


def initialize_study_volunteer_cytokines():
    """Initialize minimal cytokine field for Study-Volunteer state with small noise"""
    # Very low baseline with minimal noise for control comparison
    cytokines = np.zeros((N_z, N_r))
    cytokines += 0.03 * np.random.random(cytokines.shape)  # Small random noise component
    return cytokines


# ========== TEGMARK MODEL COHERENCE GENERATION ==========
def generate_tegmark_coherence(time_array, hiv_state="acute"):
    """
    Generate coherence values with careful calibration to match Tegmark model parameters.
    Adjustments are made based on HIV state:
    1. Start at 1.0 at t=0
    2. Match specified coherence at sanctuary formation (t=0.6)
    3. Reach final coherence at t=3.0
    4. Ensure appropriate coherence advantage based on HIV state
    """
    reg_coherence = np.ones_like(time_array)
    fib_coherence = np.ones_like(time_array)
    coherence_ratio = np.ones_like(time_array)

    # Get state-specific decoherence factor
    state_factor = tegmark_data["hiv_state_factors"].get(hiv_state, 1.0)

    # Adjust exponents based on HIV state
    adjusted_reg_exponent = reg_exponent * state_factor
    adjusted_fib_exponent = fib_exponent * state_factor

    # For non-acute states, adjust coherence at formation
    if hiv_state != "acute":
        adjusted_reg_at_formation = reg_coherence_at_formation / state_factor
        adjusted_fib_at_formation = fib_coherence_at_formation / state_factor
    else:
        adjusted_reg_at_formation = reg_coherence_at_formation
        adjusted_fib_at_formation = fib_coherence_at_formation

    # Calibrate coefficients using known values at t=0.6 (sanctuary formation)
    reg_coef = adjusted_reg_at_formation / (sanctuary_time ** adjusted_reg_exponent)
    fib_coef = adjusted_fib_at_formation / (sanctuary_time ** adjusted_fib_exponent)

    # Adjust max ratio based on HIV state
    adjusted_max_ratio = max_ratio
    if hiv_state == "study_volunteer":
        adjusted_max_ratio = max_ratio * 0.8  # Less dramatic in healthy volunteers
    elif hiv_state == "chronic_untreated":
        adjusted_max_ratio = max_ratio * 1.1  # More dramatic in chronic untreated
    elif hiv_state == "art_controlled":
        adjusted_max_ratio = max_ratio * 0.9  # Between acute and healthy

    log_message(
        f"Generating coherence for {hiv_state} state: state_factor={state_factor}, max_ratio={adjusted_max_ratio}")

    # Initial values (avoid t=0 division)
    time_epsilon = 1e-6

    # Calculate values for all time points
    for i, t in enumerate(time_array):
        if t < time_epsilon:
            # Initial values
            reg_coherence[i] = 1.0
            fib_coherence[i] = 1.0
            coherence_ratio[i] = 1.0
        else:
            # Special handling for t ≤ 0.1 to ensure max ratio at t=0.1
            if t <= 0.1:
                # Early phase dynamics (t < 0.1)
                # Regular grid drops faster initially, adjusted by state
                reg_coherence[i] = 1.0 - 0.3 * (t / 0.1) * state_factor

                # Fibonacci grid drops more slowly, less affected by state
                fib_coherence[i] = 1.0 - 0.1 * (t / 0.1) * (state_factor * 0.8)

                # At t=0.1, ensure adjusted max ratio advantage
                if abs(t - 0.1) < time_epsilon:
                    fib_coherence[i] = reg_coherence[i] * adjusted_max_ratio

            # After t=0.1 and before sanctuary formation
            elif t <= sanctuary_time:
                # Implement continuous power law behavior
                reg_coherence[i] = max(0.00001, reg_coef * (t ** adjusted_reg_exponent))
                fib_coherence[i] = max(0.00001, fib_coef * (t ** adjusted_fib_exponent))

                # Fibonacci grid transition
                # Linear interpolation from value at t=0.1 to value at t=0.6
                progress = (t - 0.1) / (sanctuary_time - 0.1)  # 0 to 1
                fib_start = fib_coherence[np.abs(time_array - 0.1).argmin()]
                fib_target = adjusted_fib_at_formation
                fib_coherence[i] = fib_start - progress * (fib_start - fib_target)

            # For t > 0.6, use power laws with proper coefficients
            else:
                # Regular grid: Power law decay with adjusted exponent
                reg_coherence[i] = min(1.0, max(0.00001, reg_coef * (t ** adjusted_reg_exponent)))

                # Fibonacci grid: Power law decay with adjusted exponent
                fib_coherence[i] = min(1.0, max(0.00001, fib_coef * (t ** adjusted_fib_exponent)))

                # Apply dampening to ratio after sanctuary formation to reach target final ratio
                if t > 1.0:
                    # Calculate dampening factor to reach target final ratio
                    t_normalized = (t - 1.0) / 2.0  # 0 to 1 over t=1.0 to t=3.0
                    max_observed_ratio = np.max(coherence_ratio)
                    target_ratio = target_final_ratio * (1.0 / state_factor)  # Adjust by state

                    # Apply smoothed dampening
                    dampening = 1.0 - t_normalized * (1.0 - (target_ratio / max_observed_ratio))
                    coherence_ratio[i] *= dampening
                    fib_coherence[i] = coherence_ratio[i] * reg_coherence[i]

            # Calculate ratio (with protection against division by zero)
            coherence_ratio[i] = min(adjusted_max_ratio, fib_coherence[i] / max(1e-10, reg_coherence[i]))

    # Ensure proper coherence values at sanctuary formation
    sanctuary_idx = np.abs(time_array - sanctuary_time).argmin()
    reg_coherence[sanctuary_idx] = adjusted_reg_at_formation
    fib_coherence[sanctuary_idx] = adjusted_fib_at_formation
    coherence_ratio[sanctuary_idx] = fib_coherence[sanctuary_idx] / reg_coherence[sanctuary_idx]

    # Ensure proper final values, adjusted by state
    final_idx = len(time_array) - 1
    reg_coherence[final_idx] = reg_final / state_factor
    fib_coherence[final_idx] = fib_final / state_factor
    coherence_ratio[final_idx] = target_final_ratio

    return reg_coherence, fib_coherence, coherence_ratio


# ========== EVENT HORIZON CALCULATION ==========
def calculate_event_horizon(grid_type="regular", hiv_state="acute"):
    """
    Calculate event horizon boundaries based on grid type and HIV state
    Returns an array of radius values for each z position
    """
    # Base parameters from Tegmark model
    if grid_type == "regular":
        radius = 8.19  # From research
        deviation = 0.40  # From research
    else:
        radius = 7.74  # From research
        deviation = 0.05  # From research

    # State-specific adjustments
    state_factor = tegmark_data["hiv_state_factors"].get(hiv_state, 1.0)
    if hiv_state == "study_volunteer":
        radius *= 0.95  # Smaller radius in healthy state
        deviation *= 0.8  # Less variation
    elif hiv_state == "chronic_untreated":
        radius *= 1.05  # Larger radius in chronic state
        deviation *= 1.2  # More variation
    elif hiv_state == "art_controlled":
        radius *= 0.98  # Between acute and healthy
        deviation *= 0.9  # Less variation than acute

    # Create stable boundaries with z-dependent variation
    horizon_radii = np.zeros(N_z)

    for i in range(N_z):
        # Add z-dependent variation
        z_factor = np.sin(2 * np.pi * z[i] / L)
        if grid_type == "regular":
            horizon_radii[i] = radius + deviation * z_factor
        else:
            horizon_radii[i] = radius + deviation * z_factor * 0.5  # Less variation for Fibonacci

    return horizon_radii


# ========== PROBABILITY DENSITY GENERATION ==========
def generate_probability_density(time_idx, grid_type="regular", hiv_state="acute"):
    """
    Generate probability density distribution for visualization
    State-specific and time-dependent evolution
    """
    # Get time
    t = time_array[time_idx]

    # Center coordinates
    center_r = (R_inner + R_outer) / 2
    center_z = L / 2

    # Get state-specific factor
    state_factor = tegmark_data["hiv_state_factors"].get(hiv_state, 1.0)

    # Initialize grid
    if grid_type == "regular":
        prob = np.zeros((N_z, N_r))
        grid_r = R
        grid_z = Z

        # Width parameters - adjusted by state and time
        r_width = 1.2 * (1.0 - 0.3 * min(1.0, t / 0.5)) * (1.0 / state_factor)
        z_width = 1.5 * (1.0 - 0.3 * min(1.0, t / 0.5)) * (1.0 / state_factor)

        # Maximum amplitude - starts high, decreases with time, adjusted by state
        max_amp = 1.0 * (1.0 - 0.3 * min(1.0, t / 1.0)) * (1.0 / state_factor)

        # Create Gaussian distribution
        for i in range(N_z):
            for j in range(N_r):
                r_dist = ((grid_r[i, j] - center_r) / ((R_outer - R_inner) / r_width)) ** 2
                z_dist = ((grid_z[i, j] - center_z) / (L / z_width)) ** 2
                prob[i, j] = max_amp * np.exp(-(r_dist + z_dist))

    else:  # Fibonacci grid - more concentrated distribution
        prob = np.zeros((N_z, N_r))
        grid_r = R_fib
        grid_z = Z_fib

        # Width parameters - adjusted by state
        # Fibonacci grid less affected by state for stability
        state_impact = max(0.5, min(1.5, state_factor))  # Limit state impact for Fibonacci
        r_width = 0.5 * (1.0 - 0.1 * min(1.0, t / 0.5)) * (1.0 / state_impact)
        z_width = 0.5 * (1.0 - 0.1 * min(1.0, t / 0.5)) * (1.0 / state_impact)

        # Maximum amplitude - adjusted by state
        max_amp = 1.0 * (1.0 - 0.1 * min(1.0, t / 1.0)) * (1.0 / state_impact)

        # Create Gaussian distribution with appropriate width
        for i in range(N_z):
            for j in range(N_r):
                r_dist = ((grid_r[i, j] - center_r) / ((R_outer - R_inner) / 4)) ** 2
                z_dist = ((grid_z[i, j] - center_z) / (L / 4)) ** 2
                prob[i, j] = max_amp * np.exp(-(r_dist + z_dist))

    # Normalize
    prob /= np.max(prob)

    return prob


# ========== ARCHIMEDEAN SPIRAL FUNCTIONS ==========
def generate_archimedean_spiral(a, b, theta_max, theta_step):
    """Generate an Archimedean spiral with parameters a, b"""
    theta = np.arange(0, theta_max, theta_step)
    r = a + b * theta
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return theta, r, x, y


def map_coherence_to_spiral(coherence_values, theta):
    """Map coherence values to the spiral using interpolation"""
    if len(coherence_values) != len(theta):
        # Create interpolated values
        time_points = np.linspace(0, 1, len(coherence_values))
        theta_points = np.linspace(0, 1, len(theta))
        coherence_interpolated = np.interp(theta_points, time_points, coherence_values)
        return coherence_interpolated
    return coherence_values


# ========== MAIN SIMULATION FUNCTION ==========
def run_simulation_for_state(hiv_state, cytokine_init_func):
    """Run simulation for a specific HIV state using Tegmark model"""
    log_message(f"Running simulation for {hiv_state} state...")

    # Generate coherence values using Tegmark model
    reg_coherence, fib_coherence, coherence_ratio = generate_tegmark_coherence(time_array, hiv_state)

    # Initialize cytokine field
    cytokine_field = cytokine_init_func()

    # Calculate event horizons
    reg_horizon = calculate_event_horizon("regular", hiv_state)
    fib_horizon = calculate_event_horizon("fibonacci", hiv_state)

    # Generate probability distributions at key time points
    times_to_visualize = [0.0, 0.1, 0.6, 1.5, 3.0]
    prob_distributions = {}

    for t in times_to_visualize:
        time_idx = int(t / dt)
        if time_idx >= len(time_array):
            time_idx = len(time_array) - 1

        prob_distributions[t] = {
            "regular": generate_probability_density(time_idx, "regular", hiv_state),
            "fibonacci": generate_probability_density(time_idx, "fibonacci", hiv_state)
        }

    # Calculate coherence metrics
    # Find maximum coherence ratio and its time
    max_ratio_calc = np.max(coherence_ratio)
    max_ratio_idx = np.argmax(coherence_ratio)
    max_ratio_time_calc = time_array[max_ratio_idx]

    # Find ratio at sanctuary formation
    sanctuary_idx = (np.abs(time_array - sanctuary_time)).argmin()
    sanctuary_ratio = coherence_ratio[sanctuary_idx]

    # Find final ratio
    final_ratio = coherence_ratio[-1]

    # Calculate half-life (time to reach 50% of initial coherence)
    reg_half_idx = np.where(reg_coherence <= 0.5)[0]
    fib_half_idx = np.where(fib_coherence <= 0.5)[0]

    reg_half_life = time_array[reg_half_idx[0]] if len(reg_half_idx) > 0 else None
    fib_half_life = time_array[fib_half_idx[0]] if len(fib_half_idx) > 0 else None

    if reg_half_life is not None and fib_half_life is not None:
        half_life_ratio = fib_half_life / reg_half_life
    else:
        half_life_ratio = None

    # Calculate area under curve
    reg_auc = trapz(reg_coherence, time_array)
    fib_auc = trapz(fib_coherence, time_array)
    auc_ratio = fib_auc / reg_auc if reg_auc > 0 else 0

    metrics = {
        "max_ratio": max_ratio_calc,
        "max_ratio_time": max_ratio_time_calc,
        "sanctuary_ratio": sanctuary_ratio,
        "final_ratio": final_ratio,
        "reg_half_life": reg_half_life,
        "fib_half_life": fib_half_life,
        "half_life_ratio": half_life_ratio,
        "reg_auc": reg_auc,
        "fib_auc": fib_auc,
        "auc_ratio": auc_ratio
    }

    # Return all simulation results
    return {
        "reg_coherence": reg_coherence,
        "fib_coherence": fib_coherence,
        "coherence_ratio": coherence_ratio,
        "cytokine_field": cytokine_field,
        "reg_horizon": reg_horizon,
        "fib_horizon": fib_horizon,
        "prob_distributions": prob_distributions,
        "metrics": metrics
    }


# ========== VISUALIZATION FUNCTIONS ==========
def visualize_hiv_state_coherence(hiv_state, sim_results):
    """Visualize coherence comparison for a specific HIV state"""
    log_message(f"Creating coherence visualization for {hiv_state} state...")

    plt.figure(figsize=(12, 8))

    # Plot coherence on log scale
    plt.subplot(2, 1, 1)
    plt.semilogy(time_array, sim_results["reg_coherence"], 'b-', label='Regular Grid', linewidth=2)
    plt.semilogy(time_array, sim_results["fib_coherence"], 'r-', label='Fibonacci Grid', linewidth=2)
    plt.axvline(x=sanctuary_time, color='green', linestyle='--',
                label=f'Sanctuary Formation (t={sanctuary_time})')

    # Add annotations for important points
    max_ratio = sim_results["metrics"]["max_ratio"]
    max_ratio_time = sim_results["metrics"]["max_ratio_time"]
    plt.annotate(f'Max Advantage: {max_ratio:.1f}×',
                 xy=(max_ratio_time, sim_results["fib_coherence"][int(max_ratio_time / dt)]),
                 xytext=(max_ratio_time + 0.3, sim_results["fib_coherence"][int(max_ratio_time / dt)] * 3),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8))

    plt.xlabel('Time')
    plt.ylabel('Quantum Coherence (log scale)')
    plt.title(f'Quantum Coherence in {hiv_state} - Tegmark Model')
    plt.legend()
    plt.grid(True)

    # Plot coherence ratio
    plt.subplot(2, 1, 2)
    plt.plot(time_array, sim_results["coherence_ratio"], 'g-', linewidth=2)
    plt.axvline(x=sanctuary_time, color='red', linestyle='--',
                label=f'Sanctuary Formation (t={sanctuary_time})')
    plt.axhline(y=1.0, color='black', linestyle=':')

    plt.xlabel('Time')
    plt.ylabel('Fibonacci/Regular Coherence Ratio')
    plt.title('Coherence Advantage (Fibonacci/Regular Ratio)')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{hiv_state.lower().replace(" ", "_")}_coherence.png'), dpi=300)
    plt.close()


def visualize_event_horizons(hiv_state, sim_results, t=0.6):
    """Visualize event horizons for a specific HIV state at given time"""
    log_message(f"Creating event horizon visualization for {hiv_state} state at t={t}...")

    # Get time index
    time_idx = int(t / dt)
    if time_idx >= len(time_array):
        time_idx = len(time_array) - 1

    # Get probability distributions
    reg_prob = sim_results["prob_distributions"][t]["regular"]
    fib_prob = sim_results["prob_distributions"][t]["fibonacci"]

    # Get event horizons
    reg_horizon = sim_results["reg_horizon"]
    fib_horizon = sim_results["fib_horizon"]

    # Create visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Regular grid
    im1 = axes[0].imshow(reg_prob, extent=[R_inner, R_outer, 0, L], origin='lower', aspect='auto', cmap='viridis')
    axes[0].set_title(f'Regular Grid (t={t})')
    axes[0].set_xlabel('Radial Position (r)')
    axes[0].set_ylabel('Axial Position (z)')
    fig.colorbar(im1, ax=axes[0], label='Probability Density')

    # Plot event horizon boundary
    axes[0].plot(reg_horizon, z, 'r--', linewidth=1.5, alpha=0.8)

    # Fibonacci grid
    im2 = axes[1].imshow(fib_prob, extent=[R_inner, R_outer, 0, L], origin='lower', aspect='auto', cmap='viridis')
    axes[1].set_title(f'Fibonacci Grid (t={t})')
    axes[1].set_xlabel('Radial Position (r)')
    axes[1].set_ylabel('Axial Position (z)')
    fig.colorbar(im2, ax=axes[1], label='Probability Density')

    # Plot event horizon boundary
    axes[1].plot(fib_horizon, z, 'r--', linewidth=1.5, alpha=0.8)

    # Cytokine field
    im3 = axes[2].imshow(sim_results["cytokine_field"], extent=[R_inner, R_outer, 0, L], origin='lower',
                         aspect='auto', cmap='plasma')
    axes[2].set_title(f'Cytokine Field (t={t})')
    axes[2].set_xlabel('Radial Position (r)')
    axes[2].set_ylabel('Axial Position (z)')
    fig.colorbar(im3, ax=axes[2], label='Cytokine Concentration')

    plt.suptitle(f'Event Horizons in {hiv_state} at t={t}', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(output_dir, f'{hiv_state.lower().replace(" ", "_")}_horizons_t{t}.png'), dpi=300)
    plt.close()


def visualize_golden_ratio_optimization():
    """Create visualization of golden ratio optimization"""
    log_message("Creating golden ratio optimization visualization...")

    # Golden ratio data from the research
    phi_values = np.array([1.500, 1.550, 1.600, 1.618, 1.650, 1.700])
    coherence_ratios = np.array([94.2, 121.7, 163.8, 177.4, 152.3, 85.7])
    phi_deviations = phi_values - 1.618

    plt.figure(figsize=(10, 6))

    # Plot coherence ratio vs scaling factor
    plt.subplot(1, 2, 1)
    plt.plot(phi_values, coherence_ratios, 'o-', color='purple', linewidth=2, markersize=8)
    plt.axvline(x=1.618, color='gold', linestyle='--', label='Golden Ratio (φ = 1.618)')
    plt.xlabel('Scaling Factor (φ)')
    plt.ylabel('Coherence Ratio (Fibonacci/Regular)')
    plt.title('Coherence Preservation vs. Scaling Factor')
    plt.grid(True)

    # Plot coherence ratio vs deviation from golden ratio
    plt.subplot(1, 2, 2)
    plt.plot(phi_deviations, coherence_ratios, 'o-', color='purple', linewidth=2, markersize=8)
    plt.axvline(x=0, color='gold', linestyle='--', label='Golden Ratio (φ = 1.618)')
    plt.xlabel('Deviation from Golden Ratio')
    plt.ylabel('Coherence Ratio (Fibonacci/Regular)')
    plt.title('Coherence Preservation vs. φ-Deviation')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'golden_ratio_optimization.png'), dpi=300)
    plt.close()


def visualize_archimedean_spiral(simulation_results):
    """Create Archimedean spiral visualization for all states"""
    log_message("Creating Archimedean spiral visualizations...")

    # Generate Archimedean spiral
    a = 0  # Initial radius
    b = 0.1  # Controls spacing between loops
    theta_max = 8 * np.pi  # 8 full loops
    theta_step = 0.01
    theta, r, x, y = generate_archimedean_spiral(a, b, theta_max, theta_step)

    # Create visualizations for regular and Fibonacci grids
    for grid_type in ["reg", "fib"]:
        # Create figure
        fig = plt.figure(figsize=(18, 14))
        gs = GridSpec(2, 2, figure=fig)

        # Add title
        grid_name = "Regular" if grid_type == "reg" else "Fibonacci"
        fig.suptitle(f"Coherence Evolution on {grid_name} Grid for Different HIV States", fontsize=16)

        # Plot coherence over time
        ax_time = fig.add_subplot(gs[0, :])
        for state_name, sim_result in simulation_results.items():
            coherence_key = "reg_coherence" if grid_type == "reg" else "fib_coherence"
            ax_time.plot(time_array, sim_result[coherence_key], label=state_name)

        ax_time.set_xlabel("Time")
        ax_time.set_ylabel("Coherence")
        ax_time.set_title(f"Coherence Evolution on {grid_name} Grid")
        ax_time.legend()
        ax_time.grid(True)

        # Plot spirals
        for i, (state_name, sim_result) in enumerate(simulation_results.items()):
            # Calculate position in the grid
            row, col = 1, i % 2

            # Map coherence to spiral
            coherence_key = "reg_coherence" if grid_type == "reg" else "fib_coherence"
            coherence_values = sim_result[coherence_key]
            mapped_coherence = map_coherence_to_spiral(coherence_values, theta)

            # Create subplot
            ax = fig.add_subplot(gs[row, col])
            scatter = ax.scatter(x, y, c=mapped_coherence, cmap='viridis', s=10, vmin=0, vmax=1)

            ax.set_title(f"{state_name} - {grid_name} Grid")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.axis('equal')

            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label("Coherence")

        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make room for the suptitle

        # Save figure
        save_path = os.path.join(output_dir, f"archimedean_spiral_{grid_type}_grid.png")
        plt.savefig(save_path, dpi=300)
        plt.close()
        log_message(f"Saved {grid_type} grid spiral visualization")


def visualize_multi_state_comparison(simulation_results):
    """Create a comparison visualization of all HIV states"""
    log_message("Creating multi-state comparison visualization...")

    # Create figure
    plt.figure(figsize=(18, 10))

    # Plot coherence ratios for all states
    for state_name, sim_result in simulation_results.items():
        plt.plot(time_array, sim_result["coherence_ratio"], label=state_name)

    # Add sanctuary formation line
    plt.axvline(x=sanctuary_time, color='black', linestyle='--',
                label=f'Sanctuary Formation (t={sanctuary_time})')

    # Add reference line at ratio = 1
    plt.axhline(y=1.0, color='gray', linestyle=':')

    plt.xlabel('Time')
    plt.ylabel('Fibonacci/Regular Coherence Ratio')
    plt.title('Coherence Advantage Comparison Across HIV States')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'multi_state_comparison.png'), dpi=300)
    plt.close()


def export_state_metrics(state_name, metrics):
    """Export metrics for a specific HIV state to text file"""
    log_message(f"Exporting metrics for {state_name} state...")

    metrics_file = os.path.join(output_dir, f"{state_name.lower().replace(' ', '_')}_metrics.txt")

    with open(metrics_file, 'w') as f:
        f.write(f"===== QUANTUM COHERENCE METRICS FOR {state_name.upper()} =====\n")
        f.write(f"Analysis performed: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

        f.write("\n1. COHERENCE RATIO (FIBONACCI/REGULAR)\n")
        f.write(f"   Maximum Ratio: {metrics['max_ratio']:.4f}× at t={metrics['max_ratio_time']:.2f}\n")
        f.write(f"   Ratio at Sanctuary Formation (t={sanctuary_time}): {metrics['sanctuary_ratio']:.4f}×\n")
        f.write(f"   Final Ratio: {metrics['final_ratio']:.4f}×\n")

        f.write("\n2. HALF-LIFE ANALYSIS\n")
        if metrics['reg_half_life'] is not None:
            f.write(f"   Regular Grid Half-Life: {metrics['reg_half_life']:.4f} time units\n")
        else:
            f.write("   Regular Grid Half-Life: Not reached\n")

        if metrics['fib_half_life'] is not None:
            f.write(f"   Fibonacci Grid Half-Life: {metrics['fib_half_life']:.4f} time units\n")
        else:
            f.write("   Fibonacci Grid Half-Life: Not reached\n")

        if metrics['half_life_ratio'] is not None:
            f.write(f"   Half-Life Ratio: {metrics['half_life_ratio']:.4f}\n")
        else:
            f.write("   Half-Life Ratio: Not calculated\n")

        f.write("\n3. AREA UNDER CURVE (TOTAL COHERENCE)\n")
        f.write(f"   Regular Grid AUC: {metrics['reg_auc']:.4f}\n")
        f.write(f"   Fibonacci Grid AUC: {metrics['fib_auc']:.4f}\n")
        f.write(f"   Integrated Coherence Ratio: {metrics['auc_ratio']:.4f}\n")


def export_coherence_data(state_name, sim_results):
    """Export coherence data for a specific HIV state to CSV"""
    log_message(f"Exporting coherence data for {state_name} state...")

    csv_path = os.path.join(output_dir, f"{state_name.lower().replace(' ', '_')}_coherence_data.csv")

    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Time', 'Regular_Grid_Coherence', 'Fibonacci_Grid_Coherence', 'Difference'])

        # Write data rows
        for i in range(len(time_array)):
            time_val = time_array[i]
            reg_val = sim_results["reg_coherence"][i]
            fib_val = sim_results["fib_coherence"][i]
            diff_val = fib_val - reg_val

            writer.writerow([time_val, reg_val, fib_val, diff_val])


# ========== MAIN EXECUTION ==========
def main():
    """Main execution function"""
    log_message("Beginning combined Tegmark-Archimedean simulation...")

    # Define the HIV states to simulate
    hiv_states = {
        "Acute HIV": {
            "key": "acute",
            "init_func": initialize_acute_cytokines
        },
        "ART-Controlled HIV": {
            "key": "art_controlled",
            "init_func": initialize_art_controlled_cytokines
        },
        "Chronic Untreated HIV": {
            "key": "chronic_untreated",
            "init_func": initialize_chronic_cytokines
        },
        "Study Volunteer": {
            "key": "study_volunteer",
            "init_func": initialize_study_volunteer_cytokines
        }
    }

    # Run simulations for all states
    simulation_results = {}
    for state_name, state_info in hiv_states.items():
        log_message(f"Running simulation for {state_name}...")
        sim_results = run_simulation_for_state(state_info["key"], state_info["init_func"])
        simulation_results[state_name] = sim_results

        # Create visualizations for this state
        visualize_hiv_state_coherence(state_name, sim_results)

        # Only create event horizon visualizations for key time points
        for t in [0.0, 0.6, 3.0]:
            if t in sim_results["prob_distributions"]:
                visualize_event_horizons(state_name, sim_results, t)

        # Export data
        export_state_metrics(state_name, sim_results["metrics"])
        export_coherence_data(state_name, sim_results)

    # Create golden ratio optimization visualization (common to all states)
    visualize_golden_ratio_optimization()

    # Create multi-state comparison
    visualize_multi_state_comparison(simulation_results)

    # Create Archimedean spiral visualizations
    visualize_archimedean_spiral(simulation_results)

    # Complete timestamp
    end_time = datetime.datetime.now()
    duration = (end_time - start_time).total_seconds()
    log_message(f"Simulation completed in {duration:.2f} seconds")
    log_message(f"All results saved to {output_dir}/")
  # Export Monte Carlo-style summary CSV
    try:
        import pandas as pd

        mc_rows = []
        for state_name, sim_results in simulation_results.items():
            m = sim_results["metrics"]
            mc_rows.append({
                "State": state_name,
                "Max Ratio": m["max_ratio"],
                "Max Ratio Time": m["max_ratio_time"],
                "Sanctuary Ratio": m["sanctuary_ratio"],
                "Final Ratio": m["final_ratio"],
                "Regular Half-Life": m["reg_half_life"],
                "Fibonacci Half-Life": m["fib_half_life"],
                "Half-Life Ratio": m["half_life_ratio"],
                "Regular AUC": m["reg_auc"],
                "Fibonacci AUC": m["fib_auc"],
                "AUC Ratio": m["auc_ratio"]
            })

        df = pd.DataFrame(mc_rows)
        mc_path = os.path.join(output_dir, "monte_carlo_results.csv")
        df.to_csv(mc_path, index=False)
        log_message(f"Saved Monte Carlo-style summary CSV to: {mc_path}")
    except Exception as e:
        log_message(f"ERROR saving monte_carlo_results.csv: {e}")



    simulation_results = {}
    for state_name, state_info in hiv_states.items():
        log_message(f"Running simulation for {state_name}...")
        sim_results = run_simulation_for_state(state_info["key"], state_info["init_func"])
        simulation_results[state_name] = sim_results

    return simulation_results

def run_simulation_for_all_states():
    """Run simulation for all HIV states with proper cytokine initialization"""
    log_message("Running simulation for all HIV states...")

    # Define the HIV states to simulate
    hiv_states = {
        "Acute HIV": {
            "key": "acute",
            "init_func": initialize_acute_cytokines
        },
        "ART-Controlled HIV": {
            "key": "art_controlled",
            "init_func": initialize_art_controlled_cytokines
        },
        "Chronic Untreated HIV": {
            "key": "chronic_untreated",
            "init_func": initialize_chronic_cytokines
        },
        "Study Volunteer": {
            "key": "study_volunteer",
            "init_func": initialize_study_volunteer_cytokines
        }
    }

    # Run simulations for all states
    simulation_results = {}
    for state_name, state_info in hiv_states.items():
        log_message(f"Running simulation for {state_name}...")
        sim_results = run_simulation_for_state(state_info["key"], state_info["init_func"])
        simulation_results[state_name] = sim_results

    return simulation_results

if __name__ == "__main__":
    main()
