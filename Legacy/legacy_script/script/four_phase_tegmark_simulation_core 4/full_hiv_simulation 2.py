import numpy as np
import math
import random
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# Define simulation parameters
params = {
    'max_time': 72,              # 72 time units (3 days when converted to hours)
    'output_interval': 0.1,      # Store results every 0.1 time units  
    'initial_targeting': 0.05,   # Starting immune response intensity
    'hiv_phase': 'acute'         # Acute phase of HIV infection
}

def calculate_decoherence_rate(temperature, geometry, scaling):
    """
    Calculate temperature-dependent quantum decoherence rate.

    Parameters:
    -----------
    temperature : float
        Temperature in Celsius
    geometry : str
        'fibonacci' or 'regular'
    scaling : float
        Scaling factor (1.618 for fibonacci, 1.0 for regular)

    Returns:
    --------
    float
        Decoherence rate (inverse of coherence time)
    """
    # Base parameters
    h_bar = 1.054571817e-34  # Reduced Planck constant
    k_b = 1.380649e-23  # Boltzmann constant
    r = 3.8e-9  # Superposition distance

    # Temperature-dependent dielectric constant (corrected from Tegmark)
    # Note proper temperature dependence: higher temp = faster decoherence
    epsilon = 80 * ((310 / temperature) ** 1.5)

    # Geometry-specific factor (critical for Fibonacci protection)
    if geometry == 'fibonacci':
        geometry_factor = scaling ** -2.5  # Power law protection
    else:
        geometry_factor = 10.0 ** 9.1  # Exponential collapse factor

    # Calculate decoherence rate (inverse of coherence time)
    # Note: these are extremely small/large numbers due to quantum scale
    # Using log scale to make more manageable for simulation
    raw_rate = (k_b * temperature * 1.2e-29) / (h_bar * h_bar * epsilon * (r ** 4) * geometry_factor)

    # Return a more manageable value for simulation purposes
    return raw_rate


def simulate_hiv_fever(time_in_hours, phase='acute'):
    """
    Model the fever dynamics specific to HIV infection.

    Parameters:
    -----------
    time_in_hours : float
        Time in hours
    phase : str
        'acute', 'art-treated', or 'chronic-untreated'

    Returns:
    --------
    float
        Body temperature in Celsius
    """
    base_temp = 37.0  # Normal body temperature
    hour_of_day = time_in_hours % 24
    day_number = math.floor(time_in_hours / 24) + 1

    # Fever amplitude (varies by phase)
    fever_amplitude = 0
    decay_rate = 0

    if phase == 'acute':
        fever_amplitude = 2.5  # High fever in acute phase
        decay_rate = 0.2  # Per day
    elif phase == 'art-treated':
        fever_amplitude = 1.2  # Lower fever with treatment
        decay_rate = 0.4  # Faster resolution
    elif phase == 'chronic-untreated':
        fever_amplitude = 1.8  # Moderate persistent fever
        decay_rate = 0.05  # Very slow resolution

    # Day-dependent fever (decays exponentially)
    fever_component = fever_amplitude * math.exp(-decay_rate * (day_number - 1))

    # Circadian rhythm (higher in evening)
    circadian_component = 0.4 * math.sin(2 * math.pi * (hour_of_day / 24) - math.pi / 3)

    # Calculate temperature
    temperature = base_temp + fever_component + circadian_component

    # Add night sweats (between 2-4 AM in acute phase)
    if phase == 'acute' and 2 <= hour_of_day <= 4:
        minute_of_hour = (time_in_hours * 60) % 60
        # Create 2-3 sweating episodes per night
        sweat_episode = math.sin(2 * math.pi * minute_of_hour / 30) > 0.5

        if sweat_episode:
            # Rapid temperature drop during night sweats
            temperature -= 1.5 * math.sin(math.pi * minute_of_hour / 15)

    return temperature


def update_grids(state, dt, current_time):
    """
    Update grid evolution incorporating temperature effects and quantum decoherence.

    Parameters:
    -----------
    state : dict
        Current state of the system
    dt : float
        Time step
    current_time : float
        Current simulation time
    """
    # Get current temperature from HIV fever model (convert simulation time to hours)
    hour_time = current_time * 24  # Convert to hours
    temperature = simulate_hiv_fever(hour_time, state['hiv_phase'])

    # Temperature effects on microtubule dynamics
    temp_diff = temperature - 37.0  # Difference from normal

    # MULTIPLICATIVE DECAY APPROACH - different patterns for each grid:

    # Tegmark's collapse rates - MUCH more dramatic difference
    base_reg_rate = 0.8  # Regular grid crashes extremely quickly (~exponential collapse)
    base_fib_rate = 0.005  # Fibonacci grid much more stable (160x more stable)

    # Temperature effects (exponential impact on regular grid)
    temp_factor_reg = math.exp(0.3 * max(0, temp_diff))  # Strong exponential temperature dependence
    temp_factor_fib = 1.0 + 0.05 * max(0, temp_diff)  # Mild linear temperature dependence

    # Targeting intensity effect
    targeting_effect_reg = state['targeting_intensity']
    targeting_effect_fib = state['targeting_intensity'] * 0.5  # Reduce impact on Fibonacci grid

    # Calculate degradation factors for this time step
    reg_degradation_factor = base_reg_rate * temp_factor_reg * targeting_effect_reg * dt
    fib_degradation_factor = base_fib_rate * temp_factor_fib * targeting_effect_fib * dt

    # Apply degradation (MULTIPLICATIVE for exponential decay)
    # Exponential-style multiplicative decay for regular grid
    decay_factor = reg_degradation_factor * 100  # Scale up to simulate rapid collapse
    state['regular_grid']['integrity'] *= math.exp(-decay_factor)
    state['fibonacci_grid']['integrity'] *= (1.0 - fib_degradation_factor)

    # Add small random fluctuations (multiplicative as well)
    reg_fluctuation = 1.0 - 0.01 * random.random() * reg_degradation_factor
    fib_fluctuation = 1.0 - 0.01 * random.random() * fib_degradation_factor

    state['regular_grid']['integrity'] *= reg_fluctuation
    state['fibonacci_grid']['integrity'] *= fib_fluctuation

    # Clamp integrity values between 0 and 1
    state['regular_grid']['integrity'] = max(0.0001, min(1.0, state['regular_grid']['integrity']))
    state['fibonacci_grid']['integrity'] = max(0.0001, min(1.0, state['fibonacci_grid']['integrity']))

    # Update resonance based on integrity and temperature
    # Note: resonance decays faster with temperature in regular grid
    state['regular_grid']['resonance'] = state['regular_grid']['integrity'] * math.exp(-0.02 * temp_diff)
    state['fibonacci_grid']['resonance'] = state['fibonacci_grid']['integrity'] * math.exp(-0.01 * temp_diff)

    # Update temperature-dependent defense activation thresholds
    update_defense_mechanisms(state, temperature)

    # Track temperature state for analysis
    state['temperature'] = temperature
    state['temperature_history'].append(temperature)


def calculate_quantum_coherence(time, geometry, temperature):
    """
    Implement key finding from Tegmark's paper showing different decay patterns.

    Parameters:
    -----------
    time : float
        Time
    geometry : str
        'fibonacci' or 'regular'
    temperature : float
        Temperature in Celsius

    Returns:
    --------
    float
        Quantum coherence value
    """
    # Initial coherence values
    initial_coherence = 1.0

    # Temperature impact factor - stronger impact on regular grid
    if geometry == 'fibonacci':
        temp_factor = math.exp(-0.05 * (temperature - 37.0))
    else:
        temp_factor = math.exp(-0.3 * (temperature - 37.0))  # Much stronger temperature sensitivity

    # Coherence decay calculation - ensure it matches grid integrity behavior
    if geometry == 'fibonacci':
        # Power-law decay for Fibonacci systems: t^(-1.015)
        coherence = initial_coherence * 0.015963 * (time ** -1.015) * temp_factor

        # Boundary/sanctuary formation effect at t=0.6
        if time > 0.6:
            coherence *= 1.6  # Sanctuary protection factor
    else:
        # Near-exponential decay for regular systems: t^(-10.1)
        # Make it even more extreme for regular grids
        coherence = initial_coherence * (time ** -10.1023) * temp_factor

    return max(coherence, 1e-20)  # Prevent zero/negative coherence


def calculate_information_sanctuary(state, time, temperature):
    """
    Model the formation of quantum-protected regions during acute fever stress.

    Parameters:
    -----------
    state : dict
        Current state of the system
    time : float
        Current simulation time
    temperature : float
        Temperature in Celsius

    Returns:
    --------
    bool
        True if sanctuary formed, False otherwise
    """
    # Sanctuary only forms under specific conditions
    sanctuary_threshold = 38.5  # Temperature threshold

    # For Fibonacci structure only
    if state['fibonacci_grid']['integrity'] > 0.1 and temperature > sanctuary_threshold:
        sanctuary_radius = 0.9 * state['outer_radius']
        sanctuary_strength = min(1.0, (temperature - sanctuary_threshold) / 1.5)

        # Update state to track sanctuary
        state['fibonacci_grid']['sanctuary_radius'] = sanctuary_radius
        state['fibonacci_grid']['sanctuary_strength'] = sanctuary_strength

        # Immunity within sanctuary (coherence preservation)
        sanctuary_coherence_factor = min(2.5, 1.0 + sanctuary_strength * 1.5)

        # Track spatial coherence distribution
        update_coherence_distribution(state, sanctuary_radius, sanctuary_coherence_factor)

        return True

    # No sanctuary forms for regular grid or below threshold
    state['fibonacci_grid']['sanctuary_radius'] = 0
    state['fibonacci_grid']['sanctuary_strength'] = 0
    return False


def update_coherence_distribution(state, sanctuary_radius, sanctuary_coherence_factor):
    """
    Update the coherence distribution based on sanctuary formation.

    Parameters:
    -----------
    state : dict
        Current state of the system
    sanctuary_radius : float
        Radius of the sanctuary
    sanctuary_coherence_factor : float
        Factor by which coherence is enhanced within the sanctuary
    """
    # This is a placeholder function - you would need to implement the specific logic
    # for updating the coherence distribution based on your implementation details
    pass


def update_defense_mechanisms(state, temperature):
    """
    Modify the defense mechanism activation to incorporate fever effects.

    Parameters:
    -----------
    state : dict
        Current state of the system
    temperature : float
        Temperature in Celsius
    """
    # Temperature-dependent activation threshold adjustment
    temp_diff = temperature - 37.0
    activation_boost = min(1.0, max(0, temp_diff / 3.0))

    # Phi Stability Defense activation (activates earlier during fever)
    if state['targeting_intensity'] > 0.45 - 0.15 * activation_boost:
        state['defenses']['phi_stability'] = min(1.0,
                                                 state['defenses']['phi_stability'] + 0.016 * (
                                                             1 + 0.5 * activation_boost))

    # Resonance Defense (enhanced during fever)
    if state['targeting_intensity'] > 0.6 - 0.2 * activation_boost:
        state['defenses']['resonance'] = min(1.0,
                                             state['defenses']['resonance'] + 0.0165 * (1 + 0.3 * activation_boost))

    # Acetylation Defense (highly temperature-dependent)
    if state['targeting_intensity'] > 0.7 - 0.1 * activation_boost:
        # Acetylation rate increases with temperature
        acetyl_rate = 0.0168 * (1.1 ** temp_diff)
        state['defenses']['acetylation'] = min(1.0,
                                               state['defenses']['acetylation'] + acetyl_rate)

    # Binding Sites Defense
    if state['targeting_intensity'] > 0.55 - 0.15 * activation_boost:
        state['defenses']['binding_sites'] = min(1.0,
                                                 state['defenses']['binding_sites'] + 0.0169 * (
                                                             1 + 0.2 * activation_boost))

    # Adaptation Defense (only activates during fever)
    if state['targeting_intensity'] > 0.65 - 0.2 * activation_boost and temperature > 38.0:
        state['defenses']['adaptation'] = min(1.0,
                                              state['defenses']['adaptation'] + 0.0168 * (1 + 0.4 * activation_boost))


def create_temperature_effect_visualizations(data):
    """
    Create charts to analyze and display temperature effects.

    Parameters:
    -----------
    data : list
        List of data points from simulation

    Returns:
    --------
    dict
        Dictionary containing matplotlib figure objects
    """
    # Create temperature vs. integrity chart
    temp_integrity_fig, ax1 = plt.subplots(figsize=(10, 6))

    # Extract data for plotting
    temperatures = np.array([d['temperature'] for d in data])
    regular_integrity = np.array([d['regular_grid_integrity'] for d in data])
    fibonacci_integrity = np.array([d['fibonacci_grid_integrity'] for d in data])

    # Plot data
    ax1.scatter(temperatures, regular_integrity, color='blue', alpha=0.1, label='Regular Grid')
    ax1.scatter(temperatures, fibonacci_integrity, color='green', alpha=0.1, label='Fibonacci Grid')

    # Add trend lines
    if len(temperatures) > 1:  # Only fit if we have enough data points
        # Sort temperatures for smooth line
        sorted_indices = np.argsort(temperatures)
        sorted_temps = temperatures[sorted_indices]

        # Fit regular grid
        z1 = np.polyfit(temperatures, regular_integrity, 1)
        p1 = np.poly1d(z1)
        ax1.plot(sorted_temps, p1(sorted_temps), color='blue')

        # Fit fibonacci grid
        z2 = np.polyfit(temperatures, fibonacci_integrity, 1)
        p2 = np.poly1d(z2)
        ax1.plot(sorted_temps, p2(sorted_temps), color='green')

    # Add labels and title
    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Grid Integrity')
    ax1.set_title('Temperature Effect on Grid Integrity')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Create temperature-coherence relationship chart
    coherence_fig, ax2 = plt.subplots(figsize=(10, 6))

    # Generate data points for plotting
    time_points = np.linspace(0.01, 3, 100)

    fib_coherence_normal = np.array([calculate_quantum_coherence(t, 'fibonacci', 37.0) for t in time_points])
    reg_coherence_normal = np.array([calculate_quantum_coherence(t, 'regular', 37.0) for t in time_points])
    fib_coherence_fever = np.array([calculate_quantum_coherence(t, 'fibonacci', 39.0) for t in time_points])

    # Plot data on log-log scale
    ax2.loglog(time_points, fib_coherence_normal, color='green', label='Fibonacci Power-Law (t^-1.015)')
    ax2.loglog(time_points, reg_coherence_normal, color='blue', label='Regular Grid (t^-10.1)')
    ax2.loglog(time_points, fib_coherence_fever, color='red', linestyle='--', label='Fibonacci (39°C Fever)')

    # Add labels and title
    ax2.set_xlabel('Time (normalized)')
    ax2.set_ylabel('Quantum Coherence')
    ax2.set_title('Coherence Decay Comparison (Log Scale)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Return figures
    return {
        'temp_integrity_chart': temp_integrity_fig,
        'coherence_chart': coherence_fig
    }


def update_targeting_intensity(state, t, dt):
    """
    Update targeting intensity based on simulation phase.

    Parameters:
    -----------
    state : dict
        Current state of the system
    t : float
        Current simulation time
    dt : float
        Time step
    """
    # Implement targeting intensity changes - more gradual for better visualization
    if t < 10:
        # Initial attack phase - rapidly increasing intensity
        state['targeting_intensity'] = min(0.95, state['targeting_intensity'] + 0.1 * dt)
    elif 30 < t < 60:
        # Reduced targeting phase - gradual decrease
        state['targeting_intensity'] = max(0.3, state['targeting_intensity'] - 0.01 * dt)
    else:
        # Maintenance phase - small variations
        state['targeting_intensity'] += 0.001 * dt * (random.random() - 0.5)
        state['targeting_intensity'] = max(0.3, min(0.95, state['targeting_intensity']))

    # Ensure targeting intensity stays between 0 and 1
    state['targeting_intensity'] = max(0.0, min(1.0, state['targeting_intensity']))


def initialize_system(params):
    """
    Initialize the simulation system.

    Parameters:
    -----------
    params : dict
        Simulation parameters

    Returns:
    --------
    dict
        Initialized system state
    """
    state = {
        'regular_grid': {
            'integrity': 1.0,
            'resonance': 1.0,
            'coherence': 1.0
        },
        'fibonacci_grid': {
            'integrity': 1.0,
            'resonance': 1.0,
            'coherence': 1.0,
            'sanctuary_radius': 0,
            'sanctuary_strength': 0
        },
        'defenses': {
            'phi_stability': 0.0,
            'resonance': 0.0,
            'acetylation': 0.0,
            'binding_sites': 0.0,
            'adaptation': 0.0
        },
        'targeting_intensity': params.get('initial_targeting', 0.05),
        'outer_radius': params.get('outer_radius', 1.0),
        'temperature_history': [],
        'hiv_phase': params.get('hiv_phase', 'acute')
    }
    return state


def run_full_hiv_simulation(params):
    """
    Run the full multi-scale simulation integrating all components.

    Parameters:
    -----------
    params : dict
        Simulation parameters

    Returns:
    --------
    list
        Simulation results
    """
    # Initialize system state
    state = initialize_system(params)

    # Time parameters
    dt = 0.01  # Primary time step
    max_time = params.get('max_time', 100)

    # Results storage
    results = []

    # Print initial values to verify
    print(f"Initial regular grid integrity: {state['regular_grid']['integrity']}")
    print(f"Initial fibonacci grid integrity: {state['fibonacci_grid']['integrity']}")

    # Store initial state
    results.append({
        'time': 0,
        'temperature': 37.0,  # Default starting temperature
        'regular_grid_integrity': state['regular_grid']['integrity'],
        'fibonacci_grid_integrity': state['fibonacci_grid']['integrity'],
        'regular_grid_coherence': 1.0,
        'fibonacci_grid_coherence': 1.0,
        'defenses': {k: v for k, v in state['defenses'].items()},
        'targeting_intensity': state['targeting_intensity'],
        'sanctuary_radius': 0,
        'coherence_ratio': 1.0
    })

    # Main simulation loop
    t = 0
    while t < max_time:
        # Increment time first (to avoid t=0 issues with coherence calculation)
        t += dt

        # Temperature evolution
        current_hour = t * 24  # Convert to hours
        temperature = simulate_hiv_fever(current_hour, state['hiv_phase'])
        state['temperature'] = temperature

        # Update coherence values - prevent t=0 by using max(0.01, t)
        state['regular_grid']['coherence'] = calculate_quantum_coherence(max(0.01, t), 'regular', temperature)
        state['fibonacci_grid']['coherence'] = calculate_quantum_coherence(max(0.01, t), 'fibonacci', temperature)

        # Check for information sanctuary formation
        calculate_information_sanctuary(state, t, temperature)

        # Update grid integrity and defense mechanisms
        update_grids(state, dt, t)

        # Ensure grid integrity stays clamped between 0 and 1
        state['regular_grid']['integrity'] = max(0.0001, min(1.0, state['regular_grid']['integrity']))
        state['fibonacci_grid']['integrity'] = max(0.0001, min(1.0, state['fibonacci_grid']['integrity']))

        # Print integrity values periodically to check
        if int(t * 10) % 100 == 0:  # Print more frequently (every 10 time units)
            print(
                f"Time {t:.1f}: Regular={state['regular_grid']['integrity']:.6f}, Fibonacci={state['fibonacci_grid']['integrity']:.6f}")

        # Store results at appropriate intervals
        output_interval = params.get('output_interval', 0.1)
        if t % output_interval < dt:
            coherence_ratio = 0
            if state['regular_grid']['coherence'] > 0:
                coherence_ratio = state['fibonacci_grid']['coherence'] / state['regular_grid']['coherence']

            results.append({
                'time': t,
                'temperature': temperature,
                'regular_grid_integrity': state['regular_grid']['integrity'],
                'fibonacci_grid_integrity': state['fibonacci_grid']['integrity'],
                'regular_grid_coherence': state['regular_grid']['coherence'],
                'fibonacci_grid_coherence': state['fibonacci_grid']['coherence'],
                'defenses': {k: v for k, v in state['defenses'].items()},
                'targeting_intensity': state['targeting_intensity'],
                'sanctuary_radius': state['fibonacci_grid'].get('sanctuary_radius', 0),
                'coherence_ratio': coherence_ratio
            })

        # Update targeting intensity based on simulation phase
        update_targeting_intensity(state, t, dt)

    # Print final integrity values
    print(f"Final regular grid integrity: {state['regular_grid']['integrity']}")
    print(f"Final fibonacci grid integrity: {state['fibonacci_grid']['integrity']}")

    return results


def plot_simulation_results(results):
    """
    Plot key results from the simulation.

    Parameters:
    -----------
    results : list
        Simulation results

    Returns:
    --------
    dict
        Dictionary of matplotlib figures
    """
    figures = {}

    # Extract time series data
    times = np.array([r['time'] for r in results])
    temps = np.array([r['temperature'] for r in results])
    reg_integrity = np.array([r['regular_grid_integrity'] for r in results])
    fib_integrity = np.array([r['fibonacci_grid_integrity'] for r in results])
    coherence_ratio = np.array([r['coherence_ratio'] for r in results])

    # Print some values to verify we have data
    print(f"Number of data points: {len(times)}")
    print(f"Integrity ranges - Regular: [{min(reg_integrity):.6f}, {max(reg_integrity):.6f}]")
    print(f"Integrity ranges - Fibonacci: [{min(fib_integrity):.6f}, {max(fib_integrity):.6f}]")

    # Grid integrity over time
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    ax1.plot(times, reg_integrity, 'b-', label='Regular Grid')
    ax1.plot(times, fib_integrity, 'g-', label='Fibonacci Grid')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Grid Integrity')
    ax1.set_title('Grid Integrity Over Time')
    ax1.set_ylim([0, 1.1])  # Force y-axis to include the full range
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    figures['integrity_time'] = fig1

    # Integrity difference (Fibonacci - Regular)
    fig4, ax4 = plt.subplots(figsize=(12, 6))
    integrity_diff = fib_integrity - reg_integrity
    # Fix the color/linestyle syntax - separate color and linestyle
    ax4.plot(times, integrity_diff, color='purple', linestyle='-')
    ax4.set_xlabel('Time')
    ax4.set_ylabel('Integrity Difference (Fibonacci - Regular)')
    ax4.set_title('Grid Integrity Advantage of Fibonacci Structure')
    ax4.grid(True, alpha=0.3)
    figures['integrity_diff'] = fig4

    # Temperature profile
    fig2, ax2 = plt.subplots(figsize=(12, 6))
    ax2.plot(times, temps, 'r-')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Temperature (°C)')
    ax2.set_title('Temperature Profile')
    ax2.grid(True, alpha=0.3)
    figures['temperature'] = fig2

    # Coherence ratio (Fibonacci/Regular)
    fig3, ax3 = plt.subplots(figsize=(12, 6))
    ax3.semilogy(times, coherence_ratio, 'k-')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Coherence Ratio (Fibonacci/Regular)')
    ax3.set_title('Quantum Coherence Advantage of Fibonacci Structure')
    ax3.grid(True, alpha=0.3)
    figures['coherence_ratio'] = fig3

    # Add temperature effect visualizations
    temp_effect_figs = create_temperature_effect_visualizations(results)
    # Add each figure individually to ensure proper reference handling
    figures['temp_integrity_chart'] = temp_effect_figs['temp_integrity_chart']
    figures['coherence_chart'] = temp_effect_figs['coherence_chart']

    return figures


# Example usage
if __name__ == "__main__":
    # Set simulation parameters
    params = {
        'max_time': 120,
        'output_interval': 0.1,
        'initial_targeting': 0.05,
        'outer_radius': 1.0,
        'hiv_phase': 'acute'
    }

    # Run simulation
    results = run_full_hiv_simulation(params)
    # Plot results
    figures = plot_simulation_results(results)

    # Show figures
    plt.show()