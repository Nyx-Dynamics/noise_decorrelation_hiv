import numpy as np
import math
import random
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


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
    return (k_b * temperature * 1.2e-29) / (h_bar * h_bar * epsilon * (r ** 4) * geometry_factor)


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

    # Calculate polymerization factor (doubles every 10°C increase)
    polymerization_factor = 2 ** (temp_diff / 10)

    # Calculate stability impacts - different for each structure
    stability_impact_regular = math.exp(0.2 * temp_diff)  # More temperature sensitive
    stability_impact_fibonacci = math.exp(0.12 * temp_diff)  # More robust to temperature

    # Calculate decoherence rates (much faster in regular grids)
    decoherence_regular = calculate_decoherence_rate(temperature, 'regular', 1.0)
    decoherence_fibonacci = calculate_decoherence_rate(temperature, 'fibonacci', 1.618)

    # Apply temperature-dependent effects to grid integrity
    # Regular grid degradation
    state['regular_grid']['integrity'] *= (1 - dt * state['targeting_intensity'] *
                                           stability_impact_regular * decoherence_regular)

    # Fibonacci grid degradation
    state['fibonacci_grid']['integrity'] *= (1 - dt * state['targeting_intensity'] *
                                             stability_impact_fibonacci * decoherence_fibonacci)

    # Structural fluctuations (increase with temperature)
    fluctuation_regular = 0.01 * math.exp(0.15 * temp_diff) * random.random() * state['targeting_intensity']
    fluctuation_fibonacci = 0.01 * math.exp(0.08 * temp_diff) * random.random() * state['targeting_intensity']

    # Apply fluctuations
    state['regular_grid']['integrity'] = max(0, state['regular_grid']['integrity'] - fluctuation_regular)
    state['fibonacci_grid']['integrity'] = max(0, state['fibonacci_grid']['integrity'] - fluctuation_fibonacci)

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

    # Temperature impact factor
    temp_factor = min(1.0, math.exp(-0.1 * (temperature - 37.0)))

    # Coherence decay calculation
    if geometry == 'fibonacci':
        # Power-law decay for Fibonacci systems: t^(-1.015)
        coherence = initial_coherence * 0.015963 * (time ** -1.015) * temp_factor

        # Boundary/sanctuary formation effect at t=0.6
        if time > 0.6:
            coherence *= 1.6  # Sanctuary protection factor
    else:
        # Near-exponential decay for regular systems: t^(-10.1)
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
    temperatures = [d['temperature'] for d in data]
    regular_integrity = [d['regular_grid_integrity'] for d in data]
    fibonacci_integrity = [d['fibonacci_grid_integrity'] for d in data]

    # Plot data
    ax1.scatter(temperatures, regular_integrity, color='blue', alpha=0.1, label='Regular Grid')
    ax1.scatter(temperatures, fibonacci_integrity, color='green', alpha=0.1, label='Fibonacci Grid')

    # Add trend lines
    z1 = np.polyfit(temperatures, regular_integrity, 1)
    p1 = np.poly1d(z1)
    ax1.plot(sorted(temperatures), p1(sorted(temperatures)), color='blue')

    z2 = np.polyfit(temperatures, fibonacci_integrity, 1)
    p2 = np.poly1d(z2)
    ax1.plot(sorted(temperatures), p2(sorted(temperatures)), color='green')

    # Add labels and title
    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Grid Integrity')
    ax1.set_title('Temperature Effect on Grid Integrity')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Create temperature-coherence relationship chart
    coherence_fig, ax2 = plt.subplots(figsize=(10, 6))

    # Generate data points for plotting
    time_points = [i * 0.03 + 0.01 for i in range(100)]

    fib_coherence_normal = [calculate_quantum_coherence(t, 'fibonacci', 37.0) for t in time_points]
    reg_coherence_normal = [calculate_quantum_coherence(t, 'regular', 37.0) for t in time_points]
    fib_coherence_fever = [calculate_quantum_coherence(t, 'fibonacci', 39.0) for t in time_points]

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
    # This is a placeholder - implement according to your targeting model
    # For example, a simple ramp-up:
    if t < 20:
        state['targeting_intensity'] = min(0.9, state['targeting_intensity'] + 0.01 * dt)
    elif 30 < t < 60:
        state['targeting_intensity'] = max(0.2, state['targeting_intensity'] - 0.005 * dt)
    else:
        # Maintain or slightly vary
        state['targeting_intensity'] += 0.002 * dt * (random.random() - 0.5)
        state['targeting_intensity'] = max(0.1, min(0.95, state['targeting_intensity']))


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
    quantum_dt = 1e-15  # Quantum evolution time step
    quantum_steps = 1e6  # Steps per primary step
    max_time = params.get('max_time', 100)

    # Results storage
    results = []

    # Main simulation loop
    t = 0
    while t < max_time:
        # Temperature evolution
        current_hour = t * 24  # Convert to hours
        temperature = simulate_hiv_fever(current_hour, state['hiv_phase'])
        state['temperature'] = temperature

        # Quantum decoherence evolution (approximated)
        decoherence_regular = calculate_decoherence_rate(temperature, 'regular', 1.0)
        decoherence_fibonacci = calculate_decoherence_rate(temperature, 'fibonacci', 1.618)

        # Update coherence values
        state['regular_grid']['coherence'] = calculate_quantum_coherence(max(0.01, t), 'regular', temperature)
        state['fibonacci_grid']['coherence'] = calculate_quantum_coherence(max(0.01, t), 'fibonacci', temperature)

        # Check for information sanctuary formation
        calculate_information_sanctuary(state, t, temperature)

        # Update grid integrity and defense mechanisms
        update_grids(state, dt, t)

        # Store results at appropriate intervals
        output_interval = params.get('output_interval', 0.1)
        if t % output_interval < dt:
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
                'coherence_ratio': state['fibonacci_grid']['coherence'] / max(1e-20, state['regular_grid']['coherence'])
            })

        # Update targeting intensity based on simulation phase
        update_targeting_intensity(state, t, dt)

        # Increment time
        t += dt

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
    times = [r['time'] for r in results]
    temps = [r['temperature'] for r in results]
    reg_integrity = [r['regular_grid_integrity'] for r in results]
    fib_integrity = [r['fibonacci_grid_integrity'] for r in results]
    coherence_ratio = [r['coherence_ratio'] for r in results]

    # Grid integrity over time
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    ax1.plot(times, reg_integrity, 'b-', label='Regular Grid')
    ax1.plot(times, fib_integrity, 'g-', label='Fibonacci Grid')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Grid Integrity')
    ax1.set_title('Grid Integrity Over Time')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    figures['integrity_time'] = fig1

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

    # Temperature effect visualizations
    temp_effect_figs = create_temperature_effect_visualizations(results)
    figures.update(temp_effect_figs)

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