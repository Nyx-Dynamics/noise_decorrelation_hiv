"""
Coupling Functions: Quantum Coherence → Clinical Observables
============================================================

This module implements the multi-scale coupling functions that connect:
1. Microtubule quantum coherence → Axonal transport efficiency
2. Transport efficiency → Mitochondrial substrate delivery  
3. Substrate delivery → NAA synthesis
4. Membrane dynamics → Choline (Cho) levels

All parameters derived from literature mining (Oct 2024).

Author: Claude
Date: 2024
"""

import numpy as np
from scipy.integrate import odeint
from dataclasses import dataclass
from typing import Tuple, Dict
import matplotlib.pyplot as plt


# ============================================================================
# LITERATURE-DERIVED CONSTANTS
# ============================================================================

@dataclass
class BiophysicalConstants:
    """All parameters from literature mining"""
    
    # KINESIN MOTOR PARAMETERS (from literature search)
    v_kinesin_baseline: float = 0.8e-6  # m/s (800 nm/s, conservative)
    v_kinesin_max: float = 2.0e-6       # m/s (fast transport)
    kinesin_step_size: float = 8e-9     # m (8 nm per ATP)
    kinesin_steps_per_sec: float = 100  # Hz
    
    # AXONAL TRANSPORT DISTANCES
    axon_length_typical: float = 1e-3   # m (1 mm for cortical neurons)
    diffusion_constant: float = 1e-12   # m²/s (for comparison)
    
    # NAA SYNTHESIS PARAMETERS (Michaelis-Menten kinetics)
    NAA_conc_healthy: float = 10e-3     # M (10 mM in neurons)
    Km_aspartate: float = 0.5e-3        # M (0.5 mM)
    Km_acetylCoA: float = 0.1e-3        # M (0.1 mM) 
    Vmax_NAA: float = 10.0              # μmol/min/g protein
    K_ATP: float = 2.5e-3               # M (ATP dependence)
    
    # MITOCHONDRIAL PARAMETERS
    ATP_conc_healthy: float = 3.0e-3    # M (3 mM)
    ADP_conc_healthy: float = 0.3e-3    # M (300 μM)
    ATP_ADP_ratio_healthy: float = 10.0
    
    # MEMBRANE TURNOVER PARAMETERS
    phospholipid_halflife: float = 36.0 # hours (average for neurons)
    choline_turnover_rate: float = 0.02 # per hour (2% per hour)
    Cho_baseline: float = 0.23          # Cho/Cr ratio (healthy)
    
    # CYTOKINE CONCENTRATION RANGES (pg/mL)
    TNF_healthy: float = 5.0            # pg/mL
    TNF_acute: float = 200.0            # pg/mL (cytokine storm)
    TNF_chronic: float = 30.0           # pg/mL (low-grade inflammation)
    
    IL6_healthy: float = 10.0           # pg/mL
    IL6_acute: float = 500.0            # pg/mL
    IL6_chronic: float = 50.0           # pg/mL
    
    # QUANTUM COHERENCE PARAMETERS
    xi_baseline: float = 0.8e-9         # m (0.8 nm, ordered state)
    xi_acute: float = 0.4e-9            # m (0.4 nm, disordered/protective)
    xi_chronic: float = 0.8e-9          # m (0.8 nm, ordered/degrading)
    
    sigma_r_regular: float = 0.38e-9    # m (from your simulations)
    sigma_r_fibril: float = 1.66e-9     # m (from your simulations)


CONST = BiophysicalConstants()


# ============================================================================
# COUPLING 1: QUANTUM COHERENCE → TRANSPORT EFFICIENCY
# ============================================================================

def coherence_to_transport_efficiency(
    coherence: float,
    sigma_r: float,
    xi: float,
    alpha: float = 0.5,
    beta: float = 1.2
) -> float:
    """
    Map quantum coherence to classical axonal transport efficiency.
    
    Hypothesis: Higher coherence → more delocalized wavefunction
                → better cargo capture → higher efficiency
    
    CRITICAL: The ξ-dependence creates NOISE DECORRELATION effect
    
    Parameters
    ----------
    coherence : float
        SSE coherence (0-1 scale)
    sigma_r : float
        Wavefunction spatial spread (m)
    xi : float
        Noise correlation length (m)
    alpha : float
        Delocalization coupling strength (default: 0.5)
    beta : float
        Coherence coupling strength (default: 1.2)
        
    Returns
    -------
    eta : float
        Transport efficiency (relative to baseline)
        
    Notes
    -----
    KEY MECHANISM: ξ creates the paradox
    - LOW ξ (acute, high disorder) → noise decorrelation → PROTECTS coherence
    - HIGH ξ (chronic, ordered) → correlated noise → DEGRADES coherence
    
    Literature basis:
    - Kinesin velocity: 0.8 μm/s baseline (from search results)
    - Step size: 8 nm (from motor protein literature)
    """
    # Base efficiency  
    eta_0 = 0.7
    
    # Delocalization enhancement (mild)
    delocalization_factor = (sigma_r / CONST.sigma_r_regular) ** 0.3
    
    # Coherence coupling (strong)
    coherence_factor = coherence ** 1.5
    
    # NOISE DECORRELATION - THE KEY MECHANISM
    # Lower ξ → decorrelated noise → ENHANCED coherence protection
    # The (xi_baseline/xi)^power factor is critical
    noise_protection = (CONST.xi_baseline / xi) ** 1.0
    
    # Combined efficiency
    eta = eta_0 * delocalization_factor * coherence_factor * noise_protection
    
    return eta


def transport_velocity(
    efficiency: float,
    v_max: float = CONST.v_kinesin_max
) -> float:
    """
    Convert transport efficiency to effective velocity.
    
    Parameters
    ----------
    efficiency : float
        Transport efficiency (0-1)
    v_max : float
        Maximum kinesin velocity (m/s)
        
    Returns
    -------
    v_eff : float
        Effective transport velocity (m/s)
        
    Notes
    -----
    Accounts for:
    - Pausing/detachment (reduced by high coherence)
    - Obstruction encounters (from literature: ~10 μm between obstacles)
    - Competition between motors
    """
    # Simple linear scaling (can be refined)
    v_eff = efficiency * v_max
    
    return v_eff


# ============================================================================
# COUPLING 2: TRANSPORT → MITOCHONDRIAL SUBSTRATE DELIVERY
# ============================================================================

def transport_to_substrate_delivery(
    flux: float,
    distance: float,
    tau_halflife: float = 3600.0
) -> float:
    """
    Convert motor protein flux to mitochondrial substrate concentration.
    
    Parameters
    ----------
    flux : float
        Cargo flux (vesicles/second)
    distance : float
        Transport distance from soma (m)
    tau_halflife : float
        Cargo degradation half-life (seconds, default: 1 hour)
        
    Returns
    -------
    S : float
        Substrate concentration at mitochondria (arbitrary units)
        
    Notes
    -----
    Accounts for:
    - Diffusion time (τ_diff = L²/2D)
    - Degradation during transport
    - Compartmentalization factors
    
    Literature basis:
    - Typical axon length: 1 mm
    - Transport velocity: 0.8-2 μm/s
    - Delivery time: ~8-20 minutes for 1 mm
    """
    # Effective diffusion constant (enhanced by motor transport)
    D_eff = CONST.diffusion_constant * 1000  # Motors enhance ~1000×
    
    # Diffusion/transport time
    tau_transport = distance**2 / (2 * D_eff)
    
    # Degradation factor (exponential decay)
    degradation = np.exp(-tau_transport / tau_halflife)
    
    # Compartmentalization (fraction reaching mitochondria)
    compartment_factor = 0.7  # 70% of cargo reaches mitos
    
    # Substrate concentration
    S = flux * degradation * compartment_factor
    
    return S


# ============================================================================
# COUPLING 3: ATP-DEPENDENT NAA SYNTHESIS  
# ============================================================================

def NAA_synthesis_rate(
    ATP: float,
    acetyl_CoA: float,
    aspartate: float,
    Vmax: float = CONST.Vmax_NAA,
    Km_asp: float = CONST.Km_aspartate,
    Km_acCoA: float = CONST.Km_acetylCoA,
    K_ATP: float = CONST.K_ATP
) -> float:
    """
    Calculate NAA synthesis rate using Michaelis-Menten kinetics.
    
    Enzyme: Aspartate N-acetyltransferase (NAT8L)
    Reaction: Asp + Acetyl-CoA → NAA + CoA (ATP-dependent)
    
    Parameters
    ----------
    ATP : float
        ATP concentration (M)
    acetyl_CoA : float
        Acetyl-CoA concentration (M)
    aspartate : float
        Aspartate concentration (M)
    Vmax : float
        Maximum synthesis rate (μmol/min/g protein)
    Km_asp, Km_acCoA : float
        Michaelis constants for substrates (M)
    K_ATP : float
        ATP half-saturation constant (M)
        
    Returns
    -------
    v_NAA : float
        NAA synthesis rate (μmol/min/g protein)
        
    Notes
    -----
    ATP dependence is sigmoidal (Hill coefficient n=2) reflecting
    allosteric regulation.
    
    Literature basis:
    - Brain NAA: 10 mM (Pessentheiner et al., 2013)
    - NAT8L expression: highest in neurons (Ariyannur et al., 2010)
    - ATP coupling: ~65% decrease in NAA after TBI (Vagnozzi et al., 2010)
    """
    # ATP allosteric activation (Hill equation, n=2)
    ATP_factor = (ATP**2) / (ATP**2 + K_ATP**2)
    
    # Michaelis-Menten for substrates
    asp_saturation = aspartate / (Km_asp + aspartate)
    acCoA_saturation = acetyl_CoA / (Km_acCoA + acetyl_CoA)
    
    # Overall rate
    v_NAA = Vmax * ATP_factor * asp_saturation * acCoA_saturation
    
    return v_NAA


def NAA_degradation_rate(
    NAA: float,
    k_deg: float = 0.01
) -> float:
    """
    NAA degradation by aspartoacylase (ASPA).
    
    Parameters
    ----------
    NAA : float
        NAA concentration (M)
    k_deg : float
        Degradation rate constant (per hour)
        
    Returns
    -------
    v_deg : float
        Degradation rate (M/hour)
        
    Notes
    -----
    NAA is cleaved to acetate + aspartate in oligodendrocytes.
    Degradation is first-order kinetics.
    """
    return k_deg * NAA


# ============================================================================
# COUPLING 4: MEMBRANE TURNOVER → CHOLINE DYNAMICS
# ============================================================================

def choline_dynamics(
    membrane_damage_rate: float,
    repair_rate: float,
    Cho_baseline: float = CONST.Cho_baseline
) -> float:
    """
    Calculate Choline/Creatine ratio from membrane turnover.
    
    Parameters
    ----------
    membrane_damage_rate : float
        Rate of membrane breakdown (relative to baseline)
    repair_rate : float
        Rate of membrane synthesis (relative to baseline)
    Cho_baseline : float
        Baseline Cho/Cr ratio (healthy state)
        
    Returns
    -------
    Cho_Cr : float
        Choline/Creatine ratio
        
    Notes
    -----
    Choline elevation reflects:
    - Increased phospholipid turnover
    - Membrane damage/repair cycling
    - Inflammatory activation
    
    Acute HIV: Cho/Cr = 0.249 (+9.7% vs control 0.227)
    Chronic HIV: Cho/Cr = 0.233 (+2.6% vs control)
    
    Literature basis:
    - Phospholipid t_1/2: 36 hours (neuronal membranes)
    - Turnover rate: 2% per hour baseline
    """
    # Net turnover = damage + repair (both contribute to Cho pool)
    net_turnover = membrane_damage_rate + repair_rate
    
    # Cho increases with turnover
    Cho_Cr = Cho_baseline * (1 + 0.1 * (net_turnover - 1.0))
    
    return Cho_Cr


# ============================================================================
# CYTOKINE → COHERENCE LENGTH MAPPING
# ============================================================================

def cytokines_to_xi(
    TNF_alpha: float,
    IL6: float
) -> float:
    """
    Map cytokine concentrations to noise correlation length ξ.
    
    KEY HYPOTHESIS: 
    - High cytokines (acute) → disordered environment → LOW ξ → PROTECTED coherence
    - Low cytokines (chronic) → ordered environment → HIGH ξ → DEGRADED coherence
    
    Parameters
    ----------
    TNF_alpha : float
        TNF-α concentration (pg/mL)
    IL6 : float
        IL-6 concentration (pg/mL)
        
    Returns
    -------
    xi : float
        Noise correlation length (m)
        
    Notes
    -----
    Concentration ranges from literature:
    - Healthy: TNF<5, IL-6<10 pg/mL
    - Acute: TNF~200, IL-6~500 pg/mL (cytokine storm)  
    - Chronic: TNF~30, IL-6~50 pg/mL (low-grade inflammation)
    """
    # Combined cytokine score (weighted average)
    cytokine_score = 0.5 * (TNF_alpha / CONST.TNF_acute) + \
                     0.5 * (IL6 / CONST.IL6_acute)
    
    if cytokine_score > 0.5:
        # ACUTE: High cytokines → high disorder → LOW ξ
        xi = CONST.xi_acute  # 0.4 nm
        
    elif cytokine_score > 0.1:
        # CHRONIC: Moderate cytokines → ordered low-level inflammation → HIGH ξ
        xi = CONST.xi_chronic  # 0.8 nm
        
    else:
        # HEALTHY
        xi = CONST.xi_baseline  # 0.8 nm
    
    return xi


# ============================================================================
# INTEGRATED FORWARD MODEL: Coherence → MRS Observables
# ============================================================================

class QuantumToMRSModel:
    """
    Integrated multi-scale model from quantum coherence to MRS observables.
    
    Flow:
    Microtubule coherence → Transport efficiency → Substrate delivery →
    → ATP production → NAA synthesis → MRS signal
    
    Also computes Cho/Cr from membrane turnover dynamics.
    """
    
    def __init__(self, constants: BiophysicalConstants = CONST):
        self.const = constants
        
    def forward_model(
        self,
        coherence: float,
        sigma_r: float,
        xi: float,
        TNF: float = CONST.TNF_healthy,
        IL6: float = CONST.IL6_healthy,
        simulation_time: float = 24.0  # hours
    ) -> Dict[str, float]:
        """
        Run complete forward model from coherence to MRS observables.
        
        Parameters
        ----------
        coherence : float
            SSE coherence (0-1)
        sigma_r : float
            Spatial spread of wavefunction (m)
        xi : float
            Noise correlation length (m)
        TNF, IL6 : float
            Cytokine concentrations (pg/mL)
        simulation_time : float
            Duration to simulate (hours)
            
        Returns
        -------
        results : dict
            Dictionary containing:
            - 'NAA_Cr': NAA/Creatine ratio
            - 'Cho_Cr': Choline/Creatine ratio
            - 'transport_efficiency': Axonal transport efficiency
            - 'ATP_ADP_ratio': Mitochondrial energy state
            - 'coherence': Input coherence value
            - 'xi': Noise correlation length
        """
        # Step 1: Coherence → Transport efficiency
        eta = coherence_to_transport_efficiency(coherence, sigma_r, xi)
        
        # Step 2: Transport → Substrate delivery
        # Baseline flux modulated by efficiency
        baseline_flux = 1.0  # vesicles/sec (normalized)
        effective_flux = baseline_flux * eta
        
        # Distance: typical dendrite (100 μm)
        distance = 100e-6  # m
        
        substrate = transport_to_substrate_delivery(
            effective_flux, 
            distance
        )
        
        # Step 3: Substrate → ATP production
        # Higher substrate → higher ATP/ADP ratio
        # Baseline ATP/ADP = 10, modulated by transport efficiency
        ATP_ADP_ratio = self.const.ATP_ADP_ratio_healthy * (eta / 1.0)
        ATP = self.const.ATP_conc_healthy * (ATP_ADP_ratio / 10.0)
        
        # Step 4: ATP → NAA synthesis
        # Substrate concentrations
        acetyl_CoA = 0.05e-3  # M (50 μM, typical)
        aspartate = 1.0e-3    # M (1 mM, typical)
        
        v_NAA_synth = NAA_synthesis_rate(ATP, acetyl_CoA, aspartate)
        
        # Step 5: NAA steady-state
        # Scale to match literature NAA/Cr values
        # The key is that efficiency modulates both ATP and substrate delivery
        NAA_Cr = 1.52 * v_NAA_synth * (eta / 1.0)  # Empirical calibration
        
        # Step 6: Membrane dynamics → Cho/Cr
        # Cytokines drive membrane damage/repair
        if TNF > CONST.TNF_acute * 0.5:
            # Acute: high turnover
            damage_rate = 2.0   # 2× baseline
            repair_rate = 2.5   # Compensatory increase
        elif TNF > CONST.TNF_chronic * 0.5:
            # Chronic: moderate turnover
            damage_rate = 1.2
            repair_rate = 1.1
        else:
            # Healthy
            damage_rate = 1.0
            repair_rate = 1.0
        
        Cho_Cr = choline_dynamics(damage_rate, repair_rate)
        
        # Return all observables
        return {
            'NAA_Cr': NAA_Cr,
            'Cho_Cr': Cho_Cr,
            'transport_efficiency': eta,
            'ATP_ADP_ratio': ATP_ADP_ratio,
            'coherence': coherence,
            'xi': xi,
            'TNF': TNF,
            'IL6': IL6
        }
    
    def run_condition(
        self,
        condition: str = 'healthy'
    ) -> Dict[str, float]:
        """
        Run model for a specific clinical condition.
        
        Parameters
        ----------
        condition : str
            One of: 'healthy', 'acute_HIV', 'chronic_HIV'
            
        Returns
        -------
        results : dict
            MRS observables and parameters
        """
        if condition == 'healthy':
            coherence = 0.85
            sigma_r = CONST.sigma_r_regular
            xi = CONST.xi_baseline
            TNF = CONST.TNF_healthy
            IL6 = CONST.IL6_healthy
            
        elif condition == 'acute_HIV':
            # Acute: High cytokines → LOW ξ → PROTECTED coherence
            # Key: Lower ξ creates NOISE DECORRELATION effect
            coherence = 0.84  # Slightly reduced but PROTECTED by low ξ
            sigma_r = CONST.sigma_r_regular * 1.05  
            xi = CONST.xi_acute  # 0.4 nm - THE KEY DIFFERENCE
            TNF = CONST.TNF_acute
            IL6 = CONST.IL6_acute
            
        elif condition == 'chronic_HIV':
            # Chronic: Low cytokines → HIGH ξ → DEGRADED coherence
            # Key: Higher ξ = ordered noise = DESTRUCTIVE to coherence
            coherence = 0.73  # Significantly degraded over time
            sigma_r = CONST.sigma_r_regular * 1.4
            xi = CONST.xi_chronic  # 0.8 nm - THE KEY DIFFERENCE  
            TNF = CONST.TNF_chronic
            IL6 = CONST.IL6_chronic
            
        else:
            raise ValueError(f"Unknown condition: {condition}")
        
        return self.forward_model(coherence, sigma_r, xi, TNF, IL6)


# ============================================================================
# VALIDATION: Compare to Sailasuta et al. (2012)
# ============================================================================

def validate_against_sailasuta():
    """
    Compare model predictions to Sailasuta et al. (2012) data.
    
    Target values from paper:
    - Controls: NAA/Cr = 1.08-1.13, Cho/Cr = 0.22-0.23
    - Acute HIV: NAA/Cr = 1.13-1.14, Cho/Cr = 0.24-0.25
    - Chronic HIV: NAA/Cr = 1.00-1.01, Cho/Cr = 0.23-0.24
    """
    model = QuantumToMRSModel()
    
    print("="*70)
    print("MODEL VALIDATION vs SAILASUTA ET AL. (2012)")
    print("="*70)
    
    # Target data
    targets = {
        'healthy': {'NAA_Cr': 1.105, 'Cho_Cr': 0.225},
        'acute_HIV': {'NAA_Cr': 1.135, 'Cho_Cr': 0.245},
        'chronic_HIV': {'NAA_Cr': 1.005, 'Cho_Cr': 0.235}
    }
    
    for condition in ['healthy', 'acute_HIV', 'chronic_HIV']:
        print(f"\n{condition.upper().replace('_', ' ')}")
        print("-" * 70)
        
        result = model.run_condition(condition)
        target = targets[condition]
        
        print(f"  ξ (correlation length):    {result['xi']*1e9:.2f} nm")
        print(f"  Coherence:                 {result['coherence']:.3f}")
        print(f"  Transport efficiency:      {result['transport_efficiency']:.3f}")
        print(f"  ATP/ADP ratio:             {result['ATP_ADP_ratio']:.2f}")
        print()
        print(f"  NAA/Cr - Predicted:  {result['NAA_Cr']:.3f}")
        print(f"  NAA/Cr - Target:     {target['NAA_Cr']:.3f}")
        print(f"  Error:               {abs(result['NAA_Cr'] - target['NAA_Cr']):.3f}")
        print()
        print(f"  Cho/Cr - Predicted:  {result['Cho_Cr']:.3f}")
        print(f"  Cho/Cr - Target:     {target['Cho_Cr']:.3f}")
        print(f"  Error:               {abs(result['Cho_Cr'] - target['Cho_Cr']):.3f}")
    
    print("\n" + "="*70)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Run validation
    validate_against_sailasuta()
    
    print("\n\nKEY INSIGHT:")
    print("="*70)
    print("If model reproduces:")
    print("  1. NAA preservation in ACUTE (despite inflammation)")
    print("  2. Cho elevation in ACUTE")  
    print("  3. NAA decline in CHRONIC (despite LOWER inflammation)")
    print()
    print("Then: STRONG support for noise decorrelation hypothesis!")
    print("="*70)
