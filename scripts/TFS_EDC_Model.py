#!/usr/bin/env python3
"""
Turbulent Flame Speed Calculation for Hydrogen-Air Combustion (Using simulate_flame_speed Module)
----------------------------------------------------------------------------------------------------
This script imports the laminar flame speed solver from a separate module
(laminar_flame.py) that defines:

    def simulate_flame_speed(phi, p, T, mech='h2o2', transport_model='multicomponent'):
         ...
         return flame_speed_cm  # in cm/s

For each equivalence ratio (ϕ) in a scan, the laminar flame speed is computed
(using the provided input signature). The value is then converted to m/s and used
in a flame stretch model to predict the turbulent flame speed as follows:

    S_T = Sₗ₀ · (1 – L_M * (u_prime/L)) · [1 + A * (u_prime/Sₗ₀)^n]

with typical parameters:
    L_M = 0.002 m, A = 0.966, n = 1.164

Experimental data (in cm/s) are provided and the computed turbulent flame speeds
are compared against them. The results are exported to a CSV file and a high-resolution
(4K) figure is saved.

References:
  - Magnussen, B. F. (1982). "A Simple Eddy Dissipation Concept Model for Turbulent Combustion."
  - Cantera laminar flame speed examples.
  - Experimental data from the referenced study.
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import pandas as pd
from scipy.interpolate import interp1d

# Import the laminar flame solver from the module.
# The simulate_flame_speed function has the signature:
#   simulate_flame_speed(phi, p, T, mech='h2o2', transport_model='multicomponent')
from LFS import simulate_flame_speed


# -----------------------------
# Part 1: Compute Turbulent Flame Speed with Flame Stretch Effects
# -----------------------------
def compute_turbulent_flame_speed(S_l0, u_prime, L, L_M=0.002, A=0.966, n=1.164):
    """
    Compute the turbulent flame speed using a flame stretch model.

    The model is:
       S_T = S_l0 * f_stretch * f_area,
    where:
       f_stretch = 1 - L_M * (u_prime / L)
       f_area    = 1 + A * (u_prime / S_l0)^n

    Parameters:
      S_l0    : Laminar flame speed [m/s]
      u_prime : Turbulence intensity [m/s]
      L       : Integral length scale [m]
      L_M     : Markstein length [m] (default 0.002 m)
      A       : Empirical constant for flame wrinkling (default 0.966)
      n       : Empirical exponent (default 1.164)

    Returns:
      S_T : Predicted turbulent flame speed [m/s]
    """
    f_stretch = 1 - L_M * (u_prime / L)
    f_stretch = max(f_stretch, 0.0)  # prevent negative factor
    f_area = 1 + A * (u_prime / S_l0) ** n
    S_T = S_l0 * f_stretch * f_area
    print(f"  f_stretch = {f_stretch:.3f}, f_area = {f_area:.3f}")
    print(f"  Predicted turbulent flame speed S_T = {S_T * 100:.2f} cm/s")
    return S_T


# -----------------------------
# Part 2: Experimental Data and Interpolation
# -----------------------------
def get_experimental_data():
    """
    Returns experimental data as two numpy arrays:
       phi_exp: Equivalence ratios
       S_T_exp: Experimental turbulent flame speeds [cm/s]
    """
    phi_exp = np.array([0.499538714, 0.602982821, 0.70170306, 0.905732286, 1.005399765])
    S_T_exp = np.array([178.2203716, 200.8061413, 214.2229173, 244.7067994, 260.8760003])
    return phi_exp, S_T_exp


def interp_experimental(phi_vals, phi_exp, S_T_exp):
    """
    Interpolate the experimental turbulent flame speed data to the phi values provided.
    Returns an array of experimental flame speeds (cm/s) corresponding to phi_vals.
    """
    S_T_interp = np.interp(phi_vals, phi_exp, S_T_exp)
    return S_T_interp


# -----------------------------
# Part 3: Scan Over Equivalence Ratios and Export Results
# -----------------------------
def scan_equivalence_ratios(phi_range, p, T, mech='h2o2', transport_model='multicomponent',
                            u_prime=2.0, L=0.01, L_M=0.002, A=0.966, n=1.164):
    """
    Scan over a range of equivalence ratios and compute the laminar flame speed and
    turbulent flame speed for each.

    The laminar flame speed is computed using simulate_flame_speed with signature:
         simulate_flame_speed(phi, p, T, mech, transport_model)
    which returns the flame speed in cm/s. We convert it to m/s for further computation.

    Returns:
      phi_list  : Array of equivalence ratios
      S_l0_list : Array of laminar flame speeds [cm/s]
      S_T_list  : Array of computed turbulent flame speeds [cm/s]
    """
    phi_list = []
    S_l0_list = []
    S_T_list = []

    for phi in phi_range:
        print("\n==============================================")
        print(f"Running for equivalence ratio ϕ = {phi:.3f}")
        # Call simulate_flame_speed with the required arguments.
        # Note: p is used instead of pressure and T instead of Tin.
        S_l0_cm = simulate_flame_speed(phi, p, T, mech=mech, transport_model=transport_model)
        if S_l0_cm is None:
            print(f"Simulation failed for phi = {phi:.3f}")
            continue
        # Convert laminar flame speed to m/s:
        S_l0 = S_l0_cm / 100.0
        # Compute turbulent flame speed with flame stretch correction:
        S_T = compute_turbulent_flame_speed(S_l0, u_prime, L, L_M=L_M, A=A, n=n)

        phi_list.append(phi)
        S_l0_list.append(S_l0_cm)  # in cm/s
        S_T_list.append(S_T * 100)  # in cm/s

    return np.array(phi_list), np.array(S_l0_list), np.array(S_T_list)


def export_results_csv(filename, phi, S_l0, S_T, S_T_exp):
    """
    Export the results to a CSV file with columns:
       equivalence_ratio, S_l0_cm_per_s, S_T_computed_cm_per_s, S_T_experimental_cm_per_s
    """
    data = np.column_stack((phi, S_l0, S_T, S_T_exp))
    header = "equivalence_ratio,S_l0_cm_per_s,S_T_computed_cm_per_s,S_T_experimental_cm_per_s"
    np.savetxt(filename, data, delimiter=",", header=header, comments='', fmt='%.8f')
    print(f"Results exported to {filename}")


# -----------------------------
# Part 4: Main Execution
# -----------------------------
def main():
    # Define the simulation conditions:
    p = ct.one_atm  # Pressure in Pa
    T = 300.0  # Temperature in K

    # Define the range of equivalence ratios to scan
    phi_range = np.linspace(0.5, 1.0, 11)

    # Fixed turbulence/stretch parameters:
    u_prime = 2.0  # Turbulence intensity [m/s]
    L = 0.01  # Integral length scale [m]
    L_M = 0.002  # Markstein length [m]
    A = 0.966  # Empirical area enhancement constant
    n = 1.164  # Empirical exponent

    # Run the scan
    phi_vals, S_l0_vals, S_T_vals = scan_equivalence_ratios(phi_range, p, T, mech='h2o2',
                                                            transport_model='multicomponent',
                                                            u_prime=u_prime, L=L, L_M=L_M, A=A, n=n)

    # Get experimental data and interpolate to the scanned phi values
    phi_exp, S_T_exp_data = get_experimental_data()
    S_T_exp_interp = np.interp(phi_vals, phi_exp, S_T_exp_data)

    # Export results to CSV
    export_results_csv("flame_speed_results.csv", phi_vals, S_l0_vals, S_T_vals, S_T_exp_interp)

    # Produce a high-resolution 4K figure (16x9 inches at 240 dpi ≈ 3840x2160 pixels)
    plt.figure(figsize=(16, 9), dpi=240)
    plt.plot(phi_vals, S_l0_vals, 'bo-', markersize=8, label='Laminar Flame Speed (Sₗ₀)')
    plt.plot(phi_vals, S_T_vals, 'rs-', markersize=8, label='Computed Turbulent Flame Speed (Sₜ)')
    plt.plot(phi_vals, S_T_exp_interp, 'k*', markersize=12, label='Experimental Sₜ (Interp.)')
    plt.xlabel("Equivalence Ratio (ϕ)", fontsize=16)
    plt.ylabel("Flame Speed [cm/s]", fontsize=16)
    plt.title("Flame speeds for hydrogen/air mixtures", fontsize=18)
    plt.xlim(0.5, 1.0)
    # plt.ylim(0, 400)
    plt.legend(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("flame_speed_comparison_4k.png", dpi=240)
    print("High-resolution figure saved as 'flame_speed_comparison_4k.png'.")
    plt.show()


if __name__ == '__main__':
    main()
