from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import pandas as pd
from scipy.interpolate import interp1d


def simulate_flame_speed(phi, p, T, mech='h2o2', transport_model='multicomponent'):
    """
    Simulate laminar flame speed for given conditions

    Parameters:
    -----------
    phi : float
        Equivalence ratio
    p : float
        Pressure in Pa
    T : float
        Temperature in K
    mech : str
        'h2o2' or 'gri30'
    transport_model : str
        'mixture-averaged' or 'multicomponent'

    Returns:
    --------
    float
        Flame speed in cm/s
    """
    # Calculate H2:O2 ratio based on equivalence ratio
    h2_ratio = phi * 2.0  # stoichiometric H2:O2 is 2:1
    reactants = f'H2:{h2_ratio}, O2:1, AR:5'

    width = 0.03  # m
    loglevel = 0  # amount of diagnostic output (0 to 8)
    # Select mechanism
    if mech == 'h2o2':
        gas = ct.Solution('h2o2.yaml')
    else:  # gri30
        gas = ct.Solution('gri30.yaml')

    gas.TPX = T, p, reactants

    # Set up flame object
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=2, slope=0.01, curve=0.01)
    f.transport_model = transport_model

    try:
        f.solve(loglevel=loglevel, auto=True)
        return f.velocity[0] * 100  # Convert to cm/s
    except Exception as e:
        print(f"Error at phi={phi} with {mech}: {e}")
        return None


def run_comparison_study():
    """Run simulations with both mechanisms and compare to experimental data"""

    # Load experimental data
    exp_data = pd.read_csv('LFS_reference.csv')

    # Simulation conditions
    p = ct.one_atm  # 1 atm pressure
    T = 300.0  # Temperature in K

    # Create simulation points
    phi_exp = exp_data['phi'].values
    phi_min, phi_max = phi_exp.min(), phi_exp.max()
    phi_sim = np.linspace(phi_min, phi_max, 40)

    # Dictionary to store results
    results = {
        'h2o2': [],
        'gri30': []
    }

    # Run simulations for both mechanisms
    for mech in results.keys():
        print(f"\nRunning simulations with {mech} mechanism...")
        speeds = []
        for phi in phi_sim:
            speed = simulate_flame_speed(phi, p, T, mech=mech)
            speeds.append(speed)
            print(f"φ = {phi:.2f}, Speed = {speed:.1f} cm/s")
        results[mech] = speeds

    # Create comparison plot
    plt.figure(figsize=(10, 8))

    # Plot simulation results
    plt.plot(phi_sim, results['h2o2'], 'k-', label='H2O2 Mechanism', linewidth=2)
    plt.plot(phi_sim, results['gri30'], 'r--', label='GRI 3.0 Mechanism', linewidth=2)

    # Plot experimental data
    plt.plot(exp_data['phi'], exp_data['LFS'], 'ko',
             markerfacecolor='none', label='Experimental', markersize=8)

    # Formatting
    plt.xlabel('Equivalence Ratio', fontsize=12)
    plt.ylabel('Flame speed [cm/s]', fontsize=12)
    plt.title('Laminar flame speed comparison for H₂\n(Pᵢ=1 atm, Tᵢ=300K)',
              fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)
    plt.xlim(0.0, 3.0)
    plt.ylim(0, 350)

    # Save the plot
    plt.savefig('mechanism_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Calculate statistics and errors
    print("\nMechanism Comparison Statistics:")

    for mech, speeds in results.items():
        # Create interpolation function for mechanism results
        sim_interp = interp1d(phi_sim, speeds)
        sim_at_exp_points = sim_interp(exp_data['phi'])

        # Calculate errors
        abs_error = np.abs(sim_at_exp_points - exp_data['LFS'])
        rel_error = abs_error / exp_data['LFS'] * 100

        print(f"\n{mech.upper()} Mechanism:")
        print(f"Maximum flame speed: {np.nanmax(speeds):.1f} cm/s")
        print(f"at equivalence ratio: {phi_sim[np.nanargmax(speeds)]:.2f}")
        print(f"Average absolute error: {np.mean(abs_error):.1f} cm/s")
        print(f"Average relative error: {np.mean(rel_error):.1f}%")
        print(f"Maximum relative error: {np.max(rel_error):.1f}%")

        # Save mechanism results
        mech_df = pd.DataFrame({
            'phi': phi_sim,
            'LFS': speeds
        })
        mech_df.to_csv(f'{mech}_results.csv', index=False)

    # Save experimental data for reference
    exp_data.to_csv('experimental_data.csv', index=False)


if __name__ == "__main__":
    # Run the complete comparison study
    run_comparison_study()
    print("\nResults and plots have been saved to files")
