"""Demonstrate the wing model - drag."""

import os.path
from subprocess import CalledProcessError
import numpy as np
from matplotlib import pyplot as plt
import gpkit
from gpkit import units, Variable, Model
from gpkit.tools.autosweep import autosweep_1d
from gpkit.small_scripts import mag

import wing
from wing_model_compare import datcom_model
import airfoil
import flight_state



def drag_plots():
    aspect = 10
    C_L = 0.3

    flight_state_model = flight_state.FlightState()
    airfoil_model = airfoil.DummyAirfoil(flight_state_model)
    wing_model = wing.WingAeroSubsonic(flight_state_model, airfoil_model, limit_sweep_to_valid_range=False)

    wing_model.substitutions['S_ref'] = 1    # [units: meter^2]
    wing_model.substitutions['aspect'] = aspect    # [units: dimensionless]
    wing_model.substitutions['p_static'] = 30e3    # [units: pascal]
    wing_model.substitutions['C_L'] = C_L    # [units: dimensionless]

    # minimize drag
    wing_model.cost = wing_model['drag']

    # Plot drag vs Mach at various sweeps
    ax1 = plt.subplot(3, 1, 1)
    plt.subplot(3, 1, 2, sharex=ax1)
    plt.subplot(3, 1, 3, sharex=ax1)

    mach = np.linspace(0.5, 0.99, 20)
    C_D = np.zeros(len(mach))
    c_l = np.zeros(len(mach))
    alpha = np.zeros(len(mach))
    sweeps = np.deg2rad([1e-6, 30, 60])
    # sweeps = np.deg2rad([1e-6, 20, 40])
    for sweep in sweeps:
        for i in range(len(mach)):
            wing_model.substitutions['mach'] = mach[i]
            wing_model.substitutions['sweep_c2'] = sweep
            try:
                sol = wing_model.solve()
                C_D[i] = sol['variables']['C_D']
                c_l[i] = sol['variables']['c_l']
                alpha[i] = sol['variables']['alpha']
            # except CalledProcessError, RuntimeWarning:
            except:
                C_D[i] = np.nan
                c_l[i] = np.nan
                alpha[i] = np.nan
        plt.subplot(3, 1, 1)
        plt.plot(mach, C_D, label='$\\Lambda = ${:.0f} deg'.format(
            np.rad2deg(sweep)))
        plt.subplot(3, 1, 2)
        plt.plot(mach, c_l, label='$\\Lambda = ${:.0f} deg'.format(
            np.rad2deg(sweep)))
        plt.subplot(3, 1, 3)
        plt.plot(mach, np.rad2deg(alpha), label='$\\Lambda = ${:.0f} deg'.format(
            np.rad2deg(sweep)))
    plt.subplot(3, 1, 1)
    plt.ylim([0, 0.025])
    plt.ylabel('$C_D$ [-]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')
    plt.legend()
    plt.title('Drag vs Mach at constant $C_L =${:.2f}'.format(C_L))

    plt.subplot(3, 1, 2)
    plt.ylabel('$c_l$ [-]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.ylabel('Angle of attack $\\alpha$ [deg]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')
    plt.legend()

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    """
    Why does C_D (at const C_L) decrease with sweep at low mach?
    I think because of cos^3 term on pressure drag in Drela 8.183?
    Also, although Drela says that drag should increase monotonically with sweep
    below the critical Mach number, the plots in Hoerner indicate that
    C_D may actually decrease slightly with sweep even at lower (~0.6)
    Mach numbers.
    """


def main():
    drag_plots()
    plt.show()


if __name__ == '__main__':
    main()
