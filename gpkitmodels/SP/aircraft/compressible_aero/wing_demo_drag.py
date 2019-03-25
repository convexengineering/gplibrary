"""Demonstrate the wing model - drag."""

import os.path
from subprocess import CalledProcessError
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import gpkit
from gpkit import units, Variable, Model
from gpkit.tools.autosweep import autosweep_1d
from gpkit.small_scripts import mag

import wing
from wing_model_compare import datcom_model
import airfoil
import flight_state



def drag_plots():
    fig_dir = 'docs/figures/wing_drag_model'
    aspect = 4
    C_L = 0.3

    flight_state_model = flight_state.FlightState()
    airfoil_model = airfoil.DummyAirfoil(flight_state_model)
    wing_geom = wing.WingSwept()
    wing_model_cos1 = wing.WingAeroSubsonic(
        flight_state_model, wing_geom, airfoil_model, pressure_drag_cos_exponent=1)
    wing_model_cos3 = wing.WingAeroSubsonic(
        flight_state_model, wing_geom, airfoil_model, pressure_drag_cos_exponent=3)

    fig1, axes = plt.subplots(ncols=2, nrows=3, sharex='all', sharey='row',
                              figsize=(8, 6))

    for wing_model, plot_col in zip([wing_model_cos1, wing_model_cos3], [0, 1]):
        wing_model.substitutions['S_ref'] = 1    # [units: meter^2]
        wing_model.substitutions['aspect'] = aspect    # [units: dimensionless]
        wing_model.substitutions['p_static'] = 30e3    # [units: pascal]
        wing_model.substitutions['C_L'] = C_L    # [units: dimensionless]
        # wing_model.substitutions['alpha'] = np.deg2rad(2.)    # [units: radian]

        # minimize drag
        wing_model.cost = wing_model['drag']

        # Plot drag vs Mach at various sweeps
        mach = np.linspace(0.5, 0.99, 20)
        C_D = np.zeros(len(mach))
        c_l = np.zeros(len(mach))
        alpha = np.zeros(len(mach))
        sweeps = np.deg2rad([1e-6, 20, 40, 60])
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

            axes[0, plot_col].plot(mach, C_D, label='$\\Lambda = ${:.0f} deg'.format(
                np.rad2deg(sweep)))
            axes[1, plot_col].plot(mach, c_l, label='$\\Lambda = ${:.0f} deg'.format(
                np.rad2deg(sweep)))
            axes[2, plot_col].plot(mach, np.rad2deg(alpha), label='$\\Lambda = ${:.0f} deg'.format(
                np.rad2deg(sweep)))

        axes[0, plot_col].set_ylim([0, 0.025])
        axes[0, plot_col].set_ylabel('$C_D$ [-]')
        axes[0, plot_col].set_title(
            'Drag vs Mach at constant $C_L =${:.2f}'.format(C_L)
            + '\n $\\cos^{{{:d}}} \\Lambda$, $AR={:.0f}$'.format(
                wing_model.pressure_drag_cos_exponent, aspect))

        axes[1, plot_col].set_ylabel('$c_l$ [-]')
        axes[1, plot_col].legend(loc='lower right')

        axes[2, plot_col].set_ylabel('Angle of attack $\\alpha$ [deg]')
        axes[2, plot_col].set_xlabel('Freestream Mach number $M_\\infty$ [-]')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.1)
    plt.savefig(os.path.join(fig_dir, 'drag_cos_exp.png'))
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
