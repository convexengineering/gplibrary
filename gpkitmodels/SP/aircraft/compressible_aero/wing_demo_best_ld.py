"""Demonstrate the wing model - bed lift/drag ratio."""

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


def ld_plots():
    fig_dir = 'docs/figures/wing_best_ld'
    aspect = 10

    flight_state_model = flight_state.FlightState()
    airfoil_model = airfoil.DummyAirfoil(flight_state_model)
    wing_geom = wing.WingSwept()
    wing_model = wing.WingAeroSubsonic(flight_state_model, wing_geom, airfoil_model)

    wing_model.substitutions['S_ref'] = 1    # [units: meter^2]
    wing_model.substitutions['aspect'] = aspect    # [units: dimensionless]
    wing_model.substitutions['mach'] = 0.5    # [units: dimensionless]
    wing_model.substitutions['p_static'] = 30e3    # [units: pascal]
    wing_model.substitutions['sweep_c2'] = 1e-6    # [units: dimensionless]
    # wing_model.substitutions['lift'] = 1e3    # [units: newton]

    # maximize lift/drag
    # Also, add a light cost to sweep - there are structural and control
    # penalties on sweep which this model does not capture.
    # Without this cost on sweep, the model would optimize to high sweep
    # angles for low-Mach wings, because the drag model predicts a very slight
    # decrease in drag from sweep even at low Mach.
    wing_model.cost = wing_model['drag']/wing_model['lift'] + 0.005 * wing_model['sweep_c2']

    sol = wing_model.solve()
    print sol.table()

    # Plot sweep for best L/D vs mach
    del wing_model.substitutions['sweep_c2']
    mach = np.linspace(0.01, 0.99, 11)
    ld = np.zeros(len(mach))
    c_l = np.zeros(len(mach))
    sweep = np.zeros(len(mach))
    gp_plot_style = {'color':'C1', 'linestyle':'', 'marker':'+', 'markersize':12}

    for i in range(len(mach)):
        wing_model.substitutions['mach'] = mach[i]
        sol = wing_model.localsolve()
        ld[i] = sol['variables']['lift'] / sol['variables']['drag']
        c_l[i] = sol['variables']['c_l']
        sweep[i] = sol['variables']['sweep_c2']

    ax1 = plt.subplot(3, 1, 1)
    plt.title('Optimizing for best $L/D$' + 
              '\n$AR={{{:.0f}}}$'.format(aspect))
    plt.plot(mach, ld, label='GP model', **gp_plot_style)
    plt.ylabel('$L/D$')

    plt.subplot(3, 1, 2, sharex=ax1)
    plt.plot(mach, np.rad2deg(sweep), label='GP model', **gp_plot_style)
    plt.ylabel('Sweep $\\Lambda$ [deg]')

    plt.subplot(3, 1, 3, sharex=ax1)
    plt.plot(mach, c_l, label='GP model', **gp_plot_style)
    plt.ylabel('2D lift coeff. $c_l$ [-]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.savefig(os.path.join(fig_dir, 'best_ld.png'))


def main():
    ld_plots()
    plt.show()


if __name__ == '__main__':
    main()
