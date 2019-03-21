"""Demonstrate the airfoil model."""

import os.path
import numpy as np
from matplotlib import pyplot as plt
import gpkit
from gpkit import units, Variable, Model
from gpkit.tools.autosweep import autosweep_1d

import airfoil
import flight_state

def main():
    airfoil_model = airfoil.DummyAirfoil(flight_state.FlightState())

    airfoil_model.substitutions['c_l'] = 0.5
    airfoil_model.substitutions['mach_perp'] = 0.5
    airfoil_model.substitutions['mach'] = 0.5

    airfoil_model.cost = airfoil_model['c_dp']

    sol = airfoil_model.solve()
    print sol.table()

    # Replicate Figure 8.34 from Drela's FVA book.
    plt.figure()
    for c_l in [0.5, 0.6, 0.7, 0.8]:
        airfoil_model.substitutions['c_l'] = c_l
        sweep = airfoil_model.autosweep({airfoil_model['mach_perp']: (0.50, 0.80)}, tol=1e-3)
        plt.plot(sweep['mach_perp'], sweep['c_dp'] + sweep['c_df'], label='$c_l = ${:.2f}'.format(c_l))
    plt.xlabel('Mach number $M$ [-]')
    plt.ylabel('2D drag coefficient $c_d$ [-]')
    plt.legend()
    plt.xlim([0.5, 0.8])
    plt.ylim([0, 0.020])
    plt.grid()

    plt.figure()
    for mach in [0.6, 0.65, 0.7, 0.73, 0.75, 0.76]:
        airfoil_model.substitutions['mach_perp'] = mach
        sweep = airfoil_model.autosweep({airfoil_model['c_l']: (0.30, 1.1)}, tol=1e-3)
        plt.plot(sweep['c_dp'] + sweep['c_df'], sweep['c_l'],
                 label='$M_\\perp = ${:.2f}'.format(mach))
    plt.xlabel('2D drag coefficient $c_d$ [-]')
    plt.ylabel('2D lift coefficient $c_l$ [-]')
    plt.legend()
    plt.xlim([0, 0.020])
    plt.ylim([0, 1.5])
    plt.grid()

    plt.show()


if __name__ == '__main__':
    main()
