""" CDRfix_output.py """
import numpy as np
from gas_male_fixCDR import GasMALEFixedEngine
from plotting import fix_vars, plot_sweep, plot_altitude_sweeps

if __name__ == "__main__":

    vars_to_fix = {"S":0.0, "b":0.0, "Vol_{fuse}":0.0001}
    deltas = []
    M = GasMALEFixedEngine()
    sol = M.solve("mosek")

    # Fix geometric paramters and sweep for performance curves
    fix_vars(M, sol, vars_to_fix)
    sol = M.solve("mosek") # check for solving errors

    # set objective to time on station after fixing variables
    del M.substitutions["t_{station}"]
    M.cost = 1/M["t_{station}"]

    # payload power vs time on station
    fig, ax = plot_sweep(M, "P_{pay}", np.linspace(10, 200, 15),
                         "t_{station}")
    fig.savefig("tvsP_pay.pdf")

    # payload weight vs time on station
    fig, ax = plot_sweep(M, "W_{pay}", np.linspace(5, 40, 15), "t_{station}")
    fig.savefig("tvsW_pay.pdf")

    # wind speed vs time on station
    M = GasMALEFixedEngine(wind=True)
    fix_vars(M, sol, vars_to_fix)
    del M.substitutions["t_{station}"]
    M.cost = 1/M["t_{station}"]
    fig, ax = plot_sweep(M, "V_{wind}", np.linspace(5, 40, 15), "t_{station}")
    fig.savefig("tvsV_wind.pdf")

    # altitude vs time on station
    plot_altitude_sweeps(np.linspace(14000, 23000, 20), {"t_{station}"},
                         vars_to_fix)
