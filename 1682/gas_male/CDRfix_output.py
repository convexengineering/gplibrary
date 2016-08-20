from gas_male_fixCDR import GasMALEFixedEngine
from gpkit.small_scripts import unitstr
import numpy as np
import matplotlib.pyplot as plt
import pdb

def fix_vars(model, solution, var_names):
    """
    Fixes variables to values found in initial solution

    Arguments
    ---------
    model: Model
    solution: soluation - initial solution with values of vars to be fixed
    var_names: String - variable names of vars to be fixed
    """
    for name in var_names:
        value = solution(name).magnitude
        model.substitutions.update({name: value + var_names[name]})

def plot_sweep(model, xvarname, xsweep, yvarname, ylim=[0, 10]):
    """
    Takes model with sweep input and returns figure with desired output var

    Arguments
    ---------
    model: Model
    xvarname: String - variable name of sweep var
    xsweep: np.array - sweep values
    yvarname: String - variable name of desired output
    ylim: 2D array - plotting limits on y axis; x axis defaults are sweep min
                     and max

    Returns
    -------
    fig, ax: figure and axis of subplot
    """

    oldsub = model.substitutions[xvarname]

    model.substitutions.update({xvarname: ("sweep", xsweep)})
    sol = model.solve("mosek", skipsweepfailures=True)

    fig, ax = plt.subplots()
    ax.plot(sol(xvarname), sol(yvarname))
    ax.set_xlabel("%s [%s]" % (model[xvarname].descr["label"],
                               unitstr(model[xvarname].units)))
    ax.set_ylabel("%s [%s]" % (model[yvarname].descr["label"],
                               unitstr(model[yvarname].units)))
    ax.set_title("CDR " + yvarname + " vs " + xvarname)
    ax.set_ylim((ylim[0], ylim[1]))
    plt.grid()

    model.substitutions.update({xvarname: oldsub})

    return fig, ax

def plot_altitude_sweeps(hvals, yvarnames, vars_to_fix):
    """
    Plots sweeps of yvarnames vs altitude. Only runs GasMALEFixedEngine().
    Arguments
    ---------
    hvals     : array, desired altitude sweep
    yvarnames : dict {variable name}, desired y variable for plotting

    Output
    ------
    Saves plot in pdf format.  Number of plots equal to number of values in
    yvarnames.  Plot names = "altitude_vs_%s.pdf"
    """

    vals = np.zeros([len(hvals), len(yvarnames)])
    M_fix = GasMALEFixedEngine()
    sol_fix = M_fix.solve("mosek", verbosity=0)

    for i, h in enumerate(hvals):
        M = GasMALEFixedEngine(h)
        fix_vars(M, sol_fix, vars_to_fix)
        del M.substitutions["t_{station}"]
        M.cost = 1/M["t_{station}"]
        sol = M.solve("mosek", verbosity=0)
        for j, yvarname in enumerate(yvarnames):
            vals[i,j] = sol(yvarname).magnitude

    for j, yvarname in enumerate(yvarnames):
        fig, ax = plt.subplots()
        ax.plot(hvals, vals[:, j])
        ax.set_xlabel("%s [%s]" % (M_fix["h_{station}"].descr["label"],
                                   unitstr(M_fix["h_{station}"].units)))
        ax.set_ylabel("%s [%s]" % (M_fix[yvarname].descr["label"],
                                   unitstr(M_fix[yvarname].units)))
        ax.set_title("CRD " + yvarname + " vs h_{station}")
        plt.grid()
        fig.savefig("altitude_vs_%s.pdf" %
                    M[yvarname].descr["label"].replace(" ", ""))

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
    M = GasMALEFixedEngine(WIND=True)
    fix_vars(M, sol, vars_to_fix)
    del M.substitutions["t_{station}"]
    M.cost = 1/M["t_{station}"]
    fig, ax = plot_sweep(M, "V_{wind}", np.linspace(5, 40, 15), "t_{station}")
    fig.savefig("tvsV_wind.pdf")

    # altitude vs time on station
    plot_altitude_sweeps(np.linspace(14000, 23000, 20), {"t_{station}"},
                         vars_to_fix)
