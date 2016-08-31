from gas_male_fixCDR import GasMALEFixedEngine
from gas_male_discritized_climb import GasMALEDiscritizedClimb
from gpkit.small_scripts import unitstr
import numpy as np
import matplotlib.pyplot as plt

def fix_vars(model, solution, var_names):
    """
    Fixes variables to values found in initial solution

    Arguments
    ---------
    model: Model
    solution: soluation - initial solution with values of vars to be fixed
    var_names: Dict - variable names of vars to be fixed and value of added
                      tolerance
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

def plot_altitude_sweeps(hvals, yvarnames, vars_to_fix, CLIMB=False):
    """
    Plots sweeps of yvarnames vs altitude. Only runs GasMALEFixedEngine().
    Arguments
    ---------
    hvals     : array, desired altitude sweep
    yvarnames : dict {variable name}, desired y variable for plotting
    vars_to_fix: dict - name of variable to fix and tolerance fixing
    CLIMB: boolean - True if using gasmale.py,
                     False if using gas_male_fixCDR.py

    Output
    ------
    Saves plot in pdf format.  Number of plots equal to number of values in
    yvarnames.  Plot names = "altitude_vs_%s.pdf"
    """

    vals = np.zeros([len(hvals), len(yvarnames)])
    if CLIMB:
        M_fix = GasMALEDiscritizedClimb()
    else:
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
        if CLIMB:
            plot_name = "altitude_vs_%s_Climb.pdf" % \
                       M[yvarname].descr["label"].replace(" ", "")
        else:
            plot_name = "altitude_vs_%s_CDR.pdf" % \
                       M[yvarname].descr["label"].replace(" ", "")

        fig.savefig(plot_name)

def plot_mission_var(model, yvarname, ylim, yaxis_name=None):
    """
    Plots a mission varible against mission time.

    Arguments
    ---------
    model: must be GasPoweredMALE from gasmale.py
    yvarname: String - variable string name

    Returns:
    fig, ax: matplotlib figure and axis - time not to scale, labeled
             by mission profile inhereint in model.
    """

    # solve model
    sol = model.solve("mosek")

    y = sol(yvarname).magnitude

    # define tick step on plot
    t = np.linspace(0, model.NSeg - 1, model.NSeg)

    # create plot
    fig, ax = plt.subplots()
    line, = ax.plot(t, y)

    # label time axis
    ax.xaxis.set_ticks(np.arange(0, model.NSeg-1, 1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[model.mStart + np.round(model.NClimb1/2)] = 'Climb'
    labels[model.mEndClimb + np.round(model.NCruise1/2)] = 'Cruise'
    labels[model.mEndCruise + np.round(model.NClimb2/2)] = 'Climb'
    labels[model.mEndClimb2 + np.round(model.NLoiter/2)] = 'Loiter'
    labels[model.mEndLoiter + np.round(model.NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)

    # mark mission profile changes
    ax.set_ylim([ylim[0], ylim[1]])
    if yaxis_name:
        ax.set_ylabel(yaxis_name)
    else:
        ax.set_ylabel("%s [%s]" % (model[yvarname][0].descr["label"],
                                   unitstr(model[yvarname][0].units)))
    ax.grid()
    ax.plot([model.mEndClimb - 1, model.mEndClimb - 1],
            [ylim[0], ylim[1]], '--', color='r')
    ax.plot([model.mEndCruise - 1, model.mEndCruise - 1],
            [ylim[0], ylim[1]], '--', color='r')
    ax.plot([model.mEndClimb2 - 1, model.mEndClimb2 - 1],
            [ylim[0], ylim[1]], '--', color='r')
    ax.plot([model.mEndLoiter - 1, model.mEndLoiter - 1],
            [ylim[0], ylim[1]], '--', color='r')

    return fig, ax
