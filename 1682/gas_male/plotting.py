from gpkit.small_scripts import unitstr
from gasmale import GasMALE
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

def plot_altitude_sweeps(hvals, yvarnames, vars_to_fix):
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
    M_fix = GasMALE(DF70=True)
    M_fix.substitutions.update({"t_{loiter}": 6})
    M_fix.cost = M_fix["MTOW"]
    sol_fix = M_fix.solve("mosek", verbosity=0)

    for i, h in enumerate(hvals):
        M = GasMALE(h_station=h, DF70=True)
        fix_vars(M, sol_fix, vars_to_fix)
        sol = M.solve("mosek", verbosity=0)
        for j, yvarname in enumerate(yvarnames):
            vals[i, j] = sol(yvarname).magnitude

    figures = []
    axis = []
    hvar = M_fix.variables_byname("h")[0]
    for j, yvarname in enumerate(yvarnames):
        fig, ax = plt.subplots()
        ax.plot(hvals, vals[:, j])
        ax.set_xlabel("%s [%s]" % (hvar.descr["label"], unitstr(hvar.units)))
        ax.set_ylabel("%s [%s]" % (M_fix[yvarname].descr["label"],
                                   unitstr(M_fix[yvarname].units)))
        ax.set_title("CRD " + yvarname + " vs h_{station}")
        plt.grid()

        figures.append(fig)
        axis.append(ax)

    return figures, axis


def plot_mission_var(model, sol, yvarname, ylim, yaxis_name=None):
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

    sv = sol(yvarname)
    for seg in ["Loiter", "Cruise", "Climb"]:
        if seg in sv.items()[0][0].descr["models"]:
            ind = sv.items()[0][0].descr["models"].index(seg)

    shape = 0
    modelnum = []
    for key, value in sv.items():
        shape += key.descr["shape"][0]
        modelnum.append(key.descr["modelnums"][ind])

    y = np.zeros(shape)
    for key, value in sv.items():
        if key.descr["models"][ind] == "Climb":
            if key.descr["modelnums"][ind] == max(modelnum):
                y[10:15] = value.magnitude[0:5]
            else:
                y[0:5] = value.magnitude[0:5]
        elif key.descr["models"][ind] == "Cruise":
            if key.descr["modelnums"][ind] == max(modelnum):
                y[35:40] = value.magnitude[0:5]
            else:
                y[5:10] = value.magnitude[0:5]
        else:
            y[15:35] = value.magnitude[0:20]

    # define tick step on plot
    t = np.linspace(0, shape - 1, shape)

    # create plot
    fig, ax = plt.subplots()
    line, = ax.plot(t, y)

    # label time axis
    ax.xaxis.set_ticks(np.arange(0, shape - 1, 1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[2] = 'Climb'
    labels[7] = 'Cruise'
    labels[12] = 'Climb'
    labels[25] = 'Loiter'
    labels[37] = 'Cruise'
    ax.set_xticklabels(labels)

    # mark mission profile changes
    ax.set_ylim([ylim[0], ylim[1]])
    if yaxis_name:
        ax.set_ylabel(yaxis_name)
    else:
        ax.set_ylabel("%s [%s]" %
                      (model.variables_byname(yvarname)[0].descr["label"],
                       unitstr(model.variables_byname(yvarname)[0].units)))
    ax.grid()
    ax.plot([4, 4], [ylim[0], ylim[1]], '--', color='r')
    ax.plot([9, 9], [ylim[0], ylim[1]], '--', color='r')
    ax.plot([14, 14], [ylim[0], ylim[1]], '--', color='r')
    ax.plot([34, 34], [ylim[0], ylim[1]], '--', color='r')

    return fig, ax
