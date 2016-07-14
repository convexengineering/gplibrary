import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from gpkit.small_scripts import unitstr
import pdb

matplotlib.rc('font', size=17)

def poor_mans_contour(model, xvarname, xsweep, zvarname, zsweep, yvarname,
                      ylimits, vref=None, vrefname=None, href=None,
                      hrefname=None):

    # save sold substitions substitute
    x_old = model.substitutions[xvarname]
    z_old = model.substitutions[zvarname]

    # substitute in new variables
    model.substitutions.update({xvarname: ("sweep", xsweep)})

    # create figure
    fig, ax = plt.subplots()
    lines = []

    # iterate through z variables
    for val in zsweep:
        model.substitutions.update({zvarname: val})
        sol = model.solve("mosek", skipsweepfailures=True)

        l, = ax.plot(sol(xvarname), sol(yvarname),
                     label="%s: %s" % (zvarname, val))
        lines.append(l)

    # return substitution values to what they were
    model.substitutions.update({xvarname: x_old})
    model.substitutions.update({zvarname: z_old})

    # format plot
    ax.set_xlabel("%s [%s]" % (model[xvarname].descr["label"],
                               unitstr(model[xvarname].units)))
    ax.set_ylabel("%s [%s]" % (model[yvarname].descr["label"],
                               unitstr(model[yvarname].units)))
    ax.set_ylim((ylimits[0], ylimits[1]))
    if vref:
        l, = ax.plot([vref, vref], ylimits, "r--", label=vrefname)
        lines.append(l)
    if href:
        l, = ax.plot([xsweep[0], xsweep[-1]], [href, href], "r--",
                     label=hrefname)
        lines.append(l)

    plt.legend(handles=lines)
    plt.grid()

    return fig, ax

def latitude_sweep(model, lat, xvarnames, xsweeps, zvarname, zsweep, yvarname, ylim, winddf=None):

    # save original values
    x_old = []
    for name in xvarnames:
        x_old.append(model.substitutions[name])

    z_old = model.substitutions[zvarname]

    fig, ax = plt.subplots()
    lines = []

    for i, zval in enumerate(zsweep):
        model.substitutions.update({zvarname: zval})
        y = []
        for j in range(0,xsweeps.size/len(xvarnames)):
            for k, name in enumerate(xvarnames):
                model.substitutions.update({name: xsweeps[k,j]})
            if winddf is not None:
                model.substitutions.update({"V_{wind}": np.array(winddf["%sth Percentile Winds" % zval])[j]})
                print "subbed %s for V_{wind}" % model.substitutions["V_{wind}"]
            try:
                sol = model.solve("mosek", verbosity=0)
                y.append(sol(yvarname).magnitude)
            except RuntimeWarning:
                y.append(np.nan)

        l, = ax.plot(lat, y, label="%s: %s" % (zvarname, zval), linewidth=2)
        lines.append(l)

    # format plot
    ax.set_xlabel("Latitude [deg]")
    ax.set_ylabel("%s [%s]" % (model[yvarname].descr["label"],
                               unitstr(model[yvarname].units)))
    ax.set_ylim((ylim[0], ylim[1]))

    plt.legend(handles=lines, loc="upper left")
    plt.grid()

    return fig, ax
