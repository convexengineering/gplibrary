import numpy as np
import matplotlib.pyplot as plt
from gpkit.small_scripts import unitstr

def poor_mans_contour(model, xvarname, xsweep, zvarname, zsweep, yvarname,
                      ylimits):

    # save sold substitions substitute
    x_old = model.substitutions[xvarname]
    z_old = model.substitutions[zvarname]

    # substitute in new variables
    model.substitutions.update({xvarname: ("sweep", xsweep)})

    # create figure
    fig, ax = plt.subplots()

    # iterate through z variables
    for val in zsweep:
        model.substitutions.update({zvarname: val})
        sol = model.solve("mosek", skipsweepfailures=True)

        ax.plot(sol(xvarname), sol(yvarname))

    # return substitution values to what they were
    model.substitutions.update({xvarname: x_old})
    model.substitutions.update({zvarname: z_old})

    # format plot
    ax.set_xlabel("%s [%s]" % (model[xvarname].descr["label"],
                               unitstr(model[xvarname].units)))
    ax.set_ylabel("%s [%s]" % (model[yvarname].descr["label"],
                               unitstr(model[yvarname].units)))
    ax.set_ylim((ylimits[0], ylimits[1]))

    return fig, ax
