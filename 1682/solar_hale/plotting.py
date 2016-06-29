import numpy as np
import matplotlib.pyplot as plt
from gpkit.small_scripts import unitstr

def poor_mans_contour(M, xvarname, xsweep, zvarname, zsweep, yvarname,
                      ylimits):

    # save sold substitions substitute
    x_old = M.substitutions[xvarname]
    z_old = M.substitutions[zvarname]

    # substitute in new variables
    M.substitutions.update({xvarname: ("sweep", xsweep)})

    # create figure
    fig, ax = plt.subplots()

    # iterate through z variables
    for val in zsweep:
        M.substitutions.update({zvarname: val})
        sol = M.solve("mosek", skipsweepfailures=True)

        ax.plot(sol(xvarname), sol(yvarname))

    # format plot
    ax.set_xlabel("%s [%s]" % (M[xvarname].descr["label"],
                               unitstr(M[xvarname].units)))
    ax.set_ylabel("%s [%s]" % (M[yvarname].descr["label"],
                               unitstr(M[yvarname].units)))
    ax.set_ylim((ylimits[0], ylimits[1]))

    return fig, ax
