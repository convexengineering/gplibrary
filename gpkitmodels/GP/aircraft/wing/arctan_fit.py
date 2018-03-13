from gpfit.fit import fit
import gpkitmodels.GP.aircraft.wing.gustloading as GustL
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams.update({'font.size':15})
GENERATE = True

def arctanfit():
    u = np.linspace(1e-15, 0.7, 100)
    w = np.arctan(u)

    x = np.log(u)
    y = np.log(w)

    cn, err = fit(x, y, 1, "MA")
    rm = err
    print "RMS error: %.4f" % rm

    yfit = cn.evaluate(x)
    df = cn.get_dataframe()
    fig, ax = plt.subplots()
    ax.plot(u, w, lw=2)
    ax.plot(u, np.exp(yfit), "--", lw=2)
    ax.set_xlim([0, 0.7])
    ax.grid()
    ax.set_xlabel("$V_{\\mathrm{gust}}/V$")
    ax.set_ylabel("$\\alpha_{\\mathrm{gust}}$")
    ax.legend(["$\\arctan{(V_{\\mathrm{gust}}/V)}$",
               "$0.905 (V_{\\mathrm{gust}}/V)^{0.961}$"], loc=2, fontsize=15)
    return df, fig, ax

if __name__ == "__main__":
    df, fig, ax = arctanfit()
    if GENERATE:
        path = os.path.dirname(GustL.__file__)
        df.to_csv(path + os.sep + "arctan_fit.csv", index=False)
    else:
        df.to_csv("arctan_fit.csv")
    if len(sys.argv) > 1:
        path = sys.argv[1]
        fig.savefig(path + "arctanfit.pdf", bbox_inches="tight")
    else:
        fig.savefig("arctanfit.pdf", bbox_inches="tight")
