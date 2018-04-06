"jho1_polarfits.py"
import numpy as np
from numpy import exp
import pandas as pd
import matplotlib.pyplot as plt
from gpfit.fit import fit
import sys
import inspect
import os

GENERATE = True
plt.rcParams.update({'font.size':15})

def text_to_df(filename):
    "parse XFOIL polars and concatente data in DataFrame"
    lines = list(open(filename))
    for i, l in enumerate(lines):
        lines[i] = l.split("\n")[0]
        for j in 10-np.arange(9):
            if " "*j in lines[i]:
                lines[i] = lines[i].replace(" "*j, " ")
            if "---" in lines[i]:
                start = i
    data = {}
    titles = lines[start-1].split(" ")[1:]
    for t in titles:
        data[t] = []

    for l in lines[start+1:]:
        for i, v in enumerate(l.split(" ")[1:]):
            data[titles[i]].append(v)

    df = pd.DataFrame(data)
    return df

def fit_setup(re_range, tau_range):
    "set up x and y parameters for gp fitting"
    CL = []
    CD = []
    RE = []
    tau = []
    for i in range(len(tau_range["scale"])):
        for r in re_range:
            dataf = text_to_df("dai1336a.%d.ncrit09.Re%dk.pol" % (
                tau_range["scale"][i], r))
            CL.append(dataf["CL"].values.astype(np.float))
            CD.append(dataf["CD"].values.astype(np.float))
            RE.append([r*1000.0]*len(CL[-1]))
            tau.append([tau_range["tau"][i]]*len(CL[-1]))

    u1 = np.hstack(CL)
    u2 = np.hstack(RE)
    u3 = np.hstack(tau)
    w = np.hstack(CD)
    u = [u1, u2, u3]
    x = np.log(u)
    y = np.log(w)
    return x, y

def plot_fits(cnstr, x, y):
    "plot fit compared to data"

    yfit = cnstr.evaluate(x)
    x2, ind = np.unique(x[2], return_index=True)
    ind = np.append(ind, len(x[2]))
    figs = []
    for i in range(1, len(ind)):
        colors = ["#084081", "#0868ac", "#2b8cbe", "#4eb3d3", "#7bccc4"]*5
        xt, yt = x[0:2, ind[i-1]:ind[i]], y[ind[i-1]:ind[i]]
        yft = yfit[ind[i-1]:ind[i]]
        x1, ind1 = np.unique(xt[1], return_index=True)
        ind1 = np.append(ind1, len(xt[1]))
        x0 = [xt[0][ind1[j-1]:ind1[j]] for j in range(1, len(ind1))]
        y0 = [yt[ind1[j-1]:ind1[j]] for j in range(1, len(ind1))]
        yf0 = [yft[ind1[j-1]:ind1[j]] for j in range(1, len(ind1))]
        fig, ax = plt.subplots()
        c = 0
        for r, cl, cd, fi in zip(exp(x1), x0, y0, yf0):
            ax.plot(exp(cl), exp(cd), "o", mec=colors[c], mfc="none", mew=1.5)
            ax.plot(exp(cl), exp(fi), c=colors[c], label="Re = %dk" % r, lw=2)
            c += 1
        ax.set_xlabel("$C_L$")
        ax.set_ylabel("$c_{d_p}$")
        ax.set_title("$\\tau = %.2f$" % exp(x2[i-1]))
        ax.legend(loc=2)
        ax.grid()
        figs.append(fig)
    return figs

if __name__ == "__main__":
    Re_r = [125, 150, 200]
    tau_r = {"scale": [80, 90, 100, 105],
             "tau": [0.1094, 0.1231, 0.13682, 0.146]}
    X, Y = fit_setup(Re_r, tau_r) # call fit(X, Y, 4, "SMA") to get fit
    np.random.seed(0)
    cn, err = fit(X, Y, 3, "SMA")
    print "RMS error: %.5f" % err
    df = cn.get_dataframe()
    df.to_csv("../../dai1336a.csv", index=False)

    Fs = plot_fits(cn, X, Y)
    for t, F in zip(tau_r["scale"], Fs):
        F.savefig("dai1336a.%d.fits.pdf" % t, bbox_inches="tight")
