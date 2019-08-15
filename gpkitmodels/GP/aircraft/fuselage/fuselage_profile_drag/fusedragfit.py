" fuselage drag fits "
from builtins import zip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def fit_setup(filename):
    "set up fitting variables"

    df = pd.read_csv(filename)

    u1 = np.array(df["tubelr"])
    u2 = np.array(df["noselr"])
    u3 = np.array(df["taillr"])

    w = np.array(df["cd_front"])

    u1 = u1.astype(np.float)
    u2 = u2.astype(np.float)
    u3 = u3.astype(np.float)
    w = w.astype(np.float)

    u = [u1, u2, u3]
    x = np.log(u)
    y = np.log(w)

    return x, y

def return_fit(u_1, u_2, u_3):
    "fit using SMA, K = 4, RMS = 0.0479"
    w = (
        0.00243049 * (u_1)**0.033607 * (u_2)**1.21682 * (u_3)**0.306251
        + 0.00255095 * (u_1)**-0.0316887 * (u_2)**-0.585489 * (u_3)**1.15394
        + 0.0436011 * (u_1)**0.0545722 * (u_2)**0.258228 * (u_3)**-1.42664
        + 0.00970479 * (u_1)**0.8661 * (u_2)**-0.209136 * (u_3)**-0.156166) \
        ** (1/0.996232)
    return w

def plot_fits(filename):
    "plot fit against data"

    df = pd.read_csv(filename)
    u1 = np.array(df["tubelr"])
    u2 = np.array(df["noselr"])
    u3 = np.array(df["taillr"])
    body = np.unique(u1)
    nose = np.unique(u2)
    tail = np.unique(u3)

    figs = []

    colors = ["k", "m", "b", "g", "y"]
    assert len(colors) == len(body)

    for n in nose:
        fig, ax = plt.subplots()
        datan = df[(df["noselr"] == n)]
        for b, clr in zip(body, colors):
            datab = datan[(datan["tubelr"] == b)]
            ax.plot(datab["taillr"], datab["cd_front"], "o", c=clr)
            trl = np.array(datab["taillr"])
            taillr = np.linspace(trl[0], trl[-1], 20)
            cd = return_fit(b, n, taillr)
            ax.plot(taillr, cd, c=clr, label="body finess ratio = %d" % b)
        ax.legend()
        ax.grid()
        ax.set_title("nose finess ratio = %d" % n)
        ax.set_xlabel("tail finess ratio")
        ax.set_ylabel("fuselage $C_{dp}$")
        fig.savefig("fuse_drag_nose%d.pdf" % n)
        figs.append(fig)

    return figs

if __name__ == "__main__":
    X, Y = fit_setup("fusedrag.csv")
    df = pd.read_csv("fusedrag.csv")
    F = plot_fits("fusedrag.csv")
