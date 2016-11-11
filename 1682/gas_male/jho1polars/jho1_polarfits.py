"jho1_polarfits.py"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def fit_setup(Re_range):
    "set up x and y parameters for gp fitting"
    CL = []
    CD = []
    RE = []
    for r in Re_range:
        dataf = text_to_df("jh01.ncrit09.Re%dk.pol" % r)
        CL.append(dataf["CL"])
        CD.append(dataf["CD"])
        RE.append([r*1000.0]*len(dataf["CL"]))

    u1 = np.hstack(CL)
    u2 = np.hstack(RE)
    w = np.hstack(CD)
    u1 = u1.astype(np.float)
    u2 = u2.astype(np.float)
    w = w.astype(np.float)
    u = [u1, u2]
    x = np.log(u)
    y = np.log(w)
    return x, y

def return_fit(cl, re):
    "polar fit for the JHO1 airfoil"
    cd = (0.0247*cl**2.49*re**-1.11 + 2.03e-7*cl**12.7*re**-0.338 +
          6.35e10*cl**-0.243*re**-3.43 + 6.49e-6*cl**-1.9*re**-0.681)**(1/3.72)
    # SMA function, K=3, max RMS error = 0.00489
    return cd

def plot_fits(re):
    "plot fit compared to data"
    colors = ["k", "m", "b", "g", "y"]
    assert len(re) == len(colors)
    fig, ax = plt.subplots()
    cls = np.linspace(0.2, 1.3, 20)
    for r, col in zip(re, colors):
        dataf = text_to_df("jh01.ncrit09.Re%dk.pol" % r)
        ax.plot(dataf["CL"], dataf["CD"], "o", c=col)
        cd = return_fit(cls, r*1000.)
        ax.plot(cls, cd, c=col, label="Re = %dk" % r)
    ax.legend(loc=2)
    ax.set_xlabel("$C_L$")
    ax.set_ylabel("$C_D$")
    ax.grid()
    return fig, ax

if __name__ == "__main__":
    Re = range(200, 750, 50)
    X, Y = fit_setup(Re) # call fit(X, Y, 4, "SMA") to get fit
    F, A = plot_fits([300, 350, 400, 450, 500])
    F.savefig("jh01polarfit.pdf")
