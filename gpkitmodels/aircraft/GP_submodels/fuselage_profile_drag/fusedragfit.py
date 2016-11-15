" fuselage drag fits "
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
    "fit using SMA, K = 4"
    w = (0.0097 * (u_1)**0.866 * (u_2)**-0.209 * (u_3)**-0.156
         + 0.0436 * (u_1)**0.0546 * (u_2)**0.258 * (u_3)**-1.43
         + 0.00243 * (u_1)**0.0336 * (u_2)**1.22 * (u_3)**0.306
         + 0.00255 * (u_1)**-0.0317 * (u_2)**-0.585 * (u_3)**1.15)**(1/0.996)
    return w

def plot_fits(filename):
    "plot fit against data"

    df = pd.read_csv(filename)
    u1 = np.array(df["tubelr"])
    u2 = np.array(df["noselr"])
    u3 = np.array(df["taillr"])
    c = np.unique(u1)
    n = np.unique(u2)
    t = np.unique(u3)

    fig, ax = plt.subplots()

    for a in c:
        d = df[(df["tubelr"]==a)]
        for b in n:
            data = d[(df["noselr"]==b)]
            ax.plot(data["taillr"], data["cd_front"])




if __name__ == "__main__":
    X, Y = fit_setup("fusedrag.csv")
    df = plot_fits("fusedrag.csv")
