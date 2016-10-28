# Gas vs Solar Trade

# Simple model of a Gas Powered Aircraft
```python
#inPDF: skip
from gassimple import GasSimple
from solar_hale import SolarSimple
from plotting import plot_sweep
from gen_tex import gen_tex_fig
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size':19})

# M1 = SolarSimple(latitude=40)
# sol1 = M1.solve("mosek")

""" comparision plots """
# M = GasSimple()
# sol = M.solve("mosek")
# exprof = {"orion": {"endurance": 5, "gtw": 11200, "span": 132, "speed": 34, 
#                     "altitude": 15000, "payload": 1500}, 
#           "predator": {"endurance": 1, "gtw": 2250, "span": 55, "speed": 70,
#                        "altitude": 25000, "payload": 450},
#           "hermes": {"endurance": 1.5, "gtw": 2600, "span": 50, "speed": 70,
#                      "altitude": 30000, "payload": 771},
#           "globalhawk": {"endurance": 1.5, "gtw": 32250, "span": 131, 
#                          "speed": 160, "altitude": 50000, "payload": 3000},
#           "penguin": {"endurance": 2, "gtw": 47, "span": 10.8, 
#                      "speed": 22, "altitude": 10000, "payload": 22}}
# mtow = []
# ext, t = [], []
# exb, b = [], []
# for plane in exprof:
#     a = exprof[plane]
#     M = GasSimple(altitude=a["altitude"])
#     M.substitutions.update({"\\eta_{prop}": np.linspace(0.65, 0.65, 5)})
#     M.substitutions.update({"c_{d_0}": 0.01})
#     M.substitutions.update({"MTOW": a["gtw"], "W_{pay}": a["payload"], 
#                             "V_{min}": a["speed"]})
#     sol = M.solve("mosek")
#     mtow.append(a["gtw"])
#     ext.append(a["endurance"])
#     exb.append(a["span"])
#     t.append(sol("t_{flight}").magnitude)
#     b.append(sol("b").magnitude)
# 
# fig, ax = plt.subplots()
# ax.plot(mtow, t, "*", markersize=10, markerfacecolor="None")
# ax.plot(mtow, ext, "o", markersize=10, markerfacecolor='None')
# ax.set_ylim([0, 5.5])
# ax.set_xlabel("Max Take Off Weight [lbf]")
# ax.set_ylabel("Endurance [days]")
# ax.legend(["Prediction Value", "Actual Value"])
# ax.grid()
# fig.savefig("figs/mtowvehiclecomp.pdf", bbox_inches="tight")
# 
# fig, ax = plt.subplots()
# ax.plot(mtow, b, "*", markersize=10, markerfacecolor="None")
# ax.plot(mtow, exb, "o", markersize=10, markerfacecolor='None')
# ax.set_ylim([0, 180])
# ax.set_xlabel("Max Take Off Weight [lbf]")
# ax.legend(["Prediction Value", "Actual Value"])
# ax.set_ylabel("Span [b]")
# ax.grid()
# fig.savefig("figs/bvehiclecomp.pdf", bbox_inches="tight")

""" latitude sweep USED """
lats = np.arange(1,60,2)
avails = [80, 85, 90, 95]
fig, ax = plt.subplots()
for a in avails:
    mtow = []
    for i, l in enumerate(lats):
        M = GasSimple(latitude=l, avail=a)
        M.substitutions.update({"W_{pay}": 10})
        M.substitutions.update({"\\eta_{prop}": np.linspace(0.75, 0.75, 5)})
        M.substitutions.update({"c_{d_0}": 0.005})
        del M.substitutions["MTOW"]
        M.substitutions.update({"t_{flight}": 8})
        M.cost = M["MTOW"]
        try:
            sol = M.solve("mosek")
            mtow.append(sol("MTOW").magnitude)
        except RuntimeWarning:
            mtow.append(np.nan)
    ax.plot(lats, mtow)

ax.set_xlabel("Latitude [deg]")
ax.set_ylabel("Max take off weight [lbf]")
ax.set_ylim([0, 1000])
ax.grid()
ax.legend(["%d Percentile Winds" % a for a in avails], loc=2, fontsize=15)
fig.savefig("latvsmtow.pdf", bbox_inches="tight")

""" gas solar mtow comparison """

fig, ax = plt.subplots()
avails = [85]
colors = ["b", "r"]
labels = []
for a, c in zip(avails, colors):
    G = GasSimple(latitude=45, avail=a, altitude=16000)
    S = SolarSimple(latitude=30, avail=a, altitude=50000)
    G.substitutions.update({"W_{pay}": 10})
    G.substitutions.update({"\\eta_{prop}": np.linspace(0.75, 0.75, 5)})
    G.substitutions.update({"c_{d_0}": 0.005})
    S.substitutions.update({"\\eta_{prop}": np.linspace(0.75, 0.75, 1)})
    S.substitutions.update({"c_{d_0}": 0.002})
    G.substitutions.update({"t_{flight}": ("sweep", np.linspace(1, 15, 20))})
    del G.substitutions["MTOW"]
    G.cost = G["MTOW"]
    S.substitutions.update({"W_{pay}": 10})
    ssol = S.solve("mosek")
    gsol = G.solve("mosek", skipsweepfailures=True)
    ax.plot(gsol("t_{flight}").magnitude, gsol("MTOW").magnitude, "--")
    ax.plot([0,15], [ssol("MTOW").magnitude]*2)

ax.set_ylabel("Max take off weight [lbf]")
ax.set_xlabel("Endurance [days]")
ax.set_ylim([0, 500])
ax.set_xlim([0, 15])
ax.grid()
ax.legend(["Gas Powered", "Solar Powered"], loc=2, fontsize=15)
fig.savefig("solargasmtow.pdf", bbox_inches="tight")

""" solar latitude USED """

fig, ax = plt.subplots()
lat = np.arange(0, 60, 1)
for a in [80, 85, 90, 95]:
    mtow = []
    for l in lat:
        S = SolarSimple(latitude=l, avail=a, altitude=50000)
        S.substitutions.update({"W_{pay}": 10})
        S.substitutions.update({"\\eta_{prop}": np.linspace(0.75, 0.75, 1)})
        S.substitutions.update({"c_{d_0}": 0.002})
        try:
            sol = S.solve("mosek")
            mtow.append(sol("MTOW").magnitude)
        except RuntimeWarning:
            mtow.append(np.nan)
    ax.plot(lat, mtow)

ax.set_ylim([0, 1000])
ax.grid()
ax.set_xlabel("Latitude [deg]")
ax.set_ylabel("Max take off weight [lbf]")
ax.legend(["%d Percentile Winds" % a for a in [80, 85, 90, 95]], loc=2, fontsize=15)
fig.savefig("mtowvslatsolar.pdf", bbox_inches="tight")

""" solar contours """
for l in [30, 35, 40]:
    for av in [80, 85, 90]:
        fig, ax = plt.subplots()
        S = SolarSimple(latitude=l, avail=av, altitude=50000)
        S.substitutions.update({"f_{structures}": ("sweep", np.linspace(0.2, 0.5, 10))})
        S.substitutions.update({"h_{batt}": ("sweep", np.linspace(250, 400, 10))})
        S.substitutions.update({"W_{pay}": 10})
        S.substitutions.update({"\\eta_{prop}": np.linspace(0.75, 0.75, 1)})
        S.substitutions.update({"c_{d_0}": 0.002})
        S.cost = S["b"]
        sol = S.solve("mosek", skipsweepfailures=True)
        x = np.reshape(sol("f_{structures}"), [10, 10])
        y = np.reshape(sol("h_{batt}"), [10, 10])
        z = np.reshape(sol("b"), [10, 10])
        levels = np.array(range(50, 2000, 50)+ [2300])
        if av == 90:
            v = np.array(range(50, 700, 50)+ [2300])
        else:
            v = np.array(range(50, 400, 50)+ [2300])
        a = ax.contour(x, y, z, levels, colors="k")
        ax.clabel(a, v, inline=1, fmt="%d [ft]")
        ax.set_xlabel("Structural Fraction")
        ax.set_ylabel("Battery Energy Density [Whr/kg]")
        fig.savefig("bcontourl%da%d.pdf" % (l, av), bbox_inches="tight")

```


