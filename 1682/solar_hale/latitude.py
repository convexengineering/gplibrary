import numpy as np
import pandas as pd
from solar_hale import SolarHALE
from plotting import latitude_sweep

M = SolarHALE()
M.substitutions.update({M["h"]: 15000})
df = pd.read_csv("solar_irr_vs_lat.csv")
df = df[df!=0.0]
df = df.dropna()
xsweeps = np.array([df.Winter_Solstice, df.DayLight, 24 - df.DayLight])
xvarnames = ["(E/S)_{irr}", "t_{day}", "t_{night}"]

fig, ax = latitude_sweep(M, np.array(df.Latitude), xvarnames, xsweeps, "V_{wind}", [20,25,30,35], "b", [0, 200])
fig.savefig("b_vs_latitude15.pdf")
