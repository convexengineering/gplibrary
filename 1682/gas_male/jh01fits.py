from gpfit.fit import fit
import pandas as pd
from numpy import log

df = pd.read_csv("jh01polar.csv")

cl = df["CL"][df["CL"] >= 0.4].values
cd = df["CD"][df["CL"] >= 0.4].values

x = log(cl)
y = log(cd)

a, b = fit(x, y, 3, "SMA")
