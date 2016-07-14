import numpy as np
import pandas as pd

df = pd.read_csv("solar_irr_vs_lat.csv")
df = df.set_index("Latitude").reindex(index=np.arange(0, 90, 2)).interpolate()
df = df[df.index <= 50]
wind = pd.read_csv("wind.csv")
windcols = ["%sth Percentile Winds" % lat for lat in (80, 90, 95, 99)]
df[windcols] = wind[windcols]
df.to_csv("solar_irr_vs_lat_interpolated.csv")
