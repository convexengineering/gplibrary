import pandas as pd 
import numpy as np 
from numpy import logspace, log, log10
import matplotlib.pyplot as plt
# you have to change it to this and run the file form gpfit
from gpfit.fit import fit

# since this is not an installable code source you have to navigate to the
# ~gpfit folder.  
df = pd.read_csv('../gpkit-models/1682/DF35_maxPvh.csv')

h = df['ft']
#h = np.delete(h,[0])
P = df['kW']
#P = np.delete(P,[0])

P_MSL = np.amax(df['kW'])

Pnorm = P/P_MSL

L = 1-Pnorm

# the variables you input into fit() can't have negatives here or Inf or -Inf
logL = log(L)
logh = log(h)

plt.plot(logL,logh)
plt.xlabel('log(h)')
plt.ylabel('log(1-Pnorm)')
#plt.show()

K = 2
Type = 'SMA'

cstrt, rms_error = fit(logh,logL,K,Type)
