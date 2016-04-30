import pandas as pd 
import numpy as np 
from numpy import logspace, log, log10
import matplotlib.pyplot as plt
# you have to change it to this and run the file form gpfit
from fit import fit

# since this is not an installable code source you have to navigate to the
# ~gpfit folder.  

# Fitting BSFC vs. RPM
# df = pd.read_csv('Dataset_BSFC_kgKwh.csv')
# RPM = df['RPM']
# RPM_max = np.amax(RPM)
# BSFC = df['BSFC']
# BSFC_min = np.amin(BSFC)
# logRPM = log(RPM/RPM_max)
# logBSFC = log(BSFC/BSFC_min)
# Type = 'SMA'
# K = 2

# cstrt, rms_error = fit(logRPM,logBSFC,K,Type)

# plt.plot(logRPM,logBSFC)
# plt.plot()

# Fitting Power vs. RPM
# df = pd.read_csv('Dataset_Power_Kw.csv')
# RPM = df['RPM']
# RPM_max = np.amax(RPM)
# P = df['P']
# P_max = np.amax(P)
# logP = log(P/P_max)
# logRPM = log(RPM/RPM_max)
# Type = 'SMA'
# K = 1
# cstrt, rms_error = fit(logRPM,logP,K,Type)

# Fitting Torque vs. RPM
df = pd.read_csv('Dataset_Torque_Nm.csv')
RPM = df['RPM']
RPM_max = np.amax(RPM)
Q = df['Q']
Q_max = np.amax(Q)
logRPM = log(RPM/RPM_max)
logQ = log(Q/Q_max)
Type = 'SMA'
K = 1
cstrt, rms_error = fit(logRPM,logQ,K,Type)

#h = np.delete(h,[0])
#P = np.delete(P,[0])

# Fitting Power vs. RPM
# P = df['kW']
# P_MSL = np.amax(df['kW'])
# Pnorm = P/P_MSL
# L = 1-Pnorm

# the variables you input into fit() can't have negatives here or Inf or -Inf
# logL = log(L)
# logh = log(h)



# plt.xlabel('log(h)')
# plt.ylabel('log(1-Pnorm)')
# #plt.show()

#cstrt, rms_error = fit(logh[1:-1],logL[1:-1],K,Type)


