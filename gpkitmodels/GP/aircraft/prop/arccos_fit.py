import unittest
from numpy import log, exp, log10, vstack
from numpy import arccos,arange
from gpfit.fit import fit
from numpy.random import random_sample
i = arange(0.0001,1,.001)
j = arccos(exp(-i))
#P = Vdd**2 + 30*Vdd*exp(-(Vth-0.06*Vdd)/0.039)
#u = vstack((Vdd,Vth))
x = log(i)
y = log(j)
K = 1

cstrt, rmsErr = fit(x,y,K,"SMA")
print rmsErr