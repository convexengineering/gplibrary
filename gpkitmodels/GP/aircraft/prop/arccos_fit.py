import unittest
from numpy import log, exp, log10, vstack
from numpy import arccos,arange
from gpfit.fit import fit
from numpy.random import random_sample
i = arange(0.0001,3,.001)
j = arccos(exp(-i))
x = log(i)
y = log(j)
K = 1

cstrt, rmsErr = fit(x,y,K,"SMA")
print rmsErr