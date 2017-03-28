from gpfit.fit import fit
from numpy import logspace, log, log10, array

n = 19

u = array([6.5617,
           9.8425,
           13.1234,
           16.4042,
           19.6850,
           22.9659,
           26.2467,
           29.5276,
           32.8084,
           36.0892])

u = u*1000*0.3048

lo = array([15.0000,
            18.5714,
            22.1429,
            25.7143,
            29.2857,
            32.8571,
            36.4286,
            40.0000,
            44.0000,
            46.0000])

x = log(u)  
y = log(lo)  # log(lo) # log(mi)

K = 3
Type = "SMA"

cstrt, rms_error = fit(x,y,K,Type)

print rms_error
