from gpfit.fit import fit
from numpy import logspace, log, log10, array

n = 19

u = array([36.0892,
           39.3701,
           42.6509,
           45.9318,
           49.2126,
           52.4934,
           55.7743,
           59.0551,
           62.3360,
           65.6168])

u = u*1000*0.3048

lo = array([46.,
           43.,
           39.,
           35.,
           32.,
           30.,
           27.,
           24.,
           24.,
           23])

x = log(u)  
y = log(lo)  # log(lo) # log(mi)

K = 3
Type = "SMA"

cstrt, rms_error = fit(x,y,K,Type)

print rms_error
