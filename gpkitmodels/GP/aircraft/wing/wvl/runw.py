import numpy as np
from wvl import wvl
from gpkit import Variable, Model, VectorVariable

"""
Geometric parameters
    XLE     location of leading edge point
    yLE     spanwise location
    z       vertical height (dihedrial)
    c       chord length (along x)
    a0      airfoil zero-lift angle (deg)
    tw      secton twist angle
"""

#                 xLE, y,   z,   c,    tw,   a0
geom = np.array([[0.0, 0.0, 0.0, 0.22, 0.0, -6.0],
                 [0.0, 1.0, 0.2, 0.11, 0.0, -6.0]])

Sref = 0.33
bref = 2.0

Sspec = Sref
ztspec = 0.0
xcaxis = 0.4

alspec = 0.0
CLspec = 1.1

bespec = 0.0
pbspec = 0.0
rbspec = 0.0

Crspec = 0.0
Cnspec = 0.0

ialspec = 0
iCLspec = 1

ibespec = 1

ipbspec = 1
irbspec = 1

iCrspec = 0
iCnspec = 0

N = 10

ispace = 2

itmax = 20

toler = 1.0e-6

ns = geom.shape[0]

nflow = 1

alspec *= (np.pi/180)
bespec *= (np.pi/180)

b = 2*geom[ns-1, 1]

geom[:, 0:4] *= 2.0/bref
geom[:, 4:6] *= (np.pi/180)

S = 0
for i in range(ns-1):
    dy = geom[i+1, 1] - geom[i, 1]
    dz = geom[i+1, 2] - geom[i, 2]
    ca = 0.5*(geom[i, 3] + geom[i+1, 3])
    S += ca*dy

S *= 2
cavg = S/b
AR = b**2/S

if not Sspec == 0.0:
    cfac = Sspec/S
    geom[:, 3] *= cfac
    S *= cfac
    cavg *= cfac
    AR /= cfac
    print ("All chords scaled by factor of %6.4f to obtain specified "
           "S=%6.3f\n" % (cfac, Sspec))

if not ztspec == 0.0:
    zfac = ztspec/geom[ns-1, 2]
    geom[:, 2] *= zfac
    print ("All z scaled by factor of %6.4f to obtain specified "
           "ztip=%6.3f\n" % (zfac, ztspec))

if xcaxis >= 0.0:
    x0 = geom[0, 0]
    c0 = geom[0, 3]
    geom[:, 0] = x0 + (c0 - geom[:,3])*xcaxis
    print ("All xLEs shifted to obtain specified x/c=%6.3f line "
           "straight\n" % xcaxis)

cref = cavg

A, yv,zv,cl,ccl,vi,wi,alpha,beta,pbar,rbar,CL,CDi,Cr,Cn,Cb = wvl(geom, N, ispace, Sref, bref, cref, itmax, toler, alspec, bespec, pbspec, rbspec, CLspec, Crspec, Cnspec, ialspec, ibespec, ipbspec, irbspec, iCLspec, iCrspec, iCnspec)

# Ainv = np.linalg.inv(A)
Ainv = np.array([1])

class WVL(Model):
    def setup(self):

        Na = Ainv.shape[0]

        Aijm = VectorVariable([Na, Na], "A_{i,j}^{-1}", Ainv, "-",
                              "AIC matrix")
        G = VectorVariable(Na, "\\Gamma", "-", "vortex filament strength")
        th = VectorVariable(Na, "\\theta", "-", "twist")
        CL = Variable("C_L", "-", "coefficient of lift")
        CDi = Variable("C_{D_i}", "-", "induced drag coefficient")

        constraints = [G >= np.sum(Aijm*th, 1)]
