import numpy as np
import numpy as np
from gasmale import GasPoweredMALE
import plot

fixed = True
# PLOTS #
fixedPLOTS = False
missionPLOTS = False
wind = False
W_pay = False
P_pay = False
altitude = True

NLoiter = 20
NClimb1, NClimb2 = 10, 10
NCruise1, NCruise2 = 5, 5
NClimb = NClimb1 + NClimb2
NCruise = NCruise1 + NCruise2
NSeg = NLoiter + NClimb + NCruise
mStart = 0
mEndClimb = NClimb1
mEndCruise = NClimb1 + NCruise1
mEndClimb2 = NClimb1 + NCruise1 + NClimb2
mEndLoiter = NSeg - NCruise2
mEnd = NSeg
iClimb1 = range(0, mEndClimb)
iClimb2 = range(mEndCruise, mEndClimb2)
iClimb = range(0, mEndClimb) + range(mEndCruise, mEndClimb2)
iCruise = range(mEndClimb, mEndCruise) + range(mEndLoiter, mEnd)
iCruise1 = range(mEndClimb, mEndCruise)
iCruise2 = range(mEndLoiter, mEnd)
iLoiter = range(mEndClimb2, mEndLoiter)
t = np.linspace(0,NSeg-1,NSeg)

M = GasPoweredMALE(NLoiter=NLoiter, NClimb1=NClimb1, NClimb2=NClimb2, NCruise1=NCruise1, NCruise2=NCruise2, wind=wind)
sol = M.solve('mosek')

if fixed:
    S_fix = sol('S').magnitude
    Sfuse_fix = sol('S_{fuse}').magnitude
    b_fix = sol('b').magnitude
    lfuse_fix = sol('l_{fuse}').magnitude
    Volfuse_fix = sol('Vol_{fuse}').magnitude
    hspar_fix = sol('h_{spar}').magnitude
    Volfuel_fix = sol('Vol_{fuel}').magnitude

    #wind = True
    M = GasPoweredMALE(NLoiter=NLoiter, NClimb1=NClimb1, NClimb2=NClimb2, NCruise1=NCruise1, NCruise2=NCruise2, wind=wind)

    M.substitutions.update({'S': S_fix})
    #M.substitutions.update({'V_{wind}': 5})
    #M.substitutions.update({'S_{fuse}': Sfuse_fix})
    M.substitutions.update({'b': b_fix})
    #M.substitutions.update({'l_{fuse}': lfuse_fix})
    M.substitutions.update({'Vol_{fuse}': Volfuse_fix+0.0001})
    M.substitutions.update({'Vol_{fuel}': Volfuel_fix+0.0001})
    #M.substitutions.update({'h_{spar}': hspar_fix})
    sol = M.solve('mosek', verbosity=0)
    V = sol('V')

    if fixedPLOTS:
        plot.fixed(t, V, NSeg, NClimb1, NCruise1, NClimb2, NLoiter, NCruise2,
                   mStart, mEndClimb, mEndCruise, mEndClimb2, mEndLoiter)

if missionPLOTS:
    plot.mission(sol, t, V, NSeg, NClimb1, NCruise1, NClimb2, NLoiter, NCruise2,
                 mStart, mEndClimb, mEndCruise, mEndClimb2, mEndLoiter)

numplots = sum([wind, P_pay, W_pay, altitude])
if numplots > 1:
    raise ValueError("more than one is true, but there can only be one!")
elif numplots == 1:
    if "t_{station}" in M.substitutions:
        del M.substitutions["t_{station}"]
    M.cost = 1/M["t_{station}"]
    if wind:
        plot.wind(M)
    elif P_pay:
        plot.P_pay(M)
    elif W_pay:
        plot.W_pay(M)
    elif altitude:
        plot.altitude(M)
