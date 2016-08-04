"""Commercial Aircraft Model loosely inspired by TASOPT"""
import numpy as np
from gpkit import Variable, Model, units


rho = Variable("\\rho", 0.4, "kg/m^3", "air density")


class Wing(Model):
    def __init__(self, **kwargs):
        e = Variable("e", 0.85, "-",
                     "Oswald efficiency, should come from aircraft wing")
        S = Variable("S", 1300, "ft^2", "wing area")
        A = Variable("A", 9.45, "-", "Aspect ratio")
        W = Variable("W", 1, "lbf", "weight")
        b = Variable("b", "ft", "wing span")
        CD = Variable("C_D", "-", "drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")

        Model.__init__(self, None, [b**2 == S/A,
                                    CD >= 0.02 + CL**2/np.pi/e/A,
                                    W == W])
        # self.varkeys.add(W)


class SteadyFlight(Model):
    def __init__(self, **kwargs):
        CL = Variable("C_L", "-", "lift coefficient")
        CD = Variable("C_D", "-", "drag coefficient")
        V = Variable("V", 200, "knots", "true airspeed")
        T = Variable("T", "lbf", "thrust")
        L = Variable("L", "lbf", "lift")
        Sref = Variable("S_{ref}", "ft^2", "reference area")
        Wfuel = Variable("W_{fuel}", 36864.9, "lbf", "fuel weight")
        Wref = Variable("W_{ref}", "lbf", "dry mission weight")
        W = Variable("W", "lbf", "total weight")

        constraints = [L <= 0.5*rho*V**2*CL*Sref,
                       T >= 0.5*rho*V**2*CD*Sref,
                       W >= Wref + Wfuel,
                       W <= L,
                       ]

        Model.__init__(self, None, constraints, **kwargs)


class CommercialAircraft(Model):
    def __init__(self, **kwargs):
        W_e = Variable("W_e", 84749.8, "lbf", "empty weight")
        Wpay = Variable("W_{pay}", 38700, "lbf", "payload weight")
        W = Variable("W", "lbf", "total weight")

        wing = Wing()
        constraints = [wing]

        sf = SteadyFlight()
        sf.subinplace({sf["S_{ref}"]: wing["S"],
                       sf["C_L"]: wing["C_L"],
                       sf["C_D"]: wing["C_D"]})
        constraints.append(sf)
        constraints.append(sf["W_{ref}"] >= W_e + Wpay + wing["W"])

        Model.__init__(self, sf["T"], constraints, **kwargs)


m = CommercialAircraft()
print m.substitutions
sol = m.solve()
print sol.table()
