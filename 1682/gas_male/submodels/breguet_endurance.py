" breguet_endurance.py "
from gpkit import Model, Variable
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS

class BreguetEndurance(Model):
    "breguet endurance model"
    def __init__(self, perf, **kwargs):
        z_bre = Variable("z_{bre}", "-", "Breguet coefficient")
        t = Variable("t", "days", "Time per flight segment")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        Wfuel = Variable("W_{fuel}", "lbf", "Segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        constraints = [
            TCS([z_bre >= (perf["P_{total}"]*t*perf["BSFC"]*g
                           / (perf["W_{end}"]*perf["W_{start}"])**0.5)]),
            f_fueloil*Wfuel/perf["W_{end}"] >= te_exp_minus1(z_bre, 3),
            perf["W_{start}"] >= perf["W_{end}"] + Wfuel
            ]

        Model.__init__(self, None, constraints, **kwargs)
