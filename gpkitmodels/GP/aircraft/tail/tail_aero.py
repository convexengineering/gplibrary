" tail aerodynamics "
from gpkit import Variable, Model
import os
from gpkitmodels.tools.fit_constraintset import FitCS
import pandas as pd

class TailAero(Model):
    "horizontal tail aero model"
    def setup(self, static, state):

        name = get_lowername(static.__class__.__name__)
        Re = Variable("Re", "-", "%s reynolds number" % name)
        Cd = Variable("C_d", "-", "%s drag coefficient" % name)

        path = os.path.dirname(__file__)
        df = pd.read_csv(path + os.sep + "tail_dragfit.csv")

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["S"]/static["b"]
                   / state["\\mu"]),
            FitCS(df, Cd, [Re, static["\\tau"]], airfoil="naca 0008",
                  err_margin="RMS")
            ]

        return constraints

def get_lowername(classname):
    start = [c for c in classname if c.isupper()]
    name = [classname]
    for t in start:
        name = name[-1].split(t)

    n = " ".join([t.lower()+n for n, t in zip(name, start)])
    return n
