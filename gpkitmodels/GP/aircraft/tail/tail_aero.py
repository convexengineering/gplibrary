" tail aerodynamics "
from gpkit import Variable, Model
import os
from gpfit.fit_constraintset import XfoilFit
import pandas as pd

class TailAero(Model):
    "horizontal tail aero model"
    def setup(self, static, state):

        name = get_lowername(static.__class__.__name__)
        Re = Variable("Re", "-", "%s reynolds number" % name)
        Cd = Variable("C_d", "-", "%s drag coefficient" % name)

        path = os.path.dirname(__file__)
        fd = pd.read_csv(path + os.sep + "tail_dragfit.csv").to_dict(
            orient="records")[0]

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["S"]/static["b"]
                   / state["\\mu"]),
            XfoilFit(fd, Cd, [Re, static["\\tau"]],
                     err_margin="RMS", airfoil="naca 0008")
            ]

        return constraints

def get_lowername(classname):
    start = [c for c in classname if c.isupper()]
    name = [classname]
    for t in start:
        name = name[-1].split(t)

    n = " ".join([t.lower()+n for n, t in zip(name, start)])
    return n
