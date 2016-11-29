" tail aerodynamics "
from gpkit import Variable, Model

class TailAero(Model):
    "horizontal tail aero model"
    def __init__(self, static, state, **kwargs):

        name = get_lowername(static.__class__.__name__)
        Re = Variable("Re", "-", "%s reynolds number" % name)
        Cd = Variable("C_d", "-", "%s drag coefficient" % name)

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["S"]/static["b"]
                   / state["\\mu"]),
            Cd**70.5599 >= (7.42688e-90*(Re/1000)**-33.0637
                            * (static["\\tau"]*100)**18.0419
                            + 5.02826e-163*(Re/1000)**-18.7959
                            * (static["\\tau"]*100)**53.1879
                            + 4.22901e-77*(Re/1000)**-41.1704
                            * (static["\\tau"]*100)**28.4609)
            ]

        Model.__init__(self, None, constraints, **kwargs)

def get_lowername(classname):
    start = [c for c in classname if c.isupper()]
    name = [classname]
    for t in start:
        name = name[-1].split(t)

    n = " ".join([t.lower()+n for n, t in zip(name, start)])
    return n
