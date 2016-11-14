" empennage.py "
import numpy as np
from gpkit import Variable, Model
from horizontal_tail import HorizontalTail
from vertical_tail import VerticalTail
from tail_boom import TailBoom, TailBoomState

class Empennage(Model):
    "empennage model, consisting of vertical, horizontal and tailboom"
    def __init__(self, **kwargs):
        mfac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = Variable("W", "lbf", "empennage weight")

        self.horizontaltail = HorizontalTail()
        self.verticaltail = VerticalTail()
        self.tailboom = TailBoom()
        self.components = [self.horizontaltail, self.verticaltail,
                           self.tailboom]

        self.loading = EmpennageLoading

        constraints = [
            W/mfac >= (self.horizontaltail["W"] + self.verticaltail["W"]
                       + self.tailboom["W"]),
            self.tailboom["l"] >= self.horizontaltail["l_h"],
            self.tailboom["l"] >= self.verticaltail["l_v"],
            ]

        Model.__init__(self, None, [self.components, constraints],
                       **kwargs)


class EmpennageLoading(Model):
    "tail boom loading case"
    def __init__(self, empennage, **kwargs):
        state = TailBoomState()

        loading = [empennage.tailboom.horizontalbending(
            empennage.tailboom, empennage.horizontaltail, state)]
        loading.append(empennage.tailboom.verticalbending(
            empennage.tailboom, empennage.verticaltail, state))
        loading.append(empennage.tailboom.verticaltorsion(
            empennage.tailboom, empennage.verticaltail, state))

        Model.__init__(self, None, loading, **kwargs)
