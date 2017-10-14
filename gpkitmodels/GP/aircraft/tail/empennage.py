" empennage.py "
from gpkit import Variable, Model
from .horizontal_tail import HorizontalTail
from .vertical_tail import VerticalTail
from .tail_boom import TailBoom, TailBoomState

#pylint: disable=attribute-defined-outside-init
class Empennage(Model):
    "empennage model, consisting of vertical, horizontal and tailboom"
    def setup(self):
        mfac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = Variable("W", "lbf", "empennage weight")

        self.htail = HorizontalTail()
        self.htail.substitutions.update({"m_{fac}": 1.1})
        self.vtail = VerticalTail()
        self.vtail.substitutions.update({"m_{fac}": 1.1})
        self.tailboom = TailBoom()
        self.components = [self.htail, self.vtail, self.tailboom]

        state = TailBoomState()
        loading = [self.tailboom.horizontalbending(self.htail, state),
                   self.tailboom.verticalbending(self.vtail, state),
                   self.tailboom.verticaltorsion(self.vtail, state)]

        constraints = [
            W/mfac >= (self.htail.topvar("W") + self.vtail.topvar("W")
                       + self.tailboom["W"]),
            self.tailboom["l"] >= self.htail["l_h"],
            self.tailboom["l"] >= self.vtail["l_v"],
            ]

        return self.components, constraints, loading
