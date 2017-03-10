" empennage.py "
from gpkit import Variable, Model
from horizontal_tail import HorizontalTail
from vertical_tail import VerticalTail
from tail_boom import TailBoom, TailBoomState

class Empennage(Model):
    "empennage model, consisting of vertical, horizontal and tailboom"
    def setup(self):
        mfac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = Variable("W", "lbf", "empennage weight")

        self.horizontaltail = HorizontalTail()
        self.verticaltail = VerticalTail()
        self.tailboom = TailBoom()
        self.components = [self.horizontaltail, self.verticaltail,
                           self.tailboom]

        constraints = [
            W/mfac >= (self.horizontaltail["W"] + self.verticaltail["W"]
                       + self.tailboom["W"]),
            self.tailboom["l"] >= self.horizontaltail["l_h"],
            self.tailboom["l"] >= self.verticaltail["l_v"],
            ]

        return self.components, constraints

    def loading(self):
        return EmpennageLoading(self)

class EmpennageLoading(Model):
    "tail boom loading case"
    def setup(self, empennage):
        state = TailBoomState()

        loading = [empennage.tailboom.horizontalbending(
            empennage.horizontaltail, state)]
        loading.append(empennage.tailboom.verticalbending(
            empennage.verticaltail, state))
        loading.append(empennage.tailboom.verticaltorsion(
            empennage.verticaltail, state))

        return loading
