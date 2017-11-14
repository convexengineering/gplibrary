" empennage.py "
from gpkit import Variable, Model
from .horizontal_tail import HorizontalTail
from .vertical_tail import VerticalTail
from .tail_boom import TailBoom, TailBoomState

#pylint: disable=attribute-defined-outside-init, no-member
class Empennage(Model):
    """empennage model, consisting of vertical, horizontal and tailboom

    Upper Unbounded
    ---------------
    W, Vv, Vh

    Lower Unbounded
    ---------------
    lv, lh, Vv, Vh, bv, bh, mh

    """
    def setup(self):
        mfac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = self.W = Variable("W", "lbf", "empennage weight")

        self.htail = HorizontalTail()
        self.htail.substitutions.update({self.htail.mfac: 1.1})
        lh = self.lh = self.htail.lh
        self.Vh = self.htail.Vh
        self.bh = self.htail.b
        self.mh = self.htail.mh
        self.vtail = VerticalTail()
        self.vtail.substitutions.update({self.vtail.mfac: 1.1})
        lv = self.lv = self.vtail.lv
        self.Vv = self.vtail.Vv
        self.bv = self.vtail.b
        self.tailboom = TailBoom()
        self.components = [self.htail, self.vtail, self.tailboom]
        l = self.l = self.tailboom.l

        state = TailBoomState()
        loading = [self.tailboom.horizontalbending(self.htail, state),
                   self.tailboom.verticalbending(self.vtail, state),
                   self.tailboom.verticaltorsion(self.vtail, state)]

        constraints = [
            W/mfac >= sum(c.W for c in self.components),
            l >= lh, l >= lv,
            ]

        return self.components, constraints, loading
