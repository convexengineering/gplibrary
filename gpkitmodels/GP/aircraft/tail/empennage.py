" empennage.py "
from gpkit import Model, parse_variables
from .horizontal_tail import HorizontalTail
from .vertical_tail import VerticalTail
from .tail_boom import TailBoom, TailBoomState

#pylint: disable=attribute-defined-outside-init, no-member, exec-used
#pylint: disable=too-many-instance-attributes, invalid-name, undefined-variable
class Empennage(Model):
    """empennage model, consisting of vertical, horizontal and tailboom

    Variables
    ---------
    mfac        1.0     [-]     tail weight margin factor
    W                   [lbf]   empennage weight

    Upper Unbounded
    ---------------
    W, Vv, Vh

    Lower Unbounded
    ---------------
    lv, lh, Vv, Vh, bv, bh, mh

    LaTex Strings
    -------------
    mfac        m_{\\mathrm{fac}}

    """
    def setup(self):
        exec parse_variables(Empennage.__doc__)

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
        self.tailboom.tailLoad.__name__ = "HTailBoomBending"
        self.hbend = self.tailboom.tailLoad(self.tailboom, self.htail, state)
        self.tailboom.tailLoad.__name__ = "VTailBoomBending"
        self.vbend = self.tailboom.tailLoad(self.tailboom, self.vtail, state)
        loading = [self.hbend, self.vbend]

        constraints = [
            W/mfac >= sum(c.W for c in self.components),
            l >= lh, l >= lv,
            ]

        return self.components, constraints, loading
