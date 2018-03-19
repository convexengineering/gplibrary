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
    W, vtail.Vv, htail.Vh, tailboom.cave (if not hSparModel)

    Lower Unbounded
    ---------------
    htail.lh, htail.Vh, htail.b, htail.mh
    vtail.lv, vtail.Vv, vtail.b
    htail.Sy (if hSparModel), htail.J (if hSparModel)
    vtail.Sy (if vSparModel), vtail.J (if vSparModel)
    tailboom.Sy, tailboom.cave (if not hSparModel), tailboom.J (if hSparModel)

    LaTex Strings
    -------------
    mfac        m_{\\mathrm{fac}}

    """
    def setup(self, N=2):
        exec parse_variables(Empennage.__doc__)

        self.htail = HorizontalTail()
        self.hSparModel = self.htail.sparModel
        self.htail.substitutions.update({self.htail.mfac: 1.1})
        lh = self.lh = self.htail.lh
        self.vtail = VerticalTail()
        self.vSparModel = self.vtail.sparModel
        self.vtail.substitutions.update({self.vtail.mfac: 1.1})
        lv = self.lv = self.vtail.lv
        self.tailboom = TailBoom(N=N)
        self.components = [self.htail, self.vtail, self.tailboom]
        l = self.l = self.tailboom.l

        constraints = [
            W/mfac >= sum(c.W for c in self.components),
            l >= lh, l >= lv,
            ]

        return self.components, constraints
