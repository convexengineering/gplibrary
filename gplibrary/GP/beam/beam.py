" discretized beam model "
from gpkit import Model, Variable, Vectorize

#pylint: disable=invalid-name

class Beam(Model):
    """discretized beam bending model

    Upper Unbounded
    ---------------
    EIbar, dbar_tip

    Lower Unbounded
    ---------------
    dx, qbar (if not qbarFun)

    """
    qbarFun = None
    SbarFun = None
    MbarFun = None

    def setup(self, N):

        with Vectorize(N-1):
            EIbar = self.EIbar = Variable("\\bar{EI}", "-",
                             "normalized YM and moment of inertia")
            dx = self.dx = Variable("dx", "-", "normalized length of element")

        with Vectorize(N):
            Sbar = Variable("\\bar{S}", self.SbarFun, "-", "normalized shear")
            Mbar = Variable("\\bar{M}", self.MbarFun, "-", "normalized moment")
            th = Variable("\\theta", "-", "deflection slope")
            dbar = Variable("\\bar{\\delta}", "-", "normalized displacement")
            self.dbar_tip = dbar[-1]


        throot = Variable("\\theta_{root}", 1e-10, "-", "Base angle")
        dbarroot = Variable("\\bar{\\delta}_{root}", 1e-10, "-",
                            "Base deflection")

        constraints = []
        if self.SbarFun is None:
            with Vectorize(N):
                qbar = self.qbar = Variable("\\bar{q}", self.qbarFun, "-",
                                            "normalized loading")
            Sbartip = Variable("\\bar{S}_{tip}", 1e-10, "-", "Tip loading")
            constraints.extend([
                Sbar[:-1] >= Sbar[1:] + 0.5*dx*(qbar[:-1] + qbar[1:]),
                Sbar[-1] >= Sbartip])

        if self.MbarFun is None:
            Mbartip = Variable("\\bar{M}_{tip}", 1e-10, "-", "Tip moment")
            constraints.extend([
                Mbar[:-1] >= Mbar[1:] + 0.5*dx*(Sbar[:-1] + Sbar[1:]),
                Mbar[-1] >= Mbartip])

        constraints.extend([
            th[0] >= throot,
            th[1:] >= th[:-1] + 0.5*dx*(Mbar[1:] + Mbar[:-1])/EIbar,
            dbar[0] >= dbarroot,
            dbar[1:] >= dbar[:-1] + 0.5*dx*(th[1:] + th[:-1]),
            ])

        return constraints
