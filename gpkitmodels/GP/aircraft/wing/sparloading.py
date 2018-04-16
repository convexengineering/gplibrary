" spar loading "
from gpkit import Model, parse_variables
from gpkitmodels.GP.beam.beam import Beam
from gpkitmodels.GP.aircraft.tail.tail_boom import TailBoomState
from numpy import pi

#pylint: disable=no-member, unused-argument, exec-used, invalid-name
#pylint: disable=undefined-variable, attribute-defined-outside-init

class SparLoading(Model):
    """ Spar Loading Model

    Variables
    ---------
    Nmax            5              [-]     max loading
    Nsafety         1.0            [-]     safety load factor
    kappa           0.2            [-]     max tip deflection ratio
    W                              [lbf]   loading weight
    twmax           30.*pi/180     [-]     max tip twist
    Stip            1e-10          [N]     tip loading
    Mtip            1e-10          [N*m]   tip moment
    throot          1e-10          [-]     root deflection angle
    wroot           1e-10          [m]     root deflection

    Variables of length wing.N
    --------------------------
    q     self.new_qbarFun                   [-]   distributed wing loading
    S                        [-]     shear along wing
    Mbar                        [-]   wing section root moment
    th                          [-]     deflection angle
    w                           [-]     wing deflection

    Variables of length wing.N-1
    ----------------------------
    Mtw                     [N*m]   local moment due to twisting
    M                        [N*m]   wing section root moment
    theta                   [-]     twist deflection
    EIbar                   [-]    EIbar

    Upper Unbounded
    ---------------
    I, Sy, J (if wingSparJ)
    theta (if not wingSparJ), M (if not wingSparJ)

    Lower Unbounded
    ---------------
    b, W, cave
    wing.planform.CM (if wingSparJ), qne (if wingSparJ)
    theta (if not wingSparJ), M (if not wingSparJ)

    LaTex Strings
    -------------
    Nmax                N_{\\mathrm{max}}
    kappa               \\kappa
    Mr                  M_r

    """
    def new_qbarFun(self, c):
        " define qbar model for chord loading "
        barc = self.wing.planform.cbar
        return [f(c) for f in self.wing.substitutions[barc]]

    new_SbarFun = None

    def setup(self, wing, state):
        self.wing = wing
        exec parse_variables(SparLoading.__doc__)

        b = self.b = self.wing.planform.b
        I = self.I = self.wing.spar.I
        Sy = self.Sy = self.wing.spar.Sy
        cave = self.cave = self.wing.planform.cave
        cbar = self.cbar = self.wing.planform.cbar
        E = self.wing.spar.material.E
        sigma = self.wing.spar.material.sigma
        deta = self.wing.planform.deta

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            # N == Nsafety*Nmax,
            # q >= cbar,
            # S[:-1] >= S[1:] + 0.5*deta*(b/2.)*(q[:-1] + q[1:]),
            S[:-1] >= S[1:] + 0.5*deta*(q[:-1] + q[1:]),
            S[-1] >= 1e-10,
            Mbar[:-1] >= Mbar[1:] + 0.5*deta*(S[:-1] + S[1:]),
            Mbar[-1] >= 1e-10,
            M >= Mbar[:-1]*Nsafety*Nmax*W*b/4,
            th[0] >= throot,
            th[1:] >= th[:1] + 0.5*deta*(Mbar[1:] + Mbar[:-1])/EIbar,
            EIbar <= 8*E*I/Nsafety/Nmax/W/b**2,
            w[0] >= 1e-10,
            w[1:] >= w[:-1] + 0.5*deta*(th[1:] + th[:-1]),
            sigma >= M/Sy,
            w[-1] <= kappa,
            ]

        self.wingSparJ = hasattr(self.wing.spar, "J")

        if self.wingSparJ:
            qne = self.qne = state.qne
            J = self.J = self.wing.spar.J
            G = self.wing.spar.shearMaterial.G
            cm = self.wing.planform.CM
            constraints.extend([
                Mtw >= cm*cave**2*qne*deta*b/2*Nsafety,
                theta[0] >= Mtw[0]/G/J[0]*deta[0]*b/2,
                theta[1:] >= theta[:-1] + Mtw[1:]/G/J[1:]*deta[1:]*b/2,
                twmax >= theta[-1]
                ])
        return constraints
