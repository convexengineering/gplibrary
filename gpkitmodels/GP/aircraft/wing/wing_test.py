" wing test "
from gpkitmodels.GP.aircraft.wing.wing import Wing
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin
from gpkitmodels.GP.aircraft.wing.wing_core import WingCore
from gpkitmodels.GP.aircraft.wing.boxspar import BoxSpar
from gpkit import Model, parse_variables
from gpkitmodels.GP.materials.cfrp import CFRPFabric
from gpkitmodels.GP.materials.foam import Foam

#pylint: disable=no-member, exec-used

class FlightState(Model):
    """ Flight State

    Variables
    ---------
    V           50      [m/s]        airspeed
    rho         1.255   [kg/m^3]     air density
    mu          1.5e-5  [N*s/m**2]   air viscosity

    """
    def setup(self):
        exec parse_variables(FlightState.__doc__)

def wing_test():
    " test wing models "

    W = Wing()
    W.substitutions[W.W] = 50
    fs = FlightState()
    perf = W.flight_model(W, fs)
    loading = [W.spar.loading(W)]
    loading[0].substitutions["W"] = 100
    loading.append(W.spar.gustloading(W))
    loading[1].substitutions["W"] = 100

    from gpkit import settings
    if settings["default_solver"] == "cvxopt":
        for l in loading:
            for v in ["\\bar{M}_{tip}", "\\bar{S}_{tip}",
                      "\\bar{\\delta}_{root}", "\\theta_{root}"]:
                l.substitutions[v] = 1e-3

    m = Model(perf.Cd, [
        loading[1].v == fs.V,
        loading[1].cl == perf.CL,
        loading[1].Ww == W.W,
        loading[1].Ww <= 0.5*fs.rho*fs.V**2*perf.CL*W.planform.S,
        W, fs, perf, loading])
    m.solve(verbosity=0)

def box_spar():
    " test wing models "

    Wing.sparModel = BoxSpar
    W = Wing()
    W.substitutions[W.W] = 50
    fs = FlightState()
    perf = W.flight_model(W, fs)
    loading = [W.spar.loading(W)]
    loading[0].substitutions["W"] = 100
    loading.append(W.spar.gustloading(W))
    loading[1].substitutions["W"] = 100

    from gpkit import settings
    if settings["default_solver"] == "cvxopt":
        for l in loading:
            for v in ["\\bar{M}_{tip}", "\\bar{S}_{tip}",
                      "\\bar{\\delta}_{root}", "\\theta_{root}"]:
                l.substitutions[v] = 1e-3

    m = Model(perf.Cd, [
        loading[1].v == fs.V,
        loading[1].cl == perf.CL,
        loading[1].Ww == W.W,
        loading[1].Ww <= 0.5*fs.rho*fs.V**2*perf.CL*W.planform.S,
        W, fs, perf, loading])
    m.solve(verbosity=0)

def materials():
    " test wing models "

    cfrp = CFRPFabric()
    foam = Foam()
    WingCore.material = foam
    WingSkin.material = cfrp
    Wing.sparModel.shearMaterial = cfrp
    Wing.sparModel.coreMaterial = foam
    Wing.sparModel = BoxSpar
    W = Wing()
    W.substitutions[W.W] = 50
    fs = FlightState()
    perf = W.flight_model(W, fs)
    loading = [W.spar.loading(W)]
    loading[0].substitutions["W"] = 100
    loading.append(W.spar.gustloading(W))
    loading[1].substitutions["W"] = 100

    from gpkit import settings
    if settings["default_solver"] == "cvxopt":
        for l in loading:
            for v in ["\\bar{M}_{tip}", "\\bar{S}_{tip}",
                      "\\bar{\\delta}_{root}", "\\theta_{root}"]:
                l.substitutions[v] = 1e-3

    m = Model(perf.Cd, [
        loading[1].v == fs.V,
        loading[1].cl == perf.CL,
        loading[1].Ww == W.W,
        loading[1].Ww <= 0.5*fs.rho*fs.V**2*perf.CL*W.planform.S,
        W, fs, perf, loading])
    m.solve(verbosity=0)

def test():
    " tests "
    wing_test()
    box_spar()
    materials()

if __name__ == "__main__":
    test()
