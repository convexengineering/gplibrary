from gpkitmodels.GP.aircraft.wing.wing import Wing
from gpkitmodels.GP.aircraft.wing.wvl.wvlsp import WVL
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState
from gpkit import Model

W = Wing()
W.flight_model = WVL
W.substitutions[W.W] = 50
W.substitutions[W.planform.tau] = 0.115
fs = FlightState()
perf = W.flight_model(W, fs)
loading = [W.spar.loading(W, fs)]
# loading[0].substitutions["W"] = 100

m = Model(perf.Cd, [W, fs, perf, loading, loading[0].W == W.W])
m.localsolve(verbosity=0)
