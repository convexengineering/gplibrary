" wing test "
from gpkitmodels.GP.aircraft.wing.wing import Wing
from gpkit import Variable, Model

class FlightState(Model):
    " state variables "
    def setup(self):

        V = Variable("V", 50, "m/s", "airspeed")
        rho = Variable("\\rho", 1.255, "kg/m^3", "air density")
        mu = Variable("\\mu", 1.5e-5, "N*s/m**2", "air viscosity")

        constraints = [V == V, rho == rho, mu == mu]

        return constraints

def wing_test():
    " test wing models "

    Wcent = Variable("W_{cent}", 100, "lbf", "center weight")

    W = Wing()
    W.substitutions[W.topvar("W")] = 50
    fs = FlightState()
    perf = W.flight_model(W, fs)
    loading = W.loading(W, Wcent, W.topvar("W"), fs["V"], perf["C_L"])

    m = Model(perf["C_d"], [W, fs, perf, loading])
    m.solve("mosek")

if __name__ == "__main__":
    wing_test()
