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

def test():
    " test wing models "
    from gpkit import units

    Wcent = Variable("W_{cent}", 100, "lbf", "center weight")

    W = Wing()
    W.substitutions[W.topvar("W")] = 50
    fs = FlightState()
    perf = W.flight_model(W, fs)
    loading = [W.spar.loading(W)]
    loading[0].substitutions["W"] = 100
    loading.append(W.spar.gustloading(W))
    loading[1].substitutions["W"] = 100

    m = Model(perf["C_d"], [loading[1]["V"] == fs["V"],
                            loading[1]["c_l"] == perf["C_L"],
                            loading[1]["W_w"] == W.topvar("W"),
                            W, fs, perf, loading])
    m.solve(verbosity=0)

if __name__ == "__main__":
    test()
