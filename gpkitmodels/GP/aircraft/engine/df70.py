" engine_model.py "
from builtins import zip
from gpkit import Model, Variable, units

class DF70(Model):
    "engine model"
    def setup(self):

        W = Variable("W", "lbf", "Installed/Total engine weight")
        mfac = Variable("m_{fac}", 1.0, "-", "Engine weight margin factor")
        bsfc_min = Variable("BSFC_{min}", 0.3162, "kg/kW/hr", "minimum BSFC")
        Wdf70 = Variable("W_{DF70}", 7.76, "lbf",
                         "Installed/Total DF70 engine weight")
        Pslmax = Variable("P_{sl-max}", 5.17, "hp",
                          "Max shaft power at sea level")
        h = Variable("h", 12, "in", "engine height")

        constraints = [W/mfac >= Wdf70,
                       Pslmax == Pslmax]

        return constraints

    def flight_model(self, state):
        return DF70Perf(self, state)

class DF70Perf(Model):
    "engine performance model"
    def setup(self, static, state):

        Pshaft = Variable("P_{shaft}", "hp", "Shaft power")
        bsfc = Variable("BSFC", "kg/kW/hr", "Brake specific fuel consumption")
        Pavn = Variable("P_{avn}", 40, "watts", "Avionics power")
        Ptotal = Variable("P_{total}", "hp", "Total power, avionics included")
        eta_alternator = Variable("\\eta_{alternator}", 0.8, "-",
                                  "alternator efficiency")
        href = Variable("h_{ref}", 1000, "ft", "reference altitude")
        h_vals = state.substitutions("h")
        if len(href) == 1:
            h_vals = [h_vals]
        lfac = [-0.035*(v/hr.value) + 1.0
                for v, hr in zip(h_vals, href)]
        Leng = Variable("L_{eng}", lfac, "-", "shaft power loss factor")
        Pshaftmax = Variable("P_{shaft-max}",
                             "hp", "Max shaft power at altitude")
        mfac = Variable("m_{fac}", 1.0, "-", "BSFC margin factor")
        rpm = Variable("RPM", "rpm", "Engine operating RPM")
        rpm_max = Variable("RPM_{max}", 7698, "rpm", "Maximum RPM")

        constraints = [
            (bsfc/mfac/static["BSFC_{min}"])**36.2209 >= (
                2.31541*(rpm/rpm_max)**8.06517
                + 0.00103364*(rpm/rpm_max)**-38.8545),
            (Ptotal/Pshaftmax)**0.1 == 0.999495*(rpm/rpm_max)**0.294421,
            rpm <= rpm_max,
            Pshaftmax/static["P_{sl-max}"] == Leng,
            Pshaftmax >= Ptotal,
            Ptotal >= Pshaft + Pavn/eta_alternator
            ]

        return constraints
