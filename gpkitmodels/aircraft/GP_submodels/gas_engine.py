" engine_model.py "
from gpkit import Model, Variable, units

class Engine(Model):
    "engine model"
    def setup(self, DF70=False):

        self.DF70 = DF70

        W = Variable("W", "lbf", "Installed/Total engine weight")
        mfac = Variable("m_{fac}", 1.0, "-", "Engine weight margin factor")

        if DF70:
            Wdf70 = Variable("W_{DF70}", 7.1, "lbf",
                             "Installed/Total DF70 engine weight")
            Pslmax = Variable("P_{sl-max}", 5.17, "hp",
                              "Max shaft power at sea level")
            h = Variable("h", 12, "in", "engine height")

            constraints = [W/mfac >= Wdf70,
                           Pslmax == Pslmax]

        else:
            Pref = Variable("P_{ref}", 10.0, "hp", "Reference shaft power")
            Wengref = Variable("W_{eng-ref}", 10.0, "lbf",
                               "Reference engine weight")
            Weng = Variable("W_{eng}", "lbf", "engine weight")
            Pslmax = Variable("P_{sl-max}", "hp",
                              "Max shaft power at sea level")

            constraints = [
                Weng/Wengref >= 1.27847*(Pslmax/Pref)**0.772392,
                W/mfac >= 2.572*Weng**0.922*units("lbf")**0.078]

        return constraints

    def flight_model(self, state):
        return EnginePerf(self, state)

class EnginePerf(Model):
    "engine performance model"
    def setup(self, static, state):

        Pshaft = Variable("P_{shaft}", "hp", "Shaft power")
        bsfc = Variable("BSFC", "lb/hr/hp", "Brake specific fuel consumption")
        Pavn = Variable("P_{avn}", 40, "watts", "Avionics power")
        Ptotal = Variable("P_{total}", "hp", "Total power, avionics included")
        eta_alternator = Variable("\\eta_{alternator}", 0.8, "-",
                                  "alternator efficiency")
        href = Variable("h_{ref}", 1000, "ft", "reference altitude")
        # lfac = [1 - 0.906**(1/0.15)*(v.value/hr.value)**0.92
        #         for v, hr in zip(state["h"], state["h_{ref}"])]
        lfac = [-0.035*(v.value/hr.value) + 1.0
                for v, hr in zip(state["h"], href)]
        Leng = Variable("L_{eng}", lfac, "-", "shaft power loss factor")
        Pshaftmax = Variable("P_{shaft-max}",
                             "hp", "Max shaft power at altitude")
        mfac = Variable("m_{fac}", 1.0, "-", "BSFC margin factor")

        if static.DF70:
            rpm = Variable("RPM", "rpm", "Engine operating RPM")
            rpm_max = Variable("RPM_{max}", 7698, "rpm", "Maximum RPM")
            bsfc_min = Variable("BSFC_{min}", 0.3162, "kg/kW/hr",
                                "Minimum BSFC")

            constraints = [
                (bsfc/mfac/bsfc_min)**36.2209 >= (
                    2.31541*(rpm/rpm_max)**8.06517
                    + 0.00103364*(rpm/rpm_max)**-38.8545),
                (Ptotal/Pshaftmax)**0.1 == 0.999495*(rpm/rpm_max)**0.294421,
                rpm <= rpm_max
                ]
        else:
            bsfc_min = Variable("BSFC_{min}", 0.316, "kg/kW/hr",
                                "Minimum BSFC")

            constraints = [
                (bsfc/mfac/bsfc_min)**18.5563 >= (
                    0.00866321 * (Ptotal/Pshaftmax)**-7.70161
                    + 1.38628 * (Ptotal/Pshaftmax)**1.12922),
                ]

        constraints.extend([Pshaftmax/static["P_{sl-max}"] == Leng,
                            Pshaftmax >= Ptotal,
                            Ptotal >= Pshaft + Pavn/eta_alternator
                           ])

        return constraints
