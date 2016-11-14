" engine_model.py "
from gpkit import Model, Variable, units

class Engine(Model):
    "engine model"
    def __init__(self, DF70=False, **kwargs):

        self.DF70 = DF70

        W = Variable("W", "lbf", "Installed/Total engine weight")
        mfac = Variable("m_{fac}", 1.0, "-", "Engine weight margin factor")

        if DF70:
            Wdf70 = Variable("W_{DF70}", 7.1, "lbf",
                             "Installed/Total DF70 engine weight")
            Pslmax = Variable("P_{sl-max}", 5.17, "hp",
                              "Max shaft power at sea level")
            constraints = [W/mfac >= Wdf70,
                           Pslmax == Pslmax]

        else:
            Pref = Variable("P_{ref}", 2.295, "hp", "Reference shaft power")
            Wengref = Variable("W_{eng-ref}", 4.4107, "lbf",
                               "Reference engine weight")
            Weng = Variable("W_{eng}", "lbf", "engine weight")
            Pslmax = Variable("P_{sl-max}", "hp",
                              "Max shaft power at sea level")

            constraints = [
                Weng/Wengref >= 0.5538*(Pslmax/Pref)**1.075,
                W/mfac >= 2.572*Weng**0.922*units("lbf")**0.078]

        self.flight_model = EnginePerf

        Model.__init__(self, None, constraints, **kwargs)

class EnginePerf(Model):
    "engine performance model"
    def __init__(self, static, state, **kwargs):

        Pshaft = Variable("P_{shaft}", "hp", "Shaft power")
        bsfc = Variable("BSFC", "lb/hr/hp", "Brake specific fuel consumption")
        rpm = Variable("RPM", "rpm", "Engine operating RPM")
        Pavn = Variable("P_{avn}", 40, "watts", "Avionics power")
        Ppay = Variable("P_{pay}", 10, "watts", "Payload power")
        Ptotal = Variable("P_{total}", "hp", "Total power, avionics included")
        eta_alternator = Variable("\\eta_{alternator}", 0.8, "-",
                                  "alternator efficiency")
        lfac = [1 - 0.906**(1/0.15)*(v.value/hr.value)**0.92
                for v, hr in zip(state["h"], state["h_{ref}"])]
        Leng = Variable("L_{eng}", lfac, "-", "shaft power loss factor")
        Pshaftmax = Variable("P_{shaft-max}",
                             "hp", "Max shaft power at altitude")
        mfac = Variable("m_{fac}", 1.0, "-", "BSFC margin factor")

        if static.DF70:
            rpm_max = Variable("RPM_{max}", 7698, "rpm", "Maximum RPM")
            bsfc_min = Variable("BSFC_{min}", 0.3162, "kg/kW/hr",
                                "Minimum BSFC")

            constraints = [
                (bsfc/mfac/bsfc_min)**35.7 >=
                (2.29*(rpm/rpm_max)**8.02 + 0.00114*(rpm/rpm_max)**-38.3),
                (Ptotal/Pshaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.294
                ]
        else:
            bsfc_min = Variable("BSFC_{min}", 0.32, "kg/kW/hr",
                                "Minimum BSFC")
            rpm_max = Variable("RPM_{max}", 9000, "rpm", "Maximum RPM")

            constraints = [
                (bsfc/mfac/bsfc_min)**0.129 >=
                (0.972*(rpm/rpm_max)**-0.141 + 0.0268*(rpm/rpm_max)**9.62),
                (Ptotal/Pshaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.292
                ]

        constraints.extend([Pshaftmax/static["P_{sl-max}"] == Leng,
                            Pshaftmax >= Ptotal,
                            rpm <= rpm_max])

        if state.onStation:
            constraints.extend([
                Ptotal >= Pshaft + (Pavn + Ppay)/eta_alternator])
        else:
            constraints.extend([Ptotal >= Pshaft + Pavn/eta_alternator])


        Model.__init__(self, None, constraints, **kwargs)
