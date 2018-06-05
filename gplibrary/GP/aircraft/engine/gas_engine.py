" engine_model.py "
from gpkit import Model, Variable, units
import os
import pandas as pd
# from gplibrary.tools.fit_constraintset import FitCS
from gpfit.fit_constraintset import FitCS


class Engine(Model):
    "engine model"
    def setup(self, DF70=False):

        self.DF70 = DF70

        W = Variable("W", "lbf", "Installed/Total engine weight")
        mfac = Variable("m_{fac}", 1.0, "-", "Engine weight margin factor")
        bsfc_min = Variable("BSFC_{min}", 0.3162, "kg/kW/hr", "minimum BSFC")
        Pref = Variable("P_{ref}", 10.0, "hp", "Reference shaft power")
        Wengref = Variable("W_{eng-ref}", 10.0, "lbf",
                           "Reference engine weight")
        Weng = Variable("W_{eng}", "lbf", "engine weight")
        Pslmax = Variable("P_{sl-max}", "hp",
                          "Max shaft power at sea level")

        path = os.path.dirname(__file__)
        df = pd.read_csv(path + os.sep + "power_lawfit.csv").to_dict(
            orient="records")[0]

        constraints = [
            FitCS(df, Weng/Wengref, [Pslmax/Pref]),
            W/mfac >= 2.572*Weng**0.922*units("lbf")**0.078]

        return constraints

    def flight_model(self, state):
        return EnginePerf(self, state)


class EnginePerf(Model):
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

        path = os.path.dirname(__file__)
        df = pd.read_csv(path + os.sep + "powerBSFCfit.csv").to_dict(
            orient="records")[0]

        constraints = [
            FitCS(df, bsfc/mfac/static["BSFC_{min}"], [Ptotal/Pshaftmax]),
            Pshaftmax/static["P_{sl-max}"] == Leng,
            Pshaftmax >= Ptotal,
            Ptotal >= Pshaft + Pavn/eta_alternator
        ]

        return constraints
