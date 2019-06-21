import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, SignomialEquality, VarKey, units, Vectorize
from SimPleAC_mission import Mission, SimPleAC
from gpkitmodels.SP.atmosphere.atmosphere import Atmosphere

# SimPleAC with multimission design (updated 5/31/2019, by Berk Ozturk)

class Multimission(Model):
    def setup(self,aircraft,Nmissions,Nsegments):
        self.aircraft = aircraft
        self.missions = []
        for i in range(0,Nmissions):
            self.missions.append(Mission(self.aircraft,Nsegments))

        # Multimission objective variables
        W_f_mm = Variable('W_{f_{mm}}','N','multimission fuel weight')

        with Vectorize(Nmissions):
            # Mission variables
            hcruise    = Variable('h_{cruise_{mm}}', 'm', 'minimum cruise altitude')
            Range      = Variable("Range_{mm}", "km", "aircraft range")
            W_p        = Variable("W_{p_{mm}}", "N", "payload weight", pr=20.)
            rho_p      = Variable("\\rho_{p_{mm}}", 1500, "kg/m^3", "payload density", pr=10.)
            V_min      = Variable("V_{min_{mm}}", 25, "m/s", "takeoff speed", pr=20.)
            cost_index = Variable("C_{mm}", '1/hr','hourly cost index')
            TOfac      = Variable('T/O factor_{mm}', 2.,'-','takeoff thrust factor')

        constraints = []

        # Setting up the missions
        for i in range(0,Nmissions):
            constraints += [
            self.missions[i]['h_{cruise_m}'] == hcruise[i],
            self.missions[i]['Range_m']      == Range[i],
            self.missions[i]['W_{p_m}']      == W_p[i],
            self.missions[i]['\\rho_{p_m}']  == rho_p[i],
            self.missions[i]['V_{min_m}']    == V_min[i],
            self.missions[i]['C_m']          == cost_index[i],
            self.missions[i]['T/O factor_m'] == TOfac[i],
            # Upper bounding relevant variables
            W_f_mm <= 1e11*units('N'),
            ]

        # Multimission constraints
        constraints += [W_f_mm >= sum(self.missions[i]['W_{f_m}'] for i in range(0,Nmissions))]

        return constraints, self.aircraft, self.missions

def test():
    Nmissions = 2
    Nsegments = 4
    aircraft = SimPleAC()
    m = Multimission(aircraft,Nmissions,Nsegments)
    m.substitutions.update({
        'h_{cruise_{mm}}':[5000*units('m'), 5000*units('m')],
        'Range_{mm}'     :[3000*units('km'), 2000*units('km')],
        'W_{p_{mm}}'     :[6250*units('N'),   8000*units('N')],
        '\\rho_{p_{mm}}' :[1500*units('kg/m^3'), 2000*units('kg/m^3')],
        'C_{mm}'         :[120*units('1/hr'), 360*units('1/hr')],
    })

    m.cost = (m.missions[0]['W_{f_m}']*units('1/N') + m.missions[1]['C_m']*m.missions[1]['t_m'])
    sol = m.localsolve(verbosity = 0)

if __name__ == "__main__":
    Nmissions = 2
    Nsegments = 4
    aircraft = SimPleAC()
    m = Multimission(aircraft,Nmissions,Nsegments)
    m.substitutions.update({
        'h_{cruise_{mm}}':[5000*units('m'), 5000*units('m')],
        'Range_{mm}'     :[3000*units('km'), 2000*units('km')],
        'W_{p_{mm}}'     :[6250*units('N'),   8000*units('N')],
        '\\rho_{p_{mm}}' :[1500*units('kg/m^3'),   2000*units('kg/m^3')],
        'C_{mm}'         :[120*units('1/hr'), 360*units('1/hr')],
    })

    m.cost = (m.missions[0]['W_{f_m}']*units('1/N') + m.missions[1]['C_m']*m.missions[1]['t_m'])
    sol = m.localsolve(verbosity = 2)

