" test tail models "
from gpkitmodels.GP.aircraft.tail.horizontal_tail import HorizontalTail
from gpkitmodels.GP.aircraft.tail.vertical_tail import VerticalTail
from gpkitmodels.GP.aircraft.tail.empennage import Empennage
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState
from gpkit import Model, Variable, units

def test_htail():

    Sw = Variable("S_w", 50, "ft**2", "wing area")
    cmac = Variable("cmac", 15, "in", "wing MAC")
    ht = HorizontalTail()
    fs = FlightState()
    ht.substitutions.update({ht.topvar("W"): 5, "m_h": 0.01, "AR": 4,
                             "V_h": 0.5, "l_h": 10})
    perf = ht.flight_model(ht, fs)

    m = Model(perf["C_d"], [ht["V_h"] <= ht["S"]*ht["l_h"]/Sw/cmac, ht, perf])
    m.solve(verbosity=0)

def test_vtail():

    Sw = Variable("S_w", 50, "ft**2", "wing area")
    bw = Variable("b_w", 20, "ft", "wing span")
    vt = VerticalTail()
    fs = FlightState()
    vt.substitutions.update({vt.topvar("W"): 5, "AR": 3, "V_v": 0.04,
                             "l_v": 10})
    perf = vt.flight_model(vt, fs)

    m = Model(perf["C_d"], [vt["V_v"] <= vt["S"]*vt["l_v"]/Sw/bw, vt, perf])
    m.solve(verbosity=0)

def test_emp():

    Sw = Variable("S_w", 50, "ft**2", "wing area")
    bw = Variable("b_w", 20, "ft", "wing span")
    cmac = Variable("cmac", 15, "in", "wing MAC")
    emp = Empennage()
    fs = FlightState()
    emp.substitutions.update({emp.topvar("W"): 10, "l": 5,
                              emp.htail.planform["AR"]: 4,
                              emp.vtail.planform["AR"]: 4,
                              emp.vtail["V_v"]: 0.04,
                              emp.htail["V_h"]: 0.4,
                              emp.htail["m_h"]: 0.01})
    htperf = emp.htail.flight_model(emp.htail, fs)
    vtperf = emp.vtail.flight_model(emp.vtail, fs)
    tbperf = emp.tailboom.flight_model(fs)

    m = Model(htperf["C_d"] + vtperf["C_d"] + tbperf["C_f"],
              [emp["l_v"] == emp["l"], emp["l_h"] == emp["l"],
               emp["V_h"] <= emp.htail["S"]*emp["l_h"]/Sw/cmac,
               emp["V_v"] <= emp.vtail["S"]*emp["l_v"]/Sw/bw,
               emp, fs, htperf, vtperf, tbperf])
    m.solve(verbosity=0)

def test():
    test_htail()
    test_vtail()
    test_emp()

if __name__ == "__main__":
    test()
