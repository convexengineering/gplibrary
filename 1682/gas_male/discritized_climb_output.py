import numpy as np
import numpy as np
from gasmale import GasMALEDiscritizedClimb
from wpay import W_pay as Wpay
from plotting import fix_vars, plot_mission_var, plot_altitude_sweeps
from plotting import plot_sweep

fixed = False
# PLOTS #
fixedPLOTS = False
missionPLOTS = False
wind = False
W_pay = False
P_pay = False
altitude = False

M = GasMALEDiscritizedClimb()
sol_fix = M.solve('mosek')

vars_to_fix = {"S":0.0, "b":0.0, "Vol_{fuse}":0.0001}
fix_vars(M, sol_fix, vars_to_fix)

# plot mission profiles
fig, ax = plot_mission_var(M, "V", [0, 40])
fig.savefig("profile_t_vs_velocity.pdf")

fig, ax = plot_mission_var(M, "\\eta_{prop}", [0, 1])
fig.savefig("profile_t_vs_etaprop.pdf")

fig, ax = plot_mission_var(M, "BSFC", [0, 2])
fig.savefig("profile_t_vs_BSFC.pdf")

fig, ax = plot_mission_var(M, "P_{shaft-max}", [0, 5])
fig.savefig("profile_t_vs_Pshaftmax.pdf")

fig, ax = plot_mission_var(M, "P_{shaft-tot}", [0, 5])
fig.savefig("profile_t_vs_Pshafttot.pdf")

fig, ax = plot_mission_var(M, "RPM", [0, 9000])
fig.savefig("profile_t_vs_RPM.pdf")

fig, ax = plot_mission_var(M, "W_{end}", [0, 150], "aircraft weight [lbf]")
fig.savefig("profile_t_vs_weight.pdf")

# plot altitude sweeps
plot_altitude_sweeps(np.linspace(14000, 24000, 20), {"t_{station}"},
                     vars_to_fix, CLIMB=True)

# plot other sweeps
del M.substitutions["t_{station}"]
M.cost = 1/M["t_{station}"]

fig, ax = plot_sweep(M, "P_{pay}", np.linspace(5, 200, 20), "t_{station}")
fig.savefig("Ppay_vs_tstation_climb.pdf")

M_wind = GasMALEDiscritizedClimb(wind=True)
sol_windfix = M_wind.solve("mosek", verbosity=0)

fix_vars(M_wind, sol_windfix, vars_to_fix)

del M_wind.substitutions["t_{station}"]
M_wind.cost = 1/M_wind["t_{station}"]

fig, ax = plot_sweep(M_wind, "V_{wind}", np.linspace(5, 40, 20),
                     "t_{station}")
fig.savefig("Vwind_vs_tstation_climb.pdf")

# uses trim drag to calculate payload weight effects
Wpay(M)
