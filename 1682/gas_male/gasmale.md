# GAS-POWERED, MEDIUM-ALTITUDE, LONG-ENDURANCE AIRCRAFT

This paper presents the designs achieved in the 82 project. This also presents all of the models used in the design of this aircraft.

## Aerodynamic Model
\input{tex/Aerodynamics.vars.generated.tex}
\input{tex/Aerodynamics.cnstrs.generated.tex}

## Atmospheric Model
\input{tex/Atmosphere.vars.generated.tex}
\input{tex/Atmosphere.cnstrs.generated.tex}

# Avionics Model
\input{tex/Avionics.vars.generated.tex}
\input{tex/Avionics.cnstrs.generated.tex}

# Beam Model
\input{tex/Beam.vars.generated.tex}
\input{tex/Beam.cnstrs.generated.tex}

## Breguet Endurance Model
\input{tex/BreguetEndurance.vars.generated.tex}
\input{tex/BreguetEndurance.cnstrs.generated.tex}

## Cap Spar Model
\input{tex/Spar.vars.generated.tex}
\input{tex/Spar.cnstrs.generated.tex}

## Climb Model
\input{tex/Climb.vars.generated.tex}
\input{tex/Climb.cnstrs.generated.tex}

## Cruise Model
\input{tex/Cruise.vars.generated.tex}
\input{tex/Cruise.cnstrs.generated.tex}

## Empennage Model
\input{tex/Empennage.vars.generated.tex}
\input{tex/Empennage.cnstrs.generated.tex}

## Engine Performance Model
\input{tex/EnginePerformance.vars.generated.tex}
\input{tex/EnginePerformance.cnstrs.generated.tex}

## Engine Weight Model
\input{tex/EngineWeight.vars.generated.tex}
\input{tex/EngineWeight.cnstrs.generated.tex}

## Fuel Model
\input{tex/Fuel.vars.generated.tex}
\input{tex/Fuel.cnstrs.generated.tex}

## Fuselage Model
\input{tex/Fuselage.vars.generated.tex}
\input{tex/Fuselage.cnstrs.generated.tex}

## Horizontal Tail Model
\input{tex/HorizontalTail.vars.generated.tex}
\input{tex/HorizontalTail.cnstrs.generated.tex}

## Loiter Model
\input{tex/Loiter.vars.generated.tex}
\input{tex/Loiter.cnstrs.generated.tex}

## Mission Profile Model
\input{tex/Mission.vars.generated.tex}
\input{tex/Mission.cnstrs.generated.tex}

## Steady Level Flight Model
\input{tex/SteadyLevelFlight.vars.generated.tex}
\input{tex/SteadyLevelFlight.cnstrs.generated.tex}

## Take Off
\input{tex/TakeOff.vars.generated.tex}
\input{tex/TakeOff.cnstrs.generated.tex}

## Wing Model
\input{tex/Wing.vars.generated.tex}
\input{tex/Wing.cnstrs.generated.tex}

## Weight Model
\input{tex/Weight.vars.generated.tex}
\input{tex/Weight.cnstrs.generated.tex}

# Tail Boom Model
\input{tex/TailBoom.vars.generated.tex}
\input{tex/TailBoom.cnstrs.generated.tex}

## Vertical Tail Model
\input{tex/VerticalTail.vars.generated.tex}
\input{tex/VerticalTail.cnstrs.generated.tex}

## Wind Model
\input{tex/Wind.vars.generated.tex}
\input{tex/Wind.cnstrs.generated.tex}

## Wing Model
\input{tex/Wing.vars.generated.tex}
\input{tex/Wing.cnstrs.generated.tex}

## Overall Model
\input{tex/GasMALE.vars.generated.tex}
\input{tex/GasMALE.cnstrs.generated.tex}

```python
#inPDF: skip
from gasmale import GasMALE
from gen_tex import gen_model_tex, find_submodels, gen_tex_fig, gen_fixvars_tex 

M = GasMALE()
models, modelnames = find_submodels([M], [])
for m in models: 
    gen_model_tex(m, m.__class__.__name__)

```

# Sizing

This model was created and then a sweep was done to determine the MTOW required to meet 5 days. 

```python
#inPDF: replace with tex/tstation_vs_MTOW_rubber.fig.generated.tex
from plotting import plot_sweep, fix_vars, plot_altitude_sweeps
import numpy as np

M = GasMALE()
M.substitutions.update({"MTOW": 150})
fig, ax = plot_sweep(M, "MTOW", np.linspace(70, 500, 15), ["t_{loiter}"])
gen_tex_fig(fig, "tstation_vs_MTOW_rubber")
```

### CDR Aircraft Sizing

After deciding on the 150 lb aircraft to meet with a 1 day margin on the loiter requirement, a DF70 engine was chosen. The aircraft was reoptimized to meet the 6 day time on station and minimize max take off weight.  This was the aircraft chosen for the CDR. The solution is shown in the tables below. 

```python
#inPDF: replace with tex/sol.generated.tex
M = GasMALE(DF70=True)
M.substitutions.update({"t_{loiter}": 6})
M.cost = M["MTOW"]
sol = M.solve("mosek")

with open("tex/sol.generated.tex", "w") as f:
    f.write(sol.table(latex=True))
```

### DF70 Engine Model
The engine model of the DF70 is shown below. 

\input{tex/DF70Performance.vars.generated.tex}
\input{tex/DF70Performance.cnstrs.generated.tex}
\input{tex/DF70Weight.vars.generated.tex}
\input{tex/DF70Weight.cnstrs.generated.tex}

```python
#inPDF: skip
models, modelnames = find_submodels([M], [])
DF70EngineP = models[modelnames.index("EnginePerformance") + 1]
DF70EngineW = models[modelnames.index("EngineWeight") + 1]
gen_model_tex(DF70EngineP, "EnginePerformance", texname="DF70Performance")
gen_model_tex(DF70EngineW, "EngineWeight", texname="DF70Weight")
```

By fixing the following variables to their respective values we were also able to generate performance curves. 

```python
#inPDF: replace with tex/fixvars.table.generated.tex

vars_to_fix = {"S":0.0, "b": 0.0, "d":0.0, "l_{fuel}": 0.0, "W_Wing, GasMALE":0.0}
gen_fixvars_tex(M, sol, vars_to_fix)

fix_vars(M, sol, vars_to_fix)
sol = M.solve("mosek") # check for solving errors
```

## Sweeps

```python
#inPDF: skip
# set objective to time on station after fixing variables
del M.substitutions["t_{loiter}"]
M.cost = 1/M["t_{loiter}"]

M.substitutions.update({"P_{pay}": 100})
# payload power vs time on station
fig, ax = plot_sweep(M, "P_{pay}", np.linspace(10, 200, 15), ["t_{loiter}"], ylim=[0,10])
gen_tex_fig(fig, "t_vs_Ppay")

# payload weight vs time on station
fig, ax = plot_sweep(M, "W_{pay}", np.linspace(5, 20, 15), ["t_{loiter}"], ylim=[0,10])
gen_tex_fig(fig, "t_vs_Wpay")

# wind speed vs time on station
M = GasMALE(wind=True, DF70=True)
fix_vars(M, sol, vars_to_fix)
M.substitutions.update({"P_{pay}": 100})
fig, ax = plot_sweep(M, "V_{wind}_Wind, Loiter, Mission, GasMALE", np.linspace(5, 40, 15), ["t_{loiter}"], ylim=[0,10])
gen_tex_fig(fig, "t_vs_Vwind")

# altitude vs time on loiter
# altitude_vars = {"t_{loiter}"}
# figs, axs = plot_altitude_sweeps(np.linspace(14000, 20000, 10), altitude_vars,
#                      vars_to_fix)
# axs[0].set_ylim([0,10])
# 
# for var, fig in zip(altitude_vars, figs):
#     gen_tex_fig(fig, "%s_vs_altitude" % var.replace("{", "").replace("}", "").replace("_", ""))
```
\input{tex/t_vs_Ppay.fig.generated.tex}
\input{tex/t_vs_Wpay.fig.generated.tex}
\input{tex/t_vs_Vwind.fig.generated.tex}

## Flight Profile

By further discretizing the climb, cruise, and loiter mission segments the following figures were generated to follow the performance over the duration of the mission. 

```python
#inPDF: skip
from plotting import plot_mission_var

Mprof = GasMALE(DF70=True)
Mprof.substitutions.update({"t_{loiter}": 6})
Mprof.cost = Mprof["MTOW"]
sol = Mprof.solve("mosek")
fix_vars(Mprof, sol, vars_to_fix)
Mprof.substitutions.update({"P_{pay}": 100})
del Mprof.substitutions["t_{loiter}"]
Mprof.cost = 1/Mprof["t_{loiter}"]
sol = Mprof.solve("mosek")

# plot mission profiles
fig, ax = plot_mission_var(Mprof, sol, "V", [0, 40])
gen_tex_fig(fig, "profile_velocity")

fig, ax = plot_mission_var(Mprof, sol, "h", [0, 20000])
gen_tex_fig(fig, "profile_altitude")

fig, ax = plot_mission_var(Mprof, sol, "\\eta_{prop}", [0, 1])
gen_tex_fig(fig, "profile_etaprop")

fig, ax = plot_mission_var(Mprof, sol, "BSFC", [0, 2])
gen_tex_fig(fig, "profile_BSFC")

fig, ax = plot_mission_var(Mprof, sol, "P_{shaft-max}", [0, 5])
gen_tex_fig(fig, "profile_Pshaftmax")

fig, ax = plot_mission_var(Mprof, sol, "P_{shaft-tot}", [0, 5])
gen_tex_fig(fig, "profile_Pshafttot")

fig, ax = plot_mission_var(Mprof, sol, "RPM", [0, 9000])
gen_tex_fig(fig, "profile_RPM")

fig, ax = plot_mission_var(Mprof, sol, "W_{N+1}", [0, 150], "aircraft weight [lbf]")
gen_tex_fig(fig, "profile_weight")
```
\input{tex/profile_velocity.fig.generated.tex}
\input{tex/profile_etaprop.fig.generated.tex}
\input{tex/profile_BSFC.fig.generated.tex}
\input{tex/profile_Pshaftmax.fig.generated.tex}
\input{tex/profile_Pshafttot.fig.generated.tex}
\input{tex/profile_RPM.fig.generated.tex}
\input{tex/profile_weight.fig.generated.tex}
