# GAS-POWERED, MEDIUM-ALTITUDE, LONG-ENDURANCE AIRCRAFT

This paper presents the designs achieved in the 82 project. This also presents all of the models used in the design of this aircraft.

## Aerodynamic Model
\input{Aerodynamics.vars.generated.tex}
\input{Aerodynamics.cnstrs.generated.tex}

## Atmospheric Model
\input{Atmosphere.vars.generated.tex}
\input{Atmosphere.cnstrs.generated.tex}

## Breguet Endurance Model
\input{BreguetEndurance.vars.generated.tex}
\input{BreguetEndurance.cnstrs.generated.tex}

## Breguet Range Model
\input{BreguetRange.vars.generated.tex}
\input{BreguetRange.cnstrs.generated.tex}

## Climb Model
\input{Climb.vars.generated.tex}
\input{Climb.cnstrs.generated.tex}

## Cruise Model
\input{Cruise.vars.generated.tex}
\input{Cruise.cnstrs.generated.tex}

## Engine Model
\input{Engine.vars.generated.tex}
\input{Engine.cnstrs.generated.tex}

## Fuel Model
\input{Fuel.vars.generated.tex}
\input{Fuel.cnstrs.generated.tex}

## Fuselage Model
\input{Fuselage.vars.generated.tex}
\input{Fuselage.cnstrs.generated.tex}

## Loiter Model
\input{Loiter.vars.generated.tex}
\input{Loiter.cnstrs.generated.tex}

## Mission Profile Model
\input{Mission.vars.generated.tex}
\input{Mission.cnstrs.generated.tex}

## Steady Level Flight Model
\input{SteadyLevelFlight.vars.generated.tex}
\input{SteadyLevelFlight.cnstrs.generated.tex}

## Structural Model
\input{Structures.vars.generated.tex}
\input{Structures.cnstrs.generated.tex}

## Weight Model
\input{Weight.vars.generated.tex}
\input{Weight.cnstrs.generated.tex}

## Wind Model
\input{Wind.vars.generated.tex}
\input{Wind.cnstrs.generated.tex}

## Overall Model
\input{GasMALE.vars.generated.tex}
\input{GasMALE.cnstrs.generated.tex}

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
#inPDF: replace with tstation_vs_MTOW_rubber.fig.generated.tex
from plotting import plot_sweep, fix_vars, plot_altitude_sweeps
import numpy as np

M = GasMALE()
M.substitutions.update({"MTOW": 150})
fig, ax = plot_sweep(M, "MTOW", np.linspace(70, 500, 15), "t_{loiter}")
gen_tex_fig(fig, "tstation_vs_MTOW_rubber")
```

### CDR Aircraft Sizing

After deciding on the 150 lb aircraft to meet with a 1 day margin on the loiter requirement, a DF70 engine was chosen. The aircraft was reoptimized to meet the 6 day time on station and minimize max take off weight.  This was the aircraft chosen for the CDR. The solution is shown in the tables below. 

```python
#inPDF: replace with sol.generated.tex
M = GasMALE(DF70=True)
M.substitutions.update({"t_{loiter}": 6})
M.cost = M["MTOW"]
sol = M.solve("mosek")

with open("sol.generated.tex", "w") as f:
    f.write(sol.table(latex=True))
```

### DF70 Engine Model
The engine model of the DF70 is shown below. 

\input{DF70.vars.generated.tex}
\input{DF70.cnstrs.generated.tex}

```python
#inPDF: skip
models, modelnames = find_submodels([M], [])
DF70Engine = models[modelnames.index("Engine") + 1]
gen_model_tex(DF70Engine, "Engine", texname="DF70")
```

By fixing the following variables to their respective values we were also able to generate performance curves. 

```python
#inPDF: replace with fixvars.table.generated.tex


vars_to_fix = {"S":0.0, "b":0.0, "Vol_{fuse}":0.00001}
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

# payload power vs time on station
fig, ax = plot_sweep(M, "P_{pay}", np.linspace(10, 200, 15),
                     "t_{loiter}")
gen_tex_fig(fig, "t_vs_Ppay")

# payload weight vs time on station
fig, ax = plot_sweep(M, "W_{pay}", np.linspace(5, 40, 15), "t_{loiter}")
gen_tex_fig(fig, "t_vs_Wpay")

# wind speed vs time on station
M = GasMALE(wind=True, DF70=True)
fix_vars(M, sol, vars_to_fix)
fig, ax = plot_sweep(M, "V_{wind}_Wind, Loiter, Mission, GasMALE", np.linspace(5, 40, 15), "t_{loiter}")
gen_tex_fig(fig, "t_vs_Vwind")

# altitude vs time on loiter
altitude_vars = {"t_{loiter}"}
figs, axs = plot_altitude_sweeps(np.linspace(14000, 20000, 20), altitude_vars,
                     vars_to_fix)

for var, fig in zip(altitude_vars, figs):
    gen_tex_fig(fig, "%s_vs_altitude" % var.replace("{", "").replace("}", "").replace("_", ""))
```
\input{t_vs_Ppay.fig.generated.tex}
\input{t_vs_Wpay.fig.generated.tex}
\input{t_vs_Vwind.fig.generated.tex}
\input{tloiter_vs_altitude.fig.generated.tex}







