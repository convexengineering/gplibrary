# Gas HALE model
```python

from numpy import pi
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
import gpkit
import numpy as np
gpkit.settings['latex_modelname'] = False

class GasPoweredHALE(Model):
    def setup(self):
        constraints = []
```


## Flight segment definitions
```python

        # define number of segments
        NSeg = 3
        #Note: NSeg has to be an odd number
        # defining indices of different flight segments
        NLoiter = (NSeg-1)/2
        if NSeg == 3:
            NCruise = [0,2]
        elif NSeg == 7:
            Nclimb = [0,2,4,6]
            NCruise = [1,5]
```


# Fuel weight model
```python

        MTOW = Variable('MTOW', 'lbf', 'max take off weight')
        W_end = VectorVariable(NSeg, 'W_{end}', 'lbf', 'segment-end weight')
        W_fuel = VectorVariable(NSeg, 'W_{fuel}', 'lbf', 'segment-fuel weight')
        W_zfw = Variable('W_{zfw}', 'lbf', 'Zero fuel weight')
        W_pay = Variable('W_{pay}',10,'lbf', 'Payload weight')
        W_avionics = Variable('W_{avionics}', 2, 'lbf', 'Avionics weight')
        f_airframe = Variable('f_{airframe}', 0.25, '-', 'airframe weight fraction')
        W_airframe = Variable('W_{airframe}', 'lbf', 'airframe weight')
        W_begin = W_end.left # define beginning of segment weight
        W_begin[0] = MTOW
```

## Assumptions

*  Mass take off weight is greater than the end of the first segment weight plus that segment fuel weight.
* The end of each flight segment weight must be greater than the next end of flight segment weight plus the fuel weight of the next flight segment.
* The end of the last flight segment weight must be greater than the zero fuel weight
```python

        constraints.extend([MTOW >= W_end[0] + W_fuel[0],
                            W_end[:-1] >= W_end[1:] + W_fuel[1:],
                            W_end[-1] >= W_zfw,
                            W_airframe >= f_airframe*MTOW])
```


# Steady level flight model
## Assumptions
* Steady state flight
* Angle of attack zero
* Thrust equals drag
* Lift equals weight
* The weight is approximated by the average weight of the flight segment.
* The average weight can be approximated by: $W_{avg} = \sqrt{W_{begin} W_{end}}$.  Where $W_{begin}$ is the beginning of the flight segment and $W_{end}$ is the end of the flight segment.

```python

        CD = VectorVariable(NSeg, 'C_D', '-', 'Drag coefficient')
    	CL = VectorVariable(NSeg, 'C_L', '-', 'Lift coefficient')
        V = VectorVariable(NSeg, 'V', 'm/s','cruise speed')
        rho = VectorVariable(NSeg, r'\rho', 'kg/m^3', 'air density')
        S = Variable('S', 16, 'ft^2', 'wing area')
        eta_prop = VectorVariable(NSeg, r'\eta_{prop}', [0.8,0.6,0.8], '-',
                                  'propulsive efficiency')
        P_shaft = VectorVariable(NSeg, 'P_{shaft}', 'hp', 'Shaft power')
```

# Climb model
```python

        h_dot = Variable(NSeg, 'h_{dot}', [200,0,0], 'ft/min', 'Climb rate')
        constraints.extend([P_shaft >= V*(W_end+W_begin)/2*CD/CL/eta_prop + W_begin*h_dot/eta_prop,
                            0.5*rho*CL*S*V**2 >= (W_end+W_begin)/2])
```

# Engine Model
```python

        W_eng = Variable('W_{eng}', 'lbf', 'Engine weight')
        W_engtot = Variable('W_{eng-tot}', 'lbf', 'Installed engine weight')
        W_engref = Variable('W_{eng-ref}', 4.4107, 'lbf', 'Reference engine weight')
        P_shaftref = Variable('P_{shaft-ref}', 2.295, 'hp', 'reference shaft power')
```

## Engine Weight constraints
```python

        constraints.extend([W_eng/W_engref >= 0.5538*(P_shaft/P_shaftref)**1.075,
                            W_engtot >= 2.572*W_eng**0.922*units('lbf')**0.078])
```

# Weight breakdown
```python

        constraints.extend([W_airframe >= f_airframe*MTOW,
                            W_zfw >= W_pay + W_avionics + W_airframe + W_engtot])
```

# Breguet Range
## Assumptions
* Constant speed during each flight section
* Constant BSFC
* The $\ln$ can be approximated using a Taylor-series expansion
```python

        z_bre = VectorVariable(NSeg, 'z_{bre}', '-', 'breguet coefficient')
        BSFC = VectorVariable(NSeg,'BSFC', [0.5,.55,0.6], 'lbf/hr/hp', 'brake specific fuel consumption')
        t = VectorVariable(NSeg, 't', 'days', 'time on station')
        R = Variable('R', 200, 'nautical_miles', 'range to station')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        constraints.extend([z_bre >= V*t*BSFC*CD/CL/eta_prop,
                            R == V[NCruise]*t[NCruise],
                            t[NLoiter] == 5*units('days'),
                            W_fuel/W_end >= te_exp_minus1(z_bre, 3)])
```


# Aerodynamics model
## Assumptions
* The wing is a box shape.
* The non-wing drag is a constant
* The stall factor is based off standard airfoil polar.
* Reference length for Reynolds number is teh chord.
* The skin friction is based off of Blasius flat plate.
* The form factor for the wing is constant.
```python


        Cd0 = Variable('C_{d0}', 0.02, '-', 'Non-wing drag coefficient')
        CLmax = Variable('C_{L-max}', 1.5, '-', 'Maximum lift coefficient')
        e = Variable('e', 0.9, '-', 'Spanwise efficiency')
        AR = Variable('AR', '-', 'Aspect ratio')
        b = Variable('b', 'ft', 'Span')
        mu = Variable(r'\mu', 1.5e-5, 'N*s/m^2', 'Dynamic viscosity')
        Re = VectorVariable(NSeg, 'Re', '-', 'Reynolds number')
        Cf = VectorVariable(NSeg, 'C_f', '-', 'wing skin friction coefficient')
        Kwing = Variable('K_{wing}', 1.3, '-', 'wing form factor')
        cl_16 = Variable('cl_{16}', 0.0001, '-', 'profile stall coefficient')

        constraints.extend([CD >= Cd0 + 2*Cf*Kwing + CL**2/(pi*e*AR) + cl_16*CL**16,
                            b**2 == S*AR,
```

In place of an actual structural model, we impose $AR \leq 20$.
```python

                            AR <= 20,
                            CL <= CLmax,
                            Re == rho*V/mu*(S/AR)**0.5,
                            Cf >= 0.074/Re**0.2])
```

# Atmosphere model
## Assumptions
* Valid only to the top of the troposphere.

## References
* [wp:Density of Air](http://en.wikipedia.org/wiki/Density_of_air#Altitude)
```python


        h = VectorVariable(NSeg, 'h', 'ft', 'Altitude')
        gamma = Variable(r'\gamma',1.4,'-', 'Heat capacity ratio of air')
        p_sl = Variable('p_{sl}', 101325, 'Pa', 'Pressure at sea level')
        T_sl = VectorVariable(NSeg, 'T_{sl}', [288.15,288.15,288.15], 'K',
                              'Temperature at sea level')
        L_atm = Variable('L_{atm}', 0.0065, 'K/m', 'Temperature lapse rate')
        T_atm = VectorVariable(NSeg, 'T_{atm}', 'K', 'Air temperature')
        a_atm = VectorVariable(NSeg, 'a_{atm}','m/s', 'Speed of sound at altitude')
        R_spec = Variable('R_{spec}', 287.058,'J/kg/K', 'Specific gas constant of air')
        TH = (g/R_spec/L_atm).value.magnitude  # dimensionless

        constraints.extend([#h <= [20000, 20000, 20000]*units.m,
                            T_sl >= T_atm + L_atm*h,  # Temp decreases w/ altitude
                            rho == p_sl*T_atm**(TH-1)/R_spec/(T_sl**TH),
                            h[NLoiter] >= 15000*units('ft'), # makes sure that the loiter occurs above minimum h
                            ])
```


# Wind speed model
```python

        V_wind = VectorVariable(NSeg, 'V_{wind}', 'm/s', 'wind speed')
        wd_cnst = Variable('wd_{cnst}', 0.0015, 'm/s/ft',
                           'wind speed constant predicted by model')
                            #0.002 is worst case, 0.0015 is mean at 45d
        wd_ln = Variable('wd_{ln}', 8.845, 'm/s',
                         'linear wind speed variable')
                        #13.009 is worst case, 8.845 is mean at 45deg
        h_min = Variable('h_{min}', 11800, 'ft', 'minimum height')
        h_max = Variable('h_{max}', 20866, 'ft', 'maximum height')

        constraints.extend([V_wind >= wd_cnst*h + wd_ln,
                            V >= V_wind,
                            h[NCruise] >= h_min])
```


# Conclusion
```python

        objective = MTOW
        return objective, constraints

if __name__ == "__main__":
    M = GasPoweredHALE()
    M.solve()
    with open("sol.tex", "w") as f:
        f.write(M.solution.table(latex=True))
```


# Solution
\input{sol.tex}
