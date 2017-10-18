SimPleAC
========

In this example, we modify the simple wing model proposed by Prof. Warren Hoburg in his thesis to create a simple aircraft model incorporating signomial programming (SP). 

```python
from gpkit import Model, Variable, SignomialsEnabled, VarKey, units
import numpy as np
import matplotlib.pyplot as plt

class SimPleAC(Model):
    def setup(self):
        # Variable definitions...
```


Weight and Lift model
-----------------

The simple aircraft borrows a majority of its aerodynamic model from Prof. Warren Hoburg's thesis, with some minor adjustments to make the model a SP. The aircraft is assumed to be in steady-level flight so that thrust equals drag, and the lift equals the weight. 

The total weight of the aircraft is the sum of the payload weight $W_0$, the wing weight $W_w$, and the fuel weight $W_f$.

\begin{equation}
    W \geq W_0 + W_w + W_f
\end{equation}

We use the naive $W \leq L$ model below for steady state flight:

\begin{equation}
    W_0 + W_w + 0.5 W_f \leq \frac{1}{2} \rho S C_L V^2
\end{equation}

where the lift of the aircraft is equal to weight of the aircraft with half-fuel. We would also like the fully-fueled aircraft to be able to fly at a minimum speed of $V_{min}$ without stalling , so we add the following constraint:

\begin{equation}
    W \leq \frac{1}{2} \rho V_{min}^2 S C_{L_{max}}
\end{equation}

The time of flight of the aircraft, which is a useful metric of performance, is simply the range over velocity:

\begin{equation}
    T_{flight} \geq \frac{Range}{V}
\end{equation}

The lift-to-drag ratio is also defined: 

\begin{equation}
    L/D = \frac{C_L}{C_D}    
\end{equation}

```python
        # Weight and lift model
        constraints += [W >= W_0 + W_w + W_f,
                W_0 + W_w + 0.5 * W_f <= 
                        0.5 * rho * S * C_L * V ** 2,
                W <= 0.5 * rho * S * C_Lmax * V_min ** 2,
                T_flight >= Range / V,
                LoD == C_L/C_D]
```

Thrust and drag model
----------

We assume a constant thrust specific fuel consumption (TSFC) for the 'engine' of the aircraft, which is assumed to provide as much thrust as needed. Since $T \geq D$:

\begin{equation}
    W_f \geq TSFC \times T_{flight} \times D
\end{equation}

the fuel weight required is the product of the TSFC, time of flight, and the total drag on the aircraft. The drag of the aircraft is the product of dynamic pressure ($\frac{1}{2} \rho V^2$), the planform area $S$, and the coefficient of drag of the aircraft:

\begin{equation} 
    D \geq \frac{1}{2}  \rho V^2 S  C_D
\label{e:d}
\end{equation}

The drag coefficient of the aircraft is assumed to be the sum of the fuselage drag, the wing profile drag, and the wing induced drag coefficients: 

\begin{equation}
    C_D \geq C_{D_{fuse}} + C_{D_{wpar}} + C_{D_{ind}}
\label{e:cd}
\end{equation}

The fuselage drag is a function of its drag area $CDA_0$ and the planform area of the wing:

\begin{equation}
    C_{D_{fuse}} = \frac{CDA_0}{S}
\label{e:cdfuse}
\end{equation}

where the $CDA_0$ is linearly proportional to the volume of fuel in the fuselage:

\begin{equation}
    V_{f_{fuse}} \leq CDA_0 \times 10 \mathrm{[meters]}
\label{e:vffuse}
\end{equation}

Note that we correct the dimensionality of the volume here, since GPkit automatically checks units. 

The wing profile drag is the product of the form factor, the friction drag coefficient, and the wetted area ratio of the wing, 

\begin{equation}
    C_{D_{wpar}} = k C_f S_{wetratio}
\label{e:cdwpar}
\end{equation}

The Reynolds number of the aircraft wing is approximated

\begin{equation}
    Re \leq \frac{\rho}{\mu} V \sqrt{\frac{S}{AR}}
\label{e:re}
\end{equation}

and used to find the friction drag coefficient of wing. We approximate the $C_f$ by assuming a turbulent Blasius flow over a flat plate as below:

\begin{equation}
    C_f \geq \frac{0.074} {Re^{0.2}}
\end{equation}

The induced drag of the wing is calculated with a span efficiency factor e, and is a function of the $C_L$ and aspect ratio $AR$ of the wing.

\begin{equation}
    C_{D_{induced}} = \frac{C_L^2}{\pi AR e}
\label{e:cdinduced}
\end{equation}

```python
        # Thrust and drag model
        C_D_fuse = CDA0 / S
        C_D_wpar = k * C_f * S_wetratio
        C_D_ind  = C_L ** 2 / (np.pi * A * e)
        constraints += [W_f >= TSFC * T_flight * D,
                    D >= 0.5 * rho * S * C_D * V ** 2,
                    C_D >= C_D_fuse + C_D_wpar + C_D_ind,
                    V_f_fuse <= 10*units('m')*CDA0,
                    Re <= (rho / mu) * V * (S / A) ** 0.5,
                    C_f >= 0.074 / Re ** 0.2]
```

Fuel volume model
-----------------

The fuel volume model is the main difference between simpleWing and SimPleAC, and introduces the only signomial constraints in the model. Firstly we define the required fuel volume using fuel density $\rho_f$. 

\begin{equation}
    V_f = \frac{W_f } {\rho_f g}
\label{e:vf}
\end{equation}

We consider wing fuel tanks and fuselage fuel tanks. The fuselage fuel was defined in the aerodynamic model, where the fuel volume contributes to drag. 

\begin{equation}
    V_{f_{wing}}^2 \leq 0.0009 \frac{S \tau^2}{AR}
\label{e:vfwing}
\end{equation}

The signomial constraint is introduced here, where the available fuel volume is upper-bounded by the sum of the fuselage and wing fuel volumes. 

\begin{equation}
    V_{f_{avail}} \leq V_{f_{wing}} + V_{f_{fuse}}
\label{vfavail}
\end{equation}

We constrain the total fuel volume to be less than the available fuel volume. 

\begin{equation}
    V_{f_{avail}} \geq V_{f}
\label{e:vfineq}
\end{equation}

```python
        # Fuel volume model 
        with SignomialsEnabled():
            constraints +=[V_f == W_f / g / rho_f,
                V_f_wing**2 <= 0.0009*S**3/A*tau**2, 
                    # linear with b and tau, quadratic with chord
                V_f_avail <= V_f_wing + V_f_fuse, #[SP]
                V_f_avail >= V_f]
```

Wing Weight Build-Up
---------------

The wing surface weight is a function of the planform area of the wing. 

\begin{equation}
W_{w_{surf}} \geq W_{w_{coeff2}} S
\label{e:wwsurf}
\end{equation}

The wing structural weight is a complex posynomial expression that takes into account the root bending moment and shear relief due to presence of fuel in the wings. 

\begin{equation}
W_{w_{strc}}^2 \geq \frac{W_{w_{coeff1}}^2}{\tau^2} (N_{ult}^2 AR ^ 3 ((W_0+\rho_fgV_{f_{fuse}}) W S))
\label{e:wwstrc}
\end{equation}

The total wing weight is lower-bounded. 

\begin{equation}
W_w \geq W_{w_{surf}} + W_{w_{strc}}
\label{e:ww}
\end{equation}

```python
        # Wing weight model
        constraints += [W_w_surf >= W_W_coeff2 * S,
            W_w_strc**2. >= W_W_coeff1**2./ tau**2. *
            (N_ult**2. * A ** 3. * ((W_0+V_f_fuse*g*rho_f) * W * S)),
            W_w >= W_w_surf + W_w_strc]

        return constraints
```

Valid objective functions
---------

We have tested a variety of potential objectives for the SimpPleAC model, some of which are as follows:
\begin{itemize}
    \item $W_f$: Fuel weight, the default objective in the ```main``` method.
    \item $W$: Total aircraft weight. Like fuel weight, but also adding extra cost for airframe weight.
    \item $D$: Drag. 
    \item $\frac{W_f}{T_{flight}}$: Product of the fuel weight and the inverse of the time of flight. 
    \item $W_{f} + c \times T_{flight}$: A linear combination of fuel weight and time of flight. This can simulate recurring costs (fuel and labor), and yield interesting results. 
\end{itemize}
