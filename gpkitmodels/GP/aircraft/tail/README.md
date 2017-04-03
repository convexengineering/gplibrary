# Aircraft Empennage 

The empennage adds both weight and drag to each aircraft.  The Empennage model in `empennage.py` creates 3 models: TailBoom, HorizontalTail, and VerticalTail. Each component has its own separate bending and aerodynamic models. 

## TailBoom 

The tail boom has a diameter $d$, root wall thickness $t_0$, root moment of inertia $I_0$, modulus $E$, and density $\rho_{\text{cfrp}} = 1.6$ [g/cm$^3$], and length $l_{\text{h}}$. 
The total mass and root bending inertia are imposed in the optimization model as 
\begin{align}
    m &\geq \pi \rho_{\text{cfrp}} t_0 d l_{\text{h}} \left( 1 - \frac{1}{2} k\right) \\
    I_0 &\leq \pi t_0 d^3/8
\end{align}

where the index $k=0$ corresponds to a uniform wall thickness and stiffness, and $k=1$ corresponds to a linear drop-off to zero.  For both the solar-electric and gas powered aircraft $k=0.8$ is assumed.  
When the tail boom is loaded at the endpoint $x=l_{\text{h}}$, by the horizontal tail lift $L_{\text{h}}$, the end deflection angle follows from standard beam analysis

\begin{align}
    \label{e:boomdefl}
    \theta &\geq \frac{L_{\text{h}} l_{\text{h}}^2}{EI_0} \frac{1+k}{2} \\
    L_{\text{h}} &= \frac{1}{2} C_{L_{\text{h}}} \rho V^2 S_{\text{h}}.
\end{align}

## HorizontalTail

The horizontal tail can be sized to satisfy a horizontal tail volume coefficient condition, $V_{\text{h}} = 0.45$,

\begin{equation}
    V_{\text{h}} = \frac{S_{\text{h}}l_{\text{h}}}{Sc}
\end{equation}

## VerticalTail

The vertical tail is sized to meet a conservative tail volume coefficient, $V_{\text{v}}= 0.04$,

\begin{equation}
    \label{e:vtv}
    V_{\text{v}} = \frac{S_{\text{v}}}{S} \frac{l_{\text{v}}}{b}
\end{equation}

where $l_{\text{v}}$ is the vertical tail moment arm, assumed to be equal to the horizontal tail moment arm, $l_{\text{v}} = l_{\text{h}}$.

## Both Tails

Both the horizontal and vertical tails are assumed to have a carbon fiber skin and solid foam interior where their respective densities are $\rho_{A_{\text{cfrp}}} = 0.049$ [g/cm$^2$], $\rho_{\text{foam}} = 1.5$ [lbf/ft$^3$]. 
The weight of the horizontal and vertical tails is

\begin{align}
    \label{e:htweight}
    W_{\text{h}}/m_{\text{fac}} &= \rho_{\text{foam}} \frac{S_{\text{h}}^2}{b_{\text{h}}} \bar{A} + g\rho_{A_{\text{cfrp}}} S_{\text{h}} \\
    \label{e:vtweight}
    W_{\text{v}}/m_{\text{fac}} &= \rho_{\text{foam}} \frac{S_{\text{v}}^2}{b_{\text{v}}} \bar{A} + g\rho_{A_{\text{cfrp}}} S_{\text{v}}
\end{align}

where $b_{\text{h}}$ and $b_{\text{v}}$ are the spans of the horizontal and vertical tails respectively and $\bar{A}$ is the cross sectional area of the NACA 0008 airfoil. The margin factor $m_{\text{fac}}=1.1$, is included to account for control surfaces, attachment joints, actuators, etc. 

## TailBoomAero

The drag of the tail boom is calculated using a turbulent flat plate model,

\begin{align}
    \label{e:boomdrag}
    D_{\text{boom}} &\geq \frac{1}{2} C_f \rho V^2 l_{\text{h}}\pi d \\
    C_f &\geq \frac{0.445}{Re_{\text{boom}}^{0.3}} \\
    Re_{\text{boom}} &= \frac{V\rho l_{\text{h}}}{\mu}
\end{align}

## TailAero

The drag of the horizontal and vertical tails is computed using a GP-compatible fit of XFOIL data for a range of Reynolds numbers and NACA airfoil thicknesses,

\begin{align}
    D_{\text{h}} &\geq \frac{1}{2} c_{d_{\text{h}}} \rho V^2 S_{\text{h}} \\
    D_{\text{v}} &\geq \frac{1}{2} c_{d_{\text{v}}} \rho V^2 S_{\text{v}} \\
    \label{e:taildrag}
    c_{d_{\text{(v,h)}}}^{70.5599} &\geq \num{7.42688d-90} \left( \frac{Re_{\text{(v,h)}}}{\num{1d3}} \right)^{-33.0637}(100 \tau_{\text{(v,h)}})^{18.0419}  \nonumber \\
                 & + \num{5.02826d-163}\left(\frac{Re_{\text{(v,h)}}}{\num{1d3}}\right)^{-18.7959} (100\tau_{\text{(v,h)}})^{53.1879} \\
                 &+ \num{4.22901d-77}\left(\frac{Re_{\text{(v,h)}}}{\num{1d3}}\right)^{-41.1704} (100\tau_{\text{(v,h)}})^{28.4609} \nonumber \\
    Re_{\text{v}} &= \frac{V\rho S_{\text{v}}/b_{\text{v}}}{\mu} \\
    Re_{\text{h}} &= \frac{V\rho S_{\text{h}}/b_{\text{h}}}{\mu} 
\end{align}

where the selected airfoil is the NACA 0008 for both the horizontal and vertical tails (i.e. $\tau_{\text{(v,h)}} = 0.08$). 
The XFOIL data was generated for a zero angle of attack, based upon steady level flight where neither surface is generating lift.  
A comparison of the XFOIL data and Equation~\eqref{e:taildrag} is shown in Figure~\ref{f:taildragpolar}.

![Thickness and Re fit for CL=0](taildragpolar.pdf)

