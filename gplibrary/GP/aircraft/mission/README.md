# Mission

This folder contains various ''mission`` models that can be used for a given aircraft. 

## BreguetEndurance

This model predicts how much fuel is needed to fly for a certain time period. It needs as an input a performance model that has the following variables:

\begin{table}[]
\centering
\begin{tabular}{ll}
Variable    & Description                      \\
\hline
$P_{total}$ & Total power used by the aircraft \\
$BSFC$      & Break specific fuel consumption  \\
$W_{end}$   & End of flight segment weight     \\
$W_{start}$ & Start of flight segment weight  
\end{tabular}
\end{table}

The following is a derivation of the form of the Breguet Range equation that is used and a description of the Taylor expanison used to make this GP-compatible. 

\begin{equation}
    \label{e:breguetendurance}
    t = \frac{W_{\text{ave}}}{P_{\text{shaft}}\text{BSFC}g} \ln{\left( \frac{W_{\text{initial}}}{W_{\text{final}}}\right)}.
\end{equation}

The derivation begins with the differential form of Breguet Range,\cite{br2}

\begin{equation}
    \label{e:breguetdiff}
    -\frac{dW}{dt} = g\dot{m}_{\text{fuel}}.
\end{equation}

Using the definition of $\text{BSFC}$

\begin{equation}
    \label{e:brBSFC}
    \text{BSFC} = \frac{\dot{m}_{\text{fuel}}}{P_{\text{shaft}}},
\end{equation}

Equation~\eqref{e:breguetdiff} can be written as

\begin{equation}
    \label{e:brdiff2}
    -dW = g P_{\text{shaft}} \text{BSFC} dt.
\end{equation}

This version comes from assuming that $\text{BSFC}$ and the power to weight ratio, $(P_{\text{shaft}}/W)$, are constant during the considered flight segment. 
One way to obtain a constant power to weight ratio is a constant velocity and constant lift coefficient.\cite{br2}
$W_{\text{ave}}$ is assumed to be the geometric mean, defined as

Dividing by $W$,

\begin{equation}
    \label{e:brdiff2}
    -\frac{dW}{W} = \frac{g P_{\text{shaft}}\text{BSFC} }{W} dt,
\end{equation}

and integrating, the Breguet Range equation can be expressed as

\begin{equation}
    \label{e:be1}
    \ln{\left( \frac{W_{\text{initial}}}{W_{\text{final}}} \right)} = \frac{gP_{\text{shaft}}\text{BSFC}}{W_{\text{ave}}} t,
\end{equation}

where $W_{\text{ave}}$ is the average weight of the aircraft during the flight segment.  

\begin{equation}
    \label{e:gpmean}
    W_{\text{ave}} = \sqrt{W_{\text{initial}}W_{\text{final}}}.
\end{equation}

    To make Equation~\eqref{e:breguetendurance} GP compatible, a Taylor expansion is used,\cite{hoburgthesis}

\begin{align}
    \label{e:brzbre}
    z_{\text{bre}} &\geq \frac{P_{\text{shaft}}t \text{BSFC} g}{W}\\
    \label{e:brtaylor}
    \frac{W_{\text{fuel}}}{W_\text{final}} &\geq z_{\text{bre}} + \frac{z_{\text{bre}}^2}{2} + \frac{z_{\text{bre}}^3}{6} + \frac{z_{\text{bre}}^4}{24} + \dots
\end{align}

    Equations~\eqref{e:brzbre} and~\eqref{e:brtaylor} are monomial and posynomial respectively and therefore GP compatible. For long-endurance aircraft, missions can last days, causing the power to weight ratio $(P_{\text{shaft}}/W)$ to vary significantly during the course of the flight.  
    Equations~\eqref{e:brzbre},~\eqref{e:brtaylor}, and~\eqref{e:slfweight} can be discretized to account for this.

\begin{align}
    \label{e:slfweightd}
    \sqrt{W_i W_{i+1}} &= \frac{1}{2} \rho_i V_i^2 C_{L_i} S \\
    \label{e:brzbred}
    z_{bre_i} &\geq \frac{P_{\text{shaft}_i}t_i \text{BSFC} g}{\sqrt{W_i W_{i+1}}}\\
    \label{e:brtaylord}
    \frac{W_{\text{fuel}_i}}{W_{i+1}} &\geq z_{bre_i} + \frac{z_{bre_i}^2}{2} + \frac{z_{bre_i}^3}{6} + \frac{z_{bre_i}^3}{24} 
    \end{align}

For evaluation of long-endurance, gas-powered aircraft a discretization of $N=5$ was used. 
