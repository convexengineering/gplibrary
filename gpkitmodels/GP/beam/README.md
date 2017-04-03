# Beam

This model is valid for a discritized beam and calculates non-dimensional sheer stresses, moments, angles and deflections for a given distributed load.  Assumes segment lengths are equal along the beam. 

Required inputs are a non-dimensional distributed load that is a vector variable and the number of nodes

\begin{table}[]
\centering
\caption{My caption}
\label{my-label}
\begin{tabular}{lll}
Input & Description                    & Type           \\
qbar  & Distributed load with N values & VectorVariable \\
N     & Number of nodes                & int           
\end{tabular}
\end{table}

The non-dimensional equation derivation is shown below for wing bending application.

Using a standard Bernoulli-Euler discretized beam model with $n$ nodes, the shear forces and moments can be computed from the distributed loads $q(y)$, with boundary conditions of zero shear forces and moments at the wing tips.

\begin{align}
    \label{e:shear}
    \mathcal{S}_i &= \mathcal{S}_{i+1} - \frac{q_{i+1} + q_i}{2}\Delta y \\
    \label{e:moment}
    M_i &= M_{i+1} - \frac{\mathcal{S}_{i+1} + \mathcal{S}_i}{2}\Delta y \\
    \label{e:shearboundary}
    \mathcal{S}_n &= 0 \\
    \label{e:momentboundary}
    M_n &= 0
\end{align}

Similarly, the angle deflection and deflection can be calculated with boundary conditions of zero angle and deflection and the wing root.

\begin{align}
    \label{e:angle}
    \Theta_{i} &= \Theta_{i+1} + \frac{1}{2} \left(\frac{M_i}{EI_i} + \frac{M_{i-1}}{EI_{i-1}} \right) \Delta y \\
    \label{e:deflection}
    w_{i} &= w_{i+1} + \frac{1}{2} (\Theta_i + \Theta_{i-1}) \Delta y \\
    \label{e:angleboundary}
    \Theta_0 &= 0 \\
    \label{e:defboundary}
    w_0 &= 0 
\end{align}
 
Equations~\eqref{e:shear}-\eqref{e:defboundary} are GP-compatible if expressed as

\begin{align}
    \label{e:sheargp}
    \mathcal{S}_{i+1} &\geq \mathcal{S}_i + \frac{q_{i+1} + q_i}{2} \Delta y \\
    \label{e:momentgp}
    \mathcal{M}_{i+1} &\geq \mathcal{M}_i + \frac{\mathcal{S}_{i+1} + \mathcal{S}_i}{2} \Delta y \\
    \label{e:anglegp}
    \Theta_{i} &\geq \Theta_{i+1} + \frac{1}{2} \left(\frac{\mathcal{M}_i}{EI_i} + \frac{\mathcal{M}_{i-1}}{EI_{i-1}} \right) \Delta y \\
    \label{e:deflection}
    w_{i} &\geq w_{i+1} + \frac{1}{2} (\Theta_i + \Theta_{i-1}) \Delta y 
\end{align}

