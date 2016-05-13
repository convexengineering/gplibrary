This is a simple template for setting up automatic report generation with GPkit.
To generate this file, run
```
make model.pdf
```
in this directory.

# Example Model

This example model minimizes $x + 1/x$ subject to the constraint
$$x \ge x_{min}.$$
The description above is hard coded as opposed to being derived from the Python
code. This is a simple and effective way to document model equations.

The full GP, derived from the Python implementation of the model, is
```python
#inPDF: replace with gp.generated.tex
from simple_gp import SimpleGP
m = SimpleGP()
with open("gp.generated.tex", "w") as f:
    f.write("$%s$" % m.latex(excluded=["models", "units"]))
```

Sweeping over $x_{min}$ shows the relationship between $x$ and $x_{min}$.
```python
#inPDF: skip
import numpy as np
import matplotlib.pyplot as plt
xmin = np.linspace(0.01, 3, 30)
m.substitutions.update({"x_{min}": ("sweep", xmin)})
sol = m.solve(verbosity=0)
plt.plot(sol["variables"]["x_{min}"], sol["cost"])
plt.ylim((0, plt.ylim()[1]))
plt.xlabel("$x_{min}$")
plt.ylabel("cost")
plt.savefig("simple_sweep.pdf")
```
```python
#inPDF: skip
# TODO: argh, the below doesn't work
# \begin{figure}[h!]
# \label{fig:sweep}
# \begin{center}
# \includegraphics[scale=0.5]{simple_sweep.pdf}
# \caption{xmin only has an effect when greater than 1}
# \end{center}
# \end{figure}
```
![xmin only has an effect when greater than 1](simple_sweep.pdf)
