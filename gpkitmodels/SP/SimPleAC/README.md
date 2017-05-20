Introduction
-------------------------

In this example, we modify the simple wing model proposed by Prof. Warren Hoburg in his thesis to create a simple aircraft model compatible with signomial programming (SP). 

```python

from gpkit import Variable, Model, SignomialsEnabled, VarKey, units
import numpy as np


class SimPleAC():
 
```


Simple Aircraft Models 
--------------------------

Equations:


\begin{equation}
\tag{1}
    C_{D_{fuse}} = \frac{CDA_0}{S}
\label{e:cdfuse}
\end{equation}



