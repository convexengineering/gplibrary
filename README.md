# gplibrary

[![Build Status](https://acdl.mit.edu/csi/job/gpkit_commons_Push_Models/badge/icon)](https://acdl.mit.edu/csi/job/gpkit_commons_Push_Models/)

This repository includes a number of useful GP- and SP- compatible models, such as:

* **Simple Models with detailed explanations** (good for learning GPkit)
  * SimpleAC: a basic aircraft model 
  * [Economic Order Quantity](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/misc/Economic%20Order%20Quantity/eoq.pdf): tradeoff between setup and holding costs
  * [Cylindrical Beam Moment of Inertia](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/misc/Moment%20of%20Inertia%20(cylindrical%20beam)/moi.pdf): GP approximation of cylindrical beam MOI
  * [Net Present Value](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/misc/Net%20Present%20Value/npv.pdf): financial tradeoff between cash and equipment
  * [Raymer Weights](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/misc/Raymer%20Weights): rule-of-thumb weight relations for aircraft design
* **GP models**
  * Aircraft
    * Engine
    * Fuselage
    * Tail
    * Wing
    * Mission: a model that unifies the subsystems and the aircraft flight profile
  * Beam: discretized beam for distributed loads
* **SP models**
  * Aircraft
    * Tail
    * Wing
