# gplibrary

[![Build Status](https://acdl.mit.edu/csi/view/convexengineering/job/gplibrary_Push_Models/badge/icon)](https://acdl.mit.edu/csi/view/convexengineering/job/gplibrary_Push_Models/)

This repository contains those GP-/SP-compatible models that we consider well documented and general enough to be useful to multiple projects.

* **Simple models with in-depth explanations** (good for learning GPkit)
  * [SimPleAC](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/SP/SimPleAC/simpleac.pdf): a basic aircraft model that captures the fundamental design tradeoffs
  * [Economic Order Quantity](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/misc/Economic%20Order%20Quantity/eoq.pdf): tradeoff between setup and holding costs
  * [Cylindrical Beam Moment of Inertia](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/misc/Moment%20of%20Inertia%20(cylindrical%20beam)/moi.pdf): GP approximation of cylindrical beam MOI
  * [Net Present Value](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/misc/Net%20Present%20Value/npv.pdf): financial tradeoff between cash and equipment
  * [Raymer Weights](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/misc/Raymer%20Weights): rule-of-thumb weight relations for aircraft design
* **GP models**
  * Aircraft
    * [Wing Structural and Aero Models](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/GP/aircraft/wing)
    * [Empennage](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/GP/aircraft/tail): TailBoom, HorizontalTail, and VerticalTail inherit from the Wing model
    * [Mission](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/GP/aircraft/mission): models that unify subsystems and flight profiles
    * [Fuselage](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/GP/aircraft/fuselage): elliptical and cylindrical fuselage models
    * [IC Gas Engine Model](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/GP/aircraft/engine)
  * [Bending Beam](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/GP/beam): discretized beam for distributed loads
* **SP models**
  * Aircraft
    * [Tail Boom Flexibility](https://github.com/convexengineering/gplibrary/tree/master/gpkitmodels/SP/aircraft/tail/tail_boom_flex.py)
    * [Wing Spanwise Effectiveness](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/SP/aircraft/wing/wing.py)
  * Atmosphere
    * [Tony Tao's fits as (efficient) signomial equalities](https://github.com/convexengineering/gplibrary/blob/master/gpkitmodels/SP/atmosphere/atmosphere.py). Valid until 10,000m of altitude. 

