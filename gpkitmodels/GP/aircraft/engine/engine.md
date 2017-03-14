# IC Gas Engine Model

The models in the this folder aim to capture both engine weight and performance for small IC engines.  The DF70 from RCV engines is also modeled with the correct weight and specified BSFC, RPM and Power values and ranges. The general IC engine model is `gas_engine.py` and the DF70 specific model is `df70.py`.

## Weight

The weight is modeled by using a power law derived from data given by the [University of North Dakota](http://media.aero.und.edu/uasresearch.org/documents/195-197 Reference-Section Engines.pdf).  Using the fitting techniques in [gpfit](https://github.com/hoburg/gpfit) an equation was derived mapping weight to maximum power and is shown below:

![Power to weight law for IC engines](powervsweightfig.pdf)

# Performance

One common way predict engine performance is to assumed a constant $BSFC$ value, or that power is directly proportional to fuel burn.  This model accounts for changes in $BSFC$ due to throttle setting.  The lower the throttle setting the higher the $BSFC$.  Using the DF70 engine, a $BSFC$ to power setting mapping is found using gpfit. 

![Power to BSFC mapping](powertobsfcfit.pdf)
