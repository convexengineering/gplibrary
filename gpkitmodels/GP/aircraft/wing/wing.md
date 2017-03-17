# Wing Structural and Aero Models

The wing GP model is mostly based off of Mark Drela's wing bending notes from MIT's [OpenCourseWare](https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl10.pdf) site.  Basic assumptions are a constant tapered wing, no sweep, bending loads are taken by the spar, torsional loads are taken by the skin.

The usage is as follows: `Wing` in `wing.py` is created either by itself or in an aircraft model.  It has two functions that return models to be created separately from `Wing`, `WingAero` and `WingLoading` (also in `wing.py`). `WingAero` is an erodynamic model of the JHO airfoil.  $C_d$ is to be used in an overall aircraft aerodyanmic model. `WingLoading` creates loading models that govern the bending loads of the spar and the torsional loads of the skin.  Inside `Wing` 2 objects are created `CapSpar` or `TubeSpar` and `WingSkin` with the option to create `WingInterior` if `hollow=True`.  Their weights are summed to give an overall weight of the wing.  Each submodel is described in detail. 

## Wing

Constrains basic shape of wing using 

$$ b^2 = SAR$$
$$\bar{c}(y) \equiv \frac{c(y)}{S/b} = \frac{2}{1+\lambda} \left( 1 + (\lambda - 1) \frac{2y}{b} \right) $$

## WingAero

Calculates the wing drag from induced drag and wing profile drag.  Wing profile drag assumes the JHO airfoil and is calculated using a posynomial fit to XFOIL data at various reynolds numbers. 

![JHO drag polars](jho1polarfit.py)

## CapSpar

Set up constraints for size of spar. Spar thickness is assumed to be no greater than the max thickness of the airfoil.  Width is no greater than %10 of the chord.  Moment of inerita is calculated using first order approximation. 

## ChordSparL

Basic assumption is that the distributed load of lift and weight is proportional to the local chord. 

$$ q(y) \approx K_q c(y) $$
$$ K_q = \frac{N_{\text{max}}W_{\text{cent}}}{S} $$

This model takes as an input `CapSpar` and `Wcent`, or the center weight of the aircaft.  Uses discretized beam theory to predict moment and deflection.  Stress and overall wing deflection are constrainted. 

## GustSparL

The main assumption is the the loading is the same as the `ChordSparL` case but there is an additional gust term that varies as $1-\cos$ along the span. So the distributed load that is used as an input to the discretized beam model is

$$\bar{q}(y) = \frac{q(y)b}{W_{\text{cent}}N_{\text{max}}} &\geq \bar{c}(y) \left[1 + \frac{c_{l_{\alpha}}}{C_L} \alpha_{\text{gust}} (y) \left(1 + \frac{W_{\text{wing}}}{W_{\text{cent}}} \right) \right] $$

where $\alpha_{\text{gust}}$ is

$$ \alpha_{\text{gust}}(y)  = \tan^{-1}\left(\frac{V_{\text{gust}}(y)}{V}
 \right). $$

 and $V_{\text{gust}}$ is:

 $$V_{\text{gust}}(y) = V_{\text{ref}} \left(1-\cos\left(\frac{2y}{b} \frac{\pi}{2} \right) \right) $$

 $\alpha_{\text{gust}}$ is approximated by a monomial fit on the range of gust velocities $V/V_{\text{gust}} \in [0, 0.7]$.

![gust loading diagram](gustloaddiagram.pdf)


## WingSkin

The wing skin calculates the weight of the wing based off of the surface aera of the wing

$$ W_{\text{skin}} \geq 2 \rho_{\text{cfrp}}} t S g. $$

## WingSkinL

The main assumptions here are that the wing has a constant thickness skin, the torsional loads are taken by the skin, and the highest torison occurs at the root.  This model only has one constraint

$$ \tau_{CFRP} \geq \frac{C_{m_w}S\rhoV_{NE}^2}{\bar{J/t} c_{root}^2 t} $$

The $\bar{J/t}$ is taken from the JHO airfoil.

## WingInterior

If this model is enabled, the wing is assummed to be filled with foam.  The cross sectional aera of the wing is assumed to be the JHO airfoil. 
