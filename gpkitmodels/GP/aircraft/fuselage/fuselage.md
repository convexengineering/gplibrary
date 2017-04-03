# Fuselage 

Included in this folder are two separate fuselage optimization models.  One is an ellipsoid shaped fuselage in `elliptical_fuselage.py` the other is a cylindrical shaped fusealge with an ellipctical nose cone and an "o-jive" profile for the rear bulkhead.

## Elliptical Fuselage

The elliptical fuselage assumes a constant skin thickness.  Drag model is based off of skin friction drag.  No loading cases are assumed. A detailed write up can be found in `./ellipsoid_fuselage.pdf`.

## Cylindrical Fusealge

The cylindrical fuselage has 3 sections, nose, body and bulkhead, whose shapes and surface area are described in the diagram below:

![Cylindrical Fuselage](cyl_fuse_drawing.pdf)

It is assumed that all of the aircraft fuel goes into the cylindrical, middle section.  This middle section must have the volume capacity then to store all of the fuel

$$ V_{body} \geq V_{fuel} $$

The fuselage weight is comprised of the skin weight.  The fuselage skin is a separate model created when the fuselage model is created.  

The aerodynamic model of the fuselage is based off of a fit to CFD runs in Solidworks for various length to radius ratios of the 3 parts: nose, body, and bulkhead. The model is called FuselageAero.

### Fuselage Skin

The fuselage skin assumes that the body section of fuselage takes all of the tosional and bending loads.  The thickness is determined by the loads or the minimum gague thickness.  A constant thickness is assumed for all 3 sections

\begin{align*}
    S \geq S_{nose} + S_{body} + S_{bulk} \\
    W >= S*rho_{Kevlar}*t*g\\
    \end{align*}

The skin is subjected to 2 loads.  The first is a pull up load constraining the thickness in either stress or deflection.  The second is a landing load. 

