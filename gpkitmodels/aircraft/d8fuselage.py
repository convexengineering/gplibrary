""" d8fuselage.py """
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units, SignomialsEnabled
from gpkit import LinkedConstraintSet
#from gpkit.tools import BoundedConstraintSet
from gpkit.tools import te_exp_minus1

class Fuselage(Model):
    def __init__(self, **kwargs):
        constraints = []
        # Will try to stick to Philippe's naming methods as closely as possible
        # for cross-compatibility (to be able to switch models quickly)

        # Fixed variables
        SPR      = Variable('SPR', 8, '-', 'Number of seats per row')
        #npass    = Variable('n_{pass}', '-', 'Number of passengers')
        #nrows    = Variable('n_{rows}', '-', 'Number of rows')
        #Nland    = Variable('N_{land}', '-', 'Emergency landing load factor')
        #Pfloor   = Variable('P_{floor}', 'N', 'Distributed floor load')
        #Pcargofloor = Variable ('P_{cargo floor}','N','Distributed cargo floor load')
        dPover = Variable('\\delta P_{over-pressure}',18.4,'psi','Cabin overpressure (2P)')
        
        # Cross sectional parameters (free)
        #Afloor   = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
        Afuse    = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        #Ahold    = Variable('A_{hold}', 'm^2', 'Cargo hold x-sectional area')
        Askin    = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        hdb = Variable('h_{db}','m', 'Web half-height')
        Rfuse    = Variable('R_{fuse}', 'm', 'Fuselage radius') # will assume for now there is
                                                                # no under-fuselage extension deltaR
        tdb    = Variable('t_{db}', 'm', 'Web thickness')
        tshell = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin  = Variable('t_{skin}', 'm', 'Skin thickness')
        waisle = Variable('w_{aisle}',0.51, 'm', 'Aisle width')
        wdb = Variable('w_{db}','m','DB added half-width')
        wfuse  = Variable('w_{fuse}', 'm', 'Fuselage width')
        wseat    = Variable('w_{seat}',0.5,'m', 'Seat width')


        # Lengths (free)
        lfuse    = Variable('l_{fuse}', 'm', 'Fuselage length')
        lnose    = Variable('l_{nose}', 'm', 'Nose length')
        lshell   = Variable('l_{shell}', 'm', 'Shell length')
      
        # Surface areas (free)
        Sbulk    = Variable('S_{bulk}', 'm^2', 'Bulkhead surface area')
        Snose    = Variable('S_{nose}', 'm^2', 'Nose surface area')

        # Volumes (free)
        Vdb    = Variable('V_{db}', 'm^3', 'Web volume')
        Vcyl   = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')
        Vnose  = Variable('V_{nose}', 'm^3', 'Nose skin volume')
        Vbulk  = Variable('V_{bulk}', 'm^3', 'Bulkhead skin volume')
        Vcabin = Variable('V_{cabin}', 'm^3', 'Cabin volume')


        # Weights (free)
        Wbuoy    = Variable('W_{buoy}', 'N', 'Buoyancy weight')
        Wfuse    = Variable('W_{fuse}', 'N', 'Fuselage weight')

        # Weights (fixed)
        Wcargo   = Variable('W_{cargo}', 10000, 'N', 'Cargo weight')
        Wavgpass = Variable('W_{avg. pass}', 180, 'lbf', 'Average passenger weight')
        Wcarryon = Variable('W_{carry on}', 15, 'lbf', 'Ave. carry-on weight')
        Wchecked = Variable('W_{checked}', 40, 'lbf', 'Ave. checked bag weight')
        Wfix     = Variable('W_{fix}', 3000, 'lbf',
                            'Fixed weights (pilots, cockpit seats, navcom)')

        # Weight fractions (fixed, with respect to the aircraft skin weight, set from PRSEUS metrics)

        ffadd     = Variable('f_{fadd}', '-','Fractional added weight of local reinforcements')
        fstring   = Variable('f_{string}','-','Fracional stringer weight')
        fframe    = Variable('f_{frame}',0.634,'-', 'Fractional frame weight')
        ffairing  = Variable('f_{fairing}','-','  Fractional fairing weight')
        fwebcargo = Variable('f_{web}', '-','Fractional web and cargo floor weight')

        # Misc free variables
        thetadb = ('\\theta_{db}','-','DB fuselage joining angle')

        # Material properties
        sigskin  = Variable('\\sigma_{skin}', 46000,'psi',
                            'Max allowable skin stress')
        sigth    = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx     = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')

        with SignomialsEnabled():
            constraints = [

            # Fuselage joint angle relations
            thetadb == wdb/Rfuse, # first order Taylor works...
            hdb >= Rfuse*(1-thetadb**2/2),
            Askin >= (2*pi + 4*thetadb)*Rfuse*tskin, #no delta R for now
            Adb >= (2*hdb)*tdb,
            Afuse >= (pi + 2*thetadb + thetadb)*Rfuse**2,
            
            # Fuselage surface area relations
            Snose >= (2*pi + 4*thetadb)*Rfuse**2 *(1/3 + 2/3*(lnose/Rfuse)**(8/5))**(5/8),
            Sbulk >= (2*pi + 4*thetadb)*Rfuse**2,

            # Fuselage length relations
            lfuse >= lnose+lshell+lcone,

            # Fuselage width relations
            wfuse >= SPR*wseat + 2*waisle + tdb + 2*tskin,


            # Fuselage volume relations
            Vcyl == Askin*lshell,
            Vnose == Snose*tskin,
            Vbulk == Sbulk*tskin,
            Vdb == Adb*lshell,

            # Fuselage weight relations
            Wskin >= rhoskin*g*(Vcyl + Vnose + Vbulk),
            Wshell >= Wskin*(1 + fstring + ffairing + fframe + ffadd +fweb),
            
            Wfuse >= Wfix + Wskin + Wshell + Wbuoy,
            
            # Stress relations
            sigth == dPover*Rfuse/tskin, 
            sigx == dPover/2*Rfuse/tshell,


            ]

        objective == 1/wfuse;            

        Model.__init__(self, None, constraints, **kwargs)

# class Aircraft(Model):
#     """
#     Combined fuselage, tail, and landing gear model
#     """

#     def __init__(self):

#         # Free variables
#         W      = Variable('W', 'N', 'Total aircraft weight')
#         Wfuse  = Variable('W_{fuse}', 'N', 'Fuselage weight')
#         Wlg    = Variable('W_{lg}', 'N', 'Landing gear weight')
#         Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
#         xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
#         xCGfu  = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
#         xCGlg  = Variable('x_{CG_{lg}}', 'm', 'x-location of landing gear CG')
#         xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of tail CG') 

#         # Fixed variables (pulled from Philippe's aircraft model)
#         Weng    = Variable('W_{eng}', 10000, 'N', 'Engine weight')
#         Wht     = Variable('W_{ht}', 5000, 'N', 'Horizontal tail weight')
#         Wwing   = Variable('W_{wing}', 30000, 'N', 'Wing weight')
#         xCGeng  = Variable('x_{CG_{eng}}', 15, 'm', 'x-location of engine CG')
#         xCGht   = Variable('x_{CG_{ht}}', 38, 'm', 'x-location of horizontal tail CG')
#         xCGwing = Variable('x_{CG_{wing}}', 15, 'm', 'x-location of wing CG')

#         fuselage = Fuselage()

#         self.submodels = [Fuselage]

#         constraints = [];

#         lc = LinkedConstraintSet([self.submodels, constraints],
#                                  include_only=INCLUDE)

#         objective = 1/W;

#         Model.__init__(self, objective, lc, **kwargs)

if __name__ == "__main__":
    M = Fuselage()
    sol = M.solve("mosek")