""" d8fuselage.py """
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units, SignomialsEnabled
from gpkit import LinkedConstraintSet as LSC
from gpkit.constraints.bounded import BoundedConstraintSet as BCS
from gpkit.tools import te_exp_minus1
from collections import defaultdict
from gpkit.small_scripts import mag



class Fuselage(Model):
    def __init__(self, **kwargs):
        constraints = []
        g = 9.81*units('m*s**-2')

        # Will try to stick to Philippe's naming methods as closely as possible
        # for cross-compatibility (to be able to switch models quickly)

        # Fixed variables
        SPR          = Variable('SPR', 8, '-', 'Number of seats per row')
        nseats       = Variable('n_{seat}',192,'-','Number of seats')
        nrows        = Variable('n_{rows}', '-', 'Number of rows')
        pitch        = Variable('p_s',81, 'cm', 'Seat pitch')
        Nland        = Variable('N_{land}',6.,'-', 'Emergency landing load factor')
        Pfloor       = Variable('P_{floor}','N', 'Distributed floor load')
        #Pcargofloor = Variable ('P_{cargo floor}','N','Distributed cargo floor load')
        dPover       = Variable('\\delta_P_{over-pressure}',18.4,'psi','Cabin overpressure')
        Sfloor       = Variable('S_{floor}', 'N', 'Maximum shear in floor beams')
        Mfloor       = Variable('M_{floor}', 'N*m', 'Max bending moment in floor beams')
        
        # Cross sectional parameters (free)
        Adb     = Variable('A_{db}', 'm^2', 'Web cross sectional area')
        Afloor  = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
        Afuse   = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        Ahold   = Variable('A_{hold}', 'm^2', 'Cargo hold x-sectional area')
        Askin   = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        hdb     = Variable('h_{db}','m', 'Web half-height')
        hfloor   = Variable('h_{floor}', 'm', 'Floor beam height')
        Rfuse   = Variable('R_{fuse}', 'm', 'Fuselage radius') # will assume for now there: no under-fuselage extension deltaR
        tdb     = Variable('t_{db}', 'm', 'Web thickness')
        tshell = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin   = Variable('t_{skin}', 'm', 'Skin thickness')
        waisle  = Variable('w_{aisle}',0.51, 'm', 'Aisle width')
        wdb     = Variable('w_{db}','m','DB added half-width')
        wfuse   = Variable('w_{fuse}', 'm', 'Fuselage width')
        wfloor  = Variable('w_{floor}', 'm', 'Floor width')
        wseat   = Variable('w_{seat}',0.5,'m', 'Seat width')


        # Lengths (free)
        lfuse    = Variable('l_{fuse}', 'm', 'Fuselage length')
        lnose    = Variable('l_{nose}', 'm', 'Nose length')
        lshell   = Variable('l_{shell}', 'm', 'Shell length')
        lfloor = Variable('l_{floor}', 'm', 'FLoor length')
      
        # Surface areas (free)
        Sbulk    = Variable('S_{bulk}', 'm^2', 'Bulkhead surface area')
        Snose    = Variable('S_{nose}', 'm^2', 'Nose surface area')

        # Volumes (free)        
        Vbulk  = Variable('V_{bulk}', 'm^3', 'Bulkhead skin volume')
        Vcabin = Variable('V_{cabin}', 'm^3', 'Cabin volume')
        Vcyl   = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')   
        Vdb    = Variable('V_{db}', 'm^3', 'Web volume')
        Vfloor   = Variable('V_{floor}', 'm^3', 'Floor volume')
        Vnose  = Variable('V_{nose}', 'm^3', 'Nose skin volume')

        # Weights (free)
        #Wbuoy    = Variable('W_{buoy}', 'N', 'Buoyancy weight')
        Wfloor   = Variable('W_{floor}', 'N', 'Floor weight')
        Wfuse    = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Wshell = Variable('W_{shell}','N','Shell weight')
        Wskin    = Variable('W_{skin}', 'N', 'Skin weight')
        Wpay     = Variable('W_{pay}', 'N', 'Payload weight')
        Wseat    = Variable('W_{seat}', 'N', 'Seating weight')

        # Weights (fixed)
        #Wcargo   = Variable('W_{cargo}', 10000, 'N', 'Cargo weight')
        #Wavgpass = Variable('W_{avg. pass}', 180, 'lbf', 'Average passenger weight')
        #Wcarryon = Variable('W_{carry on}', 15, 'lbf', 'Ave. carry-on weight')
        #Wchecked = Variable('W_{checked}', 40, 'lbf', 'Ave. checked bag weight')
        Wfix     = Variable('W_{fix}', 3000, 'lbf',
                           'Fixed weights (pilots, cockpit seats, navcom)')

        # Weight fractions (fixed, with respect to the aircraft skin weight, set from PRSEUS metrics)

        #ffadd     = Variable('f_{fadd}', '-','Fractional added weight of local reinforcements')
        fstring   = Variable('f_{string}',0.235,'-','Fractional stringer weight')
        fframe    = Variable('f_{frame}',0.634,'-', 'Fractional frame weight')
        ffairing  = Variable('f_{fairing}',0.151,'-','  Fractional fairing weight')
        fwebcargo = Variable('f_{web}',1.030, '-','Fractional web and cargo floor weight')

        # Misc free variables
        thetadb = Variable('\\theta_{db}','-','DB fuselage joining angle')

        # Material properties
        sigfloor = Variable('\\sigma_{floor}',30000/0.000145, 'N/m^2', 'Max allowable floor stress') #TASOPT value used
        rhofloor = Variable('\\rho_{floor}',2700, 'kg/m^3', 'Floor material density') #TASOPT value used
        taufloor = Variable('\\tau_{floor}',30000/0.000145, 'N/m^2', 'Max allowable shear web stress') #TASOPT value used
        Wppfloor = Variable('W\'\'_{floor}', 60,'N/m^2', 'Floor weight/area density') #TAS
        rhoskin  = Variable('\\rho_{skin}',2,'g/cm^3', 'Skin density') # notional,based on CFRP
        sigskin  = Variable('\\sigma_{skin}', 46000,'psi',
                            'Max allowable skin stress') # again notional 
        rE       = Variable('r_E', 1,'-', 'Ratio of stringer/skin moduli')
        rhobend  = Variable('\\rho_{bend}',2700, 'kg/m^3', 'Stringer density')
        sigth    = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx     = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')

        # Bending inertias (ported from TASOPT)
        Ivshell = Variable('I_{vshell}','m^4','Shell vertical bending inertia')
        Ihshell = Variable('I_{hshell}','m^4','Shell horizontal bending inertia')


        with SignomialsEnabled():
            constraints = [
            # Passenger constraints (assuming 737-sixed aircraft)
            #Temporarily
            Wpay == 150000*units('N'),
            Wseat == 50000*units('N'),

            nrows == nseats/SPR,
            lshell == nrows*pitch,

            # Fuselage joint angle relations
            thetadb == wdb/Rfuse, # first order Taylor works...
            thetadb >= 0.05, thetadb <= 0.25,
            hdb >= Rfuse*(1.0-.5*thetadb**2), #[SP]

            # Fuselage cross-sectional relations
            Askin >= (2*pi + 4*thetadb)*Rfuse*tskin + Adb, #no delta R for now
            Adb == (2*hdb)*tdb,
            Afuse >= (pi + 2*thetadb + thetadb)*Rfuse**2,
            tshell >= tskin*(1+rE*fstring*rhoskin/rhobend),
            
            # Fuselage surface area relations
            Snose >= (2*pi + 4*thetadb)*Rfuse**2 *(1/3 + 2/3*(lnose/Rfuse)**(8/5))**(5/8),
            Sbulk >= (2*pi + 4*thetadb)*Rfuse**2,

            # Fuselage length relations
            lshell == nrows*pitch,
            #lfuse >= lnose+lshell, # forget about tailcone for now
            # Temporarily
            lnose == 0.3*lshell,
            lfuse == 1.3*lshell,

            # Fuselage width relations
            wfuse >= SPR*wseat + 2*waisle + tdb + 2*tskin,
            wfuse <= 2*(Rfuse + wdb),
            wfloor >= wdb + Rfuse, # half of the total floor width in fuselage

            # Fuselage volume relations
            Vcyl == Askin*lshell,
            Vnose == Snose*tskin,
            Vbulk == Sbulk*tskin,
            Vdb == Adb*lshell,
            Vcabin >= Afuse*(lshell + 0.67*lnose + 0.67*Rfuse),

            # Fuselage weight relations
            Wskin >= rhoskin*g*(Vcyl + Vnose + Vbulk),
            Wshell >= Wskin*(1 + fstring + ffairing + fframe +fwebcargo),
            Wfuse >= Wfix + Wshell + Wfloor,
            
            ## Stress relations
            #Pressure shell loading
            tskin == dPover*Rfuse/sigskin,
            tdb == 2*dPover*wdb/sigskin,
            sigx == dPover*Rfuse/(2*tshell),
            sigth == dPover*Rfuse/tskin,

            # Floor loading (don't understand some of these relations,
            # so might be useful to go through them with someone)
            lfloor >= lshell + 2*Rfuse,
            Pfloor >= Nland*(Wpay + Wseat),
            Mfloor >= 9./256.*Pfloor*wfloor,
            Afloor >= 2.*Mfloor/(sigfloor*hfloor) + 1.5*Sfloor/taufloor,
            Vfloor == 2*wfloor*Afloor,
            Wfloor >= rhofloor*g*Vfloor + wfloor*lfloor*Wppfloor,
            Sfloor == (5./16.)*Pfloor,
            # Added synthetic constraint on hfloor to keep it from growing too large
            hfloor <= .1*Rfuse,

            # Fuselage bending model
            Ihshell <= ((pi+4*thetadb)*Rfuse**2)*Rfuse*tshell, # [SP]
            Ivshell <= (pi*Rfuse**2 + 8*wdb*Rfuse + (2*pi+4*thetadb)*wdb**2)*Rfuse*tshell #approx needs to be improved [SP]
            ]

        objective = Wfuse + Afuse*units('N/m**2') + wfloor*units('N/m') + Pfloor + Vcabin*units('N/m^3') + hfloor*units('N/m')
        Model.__init__(self, objective, constraints, **kwargs)


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
    #M = Model(M.cost, BCS(M))
    #bounds, sol = M.determine_unbounded_variables(M, solver="mosek",verbosity=4, iteration_limit=100)
    # subs = {'R_{fuse}_Fuselage':4,'w_{fuse}_Fuselage':10}
    sol = M.localsolve("mosek",tolerance = 0.01, verbosity = 1, iteration_limit=50)
    varVals = sol['variables']
    print 'Cabin volume: ' + str(varVals['V_{cabin}_Fuselage']) +'.'
    print 'Fuselage width: ' + str(varVals['w_{fuse}_Fuselage']) + '.'
    print 'Web height: ' + str(2*varVals['h_{db}_Fuselage']) + '.'
    print 'Floor width: ' + str(varVals['w_{floor}_Fuselage']) + '.'
    print 'Floor height: ' + str(varVals['h_{floor}_Fuselage']) + '.'
    print 'Floor length: ' + str(varVals['l_{floor}_Fuselage']) + '.'
    print  'Fuselage angle: ' + str(varVals['\\theta_{db}_Fuselage']) + '.'
    print 'Fuselage radius: ' + str(varVals['R_{fuse}_Fuselage']) + '.'
    print 'Floor total loading: ' + str(varVals['P_{floor}_Fuselage']) + '.'
    print 'Floor weight: ' + str(varVals['W_{floor}_Fuselage']) + '.'
    print 'Floor volume: ' + str(varVals['V_{floor}_Fuselage']) + '.'
    print 'Floor bending moment: ' + str(varVals['M_{floor}_Fuselage']) + '.'
    print 'Shell thickness: ' + str(varVals['t_{shell}_Fuselage']) + '.'
    print 'Skin hoop stress: ' + str(varVals['\\sigma_{\\theta}_Fuselage']) + '.'
    print 'Skin axial stress: ' + str(varVals['\\sigma_x_Fuselage']) + '.'