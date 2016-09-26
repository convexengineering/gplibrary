import numpy as np
from gpkit import VectorVariable, Variable, Model, SignomialsEnabled, units
from gpkit.constraints.costed import CostedConstraintSet
from gpkit import LinkedConstraintSet, ConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from gpkit.constraints.bounded import BoundedConstraintSet as BCS
from gpkit.tools import te_exp_minus1
from collections import defaultdict
from gpkit.small_scripts import mag
pi = 3.1415926
class GasMale(Model):
	def __init__(self):
		g      = 9.81*units('m*s**-2')
        W      = Variable('W',713.5,'N','Total aircraft weight')
        WNoPay = Variable('W_{NoPay}',624.3,'N','Aircraft no-payload weight')
        ACloc  = Variable('AC_{location}',0.655,'ft','Aerodynamic center location')
        M_CG   = Variable('M_{CG}',52,'N*m','Torque around AC due to CG')
        S      = Variable('S',23.69,'ft^2','Wing area')
        b      = Variable('b',25.15,'ft','Wing span')
        cbar   = Variable('c_{bar}',1.02,'ft','Mean aerodynamic chord')
        AR     = Variable('AR',26.7,'-','Aspect ratio')
        e      = Variable('e',0.95,'-','Oswald efficiency')
        N      = Variable('N','-','Load factor')

        # Takeoff conditions
        rhoTO  = Variable('\\rho_{t/o}',1.225,'kg*m^-3','Takeoff density')
        Vstall = Variable('V_{stall}','m/s','Stall speed')
        VTO    = Variable('V_{TO}','m/s', 'Takeoff speed')
        m_TO   = Variable('m_{TO}',1.3,'-','Stall margin')
        CLmax  = Variable('C_{Lmax}',1.5,'-','Maximum lift coefficient')
        TmaxTO = Variable('T_{maxTO}',13.9,'lbf','Maximum thrust at takeoff (static)') #based off DF70                                                                     # 22x8 prop
        KTO    = Variable('K_{TO}','-','Induced drag multiplier at takeoff')
        CD0TO  = Variable('CD_{0TO}',0.0250,'-','Form drag coefficient at takeoff')
        CDTO   = Variable('CD_{TO}',0.0375,'-','Drag coefficient at takeoff')

        

        # Landing conditions
        Vland = Variable('V_{land}',12,'m/s','Landing speed')
        Vwindcross = Variable('V_wind_cross}',20,'mph','Landing cross-wind speed')
class dartTail(Model):
    def __init__(self):
        

        # Airfoil properties (NACA0008)
        areaAF = Variable('A_{ratio-airfoil}',0.0548,'ft','Airfoil area/chord ratio') #for NACA0008

        # # Material properties
        # rhoCFRP     = Variable('\\rho_{CFRP}',1.6,'kg/m^3','Density of CFRP')
        rhoFoamular = Variable('\\rho_{Foamular}',1.5,'lbf/ft^3','Density of Foamular 250')
        rhoskin = Variable('\\rho_{skin}',.1,'g/cm^2','Skin density')
        # # Tail boom variables
        #rboom = Variable('r_{boom}','m','Tail boom outer radius')
        #tboom = Variable('t_{boom}','m','Tail boom thickness')
        lboom = Variable('l_{boom}',6,'ft','Tail boom length')
        #Iboom = Variable('I_{boom}','m^4','Tail boom moment of inertia')
        #Eboom = Variable('E_{boom}','N/m^2','Tail boom modulus of elasticity')
        Wboom = Variable('W_{boom}','N','Tail boom weight')
        Fboom = Variable('F_{boom}','N','Tail downforce')
        yboom = Variable('y_{boom}','m','Tail deflection at max force')
        yolboom = Variable('y/l_{boom}','-','Max tolerated tail deflection factor')
        Ffacboom = Variable('F_{fac-boom}',.5,'-','Tail downforce factor')

        # # Horizontal tail variables
        CLmaxhtail = Variable('CL_{max-htail}','-','Horizontal tail maximum lift coefficient')
        Vhtail = Variable('V_{htail}','-','Horizontal tail volume coefficient') # 0.5 common for sailplane
        Shtail      = Variable('S_{htail}','ft^2','Horizontal tail area')
        ARhtail     = Variable('AR_{htail}',5,'-','Horizontal tail aspect ratio')
        lamhtail    = Variable('\\lambda_{htail}',.6,'-','Horizontal tail taper ratio')
        # tauhtail    = Variable('\\tau_{htail}',.08,'-', 'Horizontal tail thickness ratio')
        bhtail      = Variable('b_{htail}','m', 'Horizontal tail span')
        crhtail     = Variable('c_r_{htail}','m','Horizontal tail root chord')
        deltatail = Variable('\\delta_{tail}',.2,'m','Horizontal-vertical tail offset')
        Whtail = Variable('W_{htail}','lbf','Horizontal tail weight')

        # # Vertical tail variables
        Vvtail = Variable('V_{vtail}','-','Vertical tail volume coefficient') # 0.02 common for sailplanes
        # Svtail      = Variable('S_{vtail}','m^2','Vertical tail area')
        # ARvtail     = Variable('AR_{vtail}','-','Vertical tail aspect ratio')
        # lamvtail    = Variable('\\lambda_{vtail}','-','Vertical tail taper ratio')
        # tauvtail    = Variable('\\tau_{vtail}',.08,'-', 'Vertical tail thickness ratio')
        # hvtail      = Variable('b_{vtail}','m', 'Vertical tail height')
        # crvtail     = Variable('c_r_{vtail}','m','Vertical tail root chord')
        Lmaxvtail = Variable('L_{max-vtail}','N/m','Maximum vertical tail moment')

        # Other performance variables
        nmax   = Variable('n_{max}','-','Maximum load factor')
        RturnTO = Variable('R_{turn}','m','Turning radius at takeoff')
        RpullTO = Variable('R_{pull-up}','m','Pull-up radius at takeoff')

        with SignomialsEnabled():
            #a = TmaxTO/W
            #b = .5*rhoTO*VTO**2*CD0TO/(W*S)
            #print a, a.units, a.units.to_base_units()
            #print b, b.units, b.units.to_base_units()
            constraints = [
            Vstall == (2/rhoTO*W/S/CLmax)**.5,
            #ARvtail == hvtail**2/(Svtail/2),        #ARhtail == bhtail**2/Shtail,
            VTO == 1.3*Vstall,
            KTO == 1/(pi*e*AR),
            #CD0TO >= CDTO + KTO,


            # Performance metric calculations
            nmax**2 <= .5*rhoTO*VTO**2/(KTO*W/S)*(TmaxTO/W - .5*rhoTO*VTO**2*CD0TO/(W/S)),
            RturnTO**2*g**2*(nmax**2-1.0) >= VTO**4,
            (RpullTO*g*(nmax+1))**2 >= VTO**4,

            # Boom sizing
            M_CG <= 2*Ffacboom*Fboom*(lboom-deltatail-ACloc),
            #Iboom =
            #yboom/lboom <= yolboom,
            Fboom == .5*rhoTO*Vstall**2*Shtail*CLmaxhtail,
            
            # Horizontal tail relations (sized for heavy forward CG (20 lb payload))
            bhtail**2/Shtail == ARhtail,
            Shtail <= bhtail*crhtail*(1+lamhtail)/2, #[SP]
            CLmaxhtail*(1+2/ARhtail) <= CLmax*(1+2/AR),
            # Assuming solid foam-core wing with a min-gauge Kevlar skin
            Whtail**2 >= (rhoFoamular*bhtail*areaAF)**2*(crhtail**2 + (crhtail*lamhtail)**2.)+(1.1*g*rhoskin*Shtail)**2
            
            # Vertical tail relations (sized for cross-wind landing)

            ]
        objective = (nmax)**-1 + Fboom/units('N') + 1/CLmaxhtail + Whtail/units('lbf') + RturnTO/units('m') + RpullTO/units('m')

        Model.__init__(self, objective, constraints,)

if __name__ == "__main__":
    M = dartTail()
    #M = Model(M.cost, BCS(M))
    #subs=[]
    sol = M.localsolve("mosek",tolerance = 0.01, verbosity = 1, iteration_limit=50)
    varVals = sol['variables']