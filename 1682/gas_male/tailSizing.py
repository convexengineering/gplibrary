import numpy as np
from gpkit import VectorVariable, Variable, Model, SignomialsEnabled, units
from gpkit.constraints.costed import CostedConstraintSet
from gpkit import LinkedConstraintSet as LCS
from gpkit.constraints.tight import TightConstraintSet as TCS
from gpkit.constraints.bounded import BoundedConstraintSet as BCS
from gpkit.tools import te_exp_minus1
from collections import defaultdict
from gpkit.small_scripts import mag
pi = 3.1415926
g = 9.81*units('m*s**-2')

class dartTail(Model):
    def __init__(self,**kwargs):

        # Aircraft properties
        W      = Variable('W',713.5,'N','Total aircraft weight')
        WNoPay = Variable('W_{NoPay}',624.3,'N','Aircraft no-payload weight')
        ACloc  = Variable('AC_{location}',0.655,'ft','Aerodynamic center location')
        M_CG   = Variable('M_{CG}',52,'N*m','Torque around AC due to CG')
        S      = Variable('S',23.69,'ft^2','Wing area')
        AR     = Variable('AR',26.7,'-','Aspect ratio')

        # Airfoil properties (NACA0008)
        areaAF = Variable('A_{ratio-airfoil}',0.0548,'ft^2','Airfoil area/chord ratio') #for NACA0008
        crefAF = Variable('c_{ref}_{AF}',1,'ft','Reference chord for airfoil')

        # Takeoff conditions

        rhoTO  = Variable('\\rho_{t/o}','kg*m^-3','Takeoff density')
        Vstall = Variable('V_{stall}','m/s','Stall speed')
        VTO    = Variable('V_{TO}','m/s', 'Takeoff speed')
        CLmax  = Variable('C_{Lmax}','-','Maximum lift coefficient')

        # # Material properties
        rhoCFRP     = Variable('\\rho_{CFRP}',1.6,'kg/m^3','Density of CFRP')
        rhoFoamular = Variable('\\rho_{Foamular}',1.5,'lbf/ft^3','Density of Foamular 250')
        rhoskin = Variable('\\rho_{skin}',.1,'g/cm^2','Skin density')
        
        # Tail properties
        Wtail = Variable('W_{tail}','N','Total tail weight')

        # # Tail boom variables
        Wboom = Variable('W_{boom}','N','Tail boom weight')
        lboom = Variable('l_{boom}','ft','Tail boom length')
        I0boom = Variable('I_0_{boom}','m^4','Tail boom root moment of inertia')
        d0boom = Variable('d_0_{boom}','m','Tail boom root diameter')
        t0boom = Variable('t_0_{boom}','m','Tail boom root wall thickness')
        Eboom = Variable('E_{boom}','N/m^2','Tail boom modulus of elasticity')
        Fboom = Variable('F_{boom}','N','Tail downforce')
        kboom = Variable('k_{boom}','-','Tail boom index (1-.5k)') # k = 0, uniform thickness, 1, constant taper to zero
        #yboom = Variable('y_{boom}','m','Tail deflection at max force')
        #yolboom = Variable('y/l_{boom}','-','Max tolerated tail deflection factor')
        thetaboom = Variable('\\theta_{boom}','-','Tail boom deflection angle')
        Ffacboom = Variable('F_{fac-boom}',.5,'-','Tail downforce factor')

        # # Horizontal tail variables
        CLmaxhtail = Variable('CL_{max-htail}','-','Horizontal tail maximum lift coefficient')
        Vhtail = Variable('V_{htail}','-','Horizontal tail volume coefficient') # 0.5 common for sailplane
        Shtail      = Variable('S_{htail}','ft^2','Horizontal tail area')
        ARhtail     = Variable('AR_{htail}',5,'-','Horizontal tail aspect ratio')
        lamhtail    = Variable('\\lambda_{htail}','-','Horizontal tail taper ratio')
        bhtail      = Variable('b_{htail}','m', 'Horizontal tail span')
        crhtail     = Variable('c_r_{htail}','m','Horizontal tail root chord')
        deltatail = Variable('\\delta_{tail}',.2,'m','Horizontal-vertical tail offset')
        Whtail = Variable('W_{htail}','lbf','Horizontal tail weight')

        # # Vertical tail variables
        #Vvtail = Variable('V_{vtail}','-','Vertical tail volume coefficient') # 0.02 common for sailplanes
        # Svtail      = Variable('S_{vtail}','m^2','Vertical tail area')
        # ARvtail     = Variable('AR_{vtail}','-','Vertical tail aspect ratio')
        # lamvtail    = Variable('\\lambda_{vtail}','-','Vertical tail taper ratio')
        # hvtail      = Variable('b_{vtail}','m', 'Vertical tail height')
        # crvtail     = Variable('c_r_{vtail}','m','Vertical tail root chord')

        with SignomialsEnabled():

            constraints = [        
            # Boom sizing
            kboom >= 0.6, kboom <= 1, # Constraining boom inertia variable
            M_CG <= 2*Ffacboom*Fboom*(lboom),
            TCS([I0boom == pi*t0boom*d0boom**3/8]),
            Eboom == 150*10**9*units('N/m^2'),
            Wboom == pi*g*rhoCFRP*d0boom*lboom*t0boom*(kboom),
            thetaboom <= 0.01,
            thetaboom >= Fboom*lboom**2/(Eboom*I0boom)*(kboom),
            Fboom == .5*rhoTO*Vstall**2*Shtail*CLmaxhtail,
            # Horizontal tail relations (sized for heavy forward CG (20 lb payload))
            bhtail**2/Shtail == ARhtail,
            Shtail <= bhtail*crhtail*(1+lamhtail)/2,
            TCS([CLmaxhtail*(1+2/ARhtail) <= CLmax*(1+2/AR)]),

            # Boom materials constraints
            t0boom >= 0.25*units('mm'),
            d0boom <= 1*units('in'),

            # Assuming solid foam-core wing with a min-gauge Kevlar skin
            Whtail >= (rhoFoamular*bhtail*areaAF)*((crhtail/crefAF)**2 + (crhtail*lamhtail/crefAF)**2)/2+(1.1*g*rhoskin*Shtail),
            Wtail >= Wboom + Whtail
            # Vertical tail relations (sized for cross-wind landing)

            ]

        Model.__init__(self, None, constraints,**kwargs)

class GasMALE(Model):
    def __init__(self,**kwargs):
        W      = Variable('W',713.5,'N','Total aircraft weight')
        WNoPay = Variable('W_{NoPay}',624.3,'N','Aircraft no-payload weight')
        ACloc  = Variable('AC_{location}',0.655,'ft','Aerodynamic center location')
        M_CG   = Variable('M_{CG}',52,'N*m','Torque around AC due to CG')
        S      = Variable('S',23.69,'ft^2','Wing area')
        b      = Variable('b',25.15,'ft','Wing span')
        cbar   = Variable('\\bar_c',1.02,'ft','Mean aerodynamic chord')
        AR     = Variable('AR',26.7,'-','Aspect ratio')
        e      = Variable('e',0.95,'-','Oswald efficiency')
        Wtail = Variable('W_{tail}','N','Total tail weight')

        # Takeoff conditions
        rhoTO  = Variable('\\rho_{t/o}','kg*m^-3','Takeoff density')
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

        # Other performance variables
        nmax   = Variable('n_{max}','-','Maximum load factor')
        RturnTO = Variable('R_{turn}','m','Turning radius at takeoff')
        RpullTO = Variable('R_{pull-up}','m','Pull-up radius at takeoff')
        Fboom = Variable('F_{boom}','N','Tail downforce')
        CLmaxhtail = Variable('CL_{max-htail}','-','Horizontal tail maximum lift coefficient')

        tail = dartTail()
        self.submodels = [tail]

        
        constraints = [
                  # Performance metric calculations
                    TCS([Vstall**2 == (2/rhoTO*W/S/CLmax)]),
                    TCS([VTO == 1.3*Vstall]),
                    TCS([CLmax == 1.5]),
                    TCS([rhoTO == 1.225*units('kg*m^-3')])
                    #KTO == 1/(pi*e*AR),
                    #nmax**2 <= .5*rhoTO*VTO**2/(KTO*W/S)*(TmaxTO/W - .5*rhoTO*VTO**2*CD0TO/(W/S)),
                    #RturnTO**2*g**2*(nmax**2-1.0) >= VTO**4, 
                    #(RpullTO*g*(nmax+1))**2 <= VTO**4
                    ]    

        lc = LCS([self.submodels, constraints])

        objective = tail['W_{tail}']
        Model.__init__(self, objective, lc, **kwargs)


if __name__ == "__main__":
    M = GasMALE()
    #M = Model(M.cost, BCS(M))
    #subs=[]
    sol = M.localsolve("mosek")
    varVals = sol['variables']
