import numpy as np
from gpkit import VectorVariable, Variable, Model, SignomialsEnabled, units
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
        M_CG   = Variable('M_{CG}',32,'N*m','Torque around AC due to CG')
        S      = Variable('S',23.69,'ft^2','Wing area')
        AR     = Variable('AR',26.7,'-','Aspect ratio')
        mac    = Variable('mac',1.02,'ft','Wing mean aerodynamic chord')
        rhoh   = Variable('rho_[h}',.776,'kg/m^3','Density at 15,000 ft')
        qNE    = Variable('q_{NE}','Pa','Never-exceed dynamic pressure')
        VNE    = Variable('V_{NE}',46,'m/s','Never-exceed speed')

        # Airfoil properties (NACA0008)
        areaAF = Variable('A_{ratio-airfoil}',0.0548,'ft^2','Airfoil area/chord ratio') #for NACA0008
        crefAF = Variable('c_{ref}_{AF}',1,'ft','Reference chord for airfoil')

        # Takeoff conditions

        rhoTO  = Variable('\\rho_{t/o}','kg*m^-3','Takeoff density')
        Vstall = Variable('V_{stall}','m/s','Stall speed')
        VTO    = Variable('V_{TO}','m/s', 'Takeoff speed')
        CLmax  = Variable('C_{Lmax}','-','Maximum lift coefficient')

        # # Material properties
        rhoCFRP     = Variable('\\rho_{CFRP}',1.6,'g/cm^3','Density of CFRP')
        rhoFoamular = Variable('\\rho_{Foamular}',1.5,'lbf/ft^3','Density of Foamular 250')
        rhoskin     = Variable('\\rho_{skin}',.1,'g/cm^2','Skin density')
        
        # Tail properties
        Wtail       = Variable('W_{tail}','N','Total tail weight')
        FNE         = Variable('F_{NE}','-','Boom flexibility factor')
        
        # # Tail boom variables
        Wboom       = Variable('W_{boom}','lbf','Tail boom weight')
        lboom       = Variable('l_{boom}','ft','Tail boom length')
        I0boom      = Variable('I_0_{boom}','m^4','Tail boom root moment of inertia')
        d0boom      = Variable('d_0_{boom}','m','Tail boom root diameter')
        t0boom      = Variable('t_0_{boom}','m','Tail boom root wall thickness')
        Eboom       = Variable('E_{boom}','N/m^2','Tail boom modulus of elasticity')
        Fboom       = Variable('F_{boom}','N','Tail downforce')
        Fboommax    = Variable('F_{boom-max}','N','Max tail force')
        kboom       = Variable('k_{boom}','-','Tail boom index (1-.5k)') # k = 0, uniform thickness, 1, constant taper to zero
        #yboom      = Variable('y_{boom}','m','Tail deflection at max force')
        #yolboom    = Variable('y/l_{boom}','-','Max tolerated tail deflection factor')
        thetaboom   = Variable('\\theta_{boom}','-','Tail boom deflection angle')
        Ffacboom    = Variable('F_{fac-boom}',.4,'-','Tail downforce factor')

        # # Horizontal tail variables
        CLmaxhtail  = Variable('CL_{max-htail}','-','Horizontal tail maximum lift coefficient')
        Re_tiphtail = Variable('Re_{tip-htail}','-','Horizontal tail tip Reynolds number')
        Vhtail      = Variable('V_{htail}','-','Horizontal tail volume coefficient') # 0.5 common for sailplane
        Shtail      = Variable('S_{htail}','ft^2','Horizontal tail area')
        ARhtail     = Variable('AR_{htail}',5, '-','Horizontal tail aspect ratio')
        lamhtail    = Variable('\\lambda_{htail}',.8,'-','Horizontal tail taper ratio')
        bhtail      = Variable('b_{htail}','m', 'Horizontal tail span')
        crhtail     = Variable('c_r_{htail}','m','Horizontal tail root chord')
        deltatail   = Variable('\\delta_{tail}',.2,'m','Horizontal-vertical tail offset')
        Whtail      = Variable('W_{htail}','lbf','Horizontal tail weight')
        mhtail      = Variable('m_{htail}','-','Horizontal tail moment coefficient')
        Vhtail      = Variable('V_{vtail}','-','Vertical tail volume coefficient')

        # # Vertical tail variables
        #Vvtail = Variable('V_{vtail}','-','Vertical tail volume coefficient') # 0.02 common for sailplanes
        Svtail      = Variable('S_{vtail}','ft^2','Vertical tail area')
        ARvtail     = Variable('AR_{vtail}',4,'-','Vertical tail aspect ratio')
        lamvtail    = Variable('\\lambda_{vtail}',.8,'-','Vertical tail taper ratio')
        hvtail      = Variable('h_{vtail}','m', 'Vertical tail height')
        crvtail     = Variable('c_r_{vtail}','m','Vertical tail root chord')
        CLmaxvtail  = Variable('CL_{max-vtail}','-','Vertical tail maximum lift coefficient')
        Wvtail      = Variable('W_{vtail}','lbf','Vertical tail weight')


        # Crosswind landing variables
        Vland      = Variable('V_{land}',12,'m/s','Landing speed')
        Vwindcross = Variable('V_{wind_cross}',20,'mph','Landing cross-wind speed')
        CDy = Variable('C_{Dy}',.5,'-','Crosswind drag coefficient')
        Vrel = Variable('V_{rel}','m/s','Relative wind in crosswind')

        constraints = [        
        # Boom sizing
        kboom       >= 0.75, 
        kboom       <= 1, # Constraining boom inertia variable
        M_CG        <= 2*Ffacboom*Fboom*(lboom),
        TCS([I0boom == pi*t0boom*d0boom**3/8]),
        Eboom       == 150*10**9*units('N/m^2'),
        Wboom       >= pi*g*rhoCFRP*d0boom*lboom*t0boom*(kboom),
        thetaboom   <= 0.05,
        thetaboom   >= Fboom*lboom**2/(Eboom*I0boom)*(kboom),
        Fboommax    >= .5*rhoh*(46*units('m/s'))**2*Shtail*CLmaxhtail,
        Fboommax    >= .5*rhoh*(46*units('m/s'))**2*Svtail*CLmaxvtail,
        Fboom       >= .5*rhoTO*Vrel**2*Svtail*CLmaxvtail,
        Fboom       >= .5*rhoTO*Vstall**2*Shtail*CLmaxhtail,
        FNE**-1     >= 1 + mhtail*qNE*Shtail*lboom**2*kboom/(Eboom*I0boom),
        FNE         <= 1,
        TCS([mhtail*(1+2/ARhtail) <= 2*pi]),
            
        # Horizontal tail relations (sized for heavy forward CG (20 lb payload))
        bhtail**2/Shtail              == ARhtail,
        Shtail                        <= bhtail*crhtail*(1+.8)/2, #Substituted lambda so it wouldn't be SP
        TCS([CLmaxhtail*(1+2/ARhtail) <= CLmax*(1+2/27)]), #Substituted the aspect ratio of aircraft so it wouldn't be SP
        #TCS([Vhtail                   == Shtail*lboom/(S*mac)]),
        
        # Boom physical constraints
        t0boom                        >= 0.3*units('mm'),
        lboom                         <= 7*units('ft'),
            
        # Vertical tail relations (sized for cross-wind landing)
        TCS([CLmaxvtail*(1+2/ARvtail) <= CLmax*(1+2/27)]), #Substituted the aspect ratio of aircraft so it wouldn't be SP
        hvtail**2/(Svtail)              == ARvtail,
        Svtail                        == hvtail*crvtail*(1+.8)/2, ##Substituted lambda so it wouldn't be SP
        # Landing conditions
        TCS([Vrel**2 >= Vland**2 + Vwindcross**2]),
        Vrel <= 16*units('m/s'),
        TCS([Vwindcross**2*23.67*units('ft^2')*CDy == 2*Vrel**2*Svtail*CLmaxvtail]), #substituted S because of errors with BCS... too hacky
        # Assuming solid foam-core wing with a min-gauge Kevlar skin
        Whtail >= (rhoFoamular*bhtail*areaAF)*((crhtail/crefAF)**2 + (crhtail*.8/crefAF)**2)/2+(1.1*g*rhoskin*Shtail),
        Wvtail >= (rhoFoamular*hvtail*areaAF)*((crvtail/crefAF)**2 + (crvtail*.8/crefAF)**2)/2+(1.1*g*rhoskin*Svtail),
        Wtail  >= 2*(Wboom + Whtail + Wvtail),
        ]

        Model.__init__(self, None, constraints,**kwargs)

class GasMALE(Model):
    def __init__(self,**kwargs):
        W      = Variable('W',713.5,'N','Total aircraft weight')
        WNoPay = Variable('W_{NoPay}',624.3,'N','Aircraft no-payload weight')
        ACloc  = Variable('AC_{location}',0.655,'ft','Aerodynamic center location')
        M_CG   = Variable('M_{CG}',32,'N*m','Torque around AC due to CG')
        S      = Variable('S',23.69,'ft^2','Wing area')
        b      = Variable('b',25.15,'ft','Wing span')
        cbar   = Variable('\\bar_c',1.02,'ft','Mean aerodynamic chord')
        AR     = Variable('AR',26.7,'-','Aspect ratio')
        e      = Variable('e',0.95,'-','Oswald efficiency')
        qNE    = Variable('q_{NE}','Pa','Never-exceed dynamic pressure')
        VNE    = Variable('V_{NE}',46,'m/s','Never-exceed speed')
        rhoh   = Variable('rho_[h}',.776,'kg/m^3','Density at 15,000 ft')
        mu_atm = Variable("\\mu",1.8*10**-5,"N*s/m^2", "Dynamic viscosity")

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

       
        # Other performance variables
        nmax       = Variable('n_{max}','-','Maximum load factor')
        RturnTO    = Variable('R_{turn}','m','Turning radius at takeoff')
        RpullTO    = Variable('R_{pull-up}','m','Pull-up radius at takeoff')
        Fboom      = Variable('F_{boom}','N','Tail downforce')
        CLmaxhtail = Variable('CL_{max-htail}','-','Horizontal tail maximum lift coefficient')

        tail = dartTail()
        self.submodels = [tail]

        
        constraints = [
                    # Performance metric calculations
                    TCS([Vstall**2                 == (2/rhoTO*W/S/CLmax)]),
                    TCS([VTO                       == 1.3*Vstall]),
                    TCS([CLmax                     == 1.5]),
                    TCS([rhoTO                     == 1.225*units('kg*m^-3')]),
                    TCS([qNE                       == .5*rhoh*VNE**2]),
                    #TCS([Re_tiphtail              >= 100000]),
                    #TCS([Re_tiphtail              == rhoh*crhtail*lamhtail*Vstall/mu_atm])
                    #KTO                           == 1/(pi*e*AR),
                    #nmax**2                       <= .5*rhoTO*VTO**2/(KTO*W/S)*(TmaxTO/W - .5*rhoTO*VTO**2*CD0TO/(W/S)),
                    #RturnTO**2*g**2*(nmax**2-1.0) >= VTO**4, 
                    #(RpullTO*g*(nmax+1))**2       <= VTO**4
                    ]    

        lc = LCS([self.submodels, constraints])

        objective = tail['W_{tail}']
        Model.__init__(self, objective, lc, **kwargs)


if __name__ == "__main__":
    M       = GasMALE()
    #M.substitutions.update({'S':23.69*units('ft^2')})
    #M.substitutions.update({'\\lambda_{htail}':0.8})
    M       = Model(M.cost, BCS(M))
    sol     = M.solve("mosek")
    varVals = sol['variables']
    print 'Tail downforce: ' + str(varVals['F_{boom}'])
    print 'Boom weight: ' + str(varVals['W_{boom}'])
    print 'Boom thickness: ' + str(varVals['t_0_{boom}'])
    print 'Boom length: ' + str(varVals['l_{boom}'])
    print 'Boom diameter: ' + str(varVals['d_0_{boom}'])
    print 'Boom k-value: ' + str((varVals['k_{boom}'].magnitude-0.5)*2)
    #print 'Horizontal tail volume coeff: ' + str(varVals['V_{htail}'])
    print 'Htail weight: ' + str(varVals['W_{htail}'])
    print 'Htail surface area: ' + str(varVals['S_{htail}'])
    print 'Htail span: ' + str(varVals['b_{htail}'])
    print 'Htail root chord: ' + str(varVals['c_r_{htail}'])
    print 'CLmaxhtail: ' + str(varVals['CL_{max-htail}'])
    print 'qNE: ' + str(varVals['q_{NE}'])
    print 'Vtail weight: ' + str(varVals['W_{vtail}'])
    print 'Vtail surface area: ' + str(varVals['S_{vtail}'])
    print 'Vtail height: ' + str(varVals['h_{vtail}'])
    print 'Vtail root chord: ' + str(varVals['c_r_{htail}'])
    #print 'Vtail root chord: ' + str(varVals['cr_{vtail}']) 
    print 'CLmaxvtail: ' + str(varVals['CL_{max-vtail}'])
    print 'Relative wind: ' + str(varVals['V_{rel}'])