from gpkit import Variable, Model, SignomialsEnabled, VarKey, units
import numpy as np
import matplotlib.pyplot as plt

def SimPleAC():
    # Env. constants
    g          = Variable("g", 9.81, "m/s^2", "gravitational acceleration")
    mu         = Variable("\\mu", 1.775e-5, "kg/m/s", "viscosity of air", pr=4.)
    rho        = Variable("\\rho", 1.23, "kg/m^3", "density of air", pr=5.)
    rho_f      = Variable("\\rho_f", 817, "kg/m^3", "density of fuel")
    
    # Non-dimensional constants
    C_Lmax     = Variable("C_{L,max}", 1.6, "-", "max CL with flaps down", pr=5.)
    e          = Variable("e", 0.92, "-", "Oswald efficiency factor", pr=3.)
    k          = Variable("k", 1.17, "-", "form factor", pr=10.)
    N_ult      = Variable("N_{ult}", 3.3, "-", "ultimate load factor", pr=15.)
    S_wetratio = Variable("(\\frac{S}{S_{wet}})", 2.075, "-", "wetted area ratio", pr=3.)
    tau        = Variable("\\tau", 0.12, "-", "airfoil thickness to chord ratio", pr=10.)
    W_W_coeff1 = Variable("W_{W_{coeff1}}", 2e-5, "1/m",
                          "Wing Weight Coefficent 1", pr= 30.) #orig  12e-5
    W_W_coeff2 = Variable("W_{W_{coeff2}}", 60., "Pa",
                          "Wing Weight Coefficent 2", pr=10.)
    p_labor = Variable('p_{labor}',1.,'1/min','cost of labor', pr = 20.)

    # Dimensional constants
    CDA0      = Variable("(CDA0)", "m^2", "fuselage drag area") #0.035 originally
    Range     = Variable("Range",3000, "km", "aircraft range")
    toz       = Variable("toz", 1, "-", pr=15.)
    TSFC      = Variable("TSFC", 0.6, "1/hr", "thrust specific fuel consumption")
    V_min     = Variable("V_{min}", 25, "m/s", "takeoff speed", pr=20.)
    W_0       = Variable("W_0", 6250, "N", "aircraft weight excluding wing", pr=20.)
    
    # Free Variables
    LoD       = Variable('L/D','-','lift-to-drag ratio')
    D         = Variable("D", "N", "total drag force")
    V         = Variable("V", "m/s", "cruising speed")
    W         = Variable("W", "N", "total aircraft weight")
    Re        = Variable("Re", "-", "Reynold's number")
    C_D       = Variable("C_D", "-", "drag coefficient")
    C_L       = Variable("C_L", "-", "lift coefficient of wing")
    C_f       = Variable("C_f", "-", "skin friction coefficient")
    W_f       = Variable("W_f", "N", "fuel weight")
    V_f       = Variable("V_f", "m^3", "fuel volume")
    V_f_avail = Variable("V_{f_{avail}}","m^3","fuel volume available")
    T_flight  = Variable("T_{flight}", "hr", "flight time")
    
    # Free variables (fixed for performance eval.)
    A         = Variable("A", "-", "aspect ratio",fix = True)
    S         = Variable("S", "m^2", "total wing area", fix = True)
    W_w       = Variable("W_w", "N", "wing weight")#, fix = True)
    W_w_strc  = Variable('W_w_strc','N','wing structural weight', fix = True)
    W_w_surf  = Variable('W_w_surf','N','wing skin weight', fix = True)
    V_f_wing  = Variable("V_f_wing",'m^3','fuel volume in the wing', fix = True)
    V_f_fuse  = Variable('V_f_fuse','m^3','fuel volume in the fuselage', fix = True)
    constraints = []

    # Drag model
    C_D_fuse = CDA0 / S
    C_D_wpar = k * C_f * S_wetratio
    C_D_ind = C_L ** 2 / (np.pi * A * e)
    constraints += [C_D >= C_D_fuse * toz + C_D_wpar / toz + C_D_ind * toz]


    with SignomialsEnabled():
            # Wing weight model
            #NOTE: This is a signomial constraint that has been GPified. Could revert back to signomial?
        constraints += [W_w >= W_w_surf + W_w_strc,
                        # W_w_strc >= W_W_coeff1 * (N_ult * A ** 1.5 * ((W_0+V_f_fuse*g*rho_f) * W * S) ** 0.5) / tau, #[GP]
                        W_w_strc**2. >= W_W_coeff1**2. * (N_ult**2. * A ** 3. * ((W_0+V_f_fuse*g*rho_f) * W * S)) / tau**2.,
                        W_w_surf >= W_W_coeff2 * S]

    # and the rest of the models
        constraints += [LoD == C_L/C_D,
                    D >= 0.5 * rho * S * C_D * V ** 2,
                    Re <= (rho / mu) * V * (S / A) ** 0.5,
                    C_f >= 0.074 / Re ** 0.2,
                    T_flight >= Range / V,
                    W_0 + W_w + 0.5 * W_f <= 0.5 * rho * S * C_L * V ** 2,
                    W <= 0.5 * rho * S * C_Lmax * V_min ** 2,
                    W >= W_0 + W_w + W_f,
                    V_f == W_f / g / rho_f,
                    V_f_avail <= V_f_wing + V_f_fuse, #[SP]
                    V_f_wing**2 <= 0.0009*S**3/A*tau**2, #linear with b and tau, quadratic with chord!
                    V_f_fuse <= 10*units('m')*CDA0,
                        #TODO: Reduce volume available by the structure!
                    V_f_avail >= V_f,
                    W_f >= TSFC * T_flight * D]

    # return Model(W_f/LoD, constraints)
    #return Model(D, constraints)
    return Model(W_f, constraints)
    # return Model(W,constraints)
    # return Model(W_f*T_flight,constraints)
    # return Model(W_f + 1*T_flight*units('N/min'),constraints)