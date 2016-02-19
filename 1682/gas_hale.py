"""Gas-powered high-altitude-long-endurange UAV model"""
from numpy import pi
from gpkit import Variable, Model, units
from gpkit.tools import te_exp_minus1

class GasPoweredHALE(Model):
    def setup(self):
        constraints = []

        # Steady level flight relations
        CD = Variable('C_D', '-', 'Drag coefficient')
        CL = Variable('C_L', '-', 'Lift coefficient')
        P_shaft = Variable('P_{shaft}', 'W', 'Shaft power')
        T = Variable('Thrust','N','Cruise thrust')
        S = Variable('S', 'm^2', 'Wing reference area')
        V = Variable('V', 'm/s', 'Cruise velocity')
        W = Variable('W', 'lbf', 'Aircraft weight')

        # Propulsion metrics (for a 3 bladed propeller, activity factor 100, design CL = 0.5)
        AdvRatio = Variable('J_{advance}',1.7,'-','Advance ratio')
        CPower = Variable('C_{Power}',0.2,'-','Power coefficient')
        CThrust = Variable('C_{Thrust}',0.5,'-','Thrust coefficient')
        CTorque = Variable('C_{Torque}','-','Torque coefficient')
        nRot = Variable('n_{Rot}','1/s','Propeller rotation speed')
        D_Prop = Variable('D_{Prop}',0.6,'m','Propeller diameter')

        eta_prop = Variable(r'\eta_{prop}',0.9,'-', 'Propulsive efficiency')
        rho = Variable(r'\rho', 'kg/m^3')

        constraints.extend([P_shaft == V*W*CD/CL/eta_prop,   # eta*P = D*V
                            W == 0.5*rho*V**2*CL*S])

        # Aerodynamics model
        Cd0 = Variable('C_{d0}', 0.01, '-', "non-wing drag coefficient")
        CLmax = Variable('C_{L-max}', 1.5, '-', 'maximum lift coefficient')
        e = Variable('e', 0.9, '-', "spanwise efficiency")
        A = Variable('A', 20, '-', "aspect ratio")
        b = Variable('b', 'ft', 'span')
        mu = Variable(r'\mu', 1.5e-5, 'N*s/m^2', "dynamic viscosity")
        Re = Variable("Re", '-', "Reynolds number")
        Cf = Variable("C_f", "-", "wing skin friction coefficient")
        Kwing = Variable("K_{wing}", 1.3, "-", "wing form factor")
        cl_16 = Variable("cl_{16}", 0.0001, "-", "profile stall coefficient")
        
        constraints.extend([CD >= Cd0 + 2*Cf*Kwing + CL**2/(pi*e*A) + cl_16*CL**16,
                            #T == CD*1/2*rho*V**2*S,
                            #T <= P_shaft*(CThrust/CPower)/(nRot*D_Prop),
                            #eta_prop == T*V/P_shaft,
                            #AdvRatio == V/(nRot*D_Prop),
                            #AdvRatio >= 1, AdvRatio <= 2.8,
                            #AdvRatio == 1.8/0.23*CPower + 0.23,
                            #CPower == P_shaft/(rho*nRot**3*D_Prop**5),
                            #CThrust == T/(rho*nRot**2*D_Prop**4),
                            #P_shaft >= 2*pi*nRot*(CTorque*rho*nRot**2*D_Prop**5),
                            eta_prop == 1/(2*pi)*(CThrust/CTorque)*AdvRatio,
                            b**2 == S*A,
                            CL <= CLmax, 
                            Re == rho*V/mu*(S/A)**0.5,
                            Cf >= 0.074/Re**0.2])


        # Engine Weight Model
        W_eng = Variable('W_{eng}', 'lbf', 'Engine weight')
        W_engmin = Variable('W_{min}', 11, 'lbf', 'min engine weight')
        W_engmax = Variable('W_{max}', 275, 'lbf', 'max engine weight')
        eta_t = Variable('\\eta_t', 0.5, '-', 'percent throttle')
        eng_cnst = Variable('eng_{cnst}', 0.0011, '-',
                            'engine constant based off of power weight model')
        W_eng_installed = Variable('W_{eng-installed}','lbf','Installed engine weight')
        
        constraints.extend([W_eng >= W_engmin,
                            W_eng <= W_engmax,
                            W_eng >= P_shaft*eng_cnst/eta_t * units('lbf/watt'),
                            W_eng_installed >= 2.572*W_eng**0.922*units('lbf')**0.078])

        # Weight model
        W_airframe = Variable('W_{airframe}', 'lbf', 'Airframe weight')
        W_pay = Variable(r'W_{pay}', 4, 'lbf', 'Payload weight')
        W_fuel = Variable('W_{fuel}', 'lbf', 'Fuel Weight')
        W_zfw = Variable('W_{zfw}', 'lbf', 'Zero fuel weight')
        wl = Variable('wl', 'lbf/ft^2', 'wing loading')
        
        # Higher fidelity weight modeling
        # w_wing = Variable('w_{wing}','lbf','Wing weight')
        # w_tail = Variable('w_{tail}','lbf','Tail weight')
        W_avionics = Variable('w_{avionics}',1,'lbf','Avionics weight')
        # w_boom = Variable('w_{boom}','lbf','Boom weight')

        f_airframe = Variable('f_{airframe}', 0.3, '-',
                              'Airframe weight fraction')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        constraints.extend([W_airframe >= W*f_airframe,
                            W_zfw >= W_airframe + W_eng_installed + W_pay + W_avionics,
                            wl == W/S,
                            W >= W_pay + W_eng_installed + W_airframe + W_fuel + W_avionics])

        # Breguet Range
        z_bre = Variable("z_bre", "-", "breguet coefficient")
        h_fuel = Variable("h_{fuel}", 42e6, "J/kg", "heat of combustion")
        eta_0 = Variable("\\eta_0", 0.2, "-", "overall efficiency")
        t = Variable('t', 5, 'days', 'time on station')

        constraints.extend([z_bre >= g*t*V*CD/(h_fuel*eta_0*CL),
                            W_fuel/W_zfw >= te_exp_minus1(z_bre, 3)])

        # Atmosphere model
        h = Variable("h", "ft", "Altitude")
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        T_atm = Variable("T_{atm}", "K", "air temperature")
        M_atm = Variable("M_{atm}", 0.0289644, "kg/mol",
                         "Molar mass of dry air")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K")
        TH = (g*M_atm/R_atm/L_atm).value.magnitude  # dimensionless
        constraints.extend([h <= 20000*units.m,  # Model valid to top of troposphere
                            T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                            rho <= p_sl*T_atm**(TH-1)*M_atm/R_atm/(T_sl**TH)])
            # http://en.wikipedia.org/wiki/Density_of_air#Altitude

        # station keeping requirement
        footprint = Variable("d_{footprint}", 100, 'km',
                             "station keeping footprint diameter")
        lu = Variable(r"\theta_{look-up}", 5, '-', "look up angle")
        R_earth = Variable("R_{earth}", 6371, "km", "Radius of earth")
        tan_lu = lu*pi/180. + (lu*pi/180.)**3/3.  # Taylor series expansion
        # approximate earth curvature penalty as distance^2/(2*Re)
        constraints.extend([
            h >= tan_lu*0.5*footprint + footprint**2/8./R_earth])

        # wind speed model
        V_wind = Variable('V_{wind}',43, 'm/s', 'wind speed')
        wd_cnst = Variable('wd_{cnst}', 0.002, 'm/s/ft', 
                           'wind speed constant predited by model')
                            #0.002 is worst case, 0.0015 is mean at 45d
        wd_ln = Variable('wd_{ln}', 13.009, 'm/s',
                         'linear wind speed variable')
                        #13.009 is worst case, 8.845 is mean at 45deg
        h_min = Variable('h_{min}', 11800, 'ft', 'minimum height')
        h_max = Variable('h_{max}', 15000, 'ft', 'maximum height')
        constraints.extend([V_wind >= wd_cnst*h + wd_ln, 
                            V >= V_wind,
                            h >= h_min,
                            h <= h_max])
                            # model predicting worse case scenario at 45 deg latitude

        # fuel volume model
        tc = Variable('tc', 0.1, '-', 'thickness to chord ratio')
        rho_fuel = Variable('\\rho_{fuel}', 719, 'kg/m^3', 'density of gasoline')
        constraints.extend([S**1.5*A**-0.5*tc >= W_fuel/g/rho_fuel])


        objective = W

        return objective, constraints

if __name__ == "__main__":
    M = GasPoweredHALE()
    M.solve()
