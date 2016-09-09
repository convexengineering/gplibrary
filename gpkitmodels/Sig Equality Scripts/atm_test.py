from gpkit import Variable, VectorVariable, Model, units, SignomialsEnabled, SignomialEquality
from gpkit.constraints.set import ConstraintSet

class Atmosphere(Model):
    def __init__(self, **kwargs):
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        T_atm = Variable("T_{atm}", "K", "air temperature")
        M_atm = Variable("M_{atm}", 0.0289644, "kg/mol",
                         "Molar mass of dry air")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K",
                         "air specific heating value")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        TH = (g*M_atm/R_atm/L_atm).value

        rho = Variable('\rho', 'kg/m^3', 'Density of air')

        h = Variable("h", "ft", "Altitude")
        

        with SignomialsEnabled():
            constraints = [
                h==10000*units.m,
    
                #h <= 20000*units.m,  # Model valid to top of troposphere

                # Pressure-altitude relation
                (p_atm/p_sl)**(1/TH) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),
                
                SignomialEquality(T_sl, T_atm + L_atm*h)
                # T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                ]
                # http://en.wikipedia.org/wiki/Density_of_air#Altitude

        Model.__init__(self, T_atm, constraints, **kwargs)
        
if __name__ == "__main__":
    M = Atmosphere()
    sol = M.localsolve("mosek",iteration_limit=500)
    print sol.table()
