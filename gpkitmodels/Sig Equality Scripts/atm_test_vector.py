from gpkit import Variable, VectorVariable, Model, units, SignomialsEnabled, SignomialEquality
from gpkit.constraints.set import ConstraintSet

class Atmosphere(Model):
    def __init__(self, n, **kwargs):
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        T_atm = VectorVariable(n, "T_{atm}", "K", "air temperature")
        M_atm = Variable("M_{atm}", 0.0289644, "kg/mol",
                         "Molar mass of dry air")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K",
                         "air specific heating value")
        p_atm = VectorVariable(n, "P_{atm}", "Pa", "air pressure")
        TH = (g*M_atm/R_atm/L_atm).value

        rho = VectorVariable(n, '\rho', 'kg/m^3', 'Density of air')

        h = VectorVariable(n, "h", "ft", "Altitude")
        

        with SignomialsEnabled():
            constraints = [
                h[0]==10*units.m,
                h[1]==100*units.m,
                h[2]==100*units.m,
                h[3]==1000*units.m,
                #h <= 20000*units.m,  # Model valid to top of troposphere

                # Pressure-altitude relation
                (p_atm/p_sl)**(1/TH) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                # T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                ]
                # http://en.wikipedia.org/wiki/Density_of_air#Altitude

            for i in range(0, n):
                    constraints.extend([
                        SignomialEquality(T_sl, T_atm[i] + L_atm*h[i])
                        ])

        Model.__init__(self, sum(T_atm), constraints, **kwargs)
        
if __name__ == "__main__":
    M = Atmosphere(4)
    sol = M.localsolve("mosek",iteration_limit=500)
    print sol.table()
