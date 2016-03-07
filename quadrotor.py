"""implement quadrotor model"""
from gpkit import Variable, Model, units
from numpy import pi
from troposphere import Troposphere


L = 0.0065
M = 0.0289644
R = 8.31447


class Quadrotor(Model):
    """A vtol aircraft with four rotors"""
    def setup(self):    # pylint: disable=too-many-locals
        """This method should return objective and constraints"""
        constraints = []

        # Weights
        W = Variable("W", "lbf", "takoff weight")
        W_batt = Variable("W_{batt}", "lbf", "battery weight")
        W_payload = Variable("W_{payload}", 1, "lbf", "Weight of payload")
        W_struc = Variable("W_{struc}", "lbf", "Weight of vehicle structure")
        f_struc = Variable("f_{struc}", 0.25, "-",
                           "structure weight fraction")
        g = Variable("g", 9.8, "m/s^2", "gravitational constant")
        constraints.extend([W >= W_batt + W_payload + W_struc,
                            W_struc >= f_struc*W])

        # propulsion
        h_batt = Variable("h_{batt}", .554111, "MJ/kg",
                          "Battery energy density")
        t = Variable("t", "hr", "flight time requirement")
        P = Variable("P", "W", "power")
        n_prop = Variable("n_{prop}", 4, "-", "Number of propellers")
        rho = Variable("\\rho", "kg/m^3", "Density")
        A_prop = Variable("A_{prop}", "m^2", "Area of disc formed by 1 prop")
        d_prop = Variable("d_{prop}", 1, "ft", "propeller diameter")

        constraints.extend([
            W_batt >= P*t/h_batt*g,
            P >= n_prop*(W/n_prop)**(3/2.)/(2*rho*A_prop)**0.5,
            A_prop == 0.25*pi*d_prop**2])

        # atmosphere
        th = (g.value.magnitude*M)/(R*L)  # dimensionless
        h = Variable("h", 1000, "m", "Altitude")
        T = Variable("T", "K", "Temperature")
        Latm = Variable("L", L, "K/m", "Temperature lapse rate")
        M_over_r = Variable("(M/R)", M/R, "kg*K/J", "Air property")
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        constraints.extend([  # Model only valid up to top of the troposphere
            h <= 20000*units.m,
            # Temperature decreases with altitude at a rate of L
            T_sl >= T + Latm*h,
            # Pressure-altitude relation and ideal gas law
            rho <= p_sl*T**(th-1)*M_over_r/(T_sl**th)])


        return 1/t, constraints


if __name__ == "__main__":
    Q = Quadrotor()
    Q.solve()

#  # Constants
#  V_cell = Variable("V_{cell}", 3.7, "V", "Voltage per cell")
#
#  # Free Variableiables
#  i = Variable("i","A","Amps drawn from battery")
#  Vb = Variable("V_{b}","V","Battery voltage")
