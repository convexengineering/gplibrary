from gpkit import Variable, Model
import numpy as np

class Climb(Model):
    """
    Model for climbing flight

    References
    ---------
    [1] Anderson, Introduction to flight
    [2] https://en.wikipedia.org/wiki/Taylor_series
    [3] https://en.wikipedia.org/wiki/Thrust_specific_fuel_consumption

    Assumptions
    -----------
    - Assumes you want to minimize climb time
    - Assumes constant thrust
    - Assumes constant velocity
    - Assumes constant air density
    """
    def setup(self):

        H = 30000 # final altitude
        n = 3 # number of steps

        h = np.linspace(0, H, n)
        T0 = 288.15 # [K]
        L = 0.0065 # [K/m]
        p0 = 101325 # [Pa]
        R = 8.31 # [J/(mol*K)]
        g = 9.81 # [m/s^2]
        M = 0.02896 # [kg/mol]

        T = T0 - L*h # [K]
        p = p0*(1 - L*h/T0)**(g*M/(R*L))# [Pa]
        rho = p*M/(R*T) # [kg/m^3]

        # Free variables
        dt  = Variable('\\Delta t', 's', 'Time to climb segment')
#        t   = Variable('t', 's', 'Time to climb')
        D   = Variable('D', 'N', 'Drag')
        m_f = VectorVariable('m_f', n, 'kg', 'Mass of fuel burned')

        # Fixed parameters
        CD   = Variable('C_D', 0.02, '-', 'Drag coefficient')
        dh   = Variable('\\Delta h', h/n, 'ft', 'Altitude step')
#        h    = VectorVariable('h', h, 'm', 'Final altitude')
        rho  = VectorVariable('\\rho', rho, 'kg/m^3', 'Air density')
        S    = Variable('S', 130, 'm^2', 'Reference area')
        T    = Variable('T', 1E6, 'N', 'Thrust')
        TSFC = Variable('TSFC', 8.7E-4, 'g/(kN*s)', # [3]
                        'Thrust specific fuel consumption')
        V    = Variable('V', 100, 'm/s', 'Freestream velocity')
        W    = Variable('W', 600000, 'N', 'Aircraft weight')

        objective = m_f

        constraints = [# Climb time [1]
                       # uses Maclaurin series expansion of 1/(1-T/D) (ref [2])
                       dt >= (W/(V*T))*(1 + D/T + (D/T)**2 + (D/T)**3)*dh,

                       # Drag
                       D == 0.5*rho*V**2*S*CD,

                       # Fuel burn
                       m_f == TSFC*T*dt,
                      ]

        return objective, constraints

    def test(self):
        self.solve()

if __name__ == "__main__":
    Climb().test()
