from gpkit import Variable, Model

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

        h = 30000 # final altitude
        n = 1 # number of steps

        # Free variables
        t   = Variable('t', 's', 'Time to climb')
        D   = Variable('D', 'N', 'Drag')
        m_f = Variable('m_f', 'kg', 'Mass of fuel burned')

        # Fixed parameters
        CD   = Variable('C_D', 0.02, '-', 'Drag coefficient')
        dh   = Variable('\\Delta h', h/n, 'ft', 'Altitude step')
        h    = Variable('h', h, 'ft', 'Final altitude')
        rho  = Variable('\\rho', 1.225, 'kg/m^3', 'Air density')
        S    = Variable('S', 130, 'm^2', 'Reference area')
        T    = Variable('T', 1E6, 'N', 'Thrust')
        TSFC = Variable('TSFC', 8.7E-4, 'g/(kN*s)', # [3]
                        'Thrust specific fuel consumption')
        V    = Variable('V', 100, 'm/s', 'Freestream velocity')
        W    = Variable('W', 600000, 'N', 'Aircraft weight')

        objective = m_f 

        constraints = [# Climb time [1]
                       # uses Maclaurin series expansion of 1/(1-T/D) (ref [2])
                       t >= (W/(V*T))*(1 + D/T + (D/T)**2 + (D/T)**3)*dh,

                       # Drag
                       D == 0.5*rho*V**2*S*CD,

                       # Fuel burn
                       m_f == TSFC*T*t,
                      ]

        return objective, constraints

    def test(self):
        self.solve()

if __name__ == "__main__":
    Climb().test()
