"Electric motor model "
from gpkit import Model, parse_variables
from gpkit.constraints.tight import Tight as TCS
from gpkitmodels.GP.aircraft.prop.propeller import Propeller
from gpkitmodels import g

class MotorPerf(Model):
    """ Electric Motor Performance Model

    Variables
    ---------
    Pshaft                  [kW]            motor output shaft power
    Pelec                   [kW]            motor input shaft power
    etam                    [-]             motor efficiency
    Q                       [N*m]           torque
    omega                   [rpm]           propeller rotation rate
    i                       [amps]          current
    v                       [V]             woltage
    i0         4.5          [amps]          zero-load current
    Kv         29           [rpm/V]         motor voltage constant
    R          .033         [ohms]          internal resistance
    """
    def setup(self, static, state):
        exec parse_variables(MotorPerf.__doc__)

        return [Pshaft == Q*omega,
                Pelec == v*i,
                etam == Pshaft/Pelec,
                static.Qmax >= Q,
                v <= static.V_max,
                TCS([i >= Q*Kv+i0]),
                TCS([v >= omega/Kv + i*R])
               ]

class Motor(Model):
    """ Electric Motor Model

    Variables
    ---------
    Qstar       1           [kg/(N*m)]         motor specific torque
    W                       [lbf]              motor weight
    Qmax                    [N*m]              motor max. torque
    V_max       300         [V]                motor max voltage

    """

    flight_model = MotorPerf

    def setup(self):
        exec parse_variables(Motor.__doc__)

        constraints = [W >= Qstar*Qmax*g]

        return constraints

class PropulsorPerf(Model):
    """Propulsor Performance Model

    Variables
    ---------
    V_static       1      [cm/s]           inflow velocity for static thrust

    """

    def setup(self, static, state):
        exec parse_variables(PropulsorPerf.__doc__)
        self.prop = parent.prop.flight_model(static.prop, state)
        self.motor = parent.motor.flight_model(static.motor, state)

        self.components = [self.prop, self.motor]

        constraints = [self.prop.Q == self.motor.Q,
                       self.prop.omega == self.motor.omega
                      ]

        return constraints, self.components

class Propulsor(Model):
    """Propulsor model

    Variables
    ---------
    W                       [lbf]              propulsor weight

    """
    flight_model = PropulsorPerf

    def setup(self):
        exec parse_variables(Propulsor.__doc__)

        self.prop = Propeller()
        self.motor = Motor()

        components = [self.prop, self.motor]

        return [self.W >= self.prop.W + self.motor.W], components


