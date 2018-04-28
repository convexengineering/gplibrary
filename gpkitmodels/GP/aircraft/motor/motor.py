"Electric motor model "
from gpkit import Model, parse_variables
from gpkit.constraints.tight import Tight as TCS
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp
from gpkitmodels import g

class MotorPerf(Model):
    """ Electric Motor Performance Model

    Note: The last two constraints may not be tight if there is a motor energy constraint that does not size the motor.
    TCS removed to prevent unnecessary warnings - for most normal useage they will be tight. 

    Variables
    ---------
    Pshaft                  [kW]            motor output shaft power
    Pelec                   [kW]            motor input shaft power
    etam                    [-]             motor efficiency
    Q                       [N*m]           torque
    omega                   [rpm]           propeller rotation rate
    i                       [amps]          current
    v                       [V]             woltage
    """
    def setup(self, static, state):
        exec parse_variables(MotorPerf.__doc__)

        Kv = static.Kv
        R  = static.R
        i0 = static.i0
        V_max = static.V_max

        return [Pshaft == Q*omega,
                Pelec == v*i,
                etam == Pshaft/Pelec,
                static.Qmax >= Q,
                v <= V_max,
                i >= Q*Kv+i0,
                v >= omega/Kv + i*R]

class Motor(Model):
    """ Electric Motor Model

    Variables
    ---------
    Qstar       .5          [kg/(N*m)]         motor specific torque
    W                       [lbf]              motor weight
    Qmax                    [N*m]              motor max. torque
    V_max       300         [V]                motor max voltage
    Kv_min     1            [rpm/V]         min motor voltage constant
    Kv_max     1000         [rpm/V]         max motor voltage constant
    Kv                      [rpm/V]         motor voltage constant
    i0         4.5          [amps]          zero-load current 
    R          .033         [ohms]          internal resistance
    """

    flight_model = MotorPerf

    def setup(self):
        exec parse_variables(Motor.__doc__)

        constraints = [W >= Qstar*Qmax*g,
                       Kv >= Kv_min,
                       Kv <= Kv_max]

        return constraints

class PropulsorPerf(Model):
    """Propulsor Performance Model

    """

    def setup(self, static, state):
        exec parse_variables(PropulsorPerf.__doc__)
        self.prop = static.prop.flight_model(static.prop, state)
        self.motor = static.motor.flight_model(static.motor, state)

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
    prop_flight_model = ActuatorProp


    def setup(self):
        exec parse_variables(Propulsor.__doc__)

        Propeller.flight_model = self.prop_flight_model
        self.prop = Propeller()
        self.motor = Motor()

        components = [self.prop, self.motor]

        return [self.W >= self.prop.W + self.motor.W], components




