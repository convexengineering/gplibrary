"Electric motor model "
from numpy import pi
from gpkit import Model, parse_variables, SignomialsEnabled, SignomialEquality, units
from gpkitmodels.GP.aircraft.prop.propeller import Actuator_Propeller, One_Element_Propeller
from gpkitmodels import g
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState
import matplotlib.pyplot as plt
from gpkit.constraints.tight import Tight as TCS


    
class ElecMotor_Performance(Model):
    """ Electric Motor Performance Model

    Variables
    ---------
    Pshaft                  [kW]            motor output shaft power
    Pelec                   [kW]            motor input shaft power
    etam                    [-]             motor efficiency
    Q                       [N*m]           torque
    omega                   [rpm]         propeller rotation rate 
    i                       [amps]          current
    v                       [V]             woltage
    i0         4.5           [amps]          zero-load current
    Kv         29            [rpm/V]       motor voltage constant
    R          .033             [ohms]          internal resistance
    """
    def setup(self,parent,  state):
        exec parse_variables(ElecMotor_Performance.__doc__)


        constraints = [Pshaft == Q*omega,
                Pelec == v*i,
                etam == Pshaft/Pelec, 
                parent.Qmax >= Q,
                v <= parent.V_max,
                TCS([i >= Q*Kv+i0]),
                TCS([v >= omega/Kv + i*R])
                ]
        return constraints

class ElecMotor(Model):
    """ Electric Motor Model

    Variables
    ---------
    Qstar       1           [kg/(N*m)]         motor specific torque
    W                       [lbf]              motor weight
    Qmax                    [N*m]              motor max. torque
    V_max       300         [V]                motor max voltage

    """

    flight_model = ElecMotor_Performance

    def setup(self):
        exec parse_variables(ElecMotor.__doc__)

        constraints = [W >= Qstar*Qmax*g]

        return constraints

class Propulsor_Performance(Model):
    """Propulsor Performance Model
    Variables
    ---------
    V_static       1                 [cm/s]              inflow velocity for static thrust


    """

    

    def setup(self, parent,state):
        exec parse_variables(Propulsor_Performance.__doc__)
        self.prop    = parent.prop.flight_model(parent.prop,state)
        self.motor   = parent.motor.flight_model(parent.motor,state)
        #self.stat_FS = FlightState()
        #self.stat_FS.substitutions[V_static]
        #self.stat_prop = parent.prop.flight_model(stat_FS)

        self.components = [self.prop, self.motor]

        constraints = [self.prop.Q == self.motor.Q,
                        self.prop.omega == self.motor.omega
                        ]

        return constraints, self.components#, self.stat_FS

class Propulsor(Model):
    """Propulsor model

    Variables
    ---------
    W                       [lbf]              propulsor weight


    """
    flight_model = Propulsor_Performance

    def setup(self):
        exec parse_variables(Propulsor.__doc__)

        #self.prop = Actuator_Propeller()
        self.prop = One_Element_Propeller()
        self.motor = ElecMotor()

        components = [self.prop, self.motor]

        return [self.W >= self.prop.W + self.motor.W], components

class Actuator_Propulsor(Model):
    """Propulsor model

    Variables
    ---------
    W                       [lbf]              propulsor weight

    """
    flight_model = Propulsor_Performance

    def setup(self):
        exec parse_variables(Propulsor.__doc__)

        #self.prop = Actuator_Propeller()
        self.prop = Actuator_Propeller()
        self.motor = ElecMotor()

        components = [self.prop, self.motor]

        return [self.W >= self.prop.W + self.motor.W], components

    
