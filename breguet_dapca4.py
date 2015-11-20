from gpkit.shortcuts import Var, Model
from gpkit import units
from dapca4cost import DAPCA4Cost
from breguet_range import Breguet_Range
import gpkit
    
class breguet_dapca(Model):
    """
    This is a test model just to see how I can get dapca4 and breguet_range to interact
	"""

	def setup(self):
		
		#Constants
		W_e = Var('W_{e}', 9000, "lb", "Empty Weight")
		V = Var('V', 420, "knots", "Max Velocity")
		M = Var('M', 0.68, "-", "Max Mach Number")

		#Breguet Range Model
		r = Breguet_Range(W_oew=W_e, M_max=M)

		#Dapca4 Model
		d = DAPCA4Cost(W_e=W_e, V=V)

		combined = r & d
		
		return combined


if __name__ == "__main__":
	m = breguet_dapca()
	sol = m.solve()
		

