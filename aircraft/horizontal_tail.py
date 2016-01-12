import numpy as np
from gpkit.shortcuts import Var, Model
from gpkit import units

class HorizontalTail(Model):

    def setup(self):

        Cm_a = Var("Cm_{a}","-","Aircraft Pitching Moment")
        Cl_aw  = Var("Cl_{aw}",2*np.pi,"-","Wing Lift Curve Slope")
        cg_ac_diff=Var("cg_{ac}_{diff}",2,"m","Difference in CG and Aerodynamic Center Locations")
        c = Var("c","m",.5,"Mean Aerodynamic Chord")
        Cm_af = Var("Cm_{af}",0,"-","Fuselage Moment Coefficient")
        eta = Var("eta",1,"-","Tail Efficiency Factor")
        V_h=Var("V_{h}",'-',"Tail Volume Coefficient")
        Cl_at = Var("Cl_{at}","-","Tail Lift Curve Slope")
        S_H = Var("S_{H}","m^2","Tail Area")
        l_t = Var("l_{t}","m","Distance from Tail Aerodynamic Center to CG")
        S_w = Var("S_{w}",8,"m^2","Wing Area")
        AR = Var("AR",10,"-","Wing Aspect Ratio")
        c_t= Var("c_{t}",1,"m","Horizontal Tail Chord")
        S_Hmax = Var("S_{Hmax}",44.64,"ft^2","Max Horizontal Tail Area")
        l_tmax = Var("l_{tmax}",5,"m","Max Distance from Tail Aerodynamic Center to CG")
        de_da = Var("de_{da}",.6,"-","Downwash Term")

        constraints = [Cl_at*(np.pi*(S_H/(c_t)**2)+2*np.pi)==np.pi*(S_H/(c_t**2))*2*np.pi,
                        l_t<=l_tmax,
                        S_H<=S_Hmax]

        objective=Cl_aw*(cg_ac_diff/c)+Cm_af-(eta*((S_H*l_t)/(S_w*c))*Cl_at*(de_da))

        return objective, constraints

    def test(self):
        _=self.solve()

if __name__== "__main__":
    HorizontalTail().test()

		