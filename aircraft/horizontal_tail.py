import numpy as np
from gpkit.shortcuts import Var, Model
from gpkit import units

#model to find size of horizontal tail which provides greates pitch stability
#in cruise. Assumes the tail and main wing lift curve slope is 2pi
#(thin airofil theory)

class HorizontalTail(Model):

    def setup(self):
        # define variables outside class VAR
        cl_aw=2*np.pi          #wing lift curve slope
        xcg=3                  #cg location in m from nose of ac
        xac=5                  #aerodynamic center in m from nose
        chord=.5               #main wing mean aerodynamic chord in m
        cm_af=0                #fuselage moment coefficient
        Eta=1                  #tail efficienty factor
        s_w=8                  #wing area in m^2
        ar=10                  #wing aspect ratio
        chord_tmax=100              #tail chord
        shmax=100            #max allowed tail area
        ltmax=100                #max allowed distace from tail aerodynamic center to CG
        deda=.6                #downwash term
        a=-(cl_aw*((xcg/chord)-(xac/chord))+cm_af) #dummy vairable to absorb all constants
        b=Eta*(1-deda)
        btmax=100;               #max horizontal tail span in m
        #define constant variables of class VAR
        
        Cl_aw  = Var("Cl_{aw}",cl_aw,"-","Wing Lift Curve Slope")
        cg=Var("cg",xcg,"m","CG location measured from the nose")
        ac=Var("ac",xac,"m","Wing Aerodynamic Center Location Measured from Nose")
        c = Var("c","m",chord,"Mean Aerodynamic Chord")
        Cm_af = Var("Cm_{af}",cm_af,"-","Fuselage Moment Coefficient")
        eta = Var("eta",Eta,"-","Tail Efficiency Factor")
        Cl_at = Var("Cl_{at}",cl_aw,"-","Tail Lift Curve Slope")
        S_w = Var("S_{w}",s_w,"m^2","Wing Area")
        AR = Var("AR",ar,"-","Wing Aspect Ratio")
        c_tmax= Var("c_{tmax}",chord_tmax,"m","Horizontal Tail Chord Max Value")
        S_Hmax = Var("S_{Hmax}",shmax,"m^2","Max Horizontal Tail Area")
        l_tmax = Var("l_{tmax}",ltmax,"m","Max Distance from Tail Aerodynamic Center to CG")
        b_tmax = Var("b_{tmax}",btmax,"m","Max Horizontal Tail Span")
        de_da = Var("de_{da}",deda,"-","Downwash Term")
        A=Var("A",a,"-","Dummy Variable to Absorb all Constants")
        B=Var("B",b,"-","Dummy Variable to Absorb all Constants")
        
        #define variables of class VAR to be optomized
        b_t=Var("b_{t}","m","Horizontal Tail Span")
        Cm_a = Var("Cm_{a}","-","Aircraft Pitching Moment times -1")
        V_h=Var("V_{h}",'-',"Tail Volume Coefficient")
        S_H = Var("S_{H}","m^2","Tail Area")
        l_t = Var("l_{t}","m","Distance from Tail Aerodynamic Center to CG")
        c_t=Var("c_{tmax}","m","Horizontal Tail Chord")
        AR_H=Var("AR_{H}","-","Horizontal Tail Aspect Ratio")
        C=Var("C","-","pi*AR_H+2pi")
        #define constraint
        constraints = [(A+B*((S_H*l_t)/(S_w*c))*(2*np.pi*np.pi*(b_t**2/S_H))/C)/Cm_a<=1,
                        l_t/l_tmax<=1,
                        S_H/S_Hmax<=1,
                        c_t/c_tmax<=1,
                        b_t/b_tmax<=1,
                        (np.pi*(b_t**2/S_H)+2*np.pi)/C<=1]
        
        #define objective function
        objective=1/Cm_a
        
        return objective, constraints

    def test(self):
        _=self.solve()

if __name__== "__main__":
    HorizontalTail().test()

		