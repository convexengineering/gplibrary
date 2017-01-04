"Implements Fuselage model"
import numpy as np
from gpkit import Variable, Model, SignomialsEnabled, units, SignomialEquality, ConstraintSet, LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS

from gpkit.small_scripts import mag
from collections import defaultdict

class Fuselage(Model):
    """
    Fuselage model
    """
    def __init__(self, wingbox, htail, vtail, fuselage_type = 'narrowbody', version = 'Philip'):

        if fuselage_type == 'widebody' and version == 'Philip':
            raise NameError('Widebody frame and Philip version are not compatable')


        Afuse    = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        Ahold    = Variable('A_{hold}', 'm^2', 'Cargo hold x-sectional area')
        Askin    = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        Dfuse    = Variable('D_{fuse}', 'N', 'Total drag in cruise')
        Dfrict   = Variable('D_{friction}', 'N', 'Friction drag')
        Dupswp   = Variable('D_{upsweep}', 'N', 'Drag due to fuse upsweep')
        FF       = Variable('FF', '-','Fuselage form factor')
        LF       = Variable('LF', '-', 'Load factor')
        Lvmax    = Variable('L_{v_{max}}', 'N', 'Max vertical tail load')
        Nland    = Variable('N_{land}', '-', 'Emergency landing load factor')
        Qv       = Variable('Q_v', 'N*m', 'Torsion moment imparted by tail')
        R        = Variable('R', 287, 'J/(kg*K)', 'Universal gas constant')
        Rfuse    = Variable('R_{fuse}', 'm', 'Fuselage radius')
        SPR      = Variable('SPR', '-', 'Number of seats per row')
        Sbulk    = Variable('S_{bulk}', 'm^2', 'Bulkhead surface area')
        Spassfl  = Variable('S_{pass_{floor}}', 'N', 'Maximum shear in pasenger floor beams')
        Scargofl = Variable('S_{cargo_{floor}}', 'N', 'Maximum shear in cargo floor beams')
        Sfloor   = Variable('S_{floor}', 'N', 'Maximum shear in floor beams')
        Snose    = Variable('S_{nose}', 'm^2', 'Nose surface area')
        Tcabin   = Variable('T_{cabin}', 'K', 'Cabin temperature')
        Vbulk    = Variable('V_{bulk}', 'm^3', 'Bulkhead skin volume')
        Vfuse    = Variable('V_{fuse}', 'm^3', 'Fuse volume')
        Vcabin   = Variable('V_{cabin}', 'm^3', 'Cabin volume')
        Vcargo   = Variable('V_{cargo}', 'm^3', 'Cargo volume')
        Vcone    = Variable('V_{cone}', 'm^3', 'Cone skin volume')
        Vcyl     = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')
        Vfloor   = Variable('V_{floor}', 'm^3', 'Floor volume')
        Vhold    = Variable('V_{hold}', 'm^3', 'Hold volume')
        #Vinf     = Variable('V_{\\infty}', 'm/s', 'Cruise velocity')
        Vlugg    = Variable('V_{lugg}', 'm^3', 'Luggage volume')
        Vnose    = Variable('V_{nose}', 'm^3', 'Nose skin volume')
        Wapu     = Variable('W_{apu}', 'N', 'APU weight')
        Wavgpass = Variable('W_{avg. pass}', 'lbf', 'Average passenger weight')
        Wbuoy    = Variable('W_{buoy}', 'N', 'Buoyancy weight')
        Wcargo   = Variable('W_{cargo}', 'N', 'Cargo weight')
        Wcarryon = Variable('W_{carry on}', 'lbf', 'Ave. carry-on weight')
        Wchecked = Variable('W_{checked}', 'lbf', 'Ave. checked bag weight')
        Wcone    = Variable('W_{cone}', 'N', 'Cone weight')
        Wfix     = Variable('W_{fix}', 'lbf',
                            'Fixed weights (pilots, cockpit seats, navcom)')
        Wfloor   = Variable('W_{floor}', 'N', 'Floor weight')
        Wfuse    = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Winsul   = Variable('W_{insul}', 'N', 'Insulation material weight')
        Wlugg    = Variable('W_{lugg}', 'N', 'Passenger luggage weight')
        Wpadd    = Variable('W_{padd}', 'N',
                            'Misc weights (galley, toilets, doors etc.)')
        Wpass    = Variable('W_{pass}', 'N', 'Passenger weight')
        Wpay     = Variable('W_{pay}', 'N', 'Payload weight')
        Wppfloor = Variable('W\'\'_{floor}', 'N/m^2',
                            'Floor weight/area density')
        Wppinsul = Variable('W\'\'_{insul}', 'N/m^2',
                            'Weight/area density of insulation material')
        Wpseat   = Variable('W\'_{seat}', 'N', 'Weight per seat')
        Wpwindow = Variable('W\'_{window}', 'N/m',
                            'Weight/length density of windows')
        Wseat    = Variable('W_{seat}', 'N', 'Seating weight')
        Wshell   = Variable('W_{shell}', 'N', 'Shell weight')
        Wskin    = Variable('W_{skin}', 'N', 'Skin weight')
        Wwindow  = Variable('W_{window}', 'N', 'Window weight')
        bvt      = Variable('b_{vt}', 'm', 'Vertical tail span')
        cvt      = Variable('c_{vt}', 'm', 'Vertical tail root chord')
        dh       = Variable('\\Delta h', 'm',
                            'Distance from floor to widest part of fuselage')
        dp       = Variable('\\Delta p', 'Pa',
                            'Pressure difference across fuselage skin')
        f        = Variable('f', '-', 'Fineness ratio')
        fapu     = Variable('f_{apu}', '-',
                            'APU weight as fraction of payload weight')
        ffadd    = Variable('f_{fadd}', '-',
                            'Fractional added weight of local reinforcements')
        fframe   = Variable('f_{frame}', '-', 'Fractional frame weight')
        flugg1   = Variable('f_{lugg,1}', '-',
                            'Proportion of passengers with one suitcase')
        flugg2   = Variable('f_{lugg,2}', '-',
                            'Proportion of passengers with two suitcases')
        fpadd    = Variable('f_{padd}', '-',
                            'Other misc weight as fraction of payload weight')
        fseat    = Variable('f_{seat}', '-',
                            'Seat weight as fraction of payload weight')
        fstring  = Variable('f_{string}', '-','Fractional weight of stringers')
        g        = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        
        hhold    = Variable('h_{hold}', 'm', 'Height of the cargo hold')
        lamcone  = Variable('\\lambda_{cone}', '-',
                            'Tailcone radius taper ratio (xshell2->xtail)')
        lcone    = Variable('l_{cone}', 'm', 'Cone length')
        lfloor   = Variable('l_{floor}', 'm', 'Floor length')
        lfuse    = Variable('l_{fuse}', 'm', 'Fuselage length')
        lnose    = Variable('l_{nose}', 'm', 'Nose length')
        lshell   = Variable('l_{shell}', 'm', 'Shell length')
        #mu       = Variable('\\mu', 'N*s/m^2', 'Dynamic viscosity (35,000 ft)')
        npass    = Variable('n_{pass}', '-', 'Number of passengers')
        nrows    = Variable('n_{rows}', '-', 'Number of rows')
        nseat    = Variable('n_{seat}', '-',' Number of seats')
        pcabin   = Variable('p_{cabin}', 'Pa','Cabin air pressure (8,000ft)')
        phi      = Variable('\\phi', '-', 'Upsweep angle')
        pitch    = Variable('p_s', 'in', 'Seat pitch')
        plamv    = Variable('p_{\\lambda_v}', '-', '1 + 2*Tail taper ratio')
        rE       = Variable('r_E', '-', 'Ratio of stringer/skin moduli')
        rhobend  = Variable('\\rho_{bend}', 'kg/m^3', 'Stringer density')
        rhocabin = Variable('\\rho_{cabin}', 'kg/m^3', 'Air density in cabin')
        rhocargo = Variable('\\rho_{cargo}', 'kg/m^3', 'Cargo density')
        rhocone  = Variable('\\rho_{cone}', 'kg/m^3',
                            'Cone material density')
        rhofloor = Variable('\\rho_{floor}', 'kg/m^3',
                            'Floor material density')
        rhoinf   = Variable('\\rho_{\\infty}', 'kg/m^3',
                            'Air density (35,000ft)')
        rholugg  = Variable('\\rho_{lugg}', 'kg/m^3', 'Luggage density')
        rhoskin  = Variable('\\rho_{skin}', 'kg/m^3', 'Skin density')
        sigfloor = Variable('\\sigma_{floor}', 'N/m^2',
                            'Max allowable cap stress')
        sigskin  = Variable('\\sigma_{skin}', 'N/m^2',
                            'Max allowable skin stress')
        sigth    = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx     = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')
        taucone  = Variable('\\tau_{cone}', 'N/m^2', 'Shear stress in cone')
        taufloor = Variable('\\tau_{floor}', 'N/m^2',
                            'Max allowable shear web stress')
        tcone    = Variable('t_{cone}', 'm', 'Cone thickness')
        tshell   = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin    = Variable('t_{skin}', 'm', 'Skin thickness')
        waisle   = Variable('w_{aisle}', 'm', 'Aisle width')
        
        wfuse    = Variable('w_{fuse}', 'm', 'Fuselage width')
        wseat    = Variable('w_{seat}', 'm', 'Seat width')
        wsys     = Variable('w_{sys}', 'm',
                            'Width between cabin and skin for systems')
        xCGfu    = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        xVbulk   = Variable('xVbulk', 'm^4', 'Volume moment of bulkhead')
        xVcyl    = Variable('xVcyl', 'm^4', 'Volume moment of cylinder')
        xVnose   = Variable('xVnose', 'm^4', 'Volume moment of nose')
        xWapu    = Variable('xWapu', 'N*m', 'Moment of APU')
        xWcone   = Variable('xWcone', 'N*m', 'Moment of cone')
        xWfix    = Variable('xWfix', 'N*m', 'Moment of fixed weights')
        xWfloor  = Variable('xWfloor', 'N*m', 'Moment of floor weight')
        xWfuse   = Variable('xWfuse', 'N*m', 'Fuselage moment')
        xWinsul  = Variable('xWinsul', 'N*m', 'Moment of insulation material')
        xWpadd   = Variable('xWpadd', 'N*m', 'Moment of misc weights')
        xWseat   = Variable('xWseat', 'N*m', 'Moment of seats')
        xWshell  = Variable('xWshell', 'N*m', 'Mass moment of shell')
        xWskin   = Variable('xWskin', 'N*m', 'Mass moment of skin')
        xWwindow = Variable('xWwindow', 'N*m', 'Mass moment of windows')
        x_upswp  = Variable('x_{up}', 'm', 'Fuselage upsweep point')
        xapu     = Variable('xapu', 'ft', 'x-location of APU')
        xconend  = Variable('xconend', 'm', 'x-location of cone end')
        xfix     = Variable('xfix', 'm', 'x-location of fixed weight')
        xshell1  = Variable('x_{shell1}', 'm', 'Start of cylinder section')
        xshell2  = Variable('x_{shell2}', 'm', 'End of cylinder section')




        # Bending inertias (ported from TASOPT)
        A0           = Variable('A0','m^2','Horizontal bending area constant A0') #(shell inertia contribution)
        A1           = Variable('A1','m','Horizontal bending area constant A1') #(tail impact + aero loading)
        A2           = Variable('A2','-','Horizontal bending area constant A2') #(fuselage impact)
        Ahbendb      = Variable('A_{hbendb}','m^2','Horizontal bending area at rear wingbox')
        Ahbendf      = Variable('A_{hbendf}','m^2','Horizontal bending area at front wingbox')
        Avbendb      = Variable('A_{vbendb}','m^2','Vertical bending material area at rear wingbox')
        #B0           = Variable('B0','m^2','Vertical bending area constant B0') #(shell inertia contribution)
        #B1           = Variable('B1','m','Vertical bending area constant B1') #(vertical tail bending load)
        Ihshell      = Variable('I_{hshell}','m^4','Shell horizontal bending inertia')
        #Ivshell      = Variable('I_{vshell}','m^4','Shell vertical bending inertia')
        rMh          = Variable('r_{M_h}',.4,'-','Horizontal inertial relief factor') #[TAS]
        rMv          = Variable('r_{M_v}',.7,'-','Vertical inertial relief factor') #[TAS]
        sigbend      = Variable('\\sigma_{bend}','N/m^2','Bending material stress')
        sigMh        = Variable('\\sigma_{M_h}','N/m^2','Horizontal bending material stress')
        #sigMv        = Variable('\\sigma_{M_v}','N/m^2','Vertical bending material stress')
        Vhbend       = Variable('V_{hbend}','m^3','Horizontal bending material volume')
        Vhbendb      = Variable('V_{hbendb}','m^3','Horizontal bending material volume b') #back fuselage
        Vhbendc      = Variable('V_{hbendc}','m^3','Horizontal bending material volume c') #center fuselage
        Vhbendf      = Variable('V_{hbendf}','m^3','Horizontal bending material volume f') #front fuselage
        #Vvbend       = Variable('V_{vbend}','m^3','Vertical bending material volume')
        #Vvbendb      = Variable('V_{vbendb}','m^3','Vertical bending material volume b') #back fuselage
        #Vvbendc      = Variable('V_{vbendc}','m^3','Vertical bending material volume c') #center fuselage
        Whbend       = Variable('W_{hbend}','N','Horizontal bending material weight')
        #Wvbend       = Variable('W_{vbend}','N','Vertical bending material weight')
        xhbend       = Variable('x_{hbend}','m','Horizontal zero bending location')
        #xvbend       = Variable('x_{vbend}','m','Vertical zero bending location')
        # x-location variables
        xshell1      = Variable('x_{shell1}', 'm', 'Start of cylinder section')
        xshell2      = Variable('x_{shell2}', 'm', 'End of cylinder section')
        xtail        = Variable('x_{tail}','m', 'x-location of tail')

        with SignomialsEnabled():

            if version == 'Philip':

                Afloor   = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
                Mfloor   = Variable('M_{floor}', 'N*m',
                            'Max bending moment in floor beams')
                wfloor   = Variable('w_{floor}', 'm', 'Half floor width')
                Pfloor   = Variable('P_{floor}', 'N', 'Distributed floor load')
                hfloor   = Variable('h_{floor}', 'm', 'Floor beam height')

                constraints = [
                            Pfloor >= Nland*(Wpass + Wseat + Wcarryon + Wlugg + Wcargo),
                            Sfloor == 0.5*Pfloor, # without support
                            
                            
                            Afloor >= 2*Mfloor/(sigfloor*hfloor)
                                      + 1.5*Sfloor/taufloor,

                            Vfloor >= 2*wfloor*Afloor,

                            
                            (wfloor)**2 + dh**2 >= Rfuse**2, # [SP]
                            Rfuse >= dh + hfloor + hhold,
                            Wfloor >= rhofloor*g*Vfloor
                                      + 2*wfloor*lfloor*Wppfloor,
                            Wshell >= Wskin*(1 + fstring + fframe + ffadd),

                            # Cone
                            Rfuse*taucone*(1+plamv)*Vcone*(1+lamcone)/(4*lcone)
                                >= Lvmax*bvt*plamv/3, # [SP]

                            # Total fuselage weight
                            Wfuse >= Wfix + Wapu + Wpadd + Wseat + Wshell
                                   + Wwindow + Winsul + Wcone + Wfloor + Wbuoy,
                            
                            ]

            elif version == 'NextGen':

                Apassfl  = Variable('A_{pass_{floor}}', 'm^2', 'Passenger floor beam x-sectional area')
                Acargofl = Variable('A_{cargo_{floor}}', 'm^2', 'Cargo floor beam x-sectional area')
                Mpassfl  = Variable('M_{pass_{floor}}', 'N*m',
                                    'Max bending moment in passenger floor beams')
                Mcargofl = Variable('M_{cargo_{floor}}', 'N*m',
                                    'Max bending moment in cargo floor beams')
                wpassfl  = Variable('w_{pass_{floor}}', 'm', 'Passenger half floor width')
                wcargofl = Variable('w_{cargo_{floor}}', 'm', 'Cargo half floor width')
                hpassfl  = Variable('h_{pass_{floor}', 'm', 'Passenger floor beam height')
                hcargofl = Variable('h_{cargo_{floor}}', 'm', 'Cargo floor beam height')
                Ppassfl  = Variable('P_{pass_{floor}}', 'N', 'Distributed passenger floor load')
                Pcargofl = Variable('P_{cargo_{floor}}', 'N', 'Distributed cargo floor load')
                hfuse    = Variable('h_{fuse}','m','Fuselage height')
                Wtail    = Variable('W_{tail}','N','Total tail weight')

                thetadb  = Variable('\\theta_{db}', 0, '-','DB fuselage joining angle')
                hdb      = Variable('h_{db}','m', 0, 'Web half-height')
                tdb      = Variable('t_{db}', 0, 'm', 'Web thickness')

                constraints = [
                            Ppassfl >= Nland*(Wpass + Wseat + Wcarryon), # might want to remove Wcarryon
                            Pcargofl >= Nland*(Wlugg + Wcargo),

                            Spassfl == 0.5*Ppassfl, # without support
                            Scargofl == 0.5*Pcargofl, # without support

                            Mcargofl == Pcargofl*wcargofl/4.,

                            Apassfl >= 2*Mpassfl/(sigfloor*hpassfl)
                                      + 1.5*Spassfl/taufloor,
                            Acargofl >= 2*Mcargofl/(sigfloor*hcargofl)
                                      + 1.5*Scargofl/taufloor,
                            Apassfl >= hpassfl*lfloor,
                            Acargofl >= hcargofl*lfloor,

                            Vfloor >= 2*wpassfl*Apassfl + 2*wcargofl*Acargofl,

                            (wcargofl)**2 + (dh + hhold + hpassfl)**2 <= Rfuse**2, # [SP]
                            Wfloor >= rhofloor*g*Vfloor
                                      + (2*wpassfl + 2*wcargofl)*lfloor*Wppfloor,

                            (wpassfl)**2 + dh**2 >= Rfuse**2, # [SP] # passenger floor located below widest part of fuselage
                            dh + hhold + hpassfl + hcargofl <= Rfuse, #There are two florrs, one for pass and one for cargo

                            Rfuse >= dh + hpassfl + hcargofl + hhold,
                            Wfloor >= rhofloor*g*Vfloor
                                      + 2*(wcargofl + wpassfl)*lfloor*Wppfloor,

                            Wtail >= vtail['W_{vtail}'] + htail['W_{htail}'] + Wcone,
                            Wshell >= Wskin*(1 + fstring + fframe + ffadd) + Whbend,

                            SignomialEquality(xtail, lnose + lshell + .5*lcone), #Temporarily
                            hfuse    == Rfuse,
                            xhbend >= wingbox['x_{wing}'],


                            # Tail cone sizing
                            Rfuse*taucone*(1+plamv)*Vcone*(1+lamcone)/(4*lcone)
                                >= Lvmax*bvt*plamv/3, # [SP]
                            #Vcone*(1+lamcone)*(np.pi+4*thetadb)>= vtail['Q_v']/taucone*(np.pi+2*thetadb)*(lcone/Rfuse)*2,

                            # Total fuselage weight
                            Wfuse >= Wfix + Wapu + Wpadd + Wseat + Wshell
                                   + Wwindow + Winsul + Wfloor + Wbuoy + Wtail,


                            Ihshell <= ((np.pi+4*thetadb)*Rfuse**2)*Rfuse*tshell + 2/3*hdb**3*tdb, # [SP]
                            sigbend == rE*sigskin,

                            # Horizontal bending material model
                            # Calculating xbend, the location where additional bending material is required
                            SignomialEquality(A0,A2*(xshell2-xhbend)**2 + A1*(xtail-xhbend)), #[SP] #[SPEquality] 
                            A2      >=  Nland*(Wpay+Wshell+Wwindow+Winsul+Wfloor+Wseat)/(2*lshell*hfuse*sigMh), # Landing loads constant A2
                            A1      >= (Nland*Wtail + rMh*htail['L_{h_{max}}'])/(hfuse*sigMh),                                # Aero loads constant A1
                            A0      == (Ihshell/(rE*hfuse**2)),                                                # Shell inertia constant A0
                            Ahbendf >= A2*(xshell2-wingbox['x_f'])**2 + A1*(xtail-wingbox['x_f']) - A0, #[SP]                           # Bending area forward of wingbox

                            Ahbendf >= 1*10**-30*units['m^2'],

                            Ahbendb >= A2*(xshell2-wingbox['x_b'])**2 + A1*(xtail-wingbox['x_b']) - A0, #[SP]                           # Bending area behind wingbox

                            Vhbendf >= A2/3*((xshell2-wingbox['x_f'])**3 - (xshell2-xhbend)**3) \
                                            + A1/2*((xtail-wingbox['x_f'])**2 - (xtail - xhbend)**2) \
                                            + A0*(xhbend-wingbox['x_f']), #[SP]

                            Vhbendf >= 1*10**-30*units['m^3'], #Temporary

                            Vhbendb >= A2/3*((xshell2-wingbox['x_b'])**3 - (xshell2-xhbend)**3) \
                                            + A1/2*((xtail-wingbox['x_b'])**2 - (xtail - xhbend)**2) \
                                            + A0*(xhbend-wingbox['x_b']), #[SP]
                            
                            Vhbendc >= .5*(Ahbendf + Ahbendb)*wingbox['c_0']*wingbox['\\bar_w'],
                            Vhbend >= Vhbendc + Vhbendf + Vhbendb,
                            Whbend  >= g*rhobend*Vhbend,
                            


                            # Temporary wing variable substitutions
                            wingbox['c_0']       == 0.1*lshell, #Temporarily
                            wingbox['dx_{wing}']   == 0.25*wingbox['c_0'], #Temporarily
                                        
                            sigMh   <= sigbend - rE*dp/2*Rfuse/tshell, # The stress available to the bending material reduced because of pressurization
                            #sigMv   == sigMh,

                            ]

            else:
                raise NotImplementedError

            constraints = [constraints,
                            # Geometry relations
                            lnose == xshell1,

                            # Cross section relations
                            wfuse == 2*Rfuse,
                            Askin >= 2*np.pi*Rfuse*tskin, # simplified
                            SignomialEquality(tshell,tskin*(1+rE*fstring*rhoskin/rhobend)),


                            # Pressure shell loads
                            sigx == dp/2*Rfuse/tshell,  
                            sigth == dp*Rfuse/tskin, 
                            sigskin >= sigth,
                            sigskin >= sigx,
                            Snose**(8./5) >= (2*np.pi*Rfuse**2)**(8./5) *
                                             (1./3 + (2./3)*(lnose/Rfuse)
                                             **(8./5)),
                            Sbulk == 2*np.pi*Rfuse**2,
                            Vcyl == Askin*lshell,
                            Vnose == Snose*tskin,
                            Vbulk == Sbulk*tskin,
                            Wskin >= rhoskin*g*(Vcyl + Vnose + Vbulk),
                            

                            # Fuselage volume and buoyancy weight
                            rhocabin == (1/(R*Tcabin))*pcabin,  
                            Vfuse >= Afuse*(lshell + 0.67*lnose + 0.67*Rfuse),
                            TCS([Wbuoy >= (rhocabin - rhoinf)*g*Vfuse], reltol=1E-3), # [SP]  

                            # Windows and insulation
                            Wwindow == Wpwindow * lshell,
                            Winsul >= Wppinsul*(1.1*np.pi*Rfuse*lshell
                                      + 0.55*(Snose + Sbulk)),

                            # Payload-proportional weights
                            Wapu == Wpay*fapu,
                            Wseat == Wpseat*nseat,
                            Wpadd == Wpay*fpadd,

                            # Floor
                            lfloor >= lshell + 2*Rfuse,
                            lshell >= nrows*pitch,

                            # Tail cone
                            
                            plamv >= 1.6,
                            taucone == sigskin,
                            lamcone == 0.4, # TODO remove
                            lamcone == cvt/lcone,
                            #Wcone == 0.2*units('N')*rhocone*units('m')*units('m')*units('m')/units('kg'),
                            Wcone >= rhocone*g*Vcone*(1 + fstring + fframe
                                                      + ffadd),


                            # Payload weight breakdown
                            npass == nseat*LF,
                            nseat == nrows*SPR,
                            Wpass == npass*Wavgpass,
                            Wlugg >= flugg2*npass*2*Wchecked
                                     + flugg1*npass*Wchecked, #Removed the carryon
                            Wlugg == Vlugg*g*rholugg,
                            Wcargo == Vcargo*g*rhocargo,
                            Vhold >= Vcargo + Vlugg,
                            Vhold <= Ahold*lshell,

                            # [SP] Harris stocker 1998 (wolfram)
                            
                            Wpay >= Wpass + Wlugg + Wcargo, 

                            # Drag
                            # sources: Raymer (p285), kfid 325 notes (p180)
                            lfuse >= lnose+lshell+lcone,
                            f == lfuse/((4/np.pi*Afuse)**0.5), # fineness ratio
                            FF >= 1 + 60/f**3 + f/400, # form factor
                            #Dfrict >= FF * np.pi*Rfuse * mu*Vinf
                            #          * 0.074*(rhoinf*Vinf*lfuse/mu)**0.8,  # PERFORMANCE MODEL

                            # Drag due to fuselage upsweep (Raymer p286)
                            xshell2 >= xshell1 + lshell,
                            xshell2 == x_upswp,
                            x_upswp + lcone <= lfuse,
                            1.13226*phi**1.03759 == Rfuse/lcone, # monomial fit
                                                                 # of tan(phi)
                            #Dupswp >= 3.83*phi**2.5*Afuse * 0.5*rhoinf*Vinf**2,  # PERFORMANCE MODEL
                            #Dfuse >= Dfrict + Dupswp  # PERFORMANCE MODEL
                          ]




        if fuselage_type == 'narrowbody':
            with SignomialsEnabled():
                constraints = [constraints,

                            Afuse >= np.pi*Rfuse**2, # simplified
                            lnose >= 5.2*units.m, # TODO less arbitrary
                            wfuse >= SPR*wseat + waisle + 2*wsys + 2*tshell,


                                ]

                if version == 'Philip':

                    constraints = [constraints,
                                Mfloor == Pfloor*wfloor/4.,
                                TCS([Ahold <= (2./3)*2*wfloor*hhold + hhold**3/(2*2*wfloor)], reltol=1E-5),
                                ]

                elif version == 'NextGen':
                    constraints = [constraints,
                                Mpassfl == Ppassfl*wpassfl/4.,
                                #TCS([Ahold <= (2./3)*2*wpassfl*hhold + hhold**3/(2*2*wpassfl)], reltol=1E-5),
                                Ahold <= wcargofl*2*hhold,
                                ]
                else:
                    raise NotImplementedError

        elif fuselage_type == 'widebody':
            with SignomialsEnabled():
                constraints = [constraints,

                                Afuse >= np.pi*Rfuse**2,
                                lnose >= 2*5.2*units.m, # TODO UPDATE - less arbitrary
                                wfuse >= SPR*wseat + 2*waisle + 2*wsys + 2*tshell,
                                Ahold <= wcargofl*2*hhold,

                                
#                if version == 'Philip':

#                   constraints = [constraints,
#                              Mfloor == Pfloor*wfloor/4.,
#                             ]


                                Mpassfl == Ppassfl*wpassfl/4.,
                                ]
                #else:
                 #   raise NotImplementedError

        elif fuselage_type == 'D8':
            raise NotImplementedError

            constraints = [constraints,

                            Afuse >= np.pi*Rfuse**2, # UPDATE
                            lnose >= 5.2*units.m, # TODO UPDATE - less arbitrary
                            wfuse >= SPR*wseat + 2*waisle + 2*wsys, #UPDATE
                            #Mpassfl == 9./256*Ppassfl*wpassfl,
                            Mfloor == 9./256*Pfloor*wfloor,
                            
                            # Model the additional  requred structural bending resistance
                            # due to having the engines at tha rear

                            # theta_db = arcsin(w_db/Rfuse), theta_db is the angle to vertical
                            # that the intersection between the two bubles make

                            # h_db**2 = Rfuse**2 - w_db**2, h_db vertical distance from
                            # intersection to central of height of bubble
                            # w_db is the horizontal distance from bubble center to web

                            # Asking = (2*pi + 4*theta_db)*Rfuse*tskin + 2*DRfuse*tskin
                            # DRfuse is the additional radius extension for additional cargo area

                            # A_db = (2*h_db + DRfuse)*t_db, A_db is web cross sectional area
                            # t_db is the webs thickness, from wing to wing.

                            # I_hshell = ((pi + 2*theta_db + sin(2*theta_db))*Rfuse**2 + 4*cos(theta_db)*DRfuse*RFuse + 0.5*(pi + 2*theta_db)*DRfuse**2)*Rfuse*tshell + 2./3*(h_db + DRfuse/2)**3*t_db
                            # I_hshell is the horziontal bending interia of the shell - equation A.9 from TASOPT
                            # I_vshell = ((pi + 2*theta_db - sin(2*theta_db))*Rfuse**2 + 8*cos(theta_db)*w_db*Rfuse + (2*pi + 4*theta_db)*w_db**2)*Rfuse*tshell
                            # I_vshell is the vertical bending interia of the shell - equation A.10 from TASOPT

                            # t_db = 2*dp*w_db/sigskin

                            # V_db = A_db*lshell
                            # V_db is the volume of the web

                            # xV_db = 0.5*(xshell1 + xshell2)*V_db
                            # where xV_db is the resulting moment

                            # W_db = rhoskin*g*V_db, W_db is weight of web, might want to change material properties??



                            ]
            #Add in 
            #volume of fuselage
            #sizing of internal webs


        else:
            raise NameError


        CG_constraints = [ #TASOPT and Philip might differ here in definitions of xshell1!!!!!
                          xVcyl >= 0.5*(xshell1+xshell2)*Vcyl,
                          xVnose >= 0.5*(xshell1)*Vnose,
                          xVbulk >= xshell2*Vbulk, # simplified
                          xWskin >= rhoskin*g*(xVcyl + xVnose + xVbulk),
                          xWshell >= xWskin*(1 + fstring + fframe + ffadd),
                          xWwindow >= 0.5* (xshell1+xshell2)*Wwindow,
                          xWinsul >= 0.5*(xshell1 + xshell2)*Winsul,
                          xWapu == xapu*Wapu,
                          xWseat >= 0.5*(xshell1 + xshell2)*Wseat,
                          xWpadd >= 0.5*(xshell1 + xshell2)*Wpadd,
                          xWfix == xfix*Wfix,
                          xWfloor >= 0.5*(xshell1 + xshell2)*Wfloor,
                          xWcone >= 0.5*(xshell2 + lfuse) * Wcone,
                          xWfuse >= xWfix + xWapu + xWpadd + xWseat + xWshell
                                  + xWcone + xWwindow + xWinsul + xWfloor,
                          xCGfu == xWfuse/Wfuse,
                         ]


        self.CG_constraints = CG_constraints

        if version == 'Philip':
        
            Model.__init__(self, Dfuse + 0.5*Wfuse, constraints)

        elif version == 'NextGen':

            Model.__init__(self, None, constraints)


    def defaultsubs(self, fuselage_type = 'narrowbody'):

        substitutions = {
                         'LF': 0.898, # Might want to look into other values
                         'N_{land}': 6.0, # [TAS]
                         'T_{cabin}': 300,  
                         'W\'\'_{floor}': 60, # [TAS]
                         'W\'\'_{insul}': 22, # [TAS]
                         'W\'_{seat}': 150, # Boeing
                         'W\'_{window}': 145.*3, # [TAS]
                         'W_{avg. pass}': 180,
                         'W_{carry on}': 15,
                         'W_{checked}': 40,
                         'W_{fix}': 3000, # might differe depending on aircraft
                         '\\Delta p': 83000,  
                         #'\\mu': 1.4E-5,  # PERFORMANCE MODEL
                         '\\rho_{\\infty}': 0.38,  # PERFORMANCE MODEL
                         '\\rho_{bend}': 2700, # [TAS]
                         '\\rho_{cargo}': 150, # b757 freight doc
                         '\\rho_{cone}': 2700, # [TAS]
                         '\\rho_{floor}': 2700, # [TAS]
                         '\\rho_{lugg}': 100,
                         '\\rho_{skin}': 2700, # [TAS]
                         '\\sigma_{floor}': 30000/0.000145, # [TAS]
                         '\\sigma_{skin}': 15000/0.000145, # [TAS]
                         '\\tau_{floor}': 30000/0.000145, # [TAS]
                         'f_{apu}': 0.035, # [TAS]
                         'f_{fadd}': 0.20, # [TAS]
                         'f_{frame}': 0.25,
                         'f_{lugg,1}': 0.4,
                         'f_{lugg,2}': 0.1,
                         'f_{padd}': 0.4, # [TAS]
                         'f_{string}': 0.35, # [TAS]
                         'p_s': 31,
                         'p_{cabin}': 75000,  
                         'r_E': 1.0, # [TAS]
                         'w_{aisle}': 0.51, # Boeing
                         'w_{seat}': 0.5,
                         'w_{sys}': 0.10,
                        }

        if fuselage_type == 'narrowbody':

            subs737 = {
                         'L_{v_{max}}': 35000,  
                         'SPR': 6,
                         #'V_{\\infty}': 234,  # PERFORMANCE MODEL
                         'W_{cargo}': 10000,
                         'b_{vt}': 7,
                         'c_{vt}': 4,
                         '\\Delta h': 1,
                         'n_{seat}': 186,

                        }
            substitutions.update(subs737)

        elif fuselage_type == 'widebody':

            subs777 = {
                         'L_{v_{max}}': 35000, # UPDATE vertical tail loading
                         'SPR': 9,
                         #'V_{\\infty}': 248,  # PERFORMANCE MODEL
                         'W_{cargo}': 10000*10, # up for debate, 10x of 737
                         'b_{vt}': 9.24,
                         'c_{vt}': 5.78,
                         '\\Delta h': 0.5, # UPDATE?
                         'n_{seat}': 540,                                                  

                        }
            substitutions.update(subs777)

        elif fuselage_type == 'D8':
            raise NotImplementedError

            subsD8 = {
                         'L_{v_{max}}': 35000, # UPDATE vertical tail loading 
                         'SPR': 9, #UPDATE
                         #'V_{\\infty}': 248, # UPDATE  # PERFORMANCE MODEL
                         'W_{cargo}': 100000, # up for debate
                         'b_{vt}': 9.24, # UPDATE
                         'c_{vt}': 5.78, # UPDATE
                         '\\Delta h': 1, # UPDATE
                         'n_{seat}': 540, # UPDATE                                                  

                        }
            substitutions.update(subsD8)

            #Add in the subs here - need to review most of the assumptions
            #Maybe there should be two different D8 models, 
            #one for 737 size and another for 7777

        else:
            raise NameError



        return substitutions


    def dynamic(self, state):
        #Creates an performance model
        return FuselagePerformance(self, state)


    @classmethod
    def standalonefuselage(cls, fuselage_type = 'narrowbody'):
        """Create standalone instance of fuselage model"""

        ccs = cls(fuselage_type)

        substitutions = ccs.defaultsubs(fuselage_type)

        m = Model(ccs.cost, ccs, substitutions)
        return m

    @classmethod
    def coupled737(cls):
        """Creates instance of fuselage model for use in full aircraft model"""

        ccs = cls()

        constraints = ccs + ccs.CG_constraints
 
        dsubs = ccs.defaultsubs()
        linkedsubs = ['L_{v_{max}}', 'V_{\\infty}', 'b_{vt}', 'c_{vt}']
        substitutions = {key: value for key, value in dsubs.items()
                                    if key not in linkedsubs}

        m = Model(ccs.cost, constraints, substitutions, name='Fuselage')
        return m

    @classmethod
    def test(cls, fuselage_type = 'narrowbody'):
    	
        fu = cls.standalonefuselage(fuselage_type)

        return fu.localsolve(solver = None, verbosity=4)


class Atmosphere(Model):
    def __init__(self, **kwargs):
        g = Variable('g', 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable("p_{sl}", "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", "K/m", "Temperature lapse rate")
        T_atm = Variable("T_{atm}", "K", "air temperature")
        M_atm = Variable("M_{atm}", "kg/mol",
                         "Molar mass of dry air")
        R_atm = Variable("R_{atm}", "J/mol/K",
                         "air specific heating value")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        TH = (g*M_atm/R_atm/L_atm).value

        rho = Variable('\\rho_{\\infty}', 0.38, 'kg/m^3', 'Density of air')

        h = Variable("h", "ft", "Altitude")

        """
        Dynamic viscosity (mu) as a function of temperature
        References:
        http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
            atmos/atmos.html
        http://www.cfd-online.com/Wiki/Sutherland's_law
        """
        mu  = Variable('\\mu', 1.4E-5,'kg/(m*s)', 'Dynamic viscosity')

        T_s = Variable('T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')

        t_plus_ts_approx = (T_atm + T_s).mono_approximation({T_atm: 288.15,
                                                         T_s: T_s.value})

        with SignomialsEnabled():
            constraints = [
                # Pressure-altitude relation
                (p_atm/p_sl)**(1/5.257) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),
                #rho == rho,

                #temperature equation
                SignomialEquality(T_sl, T_atm + L_atm*h),

                #constraint on mu
                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
                #mu == mu,
                ]


        Model.__init__(self, T_atm, constraints, **kwargs)


class FlightState(Atmosphere):
    """
    creates atm model for each flight segment, has variables
    such as veloicty and altitude
    """
    def __init__(self, **kwargs):
        #declare variables
        Vinf     = Variable('V_{\\infty}', 245, 'm/s', 'Cruise velocity')
        #a        = Variable('a', 'm/s', 'Speed of Sound')
        #gamma    = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        #h        = Variable('h', 'm', 'Segment Altitude [meters]')
        #hft      = Variable('hft', 35000, 'feet', 'Segment Altitude [feet]')
        #Minf     = Variable('M_{\\infty}', '-', 'Cruise Mach number')

        atm      = Atmosphere()

        #make new constraints
        constraints = []

        constraints.extend([
            Vinf == Vinf, #required so velocity variable enters the model

         #   h == hft, #convert the units on altitude

            #compute the speed of sound with the state
          #  a  == (gamma * atm['R_{atm}'] * atm['T_{atm}'] / atm['M_{atm}'])**.5,

            #compute the mach number
           # Minf == Vinf / a
            ])

        #Constraint(Vinf, Vinf, '==')

        #build the model
        Model.__init__(self, None, [atm + constraints], **kwargs)


class FuselagePerformance(Model):
    "Wing drag model"
    def __init__(self, fuselage, state, **kwargs):

        Dfuse    = Variable('D_{fuse}', 'N', 'Total drag in cruise')
        Dfrict   = Variable('D_{friction}', 'N', 'Friction drag')
        Dupswp   = Variable('D_{upsweep}', 'N', 'Drag due to fuse upsweep')
        
        
        constraints = [

                # Drag
                # sources: Raymer (p285), kfid 325 notes (p180)
                Dfrict >= fuselage['FF'] * np.pi*fuselage['R_{fuse}'] * state['\\mu'] * state['V_{\\infty}']
                          * 0.074*(state['\\rho_{\\infty}']*state['V_{\\infty}']*fuselage['l_{fuse}']/state['\\mu'])**0.8,

                # Drag due to fuselage upsweep (Raymer p286)
                Dupswp >= 3.83*fuselage['\\phi']**2.5*fuselage['A_{fuse}'] * 0.5*state['\\rho_{\\infty}']*state['V_{\\infty}']**2,
                Dfuse >= Dfrict + Dupswp
            ]

        Model.__init__(self, None, constraints, **kwargs)








class HTail(Model):


    def __init__(self,ops,**kwargs):
        self.ops = ops
        Whtail       = Variable('W_{htail}',10000, 'N', 'Horizontal tail weight') #Temporarily
        Lhmax        = Variable('L_{h_{max}}',35000,'N', 'Max horizontal tail load')
        Shtail       = Variable('S_{htail}',32*0.8,'m^2','Horizontal tail area') #Temporarily
        CLhmax       = Variable('C_{L_{h_{max}}}', 2.5, '-', 'Max lift coefficient') #Temporarily
        constraints = [#Lhmax    == 0.5*self.ops['\\rho_{\\infty}']*self.ops['V_{NE}']**2*Shtail*CLhmax,
                       Lhmax    == Lhmax,
                       Whtail   == Whtail,
                       Shtail   == Shtail,
                       CLhmax   == CLhmax]
        Model.__init__(self, None, constraints, **kwargs)

class VTail(Model):


    def __init__(self,ops,**kwargs):
        #bvt          = Variable('b_{vt}',7, 'm', 'Vertical tail span')
        #Lvmax        = Variable('L_{v_{max}}',35000,'N', 'Max vertical tail load')
        Wvtail       = Variable('W_{vtail}',10000, 'N', 'Vertical tail weight') #Temporarily
        Qv           = Variable('Q_v', 'N*m', 'Torsion moment imparted by tail')

        constraints = [#bvt == bvt, 
                       #Lvmax == Lvmax,
                       Wvtail == Wvtail,
                       Qv == Qv]
        Model.__init__(self, None, constraints, **kwargs)




class Wing(Model):


    def __init__(self,ops,**kwargs):
        self.ops = ops
        constraints = [];
        Model.__init__(self,None,constraints,**kwargs)

class WingBox(Model):
    def dynamic(self,state):
        return WingBoxP(self,state)

    def __init__(self,ops,**kwargs):
        self.ops = ops
        xf           = Variable('x_f','m','x-location of front of wingbox')
        xb           = Variable('x_b','m','x-location of back of wingbox')
        c0           = Variable('c_0','m','Root chord of the wing')
        wbar         = Variable('\\bar_w',0.5,'-','Wingbox to chord ratio') #Temporarily
        xwing        = Variable('x_{wing}','m', 'x-location of wing')
        dxwing       = Variable('dx_{wing}','m','wing box offset')
        # Setting bending area integration bounds (defining wing box locations)
        with SignomialsEnabled():
            constraints  = [SignomialEquality(xf,xwing + dxwing + .5*c0*wbar), #[SP] [SPEquality]
                        #xf >= xwing + dxwing + .5*c0*wbar,
                        SignomialEquality(xb, xwing - dxwing + .5*c0*wbar), #[SP] [SPEquality]
                        #xb <= xwing - dxwing + .5*c0*wbar, #[SP]        
                        ];
        Model.__init__(self,None,constraints,**kwargs)


















class CommericalAircraft(Model):

    """
    class to link all models needed to simulate a commercial flight
    """
    def __init__(self, fuselage_type = 'narrowbody', version = 'Philip', **kwargs):
        #defining the number of segments
        Nclimb = 2
        Ncruise = 2

        #make the segment range
        Nseg = Nclimb + Ncruise

        #define all the submodels
        state = FlightState()
        wing = Wing(state)
        wingbox = WingBox(state)
        htail = HTail(state)
        vtail = VTail(state)



        fuse = Fuselage(wingbox, htail, vtail, fuselage_type, version)
        
        fusep = FuselagePerformance(fuse, state)


        substitutions = fuse.defaultsubs(fuselage_type)

        submodels = [fuse, state, fusep, wing, wingbox, htail, vtail]


        constraints = ConstraintSet([submodels])

        objective = fusep['D_{fuse}'] + 0.5*fuse['W_{fuse}']

        lc = LinkedConstraintSet(constraints)

        if version == 'Philip':
            Model.__init__(self, fuse.cost, lc + submodels, substitutions, **kwargs)
        elif version == 'NextGen':
            Model.__init__(self, objective, lc + submodels, substitutions, **kwargs)
        else:
            raise NameError('Version not yet implemented')

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            constraints.append([ub*units >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= lb*units])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m

    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, solver='mosek', verbosity=0,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver, verbosity, **kwargs)
        solhold = sol
        lam = sol["sensitivities"]["la"][1:]
        out = defaultdict(list)
        for i, varkey in enumerate(m.bound_all["varkeys"]):
            lam_gt, lam_lt = lam[2*i], lam[2*i+1]
            if abs(lam_gt) >= 1e-7:  # arbitrary threshold
                out["sensitive to upper bound"].append(varkey)
            if abs(lam_lt) >= 1e-7:  # arbitrary threshold
                out["sensitive to lower bound"].append(varkey)
            value = mag(sol["variables"][varkey])
            distance_below = np.log(value/m.bound_all["lb"])
            distance_above = np.log(m.bound_all["ub"]/value)
            if distance_below <= 3:  # arbitrary threshold
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out, solhold

if __name__ == '__main__':
    #Fuse = Fuselage(fuselage_type = 'narrowbody')
    #State = FlightState()
    #FuseP = FuselagePerformance(Fuse, State)
    #Mission = Model(Fuse.cost, [Fuse, FuseP])

    #sol, solhold = FuseP.determine_unbounded_variables(FuseP,verbosity=4, iteration_limit=100)
    #sol, solhold = Fuse.determine_unbounded_variables(Fuse, verbosity = 4, iteration_limit = 100)
    CA = CommericalAircraft(fuselage_type = 'narrowbody', version = 'NextGen')

    #sol, solhold = CA.determine_unbounded_variables(CA, verbosity = 4, iteration_limit = 100)
    solhold = CA.localsolve(solver='mosek', verbosity = 1)
    print(solhold.table())


# if __name__ == "__main__":
    #sol = Fuselage.test(fuselage_type = 'widebody')
    #print(sol.table())

    #Fuse = Fuselage(fuselage_type = 'narrowbody')
    #State = FlightState()
    #FuseP = FuselagePerformance(Fuse, State)
    #Mission = Model(Fuse.cost, [Fuse, FuseP])
    #sol = Mission.localsolve()#(solver = 'mosek')
    #print(sol.table())












