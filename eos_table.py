"""
Classes and functions for accessing and manipulating tabulated EOS data.
"""
############################################################## 
# IMPORT PYTHON MODULES
import numpy as np
#
# GLOBAL VARIABLE
__version__ = 'v1.1.5' # 7/1/2023 added temperatures to MGR and TIL-iSALE EOS; 
#                        fixed logic for iSALE Tillotson; problems with SESAME interpolation near zero pressure and tension
#__version__ = 'v1.1.4' # 6/24/2023 added load STD NOTENSION table
#__version__ = 'v1.1.3' # 2/15/2023 added write SESAME EXT table function
#__version__ = 'v1.1.1' # 2/11/2023 minor print statement removals
#__version__ = 'v1.1' # 2/6/2023 STS adding ideal gas compatibility
#__version__ = 'v1.0.1' # 2/27/2022 STS Added to classes for plotting data against the ANEOS model
#__version__ = 'v1.0' # 09/2019 STS original module created for ANEOS-Forsterite-2019
#
############################################################## 
#
# calculate the structure for one planet
# make a class to hold the PREM data
class isentrope_class:
    """Class to isentrope data extracted from EOS table."""  # this is a documentation string for this class
    def __init__(self): # self is the default name of the object for internal referencing of the variables in the class
        """A function to initialize the class object.""" # this is a documentation string for this function
        self.ND = 0 # number of radius points
        self.density     = []   
        self.pressure    = []
        self.temperature = []
        self.soundspeed  = []
        self.energy  = []
        self.partvel = []
        self.region = [] # Tillotson region flag
        # not going to use all the variables in the file
        self.units = '' # I like to keep a text note in a structure about the units
        self.label = ''
#
# make a class to hold a curve (e.g., isentrope) through EOS surface
class EOScurve:
    """Class EOS curve."""  # this is a documentation string for this class
    def __init__(self): # self is the default name of the object for internal referencing of the variables in the class
        """A function to initialize the class object.""" # this is a documentation string for this function
        self.rho     = []   
        self.P    = []
        self.P_err    = []
        self.T = []
        self.cs  = []
        self.gamma  = []
        self.U  = []
        self.U_err  = []
        self.up = []
        self.us = [] 
        # not going to use all the variables in the file
        self.units = '' # I like to keep a text note in a structure about the units
        self.label = ''
#
class EOShugoniot:
    """Class for Hugoniot array from extEOStable."""	
    def __init__(self):
        self.NH = 0
        self.rho0 = 0.
        self.rho0_err = 0.
        self.T0 = 0.
        self.rho = np.zeros(self.NH)   
        self.rho_err = np.zeros(self.NH)   
        self.T = np.zeros(self.NH)   
        self.T_err = np.zeros(self.NH)   
        self.P = np.zeros(self.NH)   
        self.P_err = np.zeros(self.NH)   
        self.U = np.zeros(self.NH)   
        self.S = np.zeros(self.NH)   
        self.S_err = np.zeros(self.NH)   
        self.up = np.zeros(self.NH)   
        self.up_err = np.zeros(self.NH)   
        self.us = np.zeros(self.NH)
        self.us_err = np.zeros(self.NH)
        self.cs = np.zeros(self.NH)
        self.ref = np.zeros(self.NH) # reflectivity
        self.ref_err = np.zeros(self.NH) # reflectivity
        self.gamma = np.zeros(self.NH)
        self.gamma_err = np.zeros(self.NH)
        self.units = ''
        self.label = ''
#
class EOSvaporcurve:
    """Class for vapor curve from ANEOS."""	
    def __init__(self):
        self.NT = 0
        self.NV = 0
        self.T = np.zeros(self.NT)  
        self.rl = np.zeros(self.NT)  
        self.rv = np.zeros(self.NT)  
        self.Pl = np.zeros(self.NT)  
        self.Pv = np.zeros(self.NT)  
        self.Ul = np.zeros(self.NT)  
        self.Uv = np.zeros(self.NT)  
        self.Sl = np.zeros(self.NT)  
        self.Sv = np.zeros(self.NT)
        self.Gl = np.zeros(self.NT)  
        self.Gv = np.zeros(self.NT)
        self.units = ''
#
class EOSmeltcurve:
    """Class for melt curve from ANEOS."""	
    def __init__(self):
        self.NT = 0
        self.NV = 0
        self.T  = np.zeros(self.NT)  
        self.Tl = np.zeros(self.NT)  
        self.Ts = np.zeros(self.NT)  
        self.rl = np.zeros(self.NT)  
        self.rs = np.zeros(self.NT)  
        self.Pl = np.zeros(self.NT)  
        self.Ps = np.zeros(self.NT)  
        self.Ul = np.zeros(self.NT)  
        self.Us = np.zeros(self.NT)  
        self.Sl = np.zeros(self.NT)  
        self.Ss = np.zeros(self.NT)
        self.units = ''
#
class EOS1barcurve:
    """Class for 1bar curve from the EOS."""	
    def __init__(self):
        self.NT    = 0
        self.S     = np.zeros(self.NT)  
        self.T     = np.zeros(self.NT)  
        self.rho   = np.zeros(self.NT)
        self.Tvap  = 0.
        self.Tmelt = 0.
        self.Sim   = 0.
        self.Scm   = 0.
        self.Siv   = 0.
        self.Scv   = 0.
        self.rhoiv   = 0.
        self.rhocv   = 0.
        self.rhocm   = 0.
        self.rhoim   = 0.
        self.units = ''
        self.label = ''
#
class EOSpoint:
    """Class for a state on the EOS."""	
    def __init__(self):
        self.P   = 0
        self.P_err   = 0
        self.S   = 0  
        self.S_err   = 0  
        self.T   = 0 
        self.T_err   = 0 
        self.rho = 0
        self.rho_err = 0
        self.U   = 0
        self.U_err   = 0
        self.units = ''
        self.label = ''
    def __str__(self):
        """ Print the values in the EOS point """
        return '\nClass EOSpoint: \n' + \
               f'   eos_table module version {__version__} \n' + \
               f'   P, rho, T: {self.P} {self.rho} {self.T} \n' + \
               f'   U, S: {self.U} {self.S}  \n' 
       
#
class EOScriticalpoint:
    """Class for critical point state from the EOS."""	
    def __init__(self):
        self.P   = 0
        self.S   = 0  
        self.T   = 0 
        self.rho = 0
        self.U   = 0
        self.units = ''
        self.label = ''
#
class EOStriplepoint:
    """Class for triple point state from the EOS."""	
    def __init__(self):
        self.P   = 0
        self.T   = 0 
        self.Sim   = 0.
        self.Scm   = 0.
        self.Siv   = 0.
        self.Scv   = 0.
        self.rhol  = 0.
        self.units = ''
        self.label = ''
#
class EOSaneoshugoniot:
    """Class for Hugoniot calculated in ANEOS."""	
    def __init__(self):
        self.ND  = 0
        self.NV  = 0
        #self.all = np.zeros((self.ND,self.NV))
        self.rho = 0
        self.T   = 0
        self.P   = 0
        self.U   = 0
        self.S   = 0
        self.us  = 0
        self.up  = 0
        self.units = ''
#
class extEOStable:
    """Class for accessing EXTENDED SESAME-STYLE EOS tables output from ANEOS"""
    #     ANEOS KPA FLAG
    #                                TABLE          ANEOS
    #     KPAQQ=STATE INDICATOR      =1, 1p    =1, 1p    (eos without melt)
    #                                =2, 2p lv =2, 2p liquid/solid plus vapor
    #                                          =4, 1p solid  (eos with melt)
    #                                          =5, 2p melt   (eos with melt)
    #                                          =6, 1p liquid (eos with melt)
    #                                =-1 bad value of temperature
    #                                =-2 bad value of density
    #                                =-3 bad value of material number
    #
    def __init__(self):
        self.ND  = 0 # integer; number of density points in grid
        self.NT  = 0 # integer; number of temperature points in grid
        self.rho = np.zeros(self.ND)          # g/cm3, density values
        self.T   = np.zeros(self.NT)          # K, temperature values
        self.P   = np.zeros(self.ND*self.NT)  # GPA, pressure(T,rho)
        self.U   = np.zeros(self.ND*self.NT)  # MJ/kg, sp. internal energy(T,rho)
        self.A   = np.zeros(self.ND*self.NT)  # MJ/kg, Helmholtz free energy(T,rho)
        self.S   = np.zeros(self.ND*self.NT)  # MJ/K/kg, sp. entropy(T,rho)
        self.cs  = np.zeros(self.ND*self.NT)  # cm/s, sound speed(T,rho)
        self.cv  = np.zeros(self.ND*self.NT)  # MJ/K/kg, sp. heat capacity(T,rho)
        self.KPA = np.zeros(self.ND*self.NT)  # integer, ANEOS KPA flag(T,rho)
        self.MDQ = np.zeros(self.ND*self.NT)  # integer, Model Development Quality Flag(T,rho)
        self.units = ''
        self.hug = EOShugoniot()
        self.hugo = EOShugoniot()
        self.vc  = EOSvaporcurve()
        self.mc  = EOSmeltcurve()
        self.cp  = EOScriticalpoint()
        self.tp  = EOStriplepoint()
        self.onebar = EOS1barcurve()
        self.anhug = EOSaneoshugoniot()
        # these are variables needed for the sesame header
        self.MATID   = 0.
        self.DATE    = 0.
        self.VERSION = 0.
        self.FMN     = 0.
        self.FMW     = 0.
        self.R0REF   = 0.
        self.K0REF   = 0.
        self.T0REF   = 0.
        self.P0REF   = 0.
        self.CS0REF  = 0.
        # variables needed for the ANEOS gamma function
        self.gamma0 = 0.
        self.theta0 = 0.
        self.C24    = 0.
        self.C60    = 0.
        self.C61    = 0.
        self.beta   = 0.
        # model name/version string
        self.MODELNAME = ''

    def loadstdsesame(self, fname, unitstxt=None):
        """Function for loading STD SESAME-STYLE EOS table output from ANEOS"""
        data = ([])
        if unitstxt is None:
            self.units = 'Units: rho g/cm3, T K, P GPa, U MJ/kg, A MJ/kg, S MJ/K/kg, cs cm/s, cv MJ/K/kg, KPA flag. 2D arrays are (NT,ND).'
        else:
            self.units = unitstxt
        sesamefile = open(fname,"r")  
        sesamedata=sesamefile.readlines()
        sesamefile.close()
        nskip = 6 # skip standard header to get to the content of the 301 table
        # num.density, num. temps
        tmp = sesamedata[nskip][0:16]
        dlen = float(tmp)
        tmp = sesamedata[nskip][16:32]
        tlen = float(tmp)
        if (np.mod((dlen*tlen*3.0+dlen+tlen+2.0),5.0) == 0):
            neos = int((dlen*tlen*3.0+dlen+tlen+2.0)/5.0) 
        else:
            neos = int((dlen*tlen*3.0+dlen+tlen+2.0)/5.0) +1
        #print(dlen,tlen,neos,len(sesamedata))
        data = np.zeros((neos,5),dtype=float)
        for j in range(nskip,neos+nskip):
            tmp3 = sesamedata[j]
            tmp4 = list(tmp3.split())
            if len(tmp4) < 5:
                lentmp4 = len(tmp4)
                data[j-nskip,0:lentmp4] = np.asarray(tmp4[0:lentmp4])
            else:
                data[j-nskip,:] = np.asarray(tmp4)
            #print(j,eosarr[j,:])
        #print(data.shape)
        data=np.resize(data,(neos*5))
        #print(data.shape)
        self.ND  = data[0].astype(int)  # now fill the extEOStable class
        self.NT  = data[1].astype(int)
        self.rho = data[2:2+self.ND]
        self.T   = data[2+self.ND : 2+self.ND+self.NT]
        self.P   = data[2+self.ND+self.NT : 2+self.ND+self.NT+self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.U   = data[2+self.ND+self.NT+self.ND*self.NT
                            : 2+self.ND+self.NT+2*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.A   = data[2+self.ND+self.NT+2*self.ND*self.NT
                            : 2+self.ND+self.NT+3*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        #self.S   = data[2+self.ND+self.NT+3*self.ND*self.NT
        #                    : 2+self.ND+self.NT+4*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.cs  = data[2+self.ND+self.NT+4*self.ND*self.NT
        #                    : 2+self.ND+self.NT+5*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.cv  = data[2+self.ND+self.NT+5*self.ND*self.NT
        #                    : 2+self.ND+self.NT+6*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.KPA = data[2+self.ND+self.NT+6*self.ND*self.NT
        #                    : 2+self.ND+self.NT+7*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
#
    def loadextsesame(self, fname, unitstxt=None):
        """Function for loading EXTENDED SESAME-STYLE EOS table output from ANEOS"""
        data = ([])
        if unitstxt is None:
            self.units = 'Units: rho g/cm3, T K, P GPa, U MJ/kg, A MJ/kg, S MJ/K/kg, cs cm/s, cv MJ/K/kg, KPA flag. 2D arrays are (NT,ND).'
        else:
            self.units = unitstxt
        sesamefile = open(fname,"r")  
        sesamedata=sesamefile.readlines()
        sesamefile.close()
        nskip = 6 # skip standard header to get to the content of the 301 table
        # num.density, num. temps
        tmp = sesamedata[nskip][0:16]
        dlen = float(tmp)
        tmp = sesamedata[nskip][16:32]
        tlen = float(tmp)
        if (np.mod((dlen*tlen*4.0+dlen+tlen+2.0),5.0) == 0):
            neos = int((dlen*tlen*4.0+dlen+tlen+2.0)/5.0)
        else:
            neos = int((dlen*tlen*4.0+dlen+tlen+2.0)/5.0) +1
        #print(dlen,tlen,neos,len(sesamedata))
        data = np.zeros((neos,5),dtype=float)
        for j in range(nskip,neos+nskip):
            tmp3 = sesamedata[j]
            tmp4 = list(tmp3.split())
            if len(tmp4) < 5:
                lentmp4 = len(tmp4)
                data[j-nskip,0:lentmp4] = np.asarray(tmp4[0:lentmp4])
            else:
                data[j-nskip,:] = np.asarray(tmp4)        
            #print(j,eosarr[j,:])
        #print(data.shape)
        data=np.resize(data,(neos*5))
        #print(data.shape)
        self.ND  = data[0].astype(int)  # now fill the extEOStable class
        self.NT  = data[1].astype(int)
        self.rho = data[2:2+self.ND]
        self.T   = data[2+self.ND : 2+self.ND+self.NT]
        #self.P   = data[2+self.ND+self.NT : 2+self.ND+self.NT+self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.U   = data[2+self.ND+self.NT+self.ND*self.NT
        #                    : 2+self.ND+self.NT+2*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        #self.A   = data[2+self.ND+self.NT+2*self.ND*self.NT
        #                    : 2+self.ND+self.NT+3*self.ND*self.NT
        #                    ].reshape(self.NT,self.ND)
        self.S   = data[2+self.ND+self.NT+0*self.ND*self.NT
                            : 2+self.ND+self.NT+1*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.cs  = data[2+self.ND+self.NT+1*self.ND*self.NT
                            : 2+self.ND+self.NT+2*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.cv  = data[2+self.ND+self.NT+2*self.ND*self.NT
                            : 2+self.ND+self.NT+3*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.KPA = data[2+self.ND+self.NT+3*self.ND*self.NT
                            : 2+self.ND+self.NT+4*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
#
    def view(self, q='P', Tlow=None, Thigh=None, rholow=None, rhohigh=None):
        """Function for printing values from EXTENDED SESAME-STYLE EOS table."""
        if Tlow is None:
            Tlow = self.T.min()
        if Thigh is None:
            Thigh = self.T.max()
        if rholow is None:
            rholow = self.rho.min()
        if rhohigh is None:
            rhohigh = self.rho.max()
        print(self.units)
        if q == 'P':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('P:', (self.P[np.logical_and(self.T >= Tlow,
                                                self.T<=Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'U':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('U:', (self.U[np.logical_and(self.T >= Tlow,
                                                self.T <= Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'A':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('A:', (self.A[np.logical_and(self.T >= Tlow,
                                                self.T <= Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'S':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho
                                   >= rholow,self.rho<=rhohigh)
                                  ])
            print('S:', (self.S[np.logical_and(self.T >= Tlow,self.T <= Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])      
        if q == 'cs':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cs:', (self.cs[np.logical_and(self.T >= Tlow,
                                                  self.T <= Thigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'cv':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cv:', (self.cv[np.logical_and(self.T >= Tlow,
                                                  self.T <= Thigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'KPA':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('KPA:', (self.KPA[np.logical_and(self.T >= Tlow,
                                                  self.T <= Thigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])

    def calchugoniot(self, r0=None, t0=None, pmax=None, writefilename=None):
        """Function for calculating a Hugoniot from EXTENDED SESAME-STYLE EOS table."""
        if r0 is None:
            return 'Must provide r0 and t0.'
        if t0 is None:
            return 'Must provide r0 and t0.'
        if pmax is None:
            pmax=1.E4 # GPa
        self.hug.rho = []
        self.hug.P = []
        self.hug.T = []
        self.hug.U = []
        self.hug.S = []
        self.hug.up = []
        self.hug.us = []
        self.hug.cs = []
  
        it0 = int(np.round(np.interp(t0,self.T,np.arange(self.NT)))) # uses nearest value if t0 not in array
        ir0 = int(np.round(np.interp(r0,self.rho,np.arange(self.ND)))) # uses nearest value if r0 not in the array
        p0  = self.P[it0,ir0] # GPa
        #print(self.P[it0,ir0])
        e0  = self.U[it0,ir0]#np.interp(p0,self.P[it0,:],self.U[it0,:])
        s0  = self.S[it0,ir0]#np.interp(p0,self.P[it0,:],self.S[it0,:])
        up0 = 0. # no initial particle velocity
        us0 = self.cs[it0,ir0]/1.e5 # cm/s->km/s use sound velocity for initial
        cs0 = self.cs[it0,ir0]/1.e5 # cm/s->km/s use sound velocity for initial
        #print(ir0,it0,r0,t0,p0,e0,up0,us0)
        self.hug.rho = np.append(self.hug.rho, self.rho[ir0])
        self.hug.P = np.append(self.hug.P, p0)
        self.hug.T = np.append(self.hug.T, self.T[it0])
        self.hug.U = np.append(self.hug.U, e0)
        self.hug.S = np.append(self.hug.S, s0)
        self.hug.up = np.append(self.hug.up, up0)
        self.hug.us = np.append(self.hug.us, us0)
        self.hug.cs = np.append(self.hug.cs, cs0)
        
        #for iir in range(ir0+1,self.ND):
        iir=ir0+1
        pnew=p0
        while pnew<pmax:
            ediff =0.5*(self.P[it0::,iir]+p0)*(1./r0-1./self.rho[iir])+e0 -(self.U[it0::,iir])  # MJ/kg
            # np.interp wants x values increasing
            pnew = np.interp(0.,np.flip(ediff),np.flip(self.P[it0::,iir]))
            tnew = np.interp(0.,np.flip(ediff),np.flip(self.T[it0::]))
            enew = np.interp(0.,np.flip(ediff),np.flip(self.U[it0::,iir]))
            snew = np.interp(0.,np.flip(ediff),np.flip(self.S[it0::,iir]))
            upnew = np.sqrt((pnew-p0)*(1./r0-1./self.rho[iir]))
            usnew = (1./r0)*np.sqrt((pnew-p0)/(1./r0-1./self.rho[iir]))
            csnew = np.interp(0.,np.flip(ediff),np.flip(self.cs[it0::,iir]))/1.E5 # km/s
            #print(self.rho[iir],tnew,pnew,enew,upnew,usnew)
            self.hug.rho = np.append(self.hug.rho, self.rho[iir])
            self.hug.P = np.append(self.hug.P, pnew)
            self.hug.T = np.append(self.hug.T, tnew)
            self.hug.U = np.append(self.hug.U, enew)
            self.hug.S = np.append(self.hug.S, snew)
            self.hug.up = np.append(self.hug.up, upnew)
            self.hug.us = np.append(self.hug.us, usnew)
            self.hug.cs = np.append(self.hug.cs, csnew) # km/s
            iir += 1
        self.hug.NH=len(self.hug.P)
        self.hug.units='units: T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg, Up km/s, Us km/s, cs km/s'

        if writefilename:
            print('Writing Hugoniot to file: ',writefilename)
            hugoniotfile = open(writefilename,"w")  
            hugoniotfile.writelines('  Hugoniot \n') 
            hugoniotfile.writelines('  Temperature,    Density,        Pressure,       Int. Energy,    Sp. Entropy,    Part. Vel.,     Shock Vel. \n') 
            hugoniotfile.writelines('  K,              g/cm3,          GPa,            MJ/kg,          MJ/K/kg,        km/s,           km/s\n') 
            for iih in range(0,self.hug.NH):
                hugoniotfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (
                    self.hug.T[iih],self.hug.rho[iih],self.hug.P[iih],self.hug.U[iih],self.hug.S[iih],self.hug.up[iih],self.hug.us[iih]))
            hugoniotfile.close() 

    def calcOffEOSHugoniot(self, r0=None, t0=None, p0=None, e0=None, r1=None, pmax=None, writefilename=None):
        """Function for calculating a Hugoniot from EXTENDED SESAME-STYLE EOS table."""
        if r0 is None:
            return 'Must provide r0, t0, p0, e0.'
        if t0 is None:
            return 'Must provide r0, t0, p0, e0.'
        if p0 is None:
            return 'Must provide r0, t0, p0, e0.'
        if e0 is None:
            return 'Must provide r0, t0, p0, e0.'
        if pmax is None:
            pmax=1.E4 # GPa
        self.hugo.rho = []
        self.hugo.P = []
        self.hugo.T = []
        self.hugo.U = []
        self.hugo.S = []
        self.hugo.up = []
        self.hugo.us = []
        self.hugo.cs = []
        print('r1=',r1)
        it0 = int(np.round(np.interp(t0,self.T,np.arange(self.NT)))) # uses nearest value if t0 not in array
        ir0 = int(np.round(np.interp(r1,self.rho,np.arange(self.ND)))) # start with full density
        s0  = 0. # don't need it to get going
        up0 = 0. # no initial particle velocity
        us0 = self.cs[it0,ir0]/1.e5 # cm/s->km/s use sound velocity for initial
        cs0 = self.cs[it0,ir0]/1.e5 # cm/s->km/s use sound velocity for initial
        print(ir0,it0,r0,t0,p0,e0,up0,us0)
        self.hugo.rho = np.append(self.hugo.rho, r0)
        self.hugo.P = np.append(self.hugo.P, p0)
        self.hugo.T = np.append(self.hugo.T, t0)
        self.hugo.U = np.append(self.hugo.U, e0)
        self.hugo.S = np.append(self.hugo.S, s0)
        self.hugo.up = np.append(self.hugo.up, up0)
        self.hugo.us = np.append(self.hugo.us, us0)
        self.hugo.cs = np.append(self.hugo.cs, cs0)
        
        #for iir in range(ir0+1,self.ND):
        iir=ir0+1
        pnew=p0
        count=0
        while pnew<pmax:
            ediff =0.5*(self.P[it0::,iir]+p0)*(1./r0-1./self.rho[iir])+e0 -(self.U[it0::,iir])  # MJ/kg
            # np.interp wants x values increasing
            pnew = np.interp(0.,np.flip(ediff),np.flip(self.P[it0::,iir]))
            tnew = np.interp(0.,np.flip(ediff),np.flip(self.T[it0::]))
            enew = np.interp(0.,np.flip(ediff),np.flip(self.U[it0::,iir]))
            snew = np.interp(0.,np.flip(ediff),np.flip(self.S[it0::,iir]))
            upnew = np.sqrt((pnew-p0)*(1./r0-1./self.rho[iir]))
            usnew = (1./r0)*np.sqrt((pnew-p0)/(1./r0-1./self.rho[iir]))
            csnew = np.interp(0.,np.flip(ediff),np.flip(self.cs[it0::,iir]))/1.E5 # km/s
            #print(self.rho[iir],tnew,pnew,enew,upnew,usnew)
            self.hugo.rho = np.append(self.hugo.rho, self.rho[iir])
            self.hugo.P = np.append(self.hugo.P, pnew)
            self.hugo.T = np.append(self.hugo.T, tnew)
            self.hugo.U = np.append(self.hugo.U, enew)
            self.hugo.S = np.append(self.hugo.S, snew)
            self.hugo.up = np.append(self.hugo.up, upnew)
            self.hugo.us = np.append(self.hugo.us, usnew)
            self.hugo.cs = np.append(self.hugo.cs, csnew) # km/s
            count=count+1
            iir += 1
        self.hugo.NH=len(self.hug.P)
        self.hugo.units='units: T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg, Up km/s, Us km/s, cs km/s'

        if writefilename:
            print('Writing Hugoniot to file: ',writefilename)
            hugoniotfile = open(writefilename,"w")  
            hugoniotfile.writelines('  Hugoniot \n') 
            hugoniotfile.writelines('  Temperature,    Density,        Pressure,       Int. Energy,    Sp. Entropy,    Part. Vel.,     Shock Vel. \n') 
            hugoniotfile.writelines('  K,              g/cm3,          GPa,            MJ/kg,          MJ/K/kg,        km/s,           km/s\n') 
            for iih in range(0,count):
                hugoniotfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (
                    self.hugo.T[iih],self.hugo.rho[iih],self.hugo.P[iih],self.hugo.U[iih],self.hugo.S[iih],self.hugo.up[iih],self.hugo.us[iih]))
            hugoniotfile.close() 
    ###########################################################################################
    # Functions added to be helpful for integration with the KO hydrocode
    #
    #
    def querypt(self,rho,t):
        """Return EOS state for input (rho,t)"""
        print('in query pt')
        pt = EOSpoint()
        # linearly interpolate the table for return the rho,t point values
        rindex = np.where(self.rho >= rho)[0][0]
        tindex = np.where(self.T >= t)[0][0]
        P11 = self.P[tindex-1,rindex-1]
        P21 = self.P[tindex,rindex-1]
        P12 = self.P[tindex-1,rindex]
        P22 = self.P[tindex,rindex]
        U11 = self.U[tindex-1,rindex-1]
        U21 = self.U[tindex,rindex-1]
        U12 = self.U[tindex-1,rindex]
        U22 = self.U[tindex,rindex]
        r1 = np.log10(self.rho[rindex-1]) # try linearly interpreting on the log of the index
        r2 = np.log10(self.rho[rindex])
        t1 = np.log10(self.T[tindex-1])
        t2 = np.log10(self.T[tindex])
        logr = np.log10(rho)
        logt = np.log10(t)
        u1 = np.log10(U11) * (r2-logr)/(r2-r1) + np.log10(U12) * (logr-r1)/(r2-r1)
        u2 = np.log10(U21) * (r2-logr)/(r2-r1) + np.log10(U22) * (logr-r1)/(r2-r1)
        p1 = np.log10(P11) * (r2-logr)/(r2-r1) + np.log10(P12) * (logr-r1)/(r2-r1)
        p2 = np.log10(P21) * (r2-logr)/(r2-r1) + np.log10(P22) * (logr-r1)/(r2-r1)
        logp = p1 * (t2-logt)/(t2-t1) + p2 * (logt-t1)/(t2-t1)
        pnew = np.power(10.,logp)
        logu = u1 * (t2-logt)/(t2-t1) + u2 * (logt-t1)/(t2-t1)
        unew = np.power(10.,logu)
        pt.T = t
        pt.rho = rho
        pt.P = pnew
        pt.U = unew
        return pt
    def queryptlin(self,rho,t):
        """Return EOS state for input (rho,t); linear EOS"""
        print('in query pt')
        pt = EOSpoint()
        # linearly interpolate the table for return the rho,t point values
        rindex = np.where(self.rho >= rho)[0][0]
        tindex = np.where(self.T >= t)[0][0]
        P11 = self.P[tindex-1,rindex-1]
        P21 = self.P[tindex,rindex-1]
        P12 = self.P[tindex-1,rindex]
        P22 = self.P[tindex,rindex]
        U11 = self.U[tindex-1,rindex-1]
        U21 = self.U[tindex,rindex-1]
        U12 = self.U[tindex-1,rindex]
        U22 = self.U[tindex,rindex]
        r1 = self.rho[rindex-1]
        r2 = self.rho[rindex]
        t1 = self.T[tindex-1]
        t2 = self.T[tindex]
        u1 = U11 * (r2-rho)/(r2-r1) + U12 * (rho-r1)/(r2-r1)
        u2 = U21 * (r2-rho)/(r2-r1) + U22 * (rho-r1)/(r2-r1)
        p1 = P11 * (r2-rho)/(r2-r1) + P12 * (rho-r1)/(r2-r1)
        p2 = P21 * (r2-rho)/(r2-r1) + P22 * (rho-r1)/(r2-r1)
        pnew = p1 * (t2-t)/(t2-t1) + p2 * (t-t1)/(t2-t1)
        unew = u1 * (t2-t)/(t2-t1) + u2 * (t-t1)/(t2-t1)
        pt.T = t
        pt.rho = rho
        pt.P = pnew
        pt.U = unew
        return pt
    #
    ###########################################################################################
    #            
    def writestdsesame(self, writestdsesfname=None):
        """Write standard Header-201-301 SESAME EOS TABULAR EOS FILE"""
        # write a standard SESAME ascii file
        #     WRITE STANDARD Header-201-301 SESAME FILE
        #     WRITE SESAME 301 TABLE CONTAINS P, E, HFE
        #sesfile = open("NEW-SESAME-STD-NOTENSION.EOSTXT","w")  
        if writestdsesfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writestdsesfname,"w")  
        #     WRITE SESAME HEADER INFORMATION: EOS matid number, number of words in section
        #     could input matid, date, version with the grid
        # these parameters are set in the cell above that sets up the grid for ANEOS
        # THEY SHOULD MATCH.......
        # These variables are needed for the standard table output
        NWDS=9
        SESNTABLES=2.0
        TABLE1 = 201.0
        TABLE2 = 301.0
        #     5 entries in 201 table
        SESNWDS1=5.0
        #     Number of entries in STANDARD 301 table: 3 variables at each T,rho point
        SESNWDS2=2.+self.ND+self.NT+self.ND*self.NT*3.
        #     HEADER SECTION
        #sesfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (antarr[iit],rnew,pnew,enew,snew,upnew,usnew))
        sesfile.write(" INDEX      MATID ={:7d}    NWDS = {:8d}\n".format(int(self.MATID), int(NWDS)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.MATID, self.DATE, self.DATE, self.VERSION, SESNTABLES))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(TABLE1, TABLE2, SESNWDS1, SESNWDS2))
        # 201 SECTION
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE1),int(SESNWDS1)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.FMN, self.FMW, self.R0REF, self.K0REF, self.T0REF))
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE2),int(SESNWDS2)))
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NT))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     temperature array K
        for j in range(0, int(self.NT)):
            sesfile.write("{:16.8e}".format(self.T[j]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  pressure array GPa P[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.P[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  specific internal energy array MJ/kg U[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.U[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # Helmholtz free energy array in MJ/kg A[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.A[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the SESAME TABLE FILE
        sesfile.close() 
        print('Done writing the STD SESAME 301 table to local directory: ',writestdsesfname)

    def writeextsesame(self, writeextsesfname=None):
        """Write standard Header-201-301 SESAME EOS TABULAR EOS FILE"""
        #     WRITE STANDARD Header-201-301 SESAME FILE
        if writeextsesfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writeextsesfname,"w")  
        #     WRITE SESAME HEADER INFORMATION: EOS matid number, number of words in section
        #     could input matid, date, version with the grid
        # these parameters are set in the cell above that sets up the grid for ANEOS
        # THEY SHOULD MATCH.......
        # These variables are needed for the standard table output
        NWDS=9
        SESNTABLES=2.0
        TABLE1 = 201.0
        TABLE2 = 301.0
        #     5 entries in 201 table
        SESNWDS1=5.0
        #     Number of entries in STANDARD 301 table: 3 variables at each T,rho point
        SESNWDS2=2.+self.ND+self.NT+self.ND*self.NT*3.
        #     HEADER SECTION
        #sesfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (antarr[iit],rnew,pnew,enew,snew,upnew,usnew))
        sesfile.write(" INDEX      MATID ={:7d}    NWDS = {:8d}\n".format(int(self.MATID), int(NWDS)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.MATID, self.DATE, self.DATE, self.VERSION, SESNTABLES))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(TABLE1, TABLE2, SESNWDS1, SESNWDS2))
        # 201 SECTION
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE1),int(SESNWDS1)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.FMN, self.FMW, self.R0REF, self.K0REF, self.T0REF))
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE2),int(SESNWDS2)))
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NT))
        STYLE=2
#NEW-SESAME-EXT.TXT: SESAME-style table with extra variables from ANEOS. Contains the standard 201 table and non-standard 301-extra-variables EOS table. The 301 table has: density grid values, temperature grid values, sp. entropy(T,rho), sound speed(T,rho), sp. heat capacity(T,rho), KPA flag(T,rho). 2-D arrays list all densities, looping over each temperature. 301 table units: g/cm$^3$, K, MJ/K/kg, cm/s, MJ/K/kg, integer flag, integer flag. The KPA flag is an ANEOS output with phase information. <br>

        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     temperature array K
        for j in range(0, int(self.NT)):
            sesfile.write("{:16.8e}".format(self.T[j]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  specific entropy array MJ/K/kg S[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.S[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  sound speed array cm/s cs[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.cs[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # sp. heat capacity array in MJ/K/kg cv[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.cv[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # phase flag [integer] KPA[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.KPA[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the SESAME TABLE FILE
        sesfile.close() 
        print('Done writing the EXT SESAME 301 table to local directory: ',writeextsesfname)

    def writegadgetinitsesame(self, writegadgetinitsesfname=None):
        """Write standard Header-201-301 SESAME EOS TABULAR EOS FILE"""
        # write a standard SESAME ascii file
        #     WRITE GADGET INIT Header-201-301 SESAME FILE
        #     WRITE SESAME 301 TABLE CONTAINS P, E, S
        #sesfile = open("NEW-SESAME-GADGETINIT-NOTENSION.EOSTXT","w")  
        if writegadgetinitsesfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writegadgetinitsesfname,"w")  
        #     WRITE SESAME HEADER INFORMATION: EOS matid number, number of words in section
        #     could input matid, date, version with the grid
        # these parameters are set in the cell above that sets up the grid for ANEOS
        # THEY SHOULD MATCH.......
        # These variables are needed for the standard table output
        NWDS=9
        SESNTABLES=2.0
        TABLE1 = 201.0
        TABLE2 = 301.0
        #     5 entries in 201 table
        SESNWDS1=5.0
        #     Number of entries in STANDARD 301 table: 3 variables at each T,rho point
        SESNWDS2=2.+self.ND+self.NT+self.ND*self.NT*3.
        #     HEADER SECTION
        #sesfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (antarr[iit],rnew,pnew,enew,snew,upnew,usnew))
        sesfile.write(" INDEX      MATID ={:7d}    NWDS = {:8d}\n".format(int(self.MATID), int(NWDS)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.MATID, self.DATE, self.DATE, self.VERSION, SESNTABLES))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(TABLE1, TABLE2, SESNWDS1, SESNWDS2))
        # 201 SECTION
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE1),int(SESNWDS1)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.FMN, self.FMW, self.R0REF, self.K0REF, self.T0REF))
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE2),int(SESNWDS2)))
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NT))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     temperature array K
        for j in range(0, int(self.NT)):
            sesfile.write("{:16.8e}".format(self.T[j]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  pressure array GPa P[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.P[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  specific internal energy array MJ/kg U[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.U[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # specific entropy array in MJ/kg/K S[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.S[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the SESAME TABLE FILE
        sesfile.close() 
        print('Done writing the GADGET INIT SESAME 301 notension table to local directory: ',writegadgetinitsesfname)

    def writemdqsesame(self, writemdqsesfname=None):
        """Function to write a sesame 301-style ascii file with the MDQ variable"""
        if writemdqsesfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writemdqsesfname,"w")  
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NT))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     temperature array K
        for j in range(0, int(self.NT)):
            sesfile.write("{:16.8e}".format(self.T[j]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  MDQ Flag[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.MDQ[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the SESAME TABLE FILE
        sesfile.close() 
        print('Done writing the MDQ Flag as a 301-style table to local directory: ',writemdqsesfname)

    def loadaneos(self, aneosinfname=None, aneosoutfname=None, silent=False):
        """Function for reading in ANEOS INPUT and OUTPUT FILE DATA into EOS structure."""
        if aneosinfname is None:
            return 'Must provide input file name.'
        if aneosoutfname is None:
            return 'Must provide output file name.'
        # function to gather data from ANEOS input and output files
        # SESAME FILE HEADER INFORMATION MUST BE LOADED INTO THE EOS STRUCTURE BEFORE CALLING THIS FUNCTION
        #
        # READ IN ANEOS INPUT FILE
        aneosinputfile = open(aneosinfname,"r")  
        testin=aneosinputfile.readlines()   # read in the whole ascii file at once because this is fatser
        aneosinputfile.close()
        # gather EOS information from the ANEOS.OUTPUT file
        aneosoutputfile = open(aneosoutfname,"r")  
        testout=aneosoutputfile.readlines() # read everything in at once because this is faster
        aneosoutputfile.close()
        if silent == False:
            print('Done loading ANEOS files.')

        # THIS CODE PARSES THE ANEOS.OUTPUT FILE INTO ARRAYS FOR USE IN PLOTTING/USING THE EOS
        if silent == False:
            print('ANEOS WAS CALLED WITH THE FOLLOWING INPUT, LOADED FROM FILE ',aneosinfname)
        # Gather parameters for the gamma function while printing the ANEOS INPUT FILE
        aneoscount=1
        for i in np.arange(len(testin)):
            if testin[i].find('ANEOS') == 0:
                if aneoscount<9:
                    if silent == False:
                        print(' '+testin[i-3],testin[i-2],testin[i-1],testin[i])
                    aneoscount=aneoscount+1
                else:
                    if silent == False:
                        print(' '+testin[i])
                if testin[i].find('ANEOS2') == 0:
                    tmp=testin[i]
                    nelem=int(tmp[10:20])
                    #print('nelem=',nelem)
                    rho0=float(tmp[30:40])
                    print('rho0=',rho0)
                    gamma0=float(tmp[70:80])
                    #print('gamma0=',gamma0)
                    theta0=float(tmp[80:90])
                if testin[i].find('ANEOS3') == 0:
                    tmp=testin[i]
                    C24=float(tmp[20:30])/3.
                    #print('C24=',C24)
                if testin[i].find('ANEOS5') == 0:
                    tmp=testin[i]
                    C60=float(tmp[60:70])
                    C61=float(tmp[70:80])
                    #print('C60=',C60)
                if testin[i].find('ANEOS7') == 0:
                    tmp=testin[i]
                    betagamma=float(tmp[70:80])

        # some checks
        if rho0 != self.R0REF:
            #print('WARNING: rho0 does not match. ')
            #assert(False) # just a way to stop the notebook
            print('WARNING: rho0 does not match. STOPPING THIS NOTEBOOK.')
            assert(False) # just a way to stop the notebook

        # GUESS A BIG ARRAY SIZE FOR THE PHASE BOUNDARIES AND HUGONIOT IN ANEOS.OUTPUT
        # the melt curve, vapor curve and Hugoniot curves are not fixed length outputs
        nleninit=300
        meltcurve = 0

        if silent == False:
            print('READING DATA FROM ANEOS OUTPUT FILE ',aneosoutfname)

        # Read in data from the ANEOS.OUTPUT FILE
        imc = -1 # flag for no melt curve in the model
        for i in np.arange(len(testout)):
            if testout[i].find('  Data for ANEOS number') == 0:
                tmp = testout[i+2][0:50]
                eosname = tmp.strip()
            if testout[i] == '  TWO-PHASE BOUNDARIES\n':
                nvc = nleninit
                ivc = i
                vcarrtmp = np.zeros((nvc,12),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    if testout[j+i+4].find(' anphas') == 0:
                        print(testout[j+i+4])
                        vcarrtmp[j,:]=vcarrtmp[j-1,:]
                        j=j+1
                    else:
                        tmp=str.replace(testout[j+i+4],'D','E')
                        tmp3 = tmp[0:157]
                        tmp4 = list(tmp3.split())
                        if (len(tmp4) >0) and (float(tmp4[3]) > 0) and (float(tmp4[4]) > 0): # stop if the pressures become negative on the vapor curve
                            tmp5 = np.asarray(tmp4)
                            vcarrtmp[j,:] = tmp5[:]
                            j=j+1
                        else:
                            flag=1
                vcarr = np.zeros((j,12),dtype=float)
                vcarr[:,:] = vcarrtmp[0:j,:]
            if testout[i] == ' LIQUID/SOLID PHASE CURVE\n':
                nmc = nleninit
                imc = i
                meltcurve=1
                mcarrtmp = np.zeros((nmc,11),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    tmp  = str.replace(testout[j+i+5],'D','E')
                    tmp3 = tmp[0:132]
                    tmp4 = list(tmp3.split())
                    if len(tmp4) > 0:
                        tmp5 = np.asarray(tmp4)
                        mcarrtmp[j,:] = tmp5[:]
                        j=j+1
                    else:
                        flag=1
                mcarr = np.zeros((j,11),dtype=float)
                mcarr[:,:] = mcarrtmp[0:j,:]
            if testout[i] == '   HUGONIOT\n':
                nhc = nleninit
                ihc = i
                hcarrtmp = np.zeros((nhc,9),dtype=float)
                flag=0
                j=0
#                while flag == 0:
                while j<22:
                    tmp=str.replace(testout[j+i+5],'D','E')
                    #print(j,i,5,tmp)
                    tmp3 = tmp[0:109]
                    tmp4 = list(tmp3.split())
                    if len(tmp4) > 0:
                        tmp4[3]='0.0' # this column often gives problems with exponential notation so don't read it
                        tmp5 = np.asarray(tmp4)
                        hcarrtmp[j,:] = tmp5[:]
                        j=j+1
                    else:
                        flag=1
                hcarr = np.zeros((j,9),dtype=float)
                hcarr[:,:] = hcarrtmp[0:j,:]

        # UPDATE THE MAIN EOS STRUCTURE WITH GATHERED INFORMATION
        # Add variables needed to calculate the ANEOS gamma function
        self.gamma0  = gamma0
        self.theta0  = theta0
        self.C24     = C24
        self.C60     = C60
        self.C61     = C61
        self.beta    = betagamma
        #
        # ANEOS.OUTPUT UNITS ARE NOT THE SAME AS THE SESAME TABLE!
        # add the vapor curve to this EOS object extracted from the ANEOS.OUTPUT
        #  TWO-PHASE BOUNDARIES
        #       T         RHOLIQ        RHOVAP        PLIQ         PVAP        ELIQ         EVAP         SLIQ         SVAP        GLIQ         GVAP         PSILIQ      PSIVAP         NTY
        #       K         kg/m**3       kg/m**3       GPa          GPa         J/kg         J/kg        J/kg-K       J/kg-K       J/kg         J/kg
        tmp = vcarr.shape
        #put vapor curve information in nicely named structure
        self.vc.NT  = tmp[0]
        self.vc.T   = vcarr[:,0] # K
        self.vc.rl  = vcarr[:,1]/1.E3 # g/cm3
        self.vc.rv  = vcarr[:,2]/1.E3 # g/cm3
        self.vc.Pl  = vcarr[:,3] # GPa
        self.vc.Pv  = vcarr[:,4] # GPa
        self.vc.Ul  = vcarr[:,5]/1.E6 # MJ/kg
        self.vc.Uv  = vcarr[:,6]/1.E6 # MJ/kg
        self.vc.Sl  = vcarr[:,7]/1.E6 # MJ/K/kg
        self.vc.Sv  = vcarr[:,8]/1.E6 # MJ/K/kg
        self.vc.Gl  = vcarr[:,9]/1.E6 # MJ/kg
        self.vc.Gv  = vcarr[:,10]/1.E6 # MJ/kg
        self.vc.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg, G MJ/kg'
        # np.interp wants increasing x values
        self.onebar.Tvap = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.T)) # extract boiling point temperature at 1 bar, K
        self.onebar.Siv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.Sl)) # extract liquid sp. entropy at 1 bar, MJ/K/kg
        self.onebar.Scv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.Sv)) # extract vapor sp. entropy at 1 bar, MJ/K/kg
        self.onebar.rhoiv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.rl)) # extract liquid density at 1 bar, g/cm3
        self.onebar.rhocv  = np.interp(1.E-4,np.flipud(self.vc.Pl),np.flipud(self.vc.rv)) # extract vapor density at 1 bar, g/cm3
        #
        # add the ANEOS Hugoniot to this EOS object extracted from the ANEOS.OUTPUT
        #      RHO          T           P          PC           E           S           V           U       RHO/RHOO  #IT  STATE
        #    kg/m**3        K          GPa        GPa          J/kg      J/kg-K       km/sec      km/sec
        #self.anhug.all = hcarr  # 2D array of Hugoniot variables
        tmp = hcarr.shape
        self.anhug.ND = tmp[0] # number of density points on the Hugoniot
        self.anhug.rho = hcarr[:,0]/1.E3 # g/cm3
        self.anhug.T   = hcarr[:,1] # K
        self.anhug.P   = hcarr[:,2] # GPa
        self.anhug.U   = hcarr[:,4]/1.E6 # MJ/kg
        self.anhug.S   = hcarr[:,5]/1.E6 # MJ/K/kg
        self.anhug.us  = hcarr[:,6] # km/s
        self.anhug.up  = hcarr[:,7] # km/s
        self.anhug.units = 'vars: rho g/cm3, T K, P GPa, U MJ/kg, S MJ/K/kg, Us km/s, Up km/s'
        #
        # Add melt curve to EOS objects if available
        # LIQUID/SOLID PHASE CURVE
        #       T         RLIQ       RSOLID      PLIQ       PSOLID      ELIQ        ESOLID       SLIQ       SOLID        GLIQ       GSOLID        #ITER
        #       K        kg/m**3     kg/m**3      GPa         GPa       J/kg         J/kg        J/kg-K     J/kg-K       J/kg        J/kg
        if meltcurve == 1:
            # put the melt curve information in nicely named structure
            tmp=mcarr.shape
            self.mc.NT  = tmp[0]
            self.mc.T   = mcarr[:,0] # K
            self.mc.rl  = mcarr[:,1]/1.E3 # g/cm3
            self.mc.rs  = mcarr[:,2]/1.E3 # g/cm3
            self.mc.Pl  = mcarr[:,3] # GPa
            self.mc.Ps  = mcarr[:,4] # GPa
            self.mc.Ul  = mcarr[:,5]/1.E6 # MJ/kg
            self.mc.Us  = mcarr[:,6]/1.E6 # MJ/kg
            self.mc.Sl  = mcarr[:,7]/1.E6 # MJ/K/kg
            self.mc.Ss  = mcarr[:,8]/1.E6 # MJ/K/kg
            self.mc.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg'
            # NOTE THAT TRIPLE POINT AND VAPOR CURVE SOLUTIONS DO NOT ALWAYS MATCH PERFECTLY AT THE TRIPLE POINT
            tmp = np.where(mcarr[:,3] > 0.)[0] # find the triple point first entry with positive pressure
            self.tp.T = mcarr[tmp[0],0] # K
            self.tp.P = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.Pv))   # this has trouble for forsterite; use the vapor size of the VC mcarr[:,3] # GPa
            self.tp.Sim  = mcarr[tmp[0],8]/1.E6 # extract solid sp. entropy at tp, MJ/K/kg
            self.tp.Scm  = mcarr[tmp[0],7]/1.E6 # extract liquid sp. entropy at tp, MJ/K/kg
            # np.interp wants x values increasing
            self.tp.Siv  = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.Sl)) # extract liquid sp. entropy at tp, MJ/K/kg
            self.tp.Scv  = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.Sv)) # extract vapor sp. entropy at tp, MJ/K/kg
            self.tp.rhol = mcarr[tmp[0],1]/1.E3 # extract liquid density at tp, g/cm3
            self.tp.rhos = mcarr[tmp[0],2]/1.E3 # extract solid density at tp, g/cm3
            self.tp.rhov = np.interp(self.tp.T,np.flipud(self.vc.T),np.flipud(self.vc.rv)) # extract vapor density at tp, g/cm3
            self.tp.units = 'T K, P GPa, S MJ/K/kg, rho g/cm3'
            # Extract melting point
            self.onebar.Tmelt = np.interp(1.E-4,self.mc.Pl,self.mc.T) # extract melting point temperature at 1 bar, K
            self.onebar.Sim   = np.interp(1.E-4,self.mc.Pl,self.mc.Ss) # extract liquid sp. entropy at 1 bar BP, MJ/K/kg
            self.onebar.Scm   = np.interp(1.E-4,self.mc.Pl,self.mc.Sl) # extract vapor sp. entropy at 1 bar BP, MJ/K/kg
            self.onebar.rhoim   = np.interp(1.E-4,self.mc.Pl[3::],self.mc.rs[3::]) # extract solid density at 1 bar MP, MJ/K/kg
            self.onebar.rhocm   = np.interp(1.E-4,self.mc.Pl[3::],self.mc.rl[3::]) # extract liquid density at 1 bar MP, MJ/K/kg
        # put the data for the critical point in the EOS structure for easy access
        self.cp.T   = vcarr[0,0] # K
        self.cp.rho = vcarr[0,1]/1.E3 # g/cm3
        self.cp.P   = vcarr[0,3] # GPa
        self.cp.U   = vcarr[0,5]/1.E6 # MJ/kg 
        self.cp.S   = vcarr[0,7]/1.E6 # MJ/K/kg
        self.cp.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg'
        #------------------------------------------------------------------

    def loadaneosidg(self, aneosinfname=None, aneosoutfname=None, silent=False):
        """Function for reading in IDEAL GAS MODEL
           ANEOS INPUT and OUTPUT FILE DATA into EOS structure."""
        if aneosinfname is None:
            return 'Must provide input file name.'
        if aneosoutfname is None:
            return 'Must provide output file name.'
        # function to gather data from ANEOS input and output files
        # SESAME FILE HEADER INFORMATION MUST BE LOADED INTO THE EOS STRUCTURE BEFORE CALLING THIS FUNCTION
        #
        # READ IN ANEOS INPUT FILE
        aneosinputfile = open(aneosinfname,"r")  
        testin=aneosinputfile.readlines()   # read in the whole ascii file at once because this is fatser
        aneosinputfile.close()
        # gather EOS information from the ANEOS.OUTPUT file
        aneosoutputfile = open(aneosoutfname,"r")  
        testout=aneosoutputfile.readlines() # read everything in at once because this is faster
        aneosoutputfile.close()
        if silent == False:
            print('Done loading ANEOS files.')

        # THIS CODE PARSES THE ANEOS.OUTPUT FILE INTO ARRAYS FOR USE IN PLOTTING/USING THE EOS
        if silent == False:
            print('ANEOS WAS CALLED WITH THE FOLLOWING INPUT, LOADED FROM FILE ',aneosinfname)
        # Gather parameters for the gamma function while printing the ANEOS INPUT FILE
        aneoscount=1
        for i in np.arange(len(testin)):
            if testin[i].find('ANEOS') == 0:
                if aneoscount<9:
                    if silent == False:
                        print(' '+testin[i-3],testin[i-2],testin[i-1],testin[i])
                    aneoscount=aneoscount+1
                else:
                    if silent == False:
                        print(' '+testin[i])
                if testin[i].find('ANEOS2') == 0:
                    tmp=testin[i]
                    nelem=int(tmp[10:20])
                    #print('nelem=',nelem)
                    eostype=int(tmp[20:30])
                    #print('eostype=',eostype)
                    rho0=float(tmp[30:40])
                    #print('rho0=',rho0)
                    gamma0=float(tmp[70:80])+1 # input was cp/cv-1
                    #print('gamma0=',gamma0)

        # some checks
        if rho0 != self.R0REF:
            print('WARNING: rho0 does not match. ')
            #assert(False) # just a way to stop the notebook
            #print('WARNING: rho0 does not match. STOPPING THIS NOTEBOOK.')
            #assert(False) # just a way to stop the notebook

        # GUESS A BIG ARRAY SIZE FOR THE PHASE BOUNDARIES AND HUGONIOT IN ANEOS.OUTPUT
        # the melt curve, vapor curve and Hugoniot curves are not fixed length outputs
        nleninit=300
        meltcurve = 0

        if silent == False:
            print('READING DATA FROM ANEOS OUTPUT FILE ',aneosoutfname)

        # Read in data from the ANEOS.OUTPUT FILE
        imc = -1 # flag for no melt curve in the model
        for i in np.arange(len(testout)):
            if testout[i].find('  Data for ANEOS number') == 0:
                tmp = testout[i+2][0:50]
                eosname = tmp.strip()
            if testout[i] == '   HUGONIOT\n':
                nhc = nleninit
                ihc = i
                hcarrtmp = np.zeros((nhc,9),dtype=float)
                flag=0
                j=0
#                while flag == 0:
                while j<22:
                    tmp=str.replace(testout[j+i+5],'D','E')
                    #print(j,i,5,tmp)
                    tmp3 = tmp[0:109]
                    tmp4 = list(tmp3.split())
                    if len(tmp4) > 0:
                        tmp4[3]='0.0' # this column often gives problems with exponential notation so don't read it
                        tmp5 = np.asarray(tmp4)
                        hcarrtmp[j,:] = tmp5[:]
                        j=j+1
                    else:
                        flag=1
                hcarr = np.zeros((j,9),dtype=float)
                hcarr[:,:] = hcarrtmp[0:j,:]

        # UPDATE THE MAIN EOS STRUCTURE WITH GATHERED INFORMATION
        # Add variables needed to calculate the ANEOS gamma function
        self.gamma0  = gamma0
        #self.theta0  = theta0
        #self.C24     = C24
        #self.C60     = C60
        #self.C61     = C61
        #self.beta    = betagamma
        #
        # add the ANEOS Hugoniot to this EOS object extracted from the ANEOS.OUTPUT
        #      RHO          T           P          PC           E           S           V           U       RHO/RHOO  #IT  STATE
        #    kg/m**3        K          GPa        GPa          J/kg      J/kg-K       km/sec      km/sec
        #self.anhug.all = hcarr  # 2D array of Hugoniot variables
        tmp = hcarr.shape
        self.anhug.ND = tmp[0] # number of density points on the Hugoniot
        self.anhug.rho = hcarr[:,0]/1.E3 # g/cm3
        self.anhug.T   = hcarr[:,1] # K
        self.anhug.P   = hcarr[:,2] # GPa
        self.anhug.U   = hcarr[:,4]/1.E6 # MJ/kg
        self.anhug.S   = hcarr[:,5]/1.E6 # MJ/K/kg
        self.anhug.us  = hcarr[:,6] # km/s
        self.anhug.up  = hcarr[:,7] # km/s
        self.anhug.units = 'vars: rho g/cm3, T K, P GPa, U MJ/kg, S MJ/K/kg, Us km/s, Up km/s'
        #
        #------------------------------------------------------------------
        
###########################################################################################
########## GADGET STYLE TABLES
###########################################################################################
class extGADtable:
    """Class for v2019 GADGET EOS tables using SESAME UNITS"""
    def __init__(self):
        """Function to initialize GADGET TABULAR EOS structure"""
        self.ND  = 0 # integer; number of density points in grid
        self.NS  = 0 # integer; number of sp. entropy points in grid
        self.rho = np.zeros(self.ND)          # g/cm3, density values
        self.S   = np.zeros(self.NS)          # MJ/K/kg, sp. entropy values
        self.P   = np.zeros(self.ND*self.NS)  # GPA, pressure(S,rho)
        self.T   = np.zeros(self.ND*self.NS)  # K, temperature(S,rho)
        self.U   = np.zeros(self.ND*self.NS)  # MJ/kg, sp. internal energy(S,rho)
        self.A   = np.zeros(self.ND*self.NS)  # MJ/kg, Helmholtz free energy(S,rho)
        self.cs  = np.zeros(self.ND*self.NS)  # cm/s, sound speed(S,rho)
        self.cv  = np.zeros(self.ND*self.NS)  # MJ/K/kg, sp. heat capacity(S,rho)
        self.KPA = np.zeros(self.ND*self.NS)  # integer, ANEOS KPA flag(S,rho)
        self.MDQ = np.zeros(self.ND*self.NS)  # integer, Model Development Quality Flag(S,rho)
        self.units = 'Units: g/cm3, MJ/K/kg, GPa, K, MJ/kg, MJ/kg, cm/s, MJ/K/kg, KPA flag. 2D arrays are (NS,ND).'
        self.MODELNAME = ''
    #
    def view(self, q='P', Slow=None, Shigh=None, rholow=None, rhohigh=None):
        """Function to print variables from GADGET EOS table"""
        if Slow is None:
            Slow = self.S.min()
        if Shigh is None:
            Shigh = self.S.max()
        if rholow is None:
            rholow = self.rho.min()
        if rhohigh is None:
            rhohigh = self.rho.max()
        if q == 'P':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('P:', (self.P[np.logical_and(self.S >= Slow,
                                                self.S<=Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'U':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('U:', (self.U[np.logical_and(self.S >= Slow,
                                                self.S <= Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'A':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('A:', (self.A[np.logical_and(self.S >= Slow,
                                                self.S <= Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'T':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho
                                   >= rholow,self.rho<=rhohigh)
                                  ])
            print('T:', (self.T[np.logical_and(self.S >= Slow,self.S <= Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])      
        if q == 'cs':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cs:', (self.cs[np.logical_and(self.S >= Slow,
                                                  self.S <= Shigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'cv':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cv:', (self.cv[np.logical_and(self.S >= Slow,
                                                  self.S <= Shigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'KPA':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('KPA:', (self.KPA[np.logical_and(self.S >= Slow,
                                                  self.S <= Shigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
    #
    def writestdgadget(self, writestdgadgetfname=None):
        """Function to write an standard gadget2 EOS ascii file in 301 format"""
        #sesfile = open("NEW-GADGET-STD-NOTENSION.TXT","w")  
        if writestdgadgetfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writestdgadgetfname,"w")  
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NS))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     sp. entropy array MJ/K/kg -> erg/K/g
        for j in range(0, int(self.NS)):
            sesfile.write("{:16.8e}".format(self.S[j]*1.E10))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  pressure array GPa -> dynes/cm2 P[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.P[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  temperature array K T[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.T[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  specific internal energy array MJ/kg -> erg/g U[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.U[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  sound speed array cm/s cs[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.cs[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the STD GADGET TABLE FILE
        sesfile.close() 
        print('Done writing the STD GADGET TABLE FILE with notension to local directory: ',writestdgadgetfname)
    #
    def writeextgadget(self, writeextgadgetfname=None):
        """Function to write an extended gadget2 EOS ascii file in 301 format"""
        #sesfile = open("NEW-GADGET-EXT-NOTENSION.TXT","w")  
        if writeextgadgetfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writeextgadgetfname,"w")  
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NS))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     sp. entropy array MJ/K/kg -> erg/K/g
        for j in range(0, int(self.NS)):
            sesfile.write("{:16.8e}".format(self.S[j]*1.E10))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        # Helmholtz free energy array in MJ/kg -> ergs/g A[tempindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.A[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # sp. heat capacity array MJ/K/kg -> erg/K/g cv[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.cv[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # KPA array flag k[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.KPA[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # MDQ array flag MDQ[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.MDQ[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the EXTENDED GADGET TABLE FILE
        sesfile.close() 
        print('Done writing the EXTENDED GADGET TABLE FILE with notension to local directory: ',writeextgadgetfname)
    #

###########################################################################################
########## TILLOTSON TABLE CLASS
###########################################################################################
class TillotsonClass:
    """Class for tabulated Tillotson EOS tables using SESAME UNITS"""
    def __init__(self):
        """Function to initialize Tillotson TABULAR EOS structure"""
        # variables stored in SESAME units
        self.rho0  = 0. # reference density g/m3
        self.E0    = 0. # E0 constant MJ/kg
        self.EIV   = 0. # specific internal energy of incipient vaporization MJ/kg
        self.ECV   = 0. # specific internal energy of complete vaporization MJ/kg
        self.AA    = 0. # Bulk modulus K0 GPa
        self.BB    = 0. # B constant GPa
        self.a     = 0. # a constant [-]
        self.b     = 0. # b constant [-]
        self.alpha = 0. # alpha constant [-]
        self.beta  = 0. # beta constant [-]    a+b=Gruneisen gamma at rho0
        self.cv    = 0. # constant specific heat capacity MJ/K/kg
        self.t0    = 298. # K
        self.params = [] # for passing the parameter list to subroutines
        self.ND    = 0  # integer; number of density points in grid
        self.NU    = 0  # integer; number of sp. internal energy points in grid
        self.rho   = np.zeros(self.ND)          # g/cm3, density values
        self.U     = np.zeros(self.NU)          # MJ/kg, sp. internal energy values
        self.P     = np.zeros(self.ND*self.NU)  # GPA, pressure(U,rho)
        self.cs    = np.zeros(self.ND*self.NU)  # km/s, sound speed(U,rho)
        self.region= np.zeros(self.ND*self.NU)  # region flag: 1-condensed, 2-interpolated, 3-expanded
        self.T     = np.zeros(self.ND*self.NU)  # K temperature
        self.units = 'Units: g/cm3, MJ/kg, GPa, cm/s, K. 2D arrays are (NU,ND).'
        self.MODELNAME = ''
        self.hug   = EOShugoniot()
        self.ecold = EOScurve() # energy as a function of density (rho>rho0) with initial condition of E=0 at rho0
    def assignparams(self,tilleos):
        # index numbers for material parameters in the tilleos array
        r0 = 0 
        E0 = 1
        EIV = 2
        ECV=3
        AA=4
        BB=5
        a=6
        b=7
        alpha=8
        beta=9
        cv=10 
        self.params = tilleos # in MKS
        self.rho0  = tilleos[r0]/1.e3 # reference density g/m3
        self.E0    = tilleos[E0]/1.e6 # E0 constant MJ/kg
        self.EIV   = tilleos[EIV]/1.e6 # specific internal energy of incipient vaporization MJ/kg
        self.ECV   = tilleos[ECV]/1.e6 # specific internal energy of complete vaporization MJ/kg
        self.AA    = tilleos[AA]/1.e9 # Bulk modulus K0 GPa
        self.BB    = tilleos[BB]/1.e9 # B constant GPa
        self.a     = tilleos[a] # a constant [-]
        self.b     = tilleos[b] # b constant [-]
        self.alpha = tilleos[alpha] # alpha constant [-]
        self.beta  = tilleos[beta] # beta constant [-]    a+b=Gruneisen gamma at rho0
        self.cv    = tilleos[cv]/1.e6 # constant specific heat capacity MJ/K/kg 
    def calcecold(self):
        """ Calculate the internal energy of the cold curve. Passed to Tillotson EOS function with params for temperature. """
        # start by compressing a factor of 3 with 1000 steps
        # Use params in mks units
        nsteps = 100
        self.ecold.rho = np.arange(nsteps)/nsteps*(2.*self.params[0])+self.params[0] # rho is rho0 to 3*rho0
        self.ecold.U = np.zeros(nsteps) # initial energy is zero by definition
        self.ecold.P = np.zeros(nsteps) # initial pressure is zero by definition
        self.ecold.T = np.zeros(nsteps)+self.t0 # initial temperature is T0 by definition
        h = self.ecold.rho[1]-self.ecold.rho[0]
        x0 = self.ecold.rho[0]
        y  = 0.0
        for i in range(1, nsteps):
            #"Apply Runge Kutta Formulas to find next value of y"
            # Rundage 2015: dE/drho = P(rho,E)/rho^2
            # Till_P returns [pout,flag,csout,tout]
            k1 = h * Till_P(x0,y,self.params,self.ecold,calctemp=False)[0]/x0/x0
            k2 = h * Till_P(x0 + 0.5 * h,y + 0.5 * k1,self.params,self.ecold,calctemp=False)[0]/(x0 + 0.5 * h)/(x0 + 0.5 * h)
            k3 = h * Till_P(x0 + 0.5 * h,y + 0.5 * k2,self.params,self.ecold,calctemp=False)[0]/(x0 + 0.5 * h)/(x0 + 0.5 * h)
            k4 = h * Till_P(x0 + h,y + k3,self.params,self.ecold,calctemp=False)[0]/(x0+h)/(x0+h)
            #k1 = h * dydx(x0, y)
            #k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1)
            #k3 = h * dydx(x0 + 0.5 * h, y + 0.5 * k2)
            #k4 = h * dydx(x0 + h, y + k3)

            # Update next value of y
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
            
            # Update next value of x
            x0 = x0 + h
            
            # update structure arrays
            self.ecold.rho[i] = x0
            self.ecold.U[i]   = y
            self.ecold.P[i]   = Till_P(x0,y,self.params,self.ecold,calctemp=False)[0]
            #print('Tillotson cold curve: ',i,x0,y,Till_P(x0,y,self.params,self.ecold,calctemp=False)[0])
        return
    def FillTable(self, matparams=None, modelname=None):
        """Function to populate table with Tillotson EOS. Assumes ND,NU,rho,and U are filled."""
        if matparams is None:
            print('Must provide a material parameter list [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta, cv].')
            exit(0)
        # matparams: [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta, cv]
        # units:     [kg/m3, J/kg, J/kg, J/kg, Pa, Pa, [-]x4, J/K/kg]
        # print(self.NU,self.ND)
        self.cs = np.zeros((self.NU,self.ND))
        self.P  = np.zeros((self.NU,self.ND))
        self.T  = np.zeros((self.NU,self.ND))
        self.region = np.zeros((self.NU,self.ND))
#        print('Using the iSALE Tillotson implementation.')
        print('Using standard Tillotson implementation (e.g., Hosono et al. 2019, Melosh 1989).')
        print('Tillotson parameters: ',matparams)
        for ie in range(0,self.NU):
            for ir in range(0,self.ND):
                flag               = 0
                eng                = self.U[ie]*1.e6 # J/kg
                dens               = self.rho[ir]*1000. # kg/m3
                tmp                = Till_P_Hosono(dens,eng,matparams) # returns [P Pa, flag]
                self.P[ie,ir]      = tmp[0]/1.e9 # Pa to GPa
                self.region[ie,ir] = tmp[1] # cold, hot, interpolated regions
                self.cs[ie,ir]     = Till_SoundSpeed(dens,eng,matparams)/1.E3 # km/s 

    def FillTable_iSALE(self, matparams=None, modelname=None):
        """Function to populate table with Tillotson EOS. Assumes ND,NU,rho,and U are filled."""
        if matparams is None:
            print('Must provide a material parameter list [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta, cv].')
            exit(0)
        # matparams: [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta, cv]
        # units:     [kg/m3, J/kg, J/kg, J/kg, Pa, Pa, [-]x4, J/K/kg]
        # print(self.NU,self.ND)
        self.cs = np.zeros((self.NU,self.ND))
        self.P  = np.zeros((self.NU,self.ND))
        self.T  = np.zeros((self.NU,self.ND))
        self.region = np.zeros((self.NU,self.ND))
        print('Using the iSALE Tillotson implementation.')
#        print('Using standard Tillotson implementation (e.g., Hosono et al. 2019, Melosh 1989).')
        print('Tillotson parameters: ',matparams)
        # Till_P is expecting mks units
        for ie in range(0,self.NU):
            for ir in range(0,self.ND):
                flag               = 0
                eng                = self.U[ie]*1.e6 # J/kg
                dens               = self.rho[ir]*1000. # kg/m3
                tmp                = Till_P(dens,eng,matparams,self.ecold,calctemp=True) # returns [P Pa, flag]
                self.P[ie,ir]      = tmp[0]/1.e9 # Pa to GPa
                self.region[ie,ir] = tmp[1] # cold, hot, interpolated regions
                self.cs[ie,ir]     = tmp[2]/1.e3 # km/s
                self.T[ie,ir]      = tmp[3] # K
                
    def FillTable_ea(self, matparams=None, modelname=None):
        """Function to populate table with Tillotson EOS. Assumes ND,NU,rho,and U are filled."""
        if matparams is None:
            print('Must provide a material parameter list [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta, cv].')
            exit(0)
        # matparams: [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta, cv]
        # units:     [kg/m3, J/kg, J/kg, J/kg, Pa, Pa, [-]x4, J/K/kg]
        # print(self.NU,self.ND)
        self.cs = np.zeros((self.NU,self.ND))
        self.P  = np.zeros((self.NU,self.ND))
        self.T  = np.zeros((self.NU,self.ND))
        self.region = np.zeros((self.NU,self.ND))
        print('Using the iSALE Asphaug implementation.')
        print('Tillotson parameters: ',matparams)
        for ie in range(0,self.NU):
            for ir in range(0,self.ND):
                flag               = 0
                eng                = self.U[ie]*1.e6 # J/kg
                dens               = self.rho[ir]*1000. # kg/m3
                tmp                = tilleos_ea(dens,eng,matparams) # returns [P Pa, flag]
                self.P[ie,ir]      = tmp[0]/1.e9 # Pa to GPa
                self.region[ie,ir] = tmp[1] # cold, hot, interpolated regions
                self.cs[ie,ir]     = tmp[2]/1.e3 # km/s

    def calchugoniot(self, r00=None, e0=None, pmax=None, writefilename=None):
        """Function for calculating a Hugoniot from Tillotson EOS table."""
        if r00 is None:
            return 'Must provide r00 in g/cm^3 and e0 in MJ/kg.'
        if e0 is None:
            return 'Must provide r0 in g/cm^3 and e0 in MJ/kg.'
        if pmax is None:
            pmax=1.E4 # GPa
        self.hug.rho = []
        self.hug.P = []
        self.hug.T = []
        self.hug.U = []
        self.hug.up = []
        self.hug.us = []
        self.hug.cs = [] 
        self.hug.region = [] # Tillotson region flag

        ie0 = int(np.round(np.interp(e0,self.U,np.arange(self.NU)))) # uses nearest value if e0 not in array
        ir0 = int(np.round(np.interp(r00,self.rho,np.arange(self.ND)))) # uses nearest value if r0 not in the array
        print('Hug ie0,ir0=',ie0,ir0)
        r0  = self.rho[ir0] # g/cm^3
        p0  = self.P[ie0,ir0] # GPa
        e0  = self.U[ie0] # MJ/kg
        t0  = self.t0
        up0 = 0. # no initial particle velocity km/s
        us0 = self.cs[ie0,ir0] # km/s use sound velocity for initial
        cs0 = self.cs[ie0,ir0] # km/s use sound velocity for initial
        #print(ir0,ie0,r0,p0,e0,up0,us0)
        self.hug.rho = np.append(self.hug.rho, self.rho[ir0])
        self.hug.P = np.append(self.hug.P, p0)
        self.hug.U = np.append(self.hug.U, e0)
        self.hug.up = np.append(self.hug.up, up0)
        self.hug.us = np.append(self.hug.us, us0)
        self.hug.cs = np.append(self.hug.cs, cs0)
        self.hug.T = np.append(self.hug.T, t0)
        #iir=ir0+1
        #pnew=p0
        #while pnew<pmax:
        #print(iir,pnew,self.P[ie0::,iir],r0,self.rho[iir],e0,self.U[ie0::])
        for iir in range(ir0+1,self.ND):
            ediff =0.5*(self.P[ie0::,iir]+p0)*(1./r0-1./self.rho[iir])+e0 -(self.U[ie0::])  # MJ/kg
            # np.interp wants x values increasing
            pnew = np.interp(0.,np.flip(ediff),np.flip(self.P[ie0::,iir]))
            enew = np.interp(0.,np.flip(ediff),np.flip(self.U[ie0::]))
            csnew = np.interp(0.,np.flip(ediff),np.flip(self.cs[ie0::,iir]))
            tnew = np.interp(0.,np.flip(ediff),np.flip(self.T[ie0::,iir]))
            upnew = np.sqrt((pnew-p0)*(1./r0-1./self.rho[iir]))
            usnew = (1./r0)*np.sqrt((pnew-p0)/(1./r0-1./self.rho[iir]))
            #print(self.rho[iir],pnew,enew,upnew,usnew)
            self.hug.rho = np.append(self.hug.rho, self.rho[iir])
            self.hug.P = np.append(self.hug.P, pnew)
            self.hug.U = np.append(self.hug.U, enew)
            self.hug.up = np.append(self.hug.up, upnew)
            self.hug.us = np.append(self.hug.us, usnew)
            self.hug.cs = np.append(self.hug.cs, csnew)
            self.hug.T = np.append(self.hug.T, tnew)
            #print(iir,pnew)
            #iir += 1
        self.hug.NH=len(self.hug.P)
        self.hug.units='units: rho g/cm3, P GPa, U MJ/kg, Up km/s, Us km/s, cs km/s'
        print('Done calculating Hugoniot with Tillotson EOS table.')
    def __str__(self):
        """ Print the Tillotson EOS parameters """
        return f'\n{self.MODELNAME} Tillotson EOS parameters [SESAME Units]: \n' + \
               f'   rhoref: {self.rho0} \n' + \
               f'   a:      {self.a} \n' + \
               f'   b:      {self.b} \n' + \
               f'   AA:     {self.AA} \n' + \
               f'   BB:     {self.BB} \n' + \
               f'   E0:     {self.E0} \n' + \
               f'   alpha:  {self.alpha} \n' + \
               f'   beta:   {self.beta} \n' + \
               f'   Eiv:    {self.EIV} \n' + \
               f'   Ecv:    {self.ECV} \n' + \
               f'   cv:     {self.cv} \n' + \
               f' params [mks]: {self.params} \n'
#
#
###########################################################################################
########## TILLOTSON EOS FUNCTIONS
###########################################################################################
# TILLOTSON FUNCTIONS AS IMPLEMENTED BY HOSONO
# dunite tillotson parameters used by Hosono et al. 2019 
# these olivine parameters are from Marinova et al. 2011 Icarus 
# parameters: [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta]
# units:    [kg/m3, J/kg, J/kg, J/kg, Pa, Pa, [-]x4]
# dunitetill = [3500.0, 550.0e+6, 4.500e+6, 14.50e+6, 131.00e+9,  49.00e+9, 0.5, 1.4, 5.0, 5.0]
# Basalt parameters from iSALE -- from where? Benz?
# basalttill = [2650.0, 4.87E8, 4.72E6, 18.2E6, 5.3E10, 5.3E10, 0.6, 0.6, 5., 5.]
#
def Till_P_co(dens,eng,tilleos):
    r0 = 0
    E0 = 1
    EIV = 2
    ECV=3
    AA=4
    BB=5
    a=6
    b=7
    alpha=8
    beta=9
    # Hosono implementation - straight up equation
    eta  = dens/tilleos[r0]
    mu   = eta - 1.0
    pco = (tilleos[a] + tilleos[b] / (eng / tilleos[E0] / eta / eta + 1.0)) * dens * eng + \
         tilleos[AA] * mu + tilleos[BB] * mu * mu

    return pco
#
def Till_P_ex(dens,eng,tilleos):
    import numpy as np
    r0 = 0
    E0 = 1
    EIV = 2
    ECV=3
    AA=4
    BB=5
    a=6
    b=7
    alpha=8
    beta=9
    # hosono equation
    eta  = dens/tilleos[r0]
    mu   = eta - 1.0
    pex = tilleos[a] * dens * eng + \
          (tilleos[b] * dens * eng / (eng / tilleos[E0] / eta / eta + 1.0) + \
           tilleos[AA] * mu * np.exp(- tilleos[beta] * (1.0 / eta - 1.0))) * \
                                       np.exp(- tilleos[alpha] * (1.0 / eta - 1.0) * (1.0 / eta - 1.0))
    return pex
#
def Till_P_Hosono(dens,eng,tilleos):
    flag=0 # flag for the region of the EOS
    r0 = 0 # index numbers for variables in the tilleos array
    E0 = 1
    EIV = 2
    ECV=3
    AA=4
    BB=5
    a=6
    b=7
    alpha=8
    beta=9

    # Hosono has 2 options
    pmin=1.e-7
    
    # option 1
    if ( (dens>tilleos[r0]) or (eng < tilleos[EIV]) ):
        # condensed region
        flag=1
        pco = Till_P_co(dens,eng,tilleos)
        #if pco < pmin:
        #    pout = pmin
        #else:
        #    pout = pco
        # allow tension
        pout = pco
        return [pout,flag]
    else:
        if ( (dens < tilleos[r0]) and (eng >= tilleos[ECV]) ):
            # expanded region
            flag=3
            pex = Till_P_ex(dens,eng,tilleos)
            #if pex < pmin:
            #    pout = pmin
            #else:
            #    pout = pex
            pout=pex
            return [pout,flag]
        else:
            # interpolated region
            flag=2
            Pex = Till_P_ex(dens,eng,tilleos)
            Pco = Till_P_co(dens,eng,tilleos)
            # nobody seem to do this step
            #if Pco < pmin:
            #    Pco = pmin
            Pint = ((eng-tilleos[EIV])*Pex+(tilleos[ECV]-eng)*Pco)/(tilleos[ECV]-tilleos[EIV])
            #if Pint < pmin:
            #    pout = pmin
            #else:
            #    pout = Pint
            pout=Pint
            return [pout,flag]
    '''    
    # option 2
    if ( (dens>=tilleos[r0]) or (eng < tilleos[EIV]) ):
        # condensed region
        flag=1
        pout = Till_P_co(dens,eng,tilleos)
        if (dens <= 0.9*tilleos[r0]):
            pout = 1.e-16
            return [pout,flag]
    else:
        if ( (dens < tilleos[r0]) and (eng > tilleos[ECV]) ):
            # expanded region
            flag=3
            pout = Till_P_ex(dens,eng,tilleos)
        else:
            # interpolated region
            flag=2
            Pex = Till_P_ex(dens,eng,tilleos)
            Pco = Till_P_co(dens,eng,tilleos)
            pout = ((eng-tilleos[EIV])*Pex+(tilleos[ECV]-eng)*Pco)/(tilleos[ECV]-tilleos[EIV])
    if (pout < pmin):
        pout = pmin
    return [pout,flag]
    '''
    # END HOSONO TILLOTSON FUNCTION
def Till_dPdrho(dens,eng,tilleos):
    drho=0.0001
    a = Till_P_Hosono(dens+drho,eng,tilleos)
    b = Till_P_Hosono(dens-drho,eng,tilleos)
    return  (a[0] - b[0])/(2.*drho)

def Till_dPdu(dens,eng,tilleos):
    ddu=0.0001
    a=Till_P_Hosono(dens,eng+ddu,tilleos)
    b=Till_P_Hosono(dens,eng-ddu,tilleos)
    return (a[0] - b[0])/(2.*ddu)

def Till_SoundSpeed(dens,eng,tilleos):
    #if dens<= 100.:
    #    cs = 1.e-8 # m/s
    #else:
    ptil = Till_P_Hosono(dens,eng,tilleos)[0]
    cs2   = ptil/(dens*dens)*Till_dPdu(dens,eng,tilleos) + Till_dPdrho(dens,eng,tilleos)
    if cs2 < 0.:
        cs2 = np.abs(ptil/dens) # better to use ideal gas approx
    cs = np.sqrt(cs2)
    return cs
# END HOSONO FUNCTIONS
# ------------------------------
# iSALE Dellon implementation of Tillotson --> now substantially modified
# added cold curve and temperature calcs 2023-06-30 STS
# region logic fixed and temp variables removed sts 7/3/2023

def Till_P(rin,ein,tilleos,ecold,calctemp=True):
    """ tilleos is the params list for the material
    ecold is an EOScurve object with the cold curve at t0
    returns [pout,flag,csout,tout]
    """
    # output variables
    flag  = 0  # flag for the region of the EOS
    pout  = 0.
    csout = 0.
    tout  = 0.
    #
    # iSALE implementation of Tillotson
    # index numbers for material parameters in the tilleos array
    # should all be in code units/self-consistent units
    r0    = 0 
    E0    = 1
    EIV   = 2
    ECV   = 3
    AA    = 4
    BB    = 5
    a     = 6
    b     = 7
    alpha = 8
    beta  = 9
    cv    = 10 # cv is converted to code units eu/K/(original volume cm3)
    #
    # if density is zero or negative return 0 pressure
    if (rin <= 0.):
        pout  = 0.
        csout = 0.
        tout = 0.
        flag  = -1
        return [pout,flag,csout,tout]
    #
    # define intermediate variables
    eta  = rin/tilleos[r0]
    mu   = eta - 1.0
    imu  = tilleos[r0]/rin - 1.0
    #temp = tilleos[E0] * np.power(eta,2.)
    erel = ein/(tilleos[E0] * np.power(eta,2.))
    #if (temp > 0.):
    #    erel = ein/temp
    #else:
    #    # avoids divide by zero
    #    erel = 0.
    ierel = 1./(erel + 1.)
    grun  = tilleos[a] + tilleos[b]*ierel
    gterm = tilleos[a] * rin * ein
    eterm = tilleos[b] * rin * ein * ierel
    #
    # initialize hot and cold terms to the pressure
    ph = 0.
    pc = 0.
    #
    #if ( (eta < 1.) & (ein >= tilleos[ECV]) ):
    # calculate ph for expanded region above ECV
    #flag = 3
    #tmp  = tilleos[alpha]*np.power(imu,2.)
    exp_a = np.exp(-tilleos[alpha]*np.power(imu,2.))
    #if (tmp < 100.): # why 100? need to make some plots on this limit
    #    exp_a = np.exp(-tmp)
    #else:
    #    exp_a = 0.
    #tmp = tilleos[beta]*imu
    exp_b = np.exp(-tilleos[beta]*imu)
    #if (tmp < 100.): # check this limit
    #    exp_b = np.exp(-tmp)
    #else:
    #    exp_b = 0.
    # compute "hot pressure"
    ph  = gterm + (eterm + tilleos[AA]*mu*exp_b)*exp_a  
    #
    # compute sound speed from derivatives
    dpre2   = tilleos[b]*ein*ierel*(1.+2.*(tilleos[alpha]*imu/eta + erel*ierel))*exp_a
    dpre3   = (tilleos[AA]/rin)*(1.+(mu/np.power(eta,2.))*(tilleos[beta]+2.*tilleos[alpha]*imu))*exp_a*exp_b
    dpreh   = tilleos[a]*ein+dpre2+dpre3
    dperh   = rin*(tilleos[a]+tilleos[b]*ierel*(1.-erel*ierel)*exp_a)
    c_dperh = (ph/np.power(rin,2.))*dperh
    cs2_h   = dpreh + c_dperh
    cs2     = cs2_h 

    #if ( (eta >= 1.) or (ein < tilleos[EIV]) ):
    # compressed region or expanded below energy of complete vaporization
    # should update this flag and calculate rho_IV for proper use of the low energy expansion
    if (eta >= 1.):
        #flag=1
        # low energy compression
        b_temp = tilleos[BB] # in compression without cubic 
    else:
        #flag=4
        # low energy expansion
        # different versions are implemented
        # Brundage 2013 uses zero
        # iSALE Dellen uses AA or suggests BB if BB<AA
        # obviously needs some comparisons to data to know what to do here
        # b_temp = tilleos[AA] # in expansion
        b_temp = 0.0 # using Brundage option for simplicity at this point
    #pcold = (tilleos[AA] * mu) + (b_temp * np.power(mu,2.))
    pcold = (tilleos[AA] * mu) + (tilleos[BB] * np.power(mu,2.))
    pc    = gterm + eterm + pcold
    # allow tension for pyKO
    #if (pc < 0.):
    #    pc=1.e-10
    #
    # compute sound speed from derivatives
    dprec   = (tilleos[AA]+2.*b_temp*mu)/tilleos[r0] + ein*(grun+2.*tilleos[b]*erel*ierel*ierel)
    dperc   = rin * (grun-tilleos[b]*erel*ierel*ierel)
    c_dperc = (pc/(np.power(rin,2.)))*dperc
    cs2_c   = dprec + c_dperc
    cs2     = cs2_c
    
    if calctemp:
        # Estimate the temperature for rho>rho0
        # extract the cold curve energy for the current density
        # ecold needs to be in same units as this function to work
        ecold_rin = np.interp(rin,ecold.rho,ecold.U)
        tout=ecold.T[0]+(ein-ecold_rin)/tilleos[cv]
        #print(rin,ein,ecold_rin,tilleos[cv],tout)

    # Figure out the right region
    if ( (eta < 1.) & (ein >= tilleos[ECV]) ):
        # expanded region above ECV
        flag = 3
        cs2  = cs2_h 
        pout = ph
    if ( (eta >= 1.) or (ein < tilleos[EIV]) ):
        cs2  = cs2_c
        pout = pc
        if (eta >= 1.):
            flag=1
        else:
            flag=4
    #
    if (ein > tilleos[EIV]) & (ein < tilleos[ECV]) & (eta < 1.):
        # interpolated region
        flag=2
        pout = ((ein-tilleos[EIV])*ph+(tilleos[ECV]-ein)*pc)/(tilleos[ECV]-tilleos[EIV])
        # transition sound speed
        dprem = ((ein-tilleos[EIV])*dpreh + (tilleos[ECV]-ein)*dprec) / (tilleos[ECV]-tilleos[EIV])
        dperm = ((ein-tilleos[EIV])*dperh + (tilleos[ECV]-ein)*dperc+ph-pc) / (tilleos[ECV]-tilleos[EIV])
        cs2_m = dprem+(pout/np.power(rin,2.))*dperm
        cs2   = cs2_m
    #
    # limiting after tillotson EOS
    # allow negative pressures for pyKO
    #if (pout <= 1.e-10):
    #    pout = 1.e-10 # pmin
    #    cs2  = tilleos[AA]/tilleos[r0]
    #
    if cs2>0:
        csout = np.sqrt(cs2)
    else:
        #csout = np.sqrt(tilleos[AA]/tilleos[r0])
        # use ideal gas approximation when low density limit?
        csout = np.sqrt(np.abs(pout/rin))
    return [pout,flag,csout,tout]
    # END OF TILLOTSON EOS ISALE IMPLEMENTATION FUNCTIONS
##------------------------------
## Erik Asphaug Tillotson Routine converted from tilleos.f 7/3/2023
## Did not have temperature implemented
## Using this to cross-check the sound speed derivates with iSALE version
def tilleos_ea(rin,ein,tilleos):
    # tilleos is the params list for the material
    # returns [pout,flag,csout]
    # initialize output variables as zeros
    # Tillotson region flag: 1-condensed, 2-interpolated, 3-expanded, 4-low energy expansion
    flag  = 0  # flag for the region of the EOS
    pout  = 0. # P(rin,ein)
    csout = 0. # cs(rin,ein)

    # Index numbers for tilleos array; cv=10 not used here
    r0    = 0 # index numbers for material parameters in the tilleos array
    E0    = 1
    EIV   = 2
    ECV   = 3
    AA    = 4
    BB    = 5
    a     = 6
    b     = 7
    alpha = 8
    beta  = 9
    
    # assign Tillotson parameters into Asphaug variable names
    rozn    = tilleos[r0]    # rozero[imat]
    ezn     = tilleos[E0]    # uzero[imat]
    capan   = tilleos[AA]    # capa[imat]
    esn     = tilleos[EIV]   # es[imat]
    espn    = tilleos[ECV]   # esp[imat]
    difesn  = espn-esn
    capbn   = tilleos[BB]    # capb[imat]
    smallaj = tilleos[a]     # smalla[imat]
    smallbj = tilleos[b]     # smallb[imat]
    alphaj  = tilleos[alpha] # alpha[imat]
    betaj   = tilleos[beta]  # beta[imat]
    
    ui    = ein
    rhoi  = rin
    rhoi2 = rin*rin
    eta   = rhoi/rozn
    xmu   = eta - 1.
    rap   = rozn/rhoi
    
    # condensed phase, compute pressure and sound speed
    p1 = smallbj/(ui/(ezn*eta*eta) + 1.)
    pressc = (smallaj + p1) * ui * rhoi + capan * xmu + capbn * xmu * xmu

    p1 = ui/(ezn*eta*eta) + 1.
    c2 = p1*p1
    c3 = p1 - 1.
    c4 = smallbj * ui * (3. * c3 + 1.) / c2
    dpdi = smallaj * rhoi + smallbj * rhoi / c2
    dpdrho = smallaj * ui + capan / rozn + 2. * capbn * xmu / rozn + c4
    cvelc = dpdrho + dpdi * pressc / rhoi2
    
    # expanded phase, compute pressure and sound speed
    # for now, keep original Asphaug lower value limits
    p1 = ui/ezn/eta**2 + 1.
    p2 = smallaj * ui * rhoi
    p3 = smallbj * ui * rhoi / p1
    p4 = capan * xmu
    vow = 1. / eta
    p6 = betaj * (vow - 1.)
    # changed value limit from 70 to 100 to match iSALE cross check
    p6 = min(p6, 100.)
    p6 = np.exp(-p6)
    p7 = alphaj * (vow - 1.) ** 2
    p7 = min(p7, 100.)
    p7 = np.exp(-p7)
    pressv = p2 + (p3 + p4 * p6) * p7

    c1 = smallaj * ui
    c2 = p1 ** 2
    c3 = p1 - 1.
    c4 = smallbj * ui * (3. * c3 + 1.) / c2
    c5 = p6 * p7 * capan
    c6 = 2. * alphaj * (vow - 1.)
    dpdrho = c1 + p7 * c4 + p7 * p3 * rozn * c6 / rhoi2 + c5 * (1. / rozn + xmu * rozn / rhoi2 * (c6 + betaj))
    dpdi = smallaj * rhoi + p7 / c2 * smallbj * rhoi
    cvelv = dpdrho + dpdi * pressv / rhoi2
    if cvelv < 0.: 
        cvelv = 0.

    # pick up the right state
    # default is condensed region
    flag    = 1
    pri     = pressc
    vsoundi = cvelc
    if rap >= 1. and ui >= espn:
        # expanded region
        flag    = 3 
        pri     = pressv
        vsoundi = cvelv
    if rap >= 1. and ui > esn and ui < espn:
        # interpolated region
        flag    = 2
        pri     = (pressv * (ui - esn) + pressc * (espn - ui)) / difesn
        vsoundi = (cvelv * (ui - esn) + cvelc * (espn - ui)) / difesn

    # do not limit the sound speed for now; cross checking with iSALE
    #cmin = .25 * capan / rozn
    #vsoundi = max(cmin, vsoundi)
    
    pout  = pri
    if vsoundi < 0.:
        vsoundi = np.abs(pout/rin)
    csout = np.sqrt(vsoundi)
    return [pout,flag,csout]
#################
####################################################################################
### END eostable.py FILE ###
####################################################################################
