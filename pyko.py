#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyKO 1D HYDROCODE
Based on Mark Wilkins, Computer Simulation of Dynamic Phenomena, Springer-Verlag, 1999
Adapted from John Borg's Fortran KO code v11

PROGRAM DOCUMENTATION
https://impactswiki.github.io/pyko/

PROGRAM REPOSITORY
https://github.com/ImpactsWiki/pyko

AUTHOR
Created on Wed Jan 25 08:05:24 2023
@author: S. T. Stewart, U. California Davis

LICENSE
GNU General Public License v3.0

VERSIONS
    6/27/23: v0.6 first public release for beta testing
"""
############################################################## 
# IMPORT PYTHON MODULES
import numpy as np
from copy import deepcopy
from numpy import sqrt as npsqrt
from numpy import power as nppower
from numpy import where as npwhere
from numpy import absolute as npabs
from numpy import sum as npsum
import re
import sys
from os import system
from os.path import exists
import eos_table as etab # Stewart group EOS table libraries for ANEOS and Tillotson
import time
import pickle
from dataclasses import dataclass
import yaml
import pint
# see https://pint.readthedocs.io/en/develop/advanced/performance.html for ways to speed up performance with pint
#ureg = pint.UnitRegistry(cache_folder=":auto:")
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
#
# Silence NEP 18 warning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Q_([])
#
# GLOBAL VARIABLE
__version__ = 'v0.6-release-2023-06-27' 
#
############################################################## 
# KO CODE UNITS FROM WILKINS BOOK
#
# eu   = energy unit = 10^12 ergs = 100 kJ
# rho  = density = g/cm3 = 1000 kg/m3
# vr   = relative volume [to initial state] = dimless
# iev0 = internal energy per original volume = eu/g*rho0 = eu/cm3 = 100 GJ/m3
# time = time = microseconds = 1e-6 seconds
# x,y  = space coordinates = cm = 1.e-2 m
# xdot, ydot = velocity = cm/microsec = 10 km/s
# temp = temperature = K
# p    = pressure = Mbar = 1e12 dynes/cm3 = 100 GPa
# cv   = specific heat at constant volume = eu*rho0/(g K) = eu/(k cm3) = 100 GJ/K/m3
# 
############################################################## 
# Main CLASSES
#
# EOS names
#           'MGR' Mie-Gruneisen EOS
#           'IDG' ideal gas EOS
#           'TIL' Tillotson EOS
#           'SES' STD+EXT SESAME tables - needs Stewart Group eos_table.py module
#
# Constitutive model names
#           'HYDRO' Hydrodynamic; shear modulus=0
#           'VM'    Von Mises Yield
#                   more strength models coming later
#            Dynamic fracture and void closure is implemented in v0.5+
#
# Boundary Condition names
#           'FREE'  for free surface p=pvoid
#           'FIXED' for fixed surface up=0
# 
# Gravity
#           Gravitational acceleration implemented in v0.5+
#
#------------------

class RunClass:
    """ top level information for the calculation in pyKO """
    def __init__(self,fin='config.yml',fout='output.dat',ftype='YAML'):
        self.inputfiletype  = ftype # YAML OR BORG for comparison to fortran
        self.inputfilename  = fin   # name of input file, string
        self.outputfilename = fout  # name of output file, string
        self.tstop     = 0.     # time in microseconds, double
        self.ncount    = 0      # count for number of time steps, integer
        self.nmat      = 0      # number of material layers, integer
        self.ieosid    = []     # material EOS flag using tag names, list of strings length nmat
        self.ieos      = []     # list to hold EOS parameter objects, list length nmat
        self.inodes    = np.zeros(0,dtype=int)   # number of nodes for each material, intarray length nmat = 2xspatial cells
        self.ilength   = np.zeros(0)     # length of each layer, cm
        self.ixstart   = np.zeros(0)     # initial left position of each layer, cm
        self.iupstart  = np.zeros(0)     # initial particle velocity of each layer, cm/us
        self.irhostart = np.zeros(0)     # initial density of each layer, g/cm3
        self.ipstart   = np.zeros(0)     # initial pressure of each layer, Mbar
        self.itempstart   = np.zeros(0)     # initial temperature of each layer, Mbar
        self.iiev0start= np.zeros(0)     # initial internal energy per original volume of each layer, 1e12 erg/cm3
        self.istrid    = []     # material constitutive model names, list of strings length nmat
        self.istr      = []     # list to hold strength parameter objects, list length nmat
        self.ifrac     = []     # list to hold fracture parameter objects, list length nmat
        self.nbc       = []     # number of initial boundary conditions, integer
        self.ibcid     = []     # initial boundary conditions, at least [left, right]; can be more interfaces; list of strings
        self.ibc       = []     # parameters for each boundary condition
        self.grav      = GravityClass() # gravity parameters
        self.dtstart   = 0.001  # initial time step for fKO comparison
        self.dtmin     = 1.e-9  # DEBUG later add a check for too small time steps to halt the calculation
        self.pvoid     = 0.     # require for input for now
        self.dgeom     = 0.     # geometry d is DOUBLE (used in math)
        self.dflag     = ''     # geometry is a string for boolean checks
                                # 'PLA' or 'CYL' or 'SPH'
        self.avcl      = 1.     # artificial viscosity C_L
        self.avc0      = 2.     # artificial viscosity C_0
        self.step_skip = -1      # number of steps to skip for output file dumps
        self.time_skip  = 0.0   # default for fKO comp; dt between data dumps, microsec
        self.next_time_dump = self.time_skip #
        self.outputsteps = np.zeros(0,dtype=int)   #
        self.outputtimes = np.zeros(0)   #
        self.outputietot = np.zeros(0)
        self.outputketot = np.zeros(0)
        self.outputmvtot = np.zeros(0)
        # Debugging and feature flags
        self.tstepscale  = 6     # reduces the time step by this factor; 6 ok for plate impacts; 1 ok for ideal gas sod test
        self.debugflag   = False # optional debugging full grid output
        if ftype == 'BORG':
            self.binoutput   = False # flag for binary output
        else:
            self.binoutput   = True # flag for binary output
        #
    def __str__(self):
        """ Print pyKO model run parameters """
        return '\npyKO '+__version__+' run parameters\n' + \
                '   All outputs are in code units \n' + \
               f'   Input file: {self.inputfilename} \n' + \
               f'   Output file: {self.outputfilename} \n' + \
               f'   Number of materials: {self.nmat} \n' + \
               f'   Number of nodes in each material: {self.inodes} \n' + \
               f'   Length of each material: {self.ilength} \n' + \
               f'   Initial left edge of each material: {self.ixstart} \n' + \
               f'   Boundary conditions: {self.ibcid}\n' + \
               f'   Material EOS:     {self.ieosid} \n' + \
               f'   Geometry:         {self.dflag} \n' + \
               f'   Gravity:          {self.grav.gravity} \n' + \
               f'   Void pressure:    {self.pvoid} \n' + \
               f'   Time step factor: {self.tstepscale} \n' + \
               f'   Stop time:        {self.tstop}'
    def binaryoutput(self):
        """ Save the run time parameters for this simulation in a pickle file. """
        with open(self.outputfilename+'.parameters',"wb") as f: 
            #print('dumping dataout class with pickle ',self.stepn)
            pickle.dump(self,f)
    def checkinput(self):
        readinput_yaml(self,verbose=True)
        return
#
class TILClass:
    """ Tillotson EOS: Material Parameters 
        and initial state of the homogeneous material layer.
        This EOS model relies upon eos_table module. """
    def __init__(self):
        self.name   = ''  # name of the material, string 
        # Variables are stored in code units using pint during input processing
        # here showing SESAME units as a guide for the variables
        # Tillotson EOS parameters
        self.rhoref= 0. # reference density g/m3
        self.E0    = 0. # E0 constant MJ/kg
        self.EIV   = 0. # specific internal energy of incipient vaporization MJ/kg
        self.ECV   = 0. # specific internal energy of complete vaporization MJ/kg
        self.AA    = 0. # Bulk modulus K0 GPa
        self.BB    = 0. # B constant GPa
        self.a     = 0. # a constant [-]
        self.b     = 0. # b constant [-]
        self.alpha = 0. # alpha constant [-]
        self.beta  = 0. # beta constant [-]    a+b=Gruneisen gamma at rho0
        self.cs    = 0. # cm/s, sound speed(U,rho)
        self.region= 0  # region flag: 1-condensed, 2-interpolated, 3-expanded
        # Initial state parameters
        self.rho0   = 0.  # initial state density, g/cm3
        self.p0     = 0.  # initial state pressure, Mbar
        self.up0    = 0.  # initial state particle velocity, cm/us=10 km/s
        self.iev0   = 0.  # initial state internal energy per original volume, 1e12 ergs/cm3 = 100 GJ/m3
        self.v0     = 0.  # initial state specific volume, cm3/g
        self.params = np.zeros(0) # numpy array of parameters to pass to Tillotson EOS functions
    def fillparams(self):
        self.params = np.asarray([self.rhoref,self.E0,self.EIV,self.ECV, \
                           self.AA,self.BB,self.a,self.b,self.alpha,self.beta])
        return
    def __str__(self):
        """ Print the Tillotson EOS parameters """
        return f'\n{self.name} Tillotson EOS parameters: \n' + \
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
                '   initial state: \n' + \
               f'      rho0:    {self.rho0} \n' + \
               f'      p0:      {self.p0} \n' + \
               f'      iev0:    {self.iev0} \n' + \
               f'      up0:     {self.up0} \n' 
#
class MGRClass:
    """ Mie-Grueneisen EOS: Material Parameters 
        and initial state of the homogeneous material layer """
    def __init__(self):
        self.name   = ''  # name of the material, string 
        self.rhoref = 0.  # reference density for EOS g/cm3
        self.c0     = 0.  # bulk sound speed, cm/us=10 km/s, intercept for linear Us=c0+s1+up
        self.s1     = 0.  # linear Hugoniot slope, dimless 
        self.s2     = 0.  # quadratic Hugoniot coefficient, [us/cm]
        self.gamma0 = 0.  # Thermodynamic grueneisen parameter, dimless 
        self.cv     = 0.  # Heat capacity at constant volume, (Mbar cm3)/(K g)
        self.rho0   = 0.  # initial state density, g/cm3
        self.p0     = 0.  # initial state pressure, Mbar
        self.up0    = 0.  # initial state particle velocity, cm/us=10 km/s
        self.iev0   = 0.  # initial state internal energy per original volume, 1e12 ergs/cm3 = 100 GJ/m3
        self.v0     = 0.  # initial state specific volume, cm3/g
        self.coefs  = np.zeros(4) # k1,k2,k3,k4 in Wilkins Section 3.7
                                  # coefs for MG EOS from Wilkins
                                  # x = 1 - vr
                                  # P = k1*x + k2*x^2 + k3*x^3 + k4*E
                                  # alternative formula in Appendix B
                                  # eta = 1/vr = rho/rho0
                                  # P = a(eta-1)+b(eta-1)^2+c(eta-1)^3+d*eta*E
    def calccoefs(self):
        # this form in Wilkins assumes gamma_0/V_0 = constant
        # right now no treatment of porosity
        if self.rho0 != self.rhoref:
            # DEBUG: later add porosity
            sys.exit("FATAL ERROR: Mie-Grueneisen model must have rho0=rhoref.")
        self.coefs[0] = self.rhoref*self.c0*self.c0 # k1
        self.coefs[1] = self.coefs[0]*(2.*self.s1-self.gamma0/2.) # k2
        self.coefs[2] = self.coefs[0]*self.s1*(3.*self.s1-self.gamma0) # k3
        self.coefs[3] = self.gamma0 # k4     
        #                   
    def __str__(self):
        """ Print the Mie Gruneisen material parameters """
        return f'\nClass Mie Grueneisen Material Model: {self.name} \n' + \
               f'   rhoref: {self.rhoref} \n' + \
               f'   c0:     {self.c0} \n' + \
               f'   s1:     {self.s1} \n' + \
               f'   s2:     {self.s2} \n' + \
               f'   gamma0: {self.gamma0} \n' + \
               f'   cv:     {self.cv} \n' + \
                '   initial state: \n' + \
               f'      rho0:    {self.rho0} \n' + \
               f'      p0:      {self.p0} \n' + \
               f'      iev0:    {self.iev0} \n' + \
               f'      up0:     {self.up0} \n' + \
               f'   k1,k2,k3,k4= {self.coefs[0]},{self.coefs[1]},{self.coefs[2]},{self.coefs[3]} \n'
        #
class IDGClass:
    """ Ideal Gas EOS: Material Parameters 
        and initial state of the homogeneous material layer """
    def __init__(self):
        self.name   = ''  # name of the material, string 
        self.gamma0 = 0.  # Thermodynamic grueneisen parameter, dimless 
        self.cv     = 0.  # Heat capacity at constant volume, (Mbar cm3)/(K g)
        self.rho0   = 0.  # initial state density, g/cm3
        self.p0     = 0.  # initial state pressure, Mbar
        self.up0    = 0.  # initial state particle velocity, cm/us=10 km/s
        self.iev0   = 0.  # initial state internal energy per original volume, 1e12 ergs/cm3 = 100 GJ/m3
        self.v0     = 0.  # initial state specific volume, cm3/g
        self.t0     = 0.  # initial state temperature, K
        #                        
    def __str__(self):
        """ Print the Ideal Gas material parameters """
        return f'\nClass Ideal Gas: {self.name} \n' + \
               f'   gamma0: {self.gamma0} \n' + \
               f'   cv: {self.cv} \n' + \
                '   initial state: \n' + \
               f'       rho0: {self.rho0} \n' + \
               f'       p0:   {self.p0} \n' + \
               f'       iev0: {self.iev0} \n' + \
               f'       up0:  {self.up0} \n' + \
               f'       t0:  {self.t0} \n'
        #
class SESClass:
    """ ANEOS SESAME Table EOS: Material Parameters 
        and initial state of the homogeneous material layer.
        Requires eos_table.py module. """
    def __init__(self):
        self.name    = '' # name of the material, string 
        self.path    = '' # path to the ANEOS file set
        self.aneosin = 'ANEOS.INPUT' # path and name for the sesame table file, string
        self.aneosout = 'ANEOS.OUTPUT' # path and name for the sesame table file, string
        self.sesstd = 'NEW-SESAME-STD.TXT' # path and name for the sesame table file, string
        self.sesext = 'NEW-SESAME-EXT.TXT' # path and name for the sesame table file, string
        self.sesstdnt = 'NEW-SESAME-STD-NOTENSION.TXT' # path and name for the sesame table file, string
        self.gamma0 = 0.  # Thermodynamic grueneisen parameter, dimless 
        self.cv     = 0.  # Heat capacity at constant volume, (Mbar cm3)/(K g)
        self.rho0   = 0.  # initial state density, g/cm3
        self.p0     = 0.  # initial state pressure, Mbar
        self.up0    = 0.  # initial state particle velocity, cm/us=10 km/s
        self.iev0   = 0.  # initial state internal energy per original volume, 1e12 ergs/cm3 = 100 GJ/m3
        self.v0     = 0.  # initial state specific volume, cm3/g
        self.t0     = 0.  # initial state temperature, K
        self.cs0    = 0.  # initial state sound velocity
        self.EOS    = etab.extEOStable() # initialize table object (eos_table module)
        self.rhoref = 0.
        #
    def readsesext(self, verbose=False):
        """
        # READ IN NEW ANEOS MODEL and fill the extEOStable class object
        # HARD CODED UNIT CONVERSION from Stewart Group ANEOS to KO code units
        # DEPRECATED: only for testing with Borg style inputs
        # yaml configuration files includes units conversions; use readtable function instead
        """
        #------------------------------------------------------------------
        self.EOS.loadextsesame(self.path+self.sesext) #'NEW-SESAME-EXT.TXT') # LOAD THE EXTENDED 301 SESAME FILE GENERATED BY STSM VERSION OF ANEOS
        self.EOS.loadstdsesame(self.path+self.sesstd) # 'NEW-SESAME-STD.TXT') # LOAD THE STANDARD 301 SESAME FILE GENERATED BY STSM VERSION OF ANEOS
        #print(NewEOS.units) # these are the default units for SESAME rho-T tables
        #'Units: g/c3, K, GPa, MJ/kg, MJ/kg, MJ/K/kg, cm/s, MJ/K/kg, KPA flag. 2D arrays are (NT,ND).'
        # CONVERT THE TABLE UNITS TO WILKINS CODE UNITS
        # T in K; rho in g/cm3
        if verbose: print('WARNING HARD CODED SESAME UNITS CONVERSION! \n' + \
                            'U by 1/100 MJ/kg -> 100 kJ/g; not normalizing to rho0')
        self.EOS.P = self.EOS.P/100. # GPa to Mbar
        self.EOS.U = self.EOS.U/100. # MJ/kg *10/1000 -> 100 kJ/g
        self.EOS.S = self.EOS.S/100. # MJ/K/kg *10/1000 -> 100 kJ/g
        #
    def readtable(self, config, verbose=False):
        """
        # Read in SESAME style EOS tables and fill the extEOStable class object
        # Unit conversion between table units and code units requires units section
        # in yaml configuration file
        #
        # NOTE: The specific energy in the table remains in specific energy units.
        # The normalization to initial volume is done in the interpolation function
        # used in the energy conservation loop because this requires the initial state
        # information for the problem.
        """
        #------------------------------------------------------------------
        # LOAD THE STANDARD 301 SESAME FILE GENERATED BY STSM ANEOS SCRIPTS NEW-SESAME-STD.TXT
        self.EOS.loadstdsesame(self.path+self.sesstd)
        # LOAD THE EXTENDED 301 SESAME FILE GENERATED BY STSM ANEOS SCRIPTS NEW-SESAME-EXT.TXT
        self.EOS.loadextsesame(self.path+self.sesext)
        # print(NewEOS.units) # these are the units for Stewart Group SESAME T-rho tables
        #'Units: g/cm3, K, GPa, MJ/kg, MJ/kg, MJ/K/kg, cm/s, MJ/K/kg, KPA flag. 2D arrays are (NT,ND).'
        # CONVERT THE TABLE UNITS TO CODE UNITS using pint and yaml configuration file information
        # T in K; rho in g/cm3
        # print('Converting EOS Table to code units.')
        rho = Q_(self.EOS.rho,config['tableunits']['density'])
        self.EOS.rho = rho.to(config['codeunits']['density']).magnitude
        temp = Q_(self.EOS.T,config['tableunits']['temperature'])
        self.EOS.T = temp.to(config['codeunits']['temperature']).magnitude
        P = Q_(self.EOS.P,config['tableunits']['pressure'])
        self.EOS.P = P.to(config['codeunits']['pressure']).magnitude
        U = Q_(self.EOS.U,config['tableunits']['sp_energy'])
        self.EOS.U = U.to(config['codeunits']['sp_energy']).magnitude
        S = Q_(self.EOS.S,config['tableunits']['sp_entropy'])
        self.EOS.S = S.to(config['codeunits']['sp_entropy']).magnitude
        cs = Q_(self.EOS.cs,config['tableunits']['sound_speed'])
        self.EOS.cs = cs.to(config['codeunits']['velocity']).magnitude
        cv = Q_(self.EOS.cv,config['tableunits']['sp_heat_cap'])
        self.EOS.cv = cv.to(config['codeunits']['sp_heat_cap']).magnitude
        A = Q_(self.EOS.A,config['tableunits']['sp_energy'])
        self.EOS.A = A.to(config['codeunits']['sp_energy']).magnitude
        # if no sound speeds tabulated, calculate a bulk sound speed
        if np.max(self.EOS.cs) == 0:
            print('NO SOUND SPEEDS. CALCULATING BULK SOUND SPEEDS.')
            for irho in range(1,self.EOS.ND-1):
                for itemp in range(1,self.EOS.NT-1):
                    btlocal = self.EOS.rho[irho]* \
                        (self.EOS.P[itemp,irho+1]-self.EOS.P[itemp,irho-1])/ \
                        (self.EOS.rho[irho+1]-self.EOS.rho[irho-1])
                    cslocal = np.sqrt(np.abs(btlocal)/self.EOS.rho[irho])
                    #print(self.EOS.rho[irho],self.EOS.T[itemp],btlocal,cslocal,self.EOS.CS0REF)
                    # in code units
                    if cslocal == 0.0:
                        cslocal = self.EOS.CS0REF
                    self.EOS.cs[itemp,irho] = cslocal
        return
        #
    def findpve(self,p,rho):
        rindex = npwhere(self.EOS.rho >= rho)[0][0]
        tindex = npwhere(self.EOS.P[:,rindex] >= p)[0][0]
        print('findpve p, rho -> rho, p, U', p, rho, self.EOS.rho[rindex],self.EOS.P[tindex,rindex],self.EOS.U[tindex,rindex])
        #
    def sesplog(self,rho,u,p,t,eps):
        """ Testing linear interpolation of log values. Only works for ideal gas.
            Cannot use with tables with a tension region.
            Input rho, u, p, t must be in the same units as the table."""
        pnew = np.copy(p)
        tnew = np.copy(t)
        # if eps = 0., then no change in state
        ichange = npwhere(eps != 0.)[0]
        #print('ichagne = ',ichange)
        for ij in ichange:
            # bilinear interpolation
            # https://x-engineer.org/bilinear-interpolation/
            rindex = npwhere(self.EOS.rho > rho[ij])[0][0]
            tindex = npwhere(self.EOS.U[:,rindex] > u[ij])[0][0]
            P11 = self.EOS.P[tindex-1,rindex-1]
            P21 = self.EOS.P[tindex,rindex-1]
            P12 = self.EOS.P[tindex-1,rindex]
            P22 = self.EOS.P[tindex,rindex]
            U11 = self.EOS.U[tindex-1,rindex-1]
            U21 = self.EOS.U[tindex,rindex-1]
            U12 = self.EOS.U[tindex-1,rindex]
            U22 = self.EOS.U[tindex,rindex]
            r1 = np.log10(self.EOS.rho[rindex-1]) # try linearly interpreting on the log of the index
            r2 = np.log10(self.EOS.rho[rindex])
            logrij = np.log10(rho[ij])
            loguij = np.log10(u[ij])
            t1 = np.log10(self.EOS.T[tindex-1])
            t2 = np.log10(self.EOS.T[tindex])
            u1 = np.log10(U11) * (r2-logrij)/(r2-r1) + np.log10(U12) * (logrij-r1)/(r2-r1)
            u2 = np.log10(U21) * (r2-logrij)/(r2-r1) + np.log10(U22) * (logrij-r1)/(r2-r1)
            p1 = np.log10(P11) * (r2-logrij)/(r2-r1) + np.log10(P12) * (logrij-r1)/(r2-r1)
            p2 = np.log10(P21) * (r2-logrij)/(r2-r1) + np.log10(P22) * (logrij-r1)/(r2-r1)
            pij = p1 * (u2-loguij)/(u2-u1) + p2 * (loguij-u1)/(u2-u1)
            pnew[ij] = nppower(10.,pij)
            tij = t1 * (u2-loguij)/(u2-u1) + t2 * (loguij-u1)/(u2-u1)
            tnew[ij] = nppower(10.,tij)
            if True:
                print('bilinear = ',rindex,tindex,'\n',\
                      nppower(10.,r1),rho[ij],nppower(10.,r2),'\n',\
                      nppower(10.,u1),u[ij],nppower(10.,u2),'\n', \
                      nppower(10.,p1),pnew[ij],nppower(10.,p2),'\n',\
                      nppower(10.,t1),tnew[ij],nppower(10.,t2))
                print()
                print('')
        return pnew,tnew
    #
    def sesp(self,rho,u,p,t,s,eps):
        """ Linear interpolation to extract table P and T from rho and U. 
            Input rho, u, p, t must be in the same units as the table.
            Returns pnew, tnew arrays. """
        pnew = np.copy(p) # This routine only calculates cells with nonzero strain.
        tnew = np.copy(t) # copy current P & T values to populate cells with no change.
        snew = np.copy(s) # This routine only calculates cells with nonzero strain.
        # if eps = 0., then no change in state
        ichange = npwhere(eps != 0.)[0]
        #print('ichange = ',ichange)
        for ij in ichange:
            # bilinear interpolation
            # https://x-engineer.org/bilinear-interpolation/
            rindex = npwhere(self.EOS.rho > rho[ij])[0][0]
            tindex = npwhere(self.EOS.U[:,rindex] > u[ij])[0][0]
            P11 = self.EOS.P[tindex-1,rindex-1]
            P21 = self.EOS.P[tindex,rindex-1]
            P12 = self.EOS.P[tindex-1,rindex]
            P22 = self.EOS.P[tindex,rindex]
            U11 = self.EOS.U[tindex-1,rindex-1]
            U21 = self.EOS.U[tindex,rindex-1]
            U12 = self.EOS.U[tindex-1,rindex]
            U22 = self.EOS.U[tindex,rindex]
            S11 = self.EOS.S[tindex-1,rindex-1]
            S21 = self.EOS.S[tindex,rindex-1]
            S12 = self.EOS.S[tindex-1,rindex]
            S22 = self.EOS.S[tindex,rindex]
            r1 = self.EOS.rho[rindex-1]
            r2 = self.EOS.rho[rindex]
            t1 = self.EOS.T[tindex-1]
            t2 = self.EOS.T[tindex]
            u1 = U11 * (r2-rho[ij])/(r2-r1) + U12 * (rho[ij]-r1)/(r2-r1)
            u2 = U21 * (r2-rho[ij])/(r2-r1) + U22 * (rho[ij]-r1)/(r2-r1)
            p1 = P11 * (r2-rho[ij])/(r2-r1) + P12 * (rho[ij]-r1)/(r2-r1)
            p2 = P21 * (r2-rho[ij])/(r2-r1) + P22 * (rho[ij]-r1)/(r2-r1)
            pij = p1 * (u2-u[ij])/(u2-u1) + p2 * (u[ij]-u1)/(u2-u1)
            pnew[ij] = pij
            tij = t1 * (u2-u[ij])/(u2-u1) + t2 * (u[ij]-u1)/(u2-u1)
            tnew[ij] = tij
            s1 = S11 * (r2-rho[ij])/(r2-r1) + S12 * (rho[ij]-r1)/(r2-r1)
            s2 = S21 * (r2-rho[ij])/(r2-r1) + S22 * (rho[ij]-r1)/(r2-r1)
            sij = s1 * (u2-u[ij])/(u2-u1) + s2 * (u[ij]-u1)/(u2-u1)
            snew[ij] = sij
            if False:
                print('bilinear = ',rindex,tindex,'\n',\
                      r1,rho[ij],r2,'\n',\
                      u1,u[ij],u2,'\n', \
                      p1,pnew[ij],p2,'\n',\
                      t1,tnew[ij],t2,'\n',\
                      s1,snew[ij],s2)
                print('')
        return pnew,tnew,snew
    def sescs(self,rho,u):
        """ Linear interpolation to extract table cs from rho and U. 
            Input rho, u must be in the same units as the table.
            Returns cs array. """
            # bilinear interpolation
            # https://x-engineer.org/bilinear-interpolation/
        csarr = np.zeros(len(rho))
        for ij in np.arange(len(rho)):
            rindex = npwhere(self.EOS.rho > rho[ij])[0][0]
            tindex = npwhere(self.EOS.U[:,rindex] > u[ij])[0][0]
            C11 = self.EOS.cs[tindex-1,rindex-1]
            C21 = self.EOS.cs[tindex,rindex-1]
            C12 = self.EOS.cs[tindex-1,rindex]
            C22 = self.EOS.cs[tindex,rindex]
            U11 = self.EOS.U[tindex-1,rindex-1]
            U21 = self.EOS.U[tindex,rindex-1]
            U12 = self.EOS.U[tindex-1,rindex]
            U22 = self.EOS.U[tindex,rindex]
            r1 = self.EOS.rho[rindex-1]
            r2 = self.EOS.rho[rindex]
            u1 = U11 * (r2-rho[ij])/(r2-r1) + U12 * (rho[ij]-r1)/(r2-r1)
            u2 = U21 * (r2-rho[ij])/(r2-r1) + U22 * (rho[ij]-r1)/(r2-r1)
            c1 = C11 * (r2-rho[ij])/(r2-r1) + C12 * (rho[ij]-r1)/(r2-r1)
            c2 = C21 * (r2-rho[ij])/(r2-r1) + C22 * (rho[ij]-r1)/(r2-r1)
            cij = c1 * (u2-u[ij])/(u2-u1) + c2 * (u[ij]-u1)/(u2-u1)
            csarr[ij] = cij
        return csarr
    def sesphase(self,rho,u):
        """ Linear interpolation to extract table phase from rho and U. 
            Input rho, u must be in the same units as the table.
            Returns phase array. """
            # bilinear interpolation
            # https://x-engineer.org/bilinear-interpolation/
        phasearr = np.zeros(len(rho))
        for ij in np.arange(len(rho)):
            rindex = npwhere(self.EOS.rho > rho[ij])[0][0]
            tindex = npwhere(self.EOS.U[:,rindex] > u[ij])[0][0]
            K11 = self.EOS.KPA[tindex-1,rindex-1]
            K21 = self.EOS.KPA[tindex,rindex-1]
            K12 = self.EOS.KPA[tindex-1,rindex]
            K22 = self.EOS.KPA[tindex,rindex]
            U11 = self.EOS.U[tindex-1,rindex-1]
            U21 = self.EOS.U[tindex,rindex-1]
            U12 = self.EOS.U[tindex-1,rindex]
            U22 = self.EOS.U[tindex,rindex]
            r1 = self.EOS.rho[rindex-1]
            r2 = self.EOS.rho[rindex]
            u1 = U11 * (r2-rho[ij])/(r2-r1) + U12 * (rho[ij]-r1)/(r2-r1)
            u2 = U21 * (r2-rho[ij])/(r2-r1) + U22 * (rho[ij]-r1)/(r2-r1)
            k1 = K11 * (r2-rho[ij])/(r2-r1) + K12 * (rho[ij]-r1)/(r2-r1)
            k2 = K21 * (r2-rho[ij])/(r2-r1) + K22 * (rho[ij]-r1)/(r2-r1)
            kij = k1 * (u2-u[ij])/(u2-u1) + k2 * (u[ij]-u1)/(u2-u1)
            phasearr[ij] = kij
        #print('in sesphase rho, phasearr = ',rho,phasearr)
        return phasearr
    def oneptc(self,rho,u):
        """ Linear interpolation to extract table one point P, T, cs from rho and U. 
            Input rho, u must be in the same units as the table.
            Returns P, T, cs """
            # bilinear interpolation
            # https://x-engineer.org/bilinear-interpolation/
        rindex = npwhere(self.EOS.rho > rho)[0][0]
        tindex = npwhere(self.EOS.U[:,rindex] > u)[0][0]
        P11 = self.EOS.P[tindex-1,rindex-1]
        P21 = self.EOS.P[tindex,rindex-1]
        P12 = self.EOS.P[tindex-1,rindex]
        P22 = self.EOS.P[tindex,rindex]
        C11 = self.EOS.cs[tindex-1,rindex-1]
        C21 = self.EOS.cs[tindex,rindex-1]
        C12 = self.EOS.cs[tindex-1,rindex]
        C22 = self.EOS.cs[tindex,rindex]
        U11 = self.EOS.U[tindex-1,rindex-1]
        U21 = self.EOS.U[tindex,rindex-1]
        U12 = self.EOS.U[tindex-1,rindex]
        U22 = self.EOS.U[tindex,rindex]
        r1 = self.EOS.rho[rindex-1]
        r2 = self.EOS.rho[rindex]
        t1 = self.EOS.T[tindex-1]
        t2 = self.EOS.T[tindex]
        u1 = U11 * (r2-rho)/(r2-r1) + U12 * (rho-r1)/(r2-r1)
        u2 = U21 * (r2-rho)/(r2-r1) + U22 * (rho-r1)/(r2-r1)
        p1 = P11 * (r2-rho)/(r2-r1) + P12 * (rho-r1)/(r2-r1)
        p2 = P21 * (r2-rho)/(r2-r1) + P22 * (rho-r1)/(r2-r1)
        p = p1 * (u2-u)/(u2-u1) + p2 * (u-u1)/(u2-u1)
        t = t1 * (u2-u)/(u2-u1) + t2 * (u-u1)/(u2-u1)
        c1 = C11 * (r2-rho)/(r2-r1) + C12 * (rho-r1)/(r2-r1)
        c2 = C21 * (r2-rho)/(r2-r1) + C22 * (rho-r1)/(r2-r1)
        c = c1 * (u2-u)/(u2-u1) + c2 * (u-u1)/(u2-u1)
        return p,t,c
        #
    def onepuc(self,rho,t):
        """ Linear interpolation to extract table one point P, U, cs from rho and T. 
            Input rho, t must be in the same units as the table, 
            which is converted to code units.
            Returns P, U, cs """
            # bilinear interpolation
            # https://x-engineer.org/bilinear-interpolation/
        rindex = npwhere(self.EOS.rho >= rho)[0][0]
        tindex = npwhere(self.EOS.T >= t)[0][0]
        P11 = self.EOS.P[tindex-1,rindex-1]
        P21 = self.EOS.P[tindex,rindex-1]
        P12 = self.EOS.P[tindex-1,rindex]
        P22 = self.EOS.P[tindex,rindex]
        C11 = self.EOS.cs[tindex-1,rindex-1]
        C21 = self.EOS.cs[tindex,rindex-1]
        C12 = self.EOS.cs[tindex-1,rindex]
        C22 = self.EOS.cs[tindex,rindex]
        U11 = self.EOS.U[tindex-1,rindex-1]
        U21 = self.EOS.U[tindex,rindex-1]
        U12 = self.EOS.U[tindex-1,rindex]
        U22 = self.EOS.U[tindex,rindex]
        r1 = self.EOS.rho[rindex-1]
        r2 = self.EOS.rho[rindex]
        t1 = self.EOS.T[tindex-1]
        t2 = self.EOS.T[tindex]
        u1 = U11 * (r2-rho)/(r2-r1) + U12 * (rho-r1)/(r2-r1)
        u2 = U21 * (r2-rho)/(r2-r1) + U22 * (rho-r1)/(r2-r1)
        p1 = P11 * (r2-rho)/(r2-r1) + P12 * (rho-r1)/(r2-r1)
        p2 = P21 * (r2-rho)/(r2-r1) + P22 * (rho-r1)/(r2-r1)
        p = p1 * (t2-t)/(t2-t1) + p2 * (t-t1)/(t2-t1)
        u = u1 * (t2-t)/(t2-t1) + u2 * (t-t1)/(t2-t1)
        c1 = C11 * (r2-rho)/(r2-r1) + C12 * (rho-r1)/(r2-r1)
        c2 = C21 * (r2-rho)/(r2-r1) + C22 * (rho-r1)/(r2-r1)
        c = c1 * (t2-t)/(t2-t1) + c2 * (t-t1)/(t2-t1)
        return p,u,c
        #
    def oneruc(self,p,t):
        """ Linear interpolation to extract table one point rho, U, cs from P and T. 
            Altenatively, use a notension table to initialize with gravity.
            Input P and T must be in the same units as the table.
            Returns rho, U, cs """
            # bilinear interpolation
            # https://x-engineer.org/bilinear-interpolation/
        tindex = npwhere(self.EOS.T >= t)[0][0]
        rtmp = npwhere(self.EOS.P[tindex-1,:] >= p)[0]
        rindex = rtmp[0]
        if self.EOS.T[tindex] == t:
            r1 = self.EOS.rho[rindex-1]
            r2 = self.EOS.rho[rindex]
            u1 = self.EOS.U[tindex,rindex-1]
            u2 = self.EOS.U[tindex,rindex]
            p1 = self.EOS.P[tindex,rindex-1]
            p2 = self.EOS.P[tindex,rindex]
            c1 = self.EOS.cs[tindex,rindex-1]
            c2 = self.EOS.cs[tindex,rindex]
            u = u1 * (p2-p)/(p2-p1) + u2 * (p-p1)/(p2-p1)
            rho = r1 * (p2-p)/(p2-p1) + r2 * (p-p1)/(p2-p1)
            c = c1 * (p2-p)/(p2-p1) + c2 * (p-p1)/(p2-p1)
            #print('ruc = ',r1,r2,p1,p2,p,rho,u,c,tindex,rindex)
            return rho,u,c
        P11 = self.EOS.P[tindex-1,rindex-1]
        P21 = self.EOS.P[tindex,rindex-1]
        P12 = self.EOS.P[tindex-1,rindex]
        P22 = self.EOS.P[tindex,rindex]
        C11 = self.EOS.cs[tindex-1,rindex-1]
        C21 = self.EOS.cs[tindex,rindex-1]
        C12 = self.EOS.cs[tindex-1,rindex]
        C22 = self.EOS.cs[tindex,rindex]
        U11 = self.EOS.U[tindex-1,rindex-1]
        U21 = self.EOS.U[tindex,rindex-1]
        U12 = self.EOS.U[tindex-1,rindex]
        U22 = self.EOS.U[tindex,rindex]
        r1 = self.EOS.rho[rindex-1]
        r2 = self.EOS.rho[rindex]
        t1 = self.EOS.T[tindex-1]
        t2 = self.EOS.T[tindex]
        u1 = U11 * (t2-t)/(t2-t1) + U21 * (t-t1)/(t2-t1)
        u2 = U12 * (t2-t)/(t2-t1) + U22 * (t-t1)/(t2-t1)
        p1 = P11 * (t2-t)/(t2-t1) + P21 * (t-t1)/(t2-t1)
        p2 = P12 * (t2-t)/(t2-t1) + P22 * (t-t1)/(t2-t1)
        c1 = C11 * (t2-t)/(t2-t1) + C21 * (t-t1)/(t2-t1)
        c2 = C12 * (t2-t)/(t2-t1) + C22 * (t-t1)/(t2-t1)
        u   = u1 * (p2-p)/(p2-p1) + u2 * (p-p1)/(p2-p1)
        rho = r1 * (p2-p)/(p2-p1) + r2 * (p-p1)/(p2-p1)
        c   = c1 * (p2-p)/(p2-p1) + c2 * (p-p1)/(p2-p1)
        return rho,u,c
        #
    def __str__(self):
        """ Print overview of SESAME Table EOS material parameters """
        return f'\nClass SESAME: {self.name} \n' + \
               f'   eos_table module version {etab.__version__} \n' + \
               f'   table path: {self.path} \n' + \
               f'   file names: {self.sesstd} {self.sesext} \n' + \
                '   initial state: \n' + \
               f'      rho0: {self.rho0} \n' + \
               f'      p0:   {self.p0} \n' + \
               f'      iev0: {self.iev0} \n' + \
               f'      up0:  {self.up0} \n' + \
               f'      t0:   {self.t0} \n' + \
               f'      CS0REF: {self.EOS.CS0REF} \n' + \
               f'      cs0:  {self.cs0} '
      #
class HydroClass:
    """ Hydrodynamic material class """
    def __init__(self):
        self.name  = '' # name of the material, string 
    def __str__(self):
        """ Print HYDRO material descriptor """
        return f'\n{self.name}: Hydrodynamic material'
    #
class VonMisesClass:
    """ Strength class for Von Mises Yield """
    def __init__(self):
        self.name  = '' # name of the material, string 
        self.gmod  = 0. # shear modulus, Mbar
        self.ys    = 0. # von Mises yield stress, Mbar
    def __str__(self):
        """ Print Von Mises material parameters """
        return f'\n{self.name} Von Mises parameters: \n' + \
            f'   Shear modulus [Mbar]: {self.gmod} \n' + \
            f'   Yield stress [Mbar]: {self.ys}' 
    #
class FractureClass:
    """ Fracture strength class. Fracture requires -P > pfrac and rho < rhomin*rhoref. """
    def __init__(self):
        self.name     = ''  # option string for fracture parameters
        self.pfrac    = 0.  # fracture pressure; fracture requires -P > pfrac
        self.nrhomin   = 0.8  # normalized maximum distension factor; fracture requires rho < (nrhomin*rhoref); usually 0.8-0.95
    def __str__(self):
        """ Fracture parameters """
        return f'\n{self.name} Fracture parameters: \n' + \
            f'   Fracture pressure: {self.pfrac} \n' + \
            f'   Fracture maximum distension (rhomin/rhoref): {self.nrhomin}'
    #        
class BCClass:
    """ Boundary Conditions """
    def __init__(self):
        self.name    = '' # name of the boundary condition, string 
        self.pbc     = 0. # pressure in ghost cell, Mbar
        self.upbc    = 0. # particle velocity in ghost cell, cm/us
        self.rrefbc  = 0. # density in ghost cell, g/cm3
        self.iebc    = 0. # internal energy per original volume in ghost cell, 1e12 ergs/cm3
        self.vrbc    = 0. # relative volume in ghost cell, dimless    
    def __str__(self):
        """ Print boundary condition parameters """
        return f'\nClass Boundary Condition: {self.name} \n' + \
            f'   Pressure [Mbar]: {self.pbc} ' # \n' + \
            #f'   Particle velocity [cm/us]: {self.upbc} \n' + \
            #f'   Reference density [g/cm3]: {self.rrefbc} \n' + \
            #f'   Internal energy [Mbar]: {self.iebc} \n' + \
            #f'   Relative volume [-]: {self.vrbc}'
#
class GravityClass:
    """ Gravity parameters """
    def __init__(self):
        self.gravity = 0. # gravitational acceleration
        self.matflag = np.zeros(0,dtype='int')  # Flag for initialization with gravity
        self.refpos  = np.zeros(0)  # Reference position for initialization with gravity
        self.refpres = np.zeros(0) # reference pressure at reference position
    def __str__(self):
        """ Print boundary condition parameters """
        return f'\nClass Gravity: \n' + \
            f'   Gravitational acc [cm/us/us]: {self.gravity} \n' + \
            f'   Reference position [cm/us]: {self.refpos} ' # \n' + \
            #f'   Reference density [g/cm3]: {self.rrefbc} \n' + \
            #f'   Internal energy [Mbar]: {self.iebc} \n' + \
            #f'   Relative volume [-]: {self.vrbc}'
#
@dataclass
class OutputClass:
    """ Problem domain data structure for output dumps """
    def __init__(self):
        """ Initialize the main arrays for the problem domain """
        # unknown number of cells initially. Initialize data types with zero length.
        # variables in space at a particular snapshot in time
        self.stepn    = np.zeros(0,dtype='int') # counter for the number of time steps; int
        self.time     = np.zeros(0) # same time for each cell; used to make pandas conversion easier
        self.mat      = np.zeros(0,dtype='int') # material id number
        self.pos      = np.zeros(0) # position in 1D geometry (x,r_cyl,r_sph)
        self.rho0     = np.zeros(0) # initial density rho0 g/cm3; DEBUG later can free up this memory by referring to imat definition
        self.rho      = np.zeros(0) # local density g/cm3
        self.up       = np.zeros(0) # particle velocity
        #self.iev0     = np.zeros(0) # internal energy per initial volume 1e12 erg/cm3
        self.ie       = np.zeros(0) # specific internal energy 
        self.pres     = np.zeros(0) # pressure Mbar
        self.mass     = np.zeros(0) # mass in node g
        self.temp     = np.zeros(0) # temperature K
        self.sigmar   = np.zeros(0) # total radial stress sigma_r, only evaluated at n=1 current time
        self.sigmao   = np.zeros(0) # total tangential stress sigma_theta, only evaluated at n=1 current time
        self.etot     = np.zeros(0) # total IE + total KE
        self.j        = np.zeros(0) # node index for comparisons to fKO
        #self.vr       = np.zeros(0) # relative volume to initial state=rho0/rho=1/eta in Wilkins [dimless]
        #self.phi      = np.zeros(0) # rho0*(pos_j+1-pos_j)/vr_j+1/2 = (pos_j+1-pos_j)/rho_j+1/2=phi in Wilkins
        #self.q        = np.zeros(0) # artificial viscosity Mbar
        #self.eps1     = np.zeros(0) # velocity strain dup/dpos, per microsec
        #self.eps2     = np.zeros(0) # velocity strain up/pos, per microsec
        # these variables are evaluated at n=3; full time step ahead
        #self.beta     = np.zeros(0) # (sigmar-sigmao)/(0.5 dpos)*rho = beta in Wilkins
        #self.s1       = np.zeros(0) # s1 stress deviator deriv = 2G(deps1-dvr/vr/3)
        #self.s2       = np.zeros(0) # s2 stress deviator deriv
        #self.s3       = np.zeros(0) # s3 stress deviator deriv; not needed in 1D; kept for equation clarity
        self.entropy  = np.zeros(0) # entropy [code units is eu/K/g]
        self.dtminj   = np.zeros(0) # calculated minimum time step
        # variables in time only - n
        #self.ibc      = np.zeros(0,dtype='int') # boundary condition id number, integer
        #self.yld      = np.zeros(0) # yield stress
        self.pfrac    = np.zeros(0) # fracture stress
        #self.deltaz   = np.zeros(0) # distortion energy Mbar
        self.alocal   = np.zeros(0) # sound speed; intermediate variable for AV
        #
    def printunits(self):
        print('pyKO output data class keys and units: \n' +\
              f'  stepn:  {self.stepn.units} \n' + \
              f'  time:   {self.time.units} \n' + \
              f'  mat:    {self.mat.units} \n' + \
              f'  pos:    {self.pos.units} \n' + \
              f'  rho0:   {self.rho0.units} \n' + \
              f'  rho:    {self.rho.units} \n' + \
              f'  up:     {self.up.units} \n' + \
              f'  ie:     {self.ie.units} \n' + \
              f'  pres:   {self.pres.units} \n' + \
              f'  mass:   {self.mass.units} \n' + \
              f'  temp:   {self.temp.units} \n' + \
              f'  sigmar: {self.sigmar.units} \n' + \
              f'  sigmao: {self.sigmao.units} \n' + \
              f'  etot:   {self.etot.units} \n' + \
               '  j:      (dimless)' + \
              f'  entropy:{self.entropy.units} \n' + \
              f'  dtminj: {self.dtminj.units} \n' + \
              f'  pfrac:  {self.pfrac.units} \n' + \
              f'  alocal: {self.alocal.units} \n'
              )
        return
    #
@dataclass
class DebugClass:
    """ Detailed domain data structure for debugging output dumps 
    """
    def __init__(self):
        """ Initialize the main arrays for the problem domain """
        # unknown number of cells initially. Initialize data types with zero length.
        # variables in space at a particular snapshot in time        
        self.stepn    = np.zeros(0,dtype='int') # counter for the number of time steps; int
        self.time     = np.zeros(0) # same time for each cell; used to make pandas conversion easier
        self.mat      = np.zeros(0,dtype='int') # material id number
        self.pos      = np.zeros(0) # position in 1D geometry (x,r_cyl,r_sph)
        self.rho0     = np.zeros(0) # initial density rho0 g/cm3; DEBUG later can free up this memory by referring to imat definition
        self.rho      = np.zeros(0) # local density g/cm3
        self.up       = np.zeros(0) # particle velocity
        #self.iev0     = np.zeros(0) # internal energy per initial volume 1e12 erg/cm3
        self.ie       = np.zeros(0) # specific internal energy 
        self.pres     = np.zeros(0) # pressure Mbar
        self.mass     = np.zeros(0) # mass in node g
        self.temp     = np.zeros(0) # temperature K
        self.sigmar   = np.zeros(0) # total radial stress sigma_r, only evaluated at n=1 current time
        self.sigmao   = np.zeros(0) # total tangential stress sigma_theta, only evaluated at n=1 current time
        self.etot     = np.zeros(0) # total IE + total KE
        self.j        = np.zeros(0) # node index for comparisons to fKO
        self.vr       = np.zeros(0) # relative volume to initial state=rho0/rho=1/eta in Wilkins [dimless]
        self.phi      = np.zeros(0) # rho0*(pos_j+1-pos_j)/vr_j+1/2 = (pos_j+1-pos_j)/rho_j+1/2=phi in Wilkins
        self.q        = np.zeros(0) # artificial viscosity Mbar
        self.eps1     = np.zeros(0) # velocity strain dup/dpos, per microsec
        self.eps2     = np.zeros(0) # velocity strain up/pos, per microsec
        # these variables are evaluated at n=3; full time step ahead
        self.beta     = np.zeros(0) # (sigmar-sigmao)/(0.5 dpos)*rho = beta in Wilkins
        self.s1       = np.zeros(0) # s1 stress deviator deriv = 2G(deps1-dvr/vr/3)
        self.s2       = np.zeros(0) # s2 stress deviator deriv
        self.s3       = np.zeros(0) # s3 stress deviator deriv; not needed in 1D; kept for equation clarity
        self.entropy  = np.zeros(0) # entropy [DEBUG check eu/cm3 units]
        # variables in time only - n
        self.ibc      = np.zeros(0,dtype='int') # boundary condition id number, integer
        self.yld      = np.zeros(0) # yield stress
        self.pfrac    = np.zeros(0) # fracture stress
        self.deltaz   = np.zeros(0) # distortion energy Mbar
        self.dtminj   = np.zeros(0) # 
        self.alocal   = np.zeros(0) # sound speed; intermediate variable for AV
        self.phase    = np.zeros(0) # phase id in sesame tables
        self.bcpres    = np.zeros(0) # phase id in sesame tables
        #
class DomainClass:
    """ Problem domain data structure """
    def __init__(self):
        """ Initialize the main arrays for the problem domain """
        # jtot nodes x ntot time steps
        # jtot is the initial number of nodes [can change with fracture and healing]
        # ntot is the number of time steps kept in memory
        #self.jtot = npsum(run.inodes)+1 # DEBUG need to figure out how to handle more boundaries
        self.jtot = npsum(10)+1 # DEBUG need to figure out how to handle more boundaries
                                         # initial guess is the minimum default; all materials touching; add a node for the end
                                         # handle on the fly with functions to insert or remove a node
        self.ntot  = 4 #  t(0) and t(1) current, t(2) and t(3) new; keeps memory very small 
        self.ngaps = 0 # start with the assumption that there are no interior open spaces
                       # if true, then initial jtot is correct; if not, need to extend the arrays
        # variables in time and space (n,j)
        n = self.ntot                 # running time step index
        j = self.jtot                 # spatial domain index
        # these variables are evaluated at n=1; current time step
        # these variables are evaluated at n=2; half time step ahead
        self.up       = np.zeros((n,j)) # particle velocity
        self.pos      = np.zeros((n,j)) # position in 1D geometry (x,r_cyl,r_sph)
        self.vr       = np.zeros((n,j)) # relative volume to initial state=rho0/rho=1/eta in Wilkins [dimless]
        self.phi      = np.zeros((n,j)) # rho0*(pos_j+1-pos_j)/vr_j+1/2 = (pos_j+1-pos_j)/rho_j+1/2=phi in Wilkins
        self.q        = np.zeros((n,j)) # artificial viscosity Mbar
        self.eps1     = np.zeros((n,j)) # velocity strain dup/dpos, per microsec
        self.eps2     = np.zeros((n,j)) # velocity strain up/pos, per microsec
        # these variables are evaluated at n=3; full time step ahead
        self.beta     = np.zeros((n,j)) # (sigmar-sigmao)/(0.5 dpos)*rho = beta in Wilkins
        self.iev0     = np.zeros((n,j)) # internal energy per initial volume 1e12 erg/cm3
        self.pres     = np.zeros((n,j)) # pressure Mbar
        self.s1       = np.zeros((n,j)) # s1 stress deviator deriv = 2G(deps1-dvr/vr/3)
        self.s2       = np.zeros((n,j)) # s2 stress deviator deriv
        self.s3       = np.zeros((n,j)) # s3 stress deviator deriv; not needed in 1D; kept for equation clarity
        self.temp     = np.zeros((n,j)) # temperature K
        self.entropy  = np.zeros((n,j)) # entropy [DEBUG check eu/cm3 units]
        # variables in time only - n
        self.deltat   = 0. #run.dtstart      # current time step; while in time loop
        self.deltat_next = 0. #run.dtstart   # next time step; calculated at end of time loop
        self.deltat_past = 0.                # if new contact, keep track of previous time step and reduce by XX factor (nominal is 2)
        self.time     = np.zeros(4)      # current 4 half time steps
        self.stepn    = int(0)                # counter for the number of time steps; int
        self.qtotal   = 0.
        self.mvtotal  = 0.
        self.ketotal  = 0.
        self.ietotal  = 0.
        self.etotal   = 0.
        self.mtotal   = 0.
        # variables in space only - j
        self.matid    = np.zeros(j) # material id number
        self.ibc      = np.zeros(j,dtype='int') # boundary condition id number, integer
        self.yld      = np.zeros(j) # yield stress
        self.pfrac    = np.zeros(j) # fracture stress
        self.rho0     = np.zeros(j) # initial density rho0 g/cm3; DEBUG later can free up this memory by referring to imat definition
        self.rho      = np.zeros(j) # local density g/cm3
        self.mass     = np.zeros(j) # mass in node g
        self.deltaz   = np.zeros(j) # distortion energy Mbar
        self.dtminj   = np.zeros(j) # array used to find next time step and debug EOS problems
        self.sigmar   = np.zeros(j) # total radial stress sigma_r, only evaluated at n=1 current time
        self.sigmao   = np.zeros(j) # total tangential stress sigma_theta, only evaluated at n=1 current time
        self.alocal   = np.zeros(j) # intermediate variable for AV
        self.rholocal = np.zeros(j) # intermediate variable for AV
        self.phase    = np.zeros(j) # phase identification
        self.bcpres   = np.zeros(j) # interior boundary pressures
    #
    def makegrid(self,run,verbose=False):
        """ Set up the initial spatial domain. 
            Determine boundary conditions/interfaces and domain array lengths. 
            Material is divided into a Lagrangian grid:
                Space between consecutive grid lines is a zone.
                The intersection of grid lines are called zone node points.
                In 1D the node points are the edges of the zones.
                
                [make sure added interfaces follow this rule]
                Position and velocity are evaluated at even-j zone node points
                PVE etc. are evaluated at odd-j zone centers
        
                ibc defines the type of each domain array entry.
                ibc = -1 or 1 = second order left stencil or second order right stencil
                                type of interface defined by ghost cells
                                FREE, FIXED, INFLOW, OUTFLOW
                ibc = -2, -3, etc. and +2, +3, etc. = paired internal boundaries (planes or fractures)
                ibc = 0 = interior j that is not at an interface
        
                Extra domain array entries are needed for material interfaces.
                Borg has a contact scheme for interfaces. CONTACT IS NOT IN WILKINS BOOK.
            Initial domain is placed in t(2) and t(3). For consistency in output indexing.
            t(2:3) -> t(0:1) update takes place at the START of the main advance time function.
            This allows all t(0:3) variables to be populated at the end of the first loop for 
            plotting and debugging.
            Output functions should use t(2:3) half steps.
        """
        ureg.disable_contexts()
        if verbose: print("INITIALIZING SIMULATION DOMAIN")
        # initializing these values here; anything that needs the input file parameters in run class
        self.jtot = npsum(run.inodes)+1
        self.deltat   = run.dtstart      # current time step; while in time loop
        self.deltat_next = run.dtstart   # next time step; calculated at end of time loop
        # variables in time and space (n,j)
        n = self.ntot                 # running time step index
        j = self.jtot                 # spatial domain index
        # these variables are evaluated at n=1; current time step
        # these variables are evaluated at n=2; half time step ahead
        self.up       = np.zeros((n,j)) # particle velocity
        self.pos      = np.zeros((n,j)) # position in 1D geometry (x,r_cyl,r_sph)
        self.vr       = np.zeros((n,j)) # relative volume to initial state=rho0/rho=1/eta in Wilkins [dimless]
        self.phi      = np.zeros((n,j)) # rho0*(pos_j+1-pos_j)/vr_j+1/2 = (pos_j+1-pos_j)/rho_j+1/2=phi in Wilkins
        self.q        = np.zeros((n,j)) # artificial viscosity Mbar
        self.eps1     = np.zeros((n,j)) # velocity strain dup/dpos, per microsec
        self.eps2     = np.zeros((n,j)) # velocity strain up/pos, per microsec
        # these variables are evaluated at n=3; full time step ahead
        self.beta     = np.zeros((n,j)) # (sigmar-sigmao)/(0.5 dpos)*rho = beta in Wilkins
        self.iev0     = np.zeros((n,j)) # internal energy per initial volume 1e12 erg/cm3
        self.pres     = np.zeros((n,j)) # pressure Mbar
        self.s1       = np.zeros((n,j)) # s1 stress deviator deriv = 2G(deps1-dvr/vr/3)
        self.s2       = np.zeros((n,j)) # s2 stress deviator deriv
        self.s3       = np.zeros((n,j)) # s3 stress deviator deriv; not needed in 1D; kept for equation clarity
        self.temp     = np.zeros((n,j)) # temperature K
        self.entropy  = np.zeros((n,j)) # entropy [DEBUG check eu/cm3 units]
        self.matid    = np.zeros(j) # material id number
        self.ibc      = np.zeros(j,dtype='int') # boundary condition id number, integer
        self.yld      = np.zeros(j) # yield stress
        self.pfrac    = np.zeros(j) # fracture stress
        self.rho0     = np.zeros(j) # initial density rho0 g/cm3; DEBUG later can free up this memory by referring to imat definition
        self.rho      = np.zeros(j) # local density g/cm3
        self.mass     = np.zeros(j) # mass in node g
        self.deltaz   = np.zeros(j) # distortion energy Mbar
        self.dtminj   = np.zeros(j) # array used to find next time step and debug EOS problems
        self.sigmar   = np.zeros(j) # total radial stress sigma_r, only evaluated at n=1 current time
        self.sigmao   = np.zeros(j) # total tangential stress sigma_theta, only evaluated at n=1 current time
        self.alocal   = np.zeros(j) # intermediate variable for AV
        self.rholocal = np.zeros(j) # intermediate variable for AV
        self.phase    = np.zeros(j) # phase information
        self.bcpres   = np.zeros(j) # phase information
        #
        # leave time array at all zeros here; updated at beginning of advance time step function
        # loop over each material layer
        # INITIAL VARIABLES: pos,ibc,matid
        nextj = 0
        for imat in np.arange(run.nmat):
            # set t(0) and t(1) to the initial state
            j=nextj # initial node for this layer of material
            #
            dr = run.ilength[imat]/run.inodes[imat]
            #if verbose: print('imat, j, dr = ',imat,j,dr)
            self.pos[0,j:j+run.inodes[imat]] = np.arange(run.inodes[imat])*dr+run.ixstart[imat]
            self.pos[1,j:j+run.inodes[imat]] = np.arange(run.inodes[imat])*dr+run.ixstart[imat]
            self.pos[2,j:j+run.inodes[imat]] = np.arange(run.inodes[imat])*dr+run.ixstart[imat]
            self.pos[3,j:j+run.inodes[imat]] = np.arange(run.inodes[imat])*dr+run.ixstart[imat]
            self.ibc[j:j+run.inodes[imat]]   = 0 # interior cells
            self.matid[j:j+run.inodes[imat]] = imat
            if imat == 0:
                # first material has ibc assignment
                if run.ibcid[imat] == 'FREE':
                    self.ibc[j] = -1 # left free surface
                elif run.ibcid[imat] == 'FIXED':
                    self.ibc[j] = -4 # left rigid surface
                #if verbose: print('imat=',imat,' left surface j=',j, ' ibc=',self.ibc[j])
            # check if the last material
            if imat == run.nmat-1:
                # last material
                self.pos[0,j+run.inodes[imat]]=self.pos[2,j+run.inodes[imat]-1]+dr # right surface
                self.pos[1,j+run.inodes[imat]]=self.pos[3,j+run.inodes[imat]-1]+dr # right surface
                self.pos[2,j+run.inodes[imat]]=self.pos[2,j+run.inodes[imat]-1]+dr # right surface
                self.pos[3,j+run.inodes[imat]]=self.pos[3,j+run.inodes[imat]-1]+dr # right surface
                self.matid[j+run.inodes[imat]]=imat
                if run.ibcid[1] == 'FREE': # DEBUG right now no interior gaps
                    self.ibc[j+run.inodes[imat]]=1 # right free surface
                elif run.ibcid[1] == 'FIXED': # DEBUG right now no interior gaps
                    self.ibc[j+run.inodes[imat]]=4 # right rigid surface
                #if verbose: print('imat=',imat,' right surface j=',j+run.inodes[imat]+1, ' ibc=',self.ibc[j+run.inodes[imat]])
            else:
                # start j for next material
                nextj=j+run.inodes[imat]
        # Check for interior boundaries and problems with the grid
        # loop over 1st through nmat-1 material layers
        ninteriorbc = 0
        for imat in np.arange(run.nmat)[0:-1]:
            # check if next layer overlaps
            if (run.ixstart[imat+1] < run.ixstart[imat]+run.ilength[imat]):
                print('FATAL ERROR: materials overlap.')
                sys.exit("FATAL ERROR: Next layer overlaps with current layer! imat="+str(imat)) 
            # check for gap
            if (run.ixstart[imat+1] > run.ixstart[imat]+run.ilength[imat]):
                print('GAP: ',imat,imat+1)
                # add a 3 nodes: even right interior boundary, odd void node, even left interior boundary
                if ninteriorbc == 0:
                    ninteriorbc = 101
                else:
                    ninteriorbc = ninteriorbc+1
                self.createinteriorboundary(run,imat,ninteriorbc)
                print('Created gap.')
                print('')
                #sys.exit("FATAL ERROR: Interior gaps are not coded yet.") 
        # ------------
        # have a spatial mesh, matid, and boundary conditions
        # PREVIOUSLY FILLED VARIABLES: pos,matid,ibc
        # now fill in material initial conditions
        for imat in np.arange(run.nmat):
            # initialize velocities; here for all domain array points to start
            self.up[:,self.matid==imat] = run.iupstart[imat]
            #self.up[:,max(self.matid==imat)] = run.iupstart[imat]/2.
            # initialize all array points with initial rho0
            self.rho0[self.matid==imat] = run.irhostart[imat] # user input
            self.rho[self.matid==imat] = run.irhostart[imat] # user input
            self.vr[:,self.matid==imat] = 1.0 # relative volume to initial state == 1.0
            self.iev0[:,self.matid==imat] = run.iiev0start[imat] # user input
            self.pres[:,self.matid==imat] = run.ipstart[imat] # user input
            # DEBUG: initialization for entropy yet 
            # DEBUG: need to propertly initialize temps for all EOS
            self.temp[:,self.matid==imat] = run.itempstart[imat] # starting temp for ideal gas E/cv
            # initialize the material properties
            # DEBUG ADD FLAGS FOR DIFFERENT STRENGTH MODELS
            # FOR NOW JUST HYDRO and VM
            if run.istrid[imat] == 'VM':
                self.yld[self.matid==imat]=run.istr[imat].ys
                #self.pfrac[self.matid==imat]=run.istr[imat].pfrac
            #DEBUG add fracture strength later
            #if run.istrid[imat] == 'HYDRO':
            #    self.pfrac[self.matid==imat]=run.istr[imat].pfrac
            # for any material enter pfrac
            self.pfrac[self.matid==imat]=run.ifrac[imat].pfrac
            ### Wilkins Appendix B.2 Eq 1 Mass zoning
            # initialize mass in the zones j=odd
            # select interior array points for this material; exclude boundaries
            tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
            ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
            self.mass[ioddimat] = run.irhostart[imat] * \
                ( nppower(self.pos[2,ioddimat+1],run.dgeom) - \
                  nppower(self.pos[2,ioddimat-1],run.dgeom) ) / run.dgeom
            # DEBUG: initialize sound speeds and temperature properly; depends on EOS model
            # add phase information for some tabular EOS
            if (run.ieosid[imat] == 'SES'):
                sesphase = run.ieos[imat].sesphase(self.rho[ioddimat],self.iev0[2,ioddimat]/self.rho0[ioddimat])
                self.phase[ioddimat] = sesphase
                #print('sesphase = ',sesphase)
        # -----------------------------------------------------------
        # GRAVITATIONAL EQUILIBRATION DEVELOPMENT
        # must use tabular EOS; ONLY WORKS FOR ISOTHERMAL INITIALIZATION AT THIS TIME
        for imat in np.arange(run.nmat):
            if run.grav.matflag[imat]:
                print('Need to initialize gravity for mat ',imat)
                self.initgravity(run,imat)
        # -----------------------------------------------------------
        # ESTIMATE FIRST TIME STEP
        # select interior cells with material parameters
        n=1 # current time step index
        tmp = npwhere(self.ibc == 0)[0] # select all interior array indices
        ioddj = tmp[tmp %2 == 1] # select only odd indices for zone values
        #
        # estimate a first time step; user has input flag to override
        # calculate the spatial steps
        dx0 = self.pos[2,ioddj]-np.roll(self.pos[2,ioddj],1) # cm
        dx0[0]=dx0[1]
        up0 = 0.5*(self.up[n+1,ioddj-1]+self.up[n+1,ioddj+1]) # cm/microsec
        iuppos0 = npwhere(up0 != 0.)[0]
        if len(iuppos0) == 0:
            # check for initial pressure
            cs0 = npsqrt(npabs(self.pres[n+2,ioddj]))/self.rho[ioddj]
            ippos0 = npwhere(self.pres[n+2,ioddj]>run.pvoid) # where initialize greater than pvoid
            dt0 = dx0[ippos0]/cs0[ippos0]
            run.dtstart = np.min(dt0)/10.
        else:
            dt0 = np.abs(dx0[iuppos0]/up0[iuppos0])
            run.dtstart = np.min(dt0)/10.
        self.deltat   = run.dtstart      # initialize the domain object
        self.deltat_next = run.dtstart   # initialize the domain object        
        if verbose: print('pyKO FIRST TIME STEP (microsec) = ',run.dtstart)
        # ------------
        # calculate initial energies
        self.qtotal  = npsum(self.q[n+1,ioddj])
        self.mvtotal = npsum(0.5*self.mass[ioddj]*(self.up[n+1,ioddj-1]+self.up[n+1,ioddj+1]))
        self.ketotal = npsum(0.5*self.mass[ioddj]*np.square(0.5*(self.up[n+1,ioddj-1]+self.up[n+1,ioddj+1])))
        self.ietotal = npsum( (self.iev0[n+2,ioddj] /self.rho0[ioddj]*self.mass[ioddj]) )
        self.etotal  = npsum( (self.iev0[n+2,ioddj]-self.iev0[n,ioddj])/self.deltat + \
                          self.pres[n+1,ioddj]*(self.vr[n+2,ioddj]-self.vr[n,ioddj])/self.deltat )
        self.mtotal  = npsum(self.mass[ioddj])
        #
    def initgravity(self,run,imat,verbose=True):
        """ Initialize material layer with gravity. 
            Currently only isothermal. """
        print('Next initialize gravity for material',imat)
        # initially only SESAME table
        tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
        ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
        print('Start grav initialization')
        # start at the refpos and integrate down into the gravity well
        lastP = deepcopy(run.grav.refpres[imat])
        lastx = deepcopy(run.grav.refpos[imat])
        Pbotlast = (self.pos[0,ioddimat[0]]-lastx)*self.rho[ioddimat[0]]*run.grav.gravity
        print('position check; P0,roughPbot,Pbot=',lastP,Pbotlast,self.pres[0,ioddimat[0]])
        if run.grav.refpos[imat] > np.max(self.pos[0,ioddimat]):
            # run twice
            for iloop in range(1):
                lastP = deepcopy(run.grav.refpres[imat])
                lastx = deepcopy(run.grav.refpos[imat])
                lastrho = deepcopy(run.irhostart[imat])
                #print('start:', lastP,lastx,lastrho)
                for jjj in np.flip(ioddimat[::]):
                    dP = run.grav.gravity*lastrho*(self.pos[0,jjj]-lastx)
                    self.pres[:,jjj] = lastP + dP
                    rho,u,c = run.ieos[imat].oneruc(self.pres[0,jjj],self.temp[0,jjj])
                    self.rho[jjj] = rho
                    self.vr[:,jjj] = self.rho0[jjj]/rho
                    self.iev0[:,jjj] = u*self.rho0[jjj]
                    self.alocal[jjj] = c
                    lastP = deepcopy(self.pres[0,jjj])
                    lastx = deepcopy(self.pos[0,jjj])
                    lastrho = deepcopy(rho)
                    self.mass[jjj] = self.rho[jjj] * \
                        ( nppower(self.pos[2,jjj+1],run.dgeom) - \
                          nppower(self.pos[2,jjj-1],run.dgeom) ) / run.dgeom
                    #print(jjj,self.pos[0,jjj],self.pres[0,jjj],self.rho[jjj],self.temp[0,jjj],self.alocal[jjj])
        # set the inner boundary condition to the pressure of the lowermost cell
        run.ibc[0].pbc = self.pres[0,ioddimat[0]]
        print('setting inner boundary pressure = ',run.ibc[0].pbc)
        return  
    def advancetime(self,run,verbose=False):
        """ Advance hydro one time step. """
        ureg.disable_contexts()
        ### In Borg code conversion from all time steps to 4 half steps
        ### led to a v11 code with n=1 fixed value for indexing
        ### ntot=4: t(0) and t(1) current, t(2) and t(3) new
        ### In Wilkins notation, time centering finite difference:
        ### delta_t^n = t^n+1/2 - t^n-1/2 used with pressure field to advance 
        ###             Lagrangian nodes from time t^n-1/2 to t^n+1/2
        ### Position of mesh points (and thermo centered values) are advanced
        ### delta_t^n+1/2 = t^n+1-t^n == true time step
        ###                 determined from the stability conditions with 
        ###                 delta_t^n = (1/2)(delta_t^n+1/2 + delta_t^n-1/2)
        ### Time centering is upset when stability conditions permit deviation
        ### from above. Mapping the notations:
        ### t(0) = n-1/2
        ### t(1) = n
        ### t(2) = n+1/2
        ### t(3) = n+1
        ### True time step = t(3)-t(1)
        ### fixing n=1 identifies relative steps in time and makes 
        ### code translation easier
        ###
        ### Reminder: odd j are mesh centers; even j are nodes
        ###  are defined on even j nodes
        ### q,sigmar,sigmao,p,s1,s2,s3 are defined on odd j mesh centers
        ###
        #if verbose: print("\n ADVANCING ONE TIME STEP \n")
        ### 0. UPDATE main arrays
        # The first step in the main time loop is to advance the time
        # half steps. This is done at the beginning of the loop so that
        # all the half step variables are populated at the end of the loop
        # for plotting and debugging and output.
        # The initialization function must populate t(2) and t(3)
        # which are then shifted to t(0) and t(2) by the shiftdata function.
        # This order requires storing current and next time steps.
        # and updating as well.
        self.shiftdata() # shifts time(n) and cell(n,j) arrays
        #if verbose: print('dt, dt_next=',self.deltat,self.deltat_next,self.time)
        #
        # 0.5. Check for contact; if yes, adjust time step
        if max(self.ibc[1:-1]) > 0:
            # there are gaps so check for contact
            # if there is, then resets self.deltat_next
            self.checkforcontact(run)
        #
        ### 1. TIME
        ### Wilkins Section 2.1.1
        n=1 # FIXED CONSTANT FOR CURRENT TIME INDEX
        self.stepn += 1 # advance the step counter
        self.deltat = np.copy(self.deltat_next) # update the time step
        #self.deltat_next = 0.0 # reset the nnext time step
        self.time[n+1] = self.time[n]+self.deltat/2. # plus half step
        self.time[n+2] = self.time[n]+self.deltat    # plus whole step
        dthalf = self.time[n+1]-self.time[n-1] # Wilkins delta_t^n DEBUG-CHECK THIS vs p28
        #if verbose: print('Time step t_n, t_n+1:',self.time[n],self.time[n+2])
        #if verbose: print('dt, dt_next=',self.deltat,self.deltat_next,self.deltat_next/self.deltat,self.time)
        ### DEBUG future code will want to add a contact check
        ### to decrease time steps for closing voids
        ###
        ### 2. CONSERVATION OF MOMENTUM
        # ------------ calculate interior even and odd indices ----------------
        # re-create these arrays every time loop in case of fractures
        # all interior odd j nodes array indices
        tmp = npwhere(self.ibc == 0)[0] # select all interior array indices
        ioddj = tmp[tmp %2 == 1] # select only odd indices for zone values
        #print(ioddj)
        # all interior even j nodes array indices
        tmp = npwhere(self.ibc == 0)[0] # select all interior array indices
        ievenj = tmp[tmp %2 == 0] # select only even indices for zone values
        #print(ievenj)
        # this is just used to update positions at every node, including bc
        allevenj = np.arange(0,self.jtot,2) # select ALL EVEN - interior and bc
        #print(allevenj)
        alloddj = np.arange(1,self.jtot,2) # select ALL ODD - interior and bc
        # ------------ calculate inner boundary condition = -1 ----------------
        # Need sigma_r and sigma_theta for boundary condition
        self.sigmar[ioddj] = -(self.pres[n,ioddj]+self.q[n-1,ioddj])+self.s1[n,ioddj]
        self.sigmao[ioddj] = -(self.pres[n,ioddj]+self.q[n-1,ioddj])+self.s2[n,ioddj]
        # ------------ calculate inner boundary condition = -1 ----------------
        # boundary condition section must set the sigmas for the type of BC
        ibcj = 0 # faster than npwhere(self.ibc == -1)[0]; change if needed
        #if run.ibcid[0] == 'FREE': # DEBUG should this be initialization section?
            # Free surface condition, pressure is set to zero in the ghost cell
            # IF INPUT PRESSURES ARE POSITIVE MUST FLIP THE SIGN TO NEGATIVE HERE
        #    ighost_p = -run.ibc[0].pbc #for BORG input
        # the following is generic and depends on the boundary condition sigmas
        # use vector math to update sigmas on all oddj interior indices
        self.phi[n,ibcj] = 0.5*self.rho0[ibcj+1]* \
            (self.pos[n,ibcj+2]-self.pos[n,ibcj])/self.vr[n,ibcj+1]
        self.beta[n,ibcj] = (self.sigmar[ibcj+1]-self.sigmao[ibcj+1])/ \
             (0.5*(self.pos[n,ibcj]+self.pos[n,ibcj+2]))* \
                (self.vr[n,ibcj+1]/self.rho0[ibcj+1])
#         self.beta[n,ibcj] = (self.sigmar[ibcj+1]-self.sigmao[ibcj+1])/ \
#             ((self.pos[n,ibcj+1]))* \
#                (self.pos[n,ibcj+1])*(self.vr[n,ibcj+1]/self.rho0[ibcj+1])
        if run.ibcid[0] == 'FIXED': # ibc=-4
            self.up[n+1,ibcj] = 0. # fix the boundary velocity to zero
            self.pres[n+1,ibcj] = self.pres[n,ibcj+1]  # DEBUG sts added in v0.4 to test symmetry boundary 
            #print(self.pres[n,ibcj+1])
        if run.ibcid[0] == 'FREE':
            # ibc=-1
            ighost_p = -run.ibc[0].pbc #for BORG input
            # inner sigmar is ghost_p
            self.up[n+1,ibcj] = self.up[n-1,ibcj] + dthalf / self.phi[n,ibcj] * \
                                  (self.sigmar[ibcj+1]-ighost_p) + \
                                  dthalf * self.beta[n,ibcj] * (run.dgeom-1.)       
        #print('ibc ',self.stepn,self.phi[n,ibcj],self.beta[n,ibcj],self.up[n+1,ibcj],self.vr[n,ibcj+1],self.sigmar[ibcj+1],self.sigmao[ibcj+1])
        #print('dthalf ',dthalf,self.rho0[ibcj+1],self.pos[n,ibcj+1])
        # ------------ calculate outer boundary condition = 1 -----------------
        # boundary condition section must set the sigmas for the type of BC
        obcj = self.jtot-1 # faster than npwhere(self.ibc == 1)[0]; change if needed
        #if run.ibcid[run.nbc-1] == 'FREE':
            # Free surface condition, pressure is set to zero in the ghost cell
            # IF INPUT PRESSURES ARE POSITIVE MUST FLIP THE SIGN TO NEGATIVE HERE
            #oghost_p = -run.ibc[run.nbc-1].pbc #for BORG input
        # the following is generic and depends on the boundary condition sigmas
        self.phi[n,obcj] = 0.5*self.rho0[obcj-1]* \
            (self.pos[n,obcj]-self.pos[n,obcj-2])/self.vr[n,obcj-1]
        self.beta[n,obcj] = (self.sigmar[obcj-1]-self.sigmao[obcj-1])/ \
             (0.5*(self.pos[n,obcj]+self.pos[n,obcj-2]))* \
                (self.vr[n,obcj-1]/self.rho0[obcj-1])
#        self.beta[n,obcj] = (self.sigmar[obcj-1]-self.sigmao[obcj-1])/ \
#            ((self.pos[n,obcj-1]))* \
        if run.ibcid[run.nbc-1] == 'FIXED':
            self.up[n+1,obcj] = 0. # fixed boundary condition ibc=4
            self.pres[n+1,obcj] = self.pres[n,obcj-1]  # DEBUG sts added in v0.4 to test symmetry boundary 
        if run.ibcid[run.nbc-1] == 'FREE':
            # free boundary ibc=1
            oghost_p = -run.ibc[run.nbc-1].pbc #for BORG input
            # outer sigmar is ghost_p; DEBUG borg has +ghost_p using Wilkins here
            self.up[n+1,obcj] = self.up[n-1,obcj] + dthalf / self.phi[n,obcj] * \
                                  (oghost_p-self.sigmar[obcj-1]) + \
                                  dthalf * self.beta[n,obcj] * (run.dgeom-1.)        
        #print('obc ',self.stepn,self.phi[n,obcj],self.beta[n,obcj],self.up[n+1,obcj],self.vr[n,obcj-1])
        #print('rho0 ',dthalf,self.rho0[obcj-1],self.pos[n,obcj-1])
        #
        # ------------ calculate interior boundaries --------------------------
        # ONLY FREE INTERIOR BOUNDARIES CODED
        # Void spaces can be found by matid=-1
        ivoid = np.where(self.matid == -1)[0]
        if len(ivoid) > 0:
            # there are interior boundaries
            for jvoid in ivoid:
                #print('jvoid ibc = ',jvoid,self.ibc[jvoid-4:jvoid+4])
                #print('jvoid mat = ',jvoid,self.matid[jvoid-4:jvoid+4])
                #for imat in np.arange(run.nmat-1):
                #jend = max(np.where(self.matid == imat)[0])
                jend = jvoid-1
                #print('Calculating interior boundaries ',imat, jend)
                #print('')
                # there is a gap
                # outer boundary = jend
                obcj = jend 
                # Free surface condition
                # phase holds the pressure in the interior boundary pressure
                # don't forget the stress is negative!
                oghost_p = 0.#-self.bcpres[obcj+1]
                #print('obcj, pres = ',obcj+1,self.bcpres[obcj+1])
                # the following is generic and depends on the boundary condition sigmas
                self.phi[n,obcj] = 0.5*self.rho0[obcj-1]* \
                    (self.pos[n,obcj]-self.pos[n,obcj-2])/self.vr[n,obcj-1]
                self.beta[n,obcj] = (self.sigmar[obcj-1]-self.sigmao[obcj-1])/ \
                        (0.5*(self.pos[n,obcj]+self.pos[n,obcj-2]))* \
                            (self.vr[n,obcj-1]/self.rho0[obcj-1])
                # outer sigmar is ghost_p; DEBUG borg has +ghost_p using Wilkins here
                # add gravitational acceleration
                self.up[n+1,obcj] = self.up[n-1,obcj] + dthalf / self.phi[n,obcj] * \
                                       (oghost_p-self.sigmar[obcj-1]) + \
                                       dthalf * self.beta[n,obcj] * (run.dgeom-1.) + \
                                       dthalf * run.grav.gravity
                # inner boundary is jend +2
                # don't forget the stress is negative!
                ibcj = jend+2
                ighost_p = 0.#-self.bcpres[ibcj+1]
                #print('ibcj, pres = ',ibcj+1,self.bcpres[ibcj+1])
                # the following is generic and depends on the boundary condition sigmas
                # use vector math to update sigmas on all oddj interior indices
                self.phi[n,ibcj] = 0.5*self.rho0[ibcj+1]* \
                    (self.pos[n,ibcj+2]-self.pos[n,ibcj])/self.vr[n,ibcj+1]
                self.beta[n,ibcj] = (self.sigmar[ibcj+1]-self.sigmao[ibcj+1])/ \
                     (0.5*(self.pos[n,ibcj]+self.pos[n,ibcj+2]))* \
                        (self.vr[n,ibcj+1]/self.rho0[ibcj+1])
                # inner sigmar is ghost_p
                # add gravitational acceleration
                self.up[n+1,ibcj] = self.up[n-1,ibcj] + dthalf / self.phi[n,ibcj] * \
                                    (self.sigmar[ibcj+1]-ighost_p) + \
                                     dthalf * self.beta[n,ibcj] * (run.dgeom-1.) + \
                                     dthalf * run.grav.gravity
        #
        # ------------ advance interior cell velocities even j ----------------
        ### Wilkins Appendix B.2 Eq 2 Equation of motion
        # use vector math to update sigmas on all oddj interior indices
        # Need sigma_r and sigma_theta for velocity field update
        #self.sigmar[ioddj] = -(self.pres[n,ioddj]+self.q[n-1,ioddj])+self.s1[n,ioddj]
        #self.sigmao[ioddj] = -(self.pres[n,ioddj]+self.q[n-1,ioddj])+self.s2[n,ioddj]
        # Calculate phi and beta for interior evenj nodes
        self.phi[n,ievenj]  = 0.5*(self.rho0[ievenj+1]/self.vr[n,ievenj+1]) * \
                                  (self.pos[n,ievenj+2]-self.pos[n,ievenj]) + \
                              0.5*(self.rho0[ievenj-1]/self.vr[n,ievenj-1]) * \
                                  (self.pos[n,ievenj]-self.pos[n,ievenj-2])  
        self.beta[n,ievenj] = 0.5*((self.vr[n,ievenj+1]/self.rho0[ievenj+1]) * \
                                   ( (self.sigmar[ievenj+1]-self.sigmao[ievenj+1]) / \
                                    (0.5*(self.pos[n,ievenj+2]+self.pos[n,ievenj])) ) + \
                                   (self.vr[n,ievenj-1]/self.rho0[ievenj-1]) * \
                                   ( (self.sigmar[ievenj-1]-self.sigmao[ievenj-1]) / \
                                    (0.5*(self.pos[n,ievenj]+self.pos[n,ievenj-2])) ) )
        # Update velocities for interior evenj nodes
        # v0.5-dev add gravitational force
        self.up[n+1,ievenj] = self.up[n-1,ievenj] + dthalf / self.phi[n,ievenj] * \
                              (self.sigmar[ievenj+1]-self.sigmar[ievenj-1]) + \
                              dthalf * self.beta[n,ievenj] * (run.dgeom-1.) + \
                              dthalf * run.grav.gravity
        # update velocities on ioddj nodes for debugging
        self.up[n+1,ioddj] = 0.5*(self.up[n+1,ioddj+1]+self.up[n+1,ioddj-1])
        #
        # ------------ advance node position and relative volumes -------------
        ### Wilkins Appendix B.2 Eq 3 Conservation of mass
        # DEBUG will need to check this with interior interfaces
        # DEBUG check if need to re-center
        # use vector math to update positions on all even indices (include bc)
        # evaluated at pos[n+2] and midpoint position between at [n+1]
        self.pos[n+2,allevenj] = self.pos[n,allevenj] + self.up[n+1,allevenj] * self.deltat
        self.pos[n+1,allevenj] = 0.5*(self.pos[n,allevenj]+self.pos[n+2,allevenj])
        # update the odd positions to help with debugging
        self.pos[n+2,ioddj] = 0.5*(self.pos[n+2,ioddj+1]+self.pos[n+2,ioddj-1])
        self.pos[n+1,ioddj] = 0.5*(self.pos[n+1,ioddj+1]+self.pos[n+1,ioddj-1])
        # update void positions to help with debugging
        ivoid = np.where(self.matid == -1)[0]
        if len(ivoid):
            self.pos[n+2,ivoid] = 0.5*(self.pos[n+2,ivoid+1]+self.pos[n+2,ivoid-1])
            self.pos[n+1,ivoid] = 0.5*(self.pos[n+1,ivoid+1]+self.pos[n+1,ivoid-1])
        # use vector math to update all interior relative volumes
        # technically should only be evaluated at vr[n+2]
        # [n+1] should represent the midpoint relative volume based on pos[n+1] above
        self.vr[n+1,ioddj] = self.rho0[ioddj]/self.mass[ioddj] * \
                             (nppower(self.pos[n+1,ioddj+1],run.dgeom) - \
                              nppower(self.pos[n+1,ioddj-1],run.dgeom)) /\
                             run.dgeom
        self.vr[n+2,ioddj] = self.rho0[ioddj]/self.mass[ioddj] * \
                             (nppower(self.pos[n+2,ioddj+1],run.dgeom) - \
                              nppower(self.pos[n+2,ioddj-1],run.dgeom)) / \
                             run.dgeom
        # check for negative volumes
        if len(npwhere(self.vr[n+1,ioddj] <= 0)[0]) > 0:
            sys.exit("FATAL ERROR: Negative or zero volume in data.vr[n+1]")
        if len(npwhere(self.vr[n+2,ioddj] <= 0)[0]) > 0:
            sys.exit("FATAL ERROR: Negative or zero volume in data.vr[n+2]")
        #
        # ------------ update velocity strains interior odd j ---------------------
        ### Wilkins Appendix B.2 Eq 4 Calculation of velocity strains
        # Velocities are calculated on the half time steps
        self.eps1[n+1,ioddj] = (self.up[n+1,ioddj+1]-self.up[n+1,ioddj-1]) / \
                               (self.pos[n+1,ioddj+1]-self.pos[n+1,ioddj-1])
        if (run.dflag == 'PLA'):  # use string for boolean checks
            self.eps2[n+1,ioddj] = 0. # zero for planar geometry
        else:
            self.eps2[n+1,ioddj] = (self.up[n+1,ioddj+1]+self.up[n+1,ioddj-1]) / \
                                    (self.pos[n+1,ioddj+1]+self.pos[n+1,ioddj-1])
        #
        # ------------ update stresses for interior odd j ---------------------
        ### Wilkins Appendix B.2 Eq 5 Calculation of stresses
        # Stress deviators s1, s2, s3; all interior odd j        
        for imat in np.arange(run.nmat):
            # select interior array points for this material
            tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
            ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
            if (run.istrid[imat] == 'VM'):
                gmod = run.istr[imat].gmod # Mbar
            if (run.istrid[imat] == 'HYDRO'):
                gmod = 0. # Mbar
                #if verbose: print('imat, shear modulus',imat,gmod)
            self.s1[n+2,ioddimat] = self.s1[n,ioddimat] + 2. * gmod * \
                             (self.eps1[n+1,ioddimat]*self.deltat - \
                              (1./3.)*(self.vr[n+2,ioddimat]-self.vr[n,ioddimat])/self.vr[n+1,ioddimat])
            self.s2[n+2,ioddimat] = self.s2[n,ioddimat] + 2. * gmod * \
                             (self.eps2[n+1,ioddimat]*self.deltat - \
                              (1./3.)*(self.vr[n+2,ioddimat]-self.vr[n,ioddimat])/self.vr[n+1,ioddimat])
            self.s3[n+2,ioddimat] = -(self.s1[n+2,ioddimat] + self.s2[n+2,ioddimat])
            #jtest = 199
            #print('S1,S2,S3,eps1,dt= ',self.s1[n+2,jtest],self.s2[n+2,jtest],self.s3[n+2,jtest],self.eps1[n+1,jtest],self.deltat )                                 
        #
        # DEBUG Borg fills in the staggered mesh with averages; not sure I need that
        # updating pressure for next time step requires an EOS; moved down in code order
        # Wilkins uses present time step pressure[n] to calculate
        # sqrt(P/rho) artificial viscosity term
        # not using the actual EOS sound speed
        #
        # ------------ check for yielding interior odd j ---------------------
        ### Wilkins Appendix B.2 Eq 6 Von Mises yield condition
        # check for strength model; can be different for each material
        # material specific code is unnecessary if every material is VM yield strength 
        #    and each material yield values are already stored in self.yld
        # keep this code for future development of strain-dependent strength models
        for imat in np.arange(run.nmat):
            # select interior array points for this material
            tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
            ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
            if (run.istrid[imat] == 'VM'):
                # calculate Von Mises yield criterion
                # must have already populated primary self.yld array
                # calculate the deviatoric strain at n+1 and compare it to the yield strength
                k = (np.square(self.s1[n+2,ioddimat]) + \
                     np.square(self.s2[n+2,ioddimat]) + \
                     np.square(self.s3[n+2,ioddimat])) - \
                    (2./3.)*np.square(self.yld[ioddimat])
                #if verbose: print('imat, yield location =',imat, npwhere(k > 0.)[0])
                iyld = ioddimat[npwhere(k > 0.)[0]] # these array points have failed
                if len(iyld)>0:
                    # material has yielded
                    #if verbose: print('Material has yielded: imat, j vals ',imat,iyld)
                    # multiply the stress deviators to reduce to max shear stresses
                    sscale = npsqrt(2./3.)*self.yld[iyld] / \
                             npsqrt(np.square(self.s1[n+2,iyld]) + \
                                     np.square(self.s2[n+2,iyld]) + \
                                     np.square(self.s3[n+2,iyld]))
                    self.s1[n+2,iyld] *= sscale 
                    self.s2[n+2,iyld] *= sscale 
                    self.s3[n+2,iyld] *= sscale
                    #if verbose: print('sscale, self.s1=',sscale,self.s1[n+2,iyld])
        #
        # ------------ update artificial viscosity on interior odd j ---------------------
        ### Wilkins Appendix B.2 Eq 7 Artificial viscosity
        # set artificial viscosity and local sound speed arrays to zero to start
        self.q[n+1,ioddj] = 0.   
        self.alocal[ioddj] = 0.
        # calculate true density for each ioddj cell
        self.rholocal[ioddj] = self.rho0[ioddj]/self.vr[n+1,ioddj] # density for n+1 time
        # ----------------------------- CALCULATE SOUND SPEED ------------------------
        # At this point only p[n,:] and u[n,:] are available; p[n+3,:] is calculated below
        # a is sound speed of an ideal gas - use this as the reference case
        self.alocal[ioddj] = npsqrt(npabs(self.pres[n,ioddj]) / self.rho[ioddj])
        # Do not let sound speed become a nonphysically tiny number
        ipsmall = npwhere(self.pres[n,ioddj]<run.pvoid)[0]
        if len(ipsmall) > 0:
            self.alocal[ipsmall]=0.
        # sound speed adjustments for non-ideal gases
        for imat in np.arange(run.nmat):
            # select interior array points for this material
            tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
            ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
            ########### MIE GRUNEISEN EOS #################
            if (run.ieosid[imat] == 'MGR'):
                # use c0 value if greater than a
                icheck = npwhere(self.alocal[ioddimat] < run.ieos[imat].c0)[0]
                self.alocal[ioddimat[icheck]] = run.ieos[imat].c0
            ########### SESAME EOS #################
            if (run.ieosid[imat] == 'SES'):
                # use SES value if greater than a
                # interpolate for cs in the table using rho and specific internal energy 
                sescs = run.ieos[imat].sescs(self.rholocal[ioddimat],self.iev0[n,ioddimat]/self.rho0[ioddimat])
#                icheck = npwhere(self.alocal[ioddimat] < sescs)[0]
#                self.alocal[ioddimat[icheck]] = sescs[icheck]
                self.alocal[ioddimat] = sescs
        # ----------------------------- END SOUND SPEED CALC ------------------
        # AV calculation on the entire domain
        # Wilkins two criteria for applying artificial viscosity
        # - velocity gradient and compression
        # find the indices that satisfy the criteria
        iavon = ioddj[npwhere(((self.up[n+1,ioddj+1] < self.up[n+1,ioddj-1]) & \
                                ( (self.vr[n+2,ioddj] - self.vr[n,ioddj]) < 0. ) ))[0]]
        if len(iavon)>0:
            # calculate sound speed
            # pressure not defined yet at n+1; using n
            # sound speed for ideal gas using P[n] and rho[n+1]
            # make sure pressure is positive in the sqrt
            self.alocal[iavon] = npsqrt(npabs(self.pres[n,iavon])/self.rholocal[iavon])
             # Borg doesn't let a become a tiny number
            ipsmall = iavon[npwhere(npabs(self.pres[n,iavon])<run.pvoid)[0]] 
            if len(ipsmall) > 0:
                  self.alocal[ipsmall]=0.
            self.q[n+1,iavon] = np.square(run.avc0)*self.rholocal[iavon]* \
                               np.square(self.up[n+1,iavon+1]-self.up[n+1,iavon-1]) + \
                                run.avcl*self.alocal[iavon]*self.rholocal[iavon]* \
                                npabs(self.up[n+1,iavon+1]-self.up[n+1,iavon-1])
        #
        # ------------ update internal energy on interior odd j ---------------------
        ### Wilkins Appendix B.2 Eq 8 Energy equations
        # Main conservation questoin: delta-E = -(P+q)*delta-V + delta-Z
        # delta-Z is the change in distortion energy
        # loop through each material to use its own EOS model
        # the EOS provides the P-V-E relationships
        # The EOS-specific sections also set the local sound speed
        # used in the artificial viscosity section below.
              ### Within each material section; update E and then P
              ### Wilkins Appendix B.2 Eq 5: Update Pressure here after energy calc
              # ------------ update pressure on interior odd j ---------------------

        for imat in np.arange(run.nmat):
            # select interior array points for this material
            tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
            ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
            ########### MIE GRUNEISEN EOS #################
            if (run.ieosid[imat] == 'MGR'):
                # use Mie-Grueneisen EOS
                # fill in staggered mesh with s1[n+1] stress deviators
                # doing it here after any yielding calculation
                self.s1[n+1,ioddimat] = 0.5*(self.s1[n,ioddimat]+self.s1[n+2,ioddimat])
                self.s2[n+1,ioddimat] = 0.5*(self.s2[n,ioddimat]+self.s2[n+2,ioddimat])
                self.s3[n+1,ioddimat] = 0.5*(self.s3[n,ioddimat]+self.s3[n+2,ioddimat])
                # delta-Z calculated for n+1 time and no other time
                self.deltaz[ioddimat] = self.vr[n+1,ioddimat]*self.deltat* \
                          (self.s1[n+1,ioddimat]*self.eps1[n+1,ioddimat] + \
                           (run.dgeom-1.)*self.s2[n+1,ioddimat]*self.eps2[n+1,ioddimat])
                # calculate averaged artificial viscosity Mbar
                qbar = 0.5*(self.q[n+1,ioddimat]+self.q[n-1,ioddimat])
                ### Wilkins Appendix B.2 Eq 5 Update Pressure here after energy calc
                # ------------ update pressure on interior odd j ---------------------
                # Wilkins formulation for Mie-Grueneisen EOS
                strain = 1.-self.vr[n+2,ioddimat]
                isneg = npwhere(strain <  0.)[0] # based on strain
                ispos = npwhere(strain >= 0.)[0] # based on strain               
                if len(isneg)>0:
                    # negative strain full step ahead; drop the strain squared term
                    mga = run.ieos[imat].coefs[0]*strain[isneg]+\
                                           run.ieos[imat].coefs[2]*nppower(strain[isneg],3.)
                    mgb = run.ieos[imat].coefs[3]
                                           # calculate total internal energy E
                                           # uses current time step pressure P[n]
                    self.iev0[n+2,ioddimat[isneg]]=\
                       ( (self.iev0[n,ioddimat[isneg]]- \
                          (0.5*(mga+self.pres[n,ioddimat[isneg]])+qbar[isneg])* \
                          (self.vr[n+2,ioddimat[isneg]]-self.vr[n,ioddimat[isneg]])+ \
                              self.deltaz[ioddimat[isneg]]) / \
                         (1.+0.5*mgb*(self.vr[n+2,ioddimat[isneg]]-self.vr[n,ioddimat[isneg]]) ) )
                    self.pres[n+2,ioddimat[isneg]] = run.ieos[imat].coefs[0]*strain[isneg] + \
                        run.ieos[imat].coefs[2]*nppower(strain[isneg],3) + \
                        run.ieos[imat].coefs[3]*self.iev0[n+2,ioddimat[isneg]]
                if len(ispos)>0:
                    # positive strain full step ahead
                    mga = run.ieos[imat].coefs[0]*strain[ispos]+\
                                           run.ieos[imat].coefs[1]*np.square(strain[ispos])+\
                                           run.ieos[imat].coefs[2]*nppower(strain[ispos],3.)
                    mgb = run.ieos[imat].coefs[3]
                                           # calculate total internal energy E
                                           # uses current time step pressure P[n]
                    self.iev0[n+2,ioddimat[ispos]]=\
                       ( (self.iev0[n,ioddimat[ispos]]- \
                          (0.5*(mga+self.pres[n,ioddimat[ispos]])+qbar[ispos])* \
                          (self.vr[n+2,ioddimat[ispos]]-self.vr[n,ioddimat[ispos]])+ \
                              self.deltaz[ioddimat[ispos]]) / \
                         (1.+0.5*mgb*(self.vr[n+2,ioddimat[ispos]]-self.vr[n,ioddimat[ispos]]) ) )
                    self.pres[n+2,ioddimat[ispos]] = run.ieos[imat].coefs[0]*strain[ispos] + \
                        run.ieos[imat].coefs[1]*np.square(strain[ispos]) + \
                        run.ieos[imat].coefs[2]*nppower(strain[ispos],3) + \
                        run.ieos[imat].coefs[3]*self.iev0[n+2,ioddimat[ispos]]
                # calculate energy average for E[n+1]
                self.iev0[n+1,ioddimat]=0.5*(self.iev0[n,ioddimat]+self.iev0[n+2,ioddimat])                    
                # half step ahead n+1; calculate the associated pressure
                if len(isneg)>0:
                    # energy is interpolated above
                    self.pres[n+1,ioddimat[isneg]] = run.ieos[imat].coefs[0]*strain[isneg] + \
                        run.ieos[imat].coefs[2]*nppower(strain[isneg],3) + \
                        run.ieos[imat].coefs[3]*self.iev0[n+1,ioddimat[isneg]]
                if len(ispos)>0:
                    # energy is interpolated above
                    self.pres[n+1,ioddimat[ispos]] = run.ieos[imat].coefs[0]*strain[ispos] + \
                        run.ieos[imat].coefs[1]*np.square(strain[ispos]) + \
                        run.ieos[imat].coefs[2]*nppower(strain[ispos],3) + \
                        run.ieos[imat].coefs[3]*self.iev0[n+1,ioddimat[ispos]]
            ########### IDEAL GAS EOS #################
            if (run.ieosid[imat] == 'IDG'):
                # use Ideal Gas EOS
                # fill in staggered mesh with s1[n+1] stress deviators
                # doing it here after any yielding calculation
                # keep this until sure not needed for mixed material boundaries
                self.s1[n+1,ioddimat] = 0.5*(self.s1[n,ioddimat]+self.s1[n+2,ioddimat])
                self.s2[n+1,ioddimat] = 0.5*(self.s2[n,ioddimat]+self.s2[n+2,ioddimat])
                self.s3[n+1,ioddimat] = 0.5*(self.s3[n,ioddimat]+self.s3[n+2,ioddimat])
                # delta-Z calculated for n+1 time and no other time
                self.deltaz[ioddimat] = self.vr[n+1,ioddimat]*self.deltat* \
                          (self.s1[n+1,ioddimat]*self.eps1[n+1,ioddimat] + \
                           (run.dgeom-1.)*self.s2[n+1,ioddimat]*self.eps2[n+1,ioddimat])
                # calculate averaged artificial viscosity Mbar
                qbar = 0.5*(self.q[n+1,ioddimat]+self.q[n-1,ioddimat])
                # calculate mga and mbg at n+2
                #x = 1.-self.vr[n+2,ioddimat]
                mga = 0.
                mgb = (run.ieos[imat].gamma0-1.)/self.vr[n+2,ioddimat]
                # calculate total internal energy E
                # uses current time step pressure P[n]
                self.iev0[n+2,ioddimat]=\
                    ( (self.iev0[n,ioddimat]- \
                       (0.5*(mga+self.pres[n,ioddimat])+qbar)* \
                       (self.vr[n+2,ioddimat]-self.vr[n,ioddimat])+self.deltaz[ioddimat]) / \
                      (1.+0.5*mgb*(self.vr[n+2,ioddimat]-self.vr[n,ioddimat]) ) )
                # calculate energy average for E[n+1]
                self.iev0[n+1,ioddimat]=0.5*(self.iev0[n,ioddimat]+self.iev0[n+2,ioddimat])
                #
                # ------------ update pressure on interior odd j ---------------------
                self.pres[n+2,ioddimat] = (run.ieos[imat].gamma0-1.)*self.iev0[n+2,ioddimat]/self.vr[n+2,ioddimat]
                self.pres[n+1,ioddimat] = (self.pres[n,ioddimat]+self.pres[n+2,ioddimat])/2.
                # ------------ update pressure on interior odd j ---------------------
                self.temp[n+2,ioddimat] = self.iev0[n+2,ioddimat]/run.ieos[imat].cv
                self.temp[n+1,ioddimat] = (self.temp[n,ioddimat]+self.temp[n+2,ioddimat])/2.
            ########### SESAME TABLE EOS #################
            if (run.ieosid[imat] == 'SES'):
                # use SESAME TABLE EOS
                # dV from conservation of momentum above
                # dE = -(P+q)dV + dZ from conservation of energy
                # The EOS provides hydrostatic P-V-E relationship to solve new E,P
                #
                # Fill in staggered mesh with s1[n+1] stress deviators
                # doing it here after any yielding calculation
                self.s1[n+1,ioddimat] = 0.5*(self.s1[n,ioddimat]+self.s1[n+2,ioddimat])
                self.s2[n+1,ioddimat] = 0.5*(self.s2[n,ioddimat]+self.s2[n+2,ioddimat])
                self.s3[n+1,ioddimat] = 0.5*(self.s3[n,ioddimat]+self.s3[n+2,ioddimat])
                # delta-Z calculated for n+1 time and no other time
                self.deltaz[ioddimat] = self.vr[n+1,ioddimat]*self.deltat* \
                          (self.s1[n+1,ioddimat]*self.eps1[n+1,ioddimat] + \
                           (run.dgeom-1.)*self.s2[n+1,ioddimat]*self.eps2[n+1,ioddimat])
                # calculate averaged artificial viscosity Mbar
                qbar = 0.5*(self.q[n+1,ioddimat]+self.q[n-1,ioddimat])
                # then calculate E[n+2](test), then P[n+2](test)
                # do a convergence loop for E[n+2],P[n+2],rho[n+2]
                # loop until new E[n+2] does not change by some fraction
                # necessary loops may depend on quality/gridding of table
                # CONVERGENCE LOOP STARTS HERE
                etest = 1. # convergence test variable initialize with large value
                ecount = 0 # count the number of iterations
                            # limit the number of EOS loops
                pnew1 = np.copy(self.pres[n,ioddimat]) # start with present P[n]
                dV = (self.vr[n+2,ioddimat]-self.vr[n,ioddimat])
                tnew1 = np.copy(self.temp[n,ioddimat])
                snew1 = np.copy(self.entropy[n,ioddimat])
                while (etest > 1.e-7) and (ecount < 5):
                    # start with P[n], rho[n+2]
                    # calculate total internal energy E
                    # uses current time step pressure P[n]
                    enew = self.iev0[n,ioddimat] - (pnew1+qbar)*dV + self.deltaz[ioddimat]
                    #
                    # ------------ update pressure on interior odd j ---------------------
                    # Table EOS returns P(E,rho) on the vector of points
                    # pnew1 is an input to keep track of the points with no strain/change in state
                    pnew2,tnew2,snew2 = run.ieos[imat].sesp(self.rho0[ioddimat]/self.vr[n+2,ioddimat], \
                                               enew/self.rho0[ioddimat], \
                                               pnew1, tnew1, snew1, \
                                               self.eps1[n+1,ioddimat])
                    # check if the table pressure satisfies the conservation equation
                    enew2 = self.iev0[n,ioddimat] - (pnew1+qbar)*dV + self.deltaz[ioddimat]
                    etest = max(npabs((enew2-enew)/enew))
                    if etest > 1.e-7:
                        pnew1 = np.copy(pnew2)
                    #if np.mod(self.stepn, 100) == 0:
                    #    print('stepn, imat, ecount, etest=',self.stepn,imat,ecount,etest,max(npabs(self.deltaz[ioddimat])))
                    ecount = ecount+1
                    # CONVERGENCE LOOP ENDS HERE
                #print('out of loop: ',imat,ecount,etest)
                # put the converged values for e and p in the main data object
                self.iev0[n+2,ioddimat]=enew2
                self.pres[n+2,ioddimat]=pnew2
                self.temp[n+2,ioddimat]=tnew2
                self.entropy[n+2,ioddimat]=snew2
                #print(etest)
                # calculate energy and pressure averages E[n+1], P[n+1]
                self.iev0[n+1,ioddimat]=0.5*(self.iev0[n,ioddimat]+self.iev0[n+2,ioddimat])
                self.pres[n+1,ioddimat]=0.5*(self.pres[n,ioddimat]+self.pres[n+2,ioddimat])
                # ------------ update temperature on interior odd j ---------------------
                # Now EOS returns T(E,rho) after P-V-E convergence
                # temperature for n+2 time using updated E[n+2]
                #self.temp[n+2,ioddimat] = self.iev0[n+2,ioddimat]/run.ieos[imat].cv
                self.temp[n+1,ioddimat] = (self.temp[n,ioddimat]+self.temp[n+2,ioddimat])/2.
                # ------------ update entropy on interior odd j ---------------------
                # Now EOS returns S(E,rho) after P-V-E convergence
                # entropy for n+2 time using updated E[n+2]
                #self.entropy[n+2,ioddimat] = 0.#self.iev0[n+2,ioddimat]/run.ieos[imat].cv
                self.entropy[n+1,ioddimat] = (self.entropy[n,ioddimat]+self.entropy[n+2,ioddimat])/2.
        #       # END SESAME TABLE OPTION
            ########### TILLOTSON EOS #################
            if (run.ieosid[imat] == 'TIL'):
                # use TILLOTSON EOS; no temperature implemented at this time
                # dV from conservation of momentum above
                # dE = -(P+q)dV + dZ from conservation of energy
                # The EOS provides hydrostatic P-V-E relationship to solve new E,P
                #
                # Fill in staggered mesh with s1[n+1] stress deviators
                # doing it here after any yielding calculation
                self.s1[n+1,ioddimat] = 0.5*(self.s1[n,ioddimat]+self.s1[n+2,ioddimat])
                self.s2[n+1,ioddimat] = 0.5*(self.s2[n,ioddimat]+self.s2[n+2,ioddimat])
                self.s3[n+1,ioddimat] = 0.5*(self.s3[n,ioddimat]+self.s3[n+2,ioddimat])
                # delta-Z calculated for n+1 time and no other time
                self.deltaz[ioddimat] = self.vr[n+1,ioddimat]*self.deltat* \
                          (self.s1[n+1,ioddimat]*self.eps1[n+1,ioddimat] + \
                           (run.dgeom-1.)*self.s2[n+1,ioddimat]*self.eps2[n+1,ioddimat])
                # calculate averaged artificial viscosity Mbar
                qbar = 0.5*(self.q[n+1,ioddimat]+self.q[n-1,ioddimat])
                # then calculate E[n+2](test), then P[n+2](test)
                # do a convergence loop for E[n+2],P[n+2],rho[n+2]
                # loop until new E[n+2] does not change by some fraction
                # necessary loops may depend on quality/gridding of table
                # CONVERGENCE LOOP STARTS HERE 
                # may not be needed for Tillotson - here to verify implementation
                etest = 1. # convergence test variable initialize with large value
                ecount = 0 # count the number of iterations
                            # limit the number of EOS loops
                # calculate Tillotson pressure from new internal energy and density
                # start with a loop over the ioddimat; think about how to vectorize
                # pnewx and qbar are length of ioddimat
                pnew1 = np.copy(self.pres[n,ioddimat])                 
                pnew2 = np.copy(self.pres[n,ioddimat])                 
                dV = (self.vr[n+2,ioddimat]-self.vr[n,ioddimat])
                enew = self.iev0[n,ioddimat] - (pnew1+qbar)*dV + self.deltaz[ioddimat]
                for iii in range(len(ioddimat)):
                    # Till_P(rho,E) inputs and outputs a single point
                    pnew1[iii] = etab.Till_P(self.rho0[ioddimat[iii]]/self.vr[n+2,ioddimat[iii]], \
                                             enew[iii],run.ieos[imat].params)[0]
                while (etest > 1.e-7) and (ecount < 5):
                    # start with P[n], rho[n+2]
                    # calculate total internal energy E
                    # uses current time step pressure P[n]
                    enew = self.iev0[n,ioddimat] - (pnew1+qbar)*dV + self.deltaz[ioddimat]
                    #
                    # ------------ update pressure on interior odd j ---------------------
                    for iii in range(len(ioddimat)):
                        # Till_P(rho,E) inputs and outputs a single point
                        # returns [pout,flag,csout]
                        Till_out = etab.Till_P(self.rho0[ioddimat[iii]]/self.vr[n+2,ioddimat[iii]], \
                                                 enew[iii],run.ieos[imat].params)
                        pnew2[iii] = Till_out[0]
                        self.phase[ioddimat[iii]] = Till_out[1]
                        self.alocal[ioddimat[iii]] = Till_out[2]
                    # check if the table pressure satisfies the conservation equation
                    enew2 = self.iev0[n,ioddimat] - (pnew1+qbar)*dV + self.deltaz[ioddimat]
                    etest = max(npabs((enew2-enew)/enew))
                    if etest > 1.e-7:
                        pnew1 = np.copy(pnew2)
                    #if np.mod(self.stepn, 100) == 0:
                    #    print('stepn, imat, ecount, etest=',self.stepn,imat,ecount,etest,max(npabs(self.deltaz[ioddimat])))
                    ecount = ecount+1
                    # CONVERGENCE LOOP ENDS HERE
                #print('out of loop: ',imat,ecount,etest)
                # put the converged values for e and p in the main data object
                self.iev0[n+2,ioddimat]=enew2
                self.pres[n+2,ioddimat]=pnew2
                self.temp[n+2,ioddimat]=0.
                #print(etest)
                # calculate energy and pressure averages E[n+1], P[n+1]
                self.iev0[n+1,ioddimat]=0.5*(self.iev0[n,ioddimat]+self.iev0[n+2,ioddimat])
                self.pres[n+1,ioddimat]=0.5*(self.pres[n,ioddimat]+self.pres[n+2,ioddimat])
                # ------------ update temperature on interior odd j ---------------------
                # Now EOS returns T(E,rho) after P-V-E convergence
                # temperature for n+2 time using updated E[n+2]
                #self.temp[n+2,ioddimat] = self.iev0[n+2,ioddimat]/run.ieos[imat].cv
                self.temp[n+1,ioddimat] = 0.
                # ------------ update entropy on interior odd j ---------------------
                # Now EOS returns S(E,rho) after P-V-E convergence
                # entropy for n+2 time using updated E[n+2]
                self.entropy[n+2,ioddimat] = 0.#self.iev0[n+2,ioddimat]/run.ieos[imat].cv
                self.entropy[n+1,ioddimat] = 0.#(self.temp[n,ioddimat]+self.temp[n+2,ioddimat])/2.
        #       # END TILLOTSON EOS OPTION
        # ------------ Update the conservation sums for this time step ---------------------
        self.qtotal  = npsum(self.q[n+1,ioddj])
        self.mvtotal = npsum(0.5*self.mass[ioddj]*(self.up[n+1,ioddj-1]+self.up[n+1,ioddj+1]))
        self.ketotal = npsum(0.5*self.mass[ioddj]*np.square(0.5*(self.up[n+1,ioddj-1]+self.up[n+1,ioddj+1])))
        self.ietotal = npsum( (self.iev0[n+2,ioddj]/self.rho0[ioddj]*self.mass[ioddj]) )
        self.etotal  = npsum( (self.iev0[n+2,ioddj]-self.iev0[n,ioddj])/self.deltat + \
                          self.pres[n+1,ioddj]*(self.vr[n+2,ioddj]-self.vr[n,ioddj])/self.deltat )
        self.mtotal = npsum(self.mass[ioddj])
        #
        # ------------ check for fracture ---------------------
        # check if any cell for imat has exceeded the spall strength
        # use values at n+2 to create the fracture for the start of the next time step
        # fracture needs more testing.
        #
        for imat in np.arange(run.nmat):
            # Do not check for fractures if the material is not Von Mises material. DEBUG stsm check fracture conditions
            if (run.istrid[imat] == 'VM'):
                # select interior array points for this material
                tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
                ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
                ispall = np.where((-self.pres[n+2,ioddimat] > self.pfrac[ioddimat]) & \
                                  (self.rholocal[ioddimat] < run.ifrac[imat].nrhomin*run.ieos[imat].rhoref))[0]
                if len(ispall) > 0:
                    #print('ispall = ',self.ibc[ioddimat[ispall]])
                    for icheck in ispall:
                        if True:
                            # check if right next to an existing fracture
                            lookindex = np.append(ioddimat[icheck],[np.min(ioddimat[icheck])-6,np.min(ioddimat[icheck])-4,np.min(ioddimat[icheck])-2,np.max(ioddimat[icheck])+2, np.max(ioddimat[icheck])+4, np.max(ioddimat[icheck])+6])
                            tmp = np.where(self.ibc[lookindex] != 0)[0]
                            if (len(tmp)==0):
                                imax = icheck
                                #print('SPALL imat, ioddimat[ispall], imax = ',imat, ioddimat[ispall], ioddimat[imax])
                                #for nspall in range(len(ispall)):
                                nextibc = max(self.ibc)+1
                                #self.binarydebugoutputpint(run)
                                self.createinteriorfracture(run,imat,ioddimat[imax],nextibc)
                                #print('created a fracture')
                                #self.stepn = self.stepn+1
                                #self.binarydebugoutputpint(run)
                                #sys.exit("stop and check debug files")
        # ------------ update time step ---------------------
        # if a fracture formed then the ioddj index is wrong
        # recalculate just in case
        tmp = npwhere(self.ibc == 0)[0] # select all interior array indices
        ioddj = tmp[tmp %2 == 1] # select only odd indices for zone values
        ### Wilkins Appendix B.2 Eq 9 Time steps
        # reset the dtminj array
        self.dtminj[:] = 1.1*dthalf # make the array the maximum number
        # Needs the local sound speed cs = sqrt(Kmod/rho)
        # For ideal gas Kmod = P
        # Wilkins uses current sqrt(P/rho) for sound speed for all materials
        # calculations on all ioddj
        vdot = ( self.vr[n+2,ioddj]-self.vr[n,ioddj] )/self.deltat # change in volume
        deltar = npabs(self.pos[n+2,ioddj+1]-self.pos[n+2,ioddj-1])
        #ibpos = npwhere((vdot/self.vr[n+1,ioddj]) >= 0.)[0]
        ibneg = npwhere((vdot/self.vr[n+1,ioddj]) < 0.)[0]
        #if verbose: print('ibpos,ibneg',ibpos,ibneg)
        b = np.zeros(len(ioddj))
        b[ibneg] = 8.*(np.square(run.avc0)+run.avcl)*deltar[ibneg]*(vdot[ibneg]/self.vr[n+1,ioddj[ibneg]])
        self.rho[ioddj] = self.rho0[ioddj]/self.vr[n+2,ioddj]
        # ----------------------------- CALCULATE SOUND SPEED ------------------------
        # sound speed now with p[n+1,:]
        # a is sound speed of an ideal gas; used in KO for all materials but may need an extra reduction factor
        # user input has reduction factor option run.tstepscale
        self.alocal[ioddj] = npsqrt(npabs(self.pres[n+2,ioddj]) / self.rho[ioddj])
        # sound speed adjustments for non-ideal gases
        for imat in np.arange(run.nmat):
            # select interior array points for this material
            tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
            ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
            ########### MIE GRUNEISEN EOS #################
            if (run.ieosid[imat] == 'MGR'):
                # use c0 value if greater than a
                icheck = npwhere(self.alocal[ioddimat] < run.ieos[imat].c0)[0]
                self.alocal[ioddimat[icheck]] = run.ieos[imat].c0
            ########### SESAME EOS #################
            if (run.ieosid[imat] == 'SES'):
                # use SES value if greater than a
                # interpolate for cs in the table using rho and specific internal energy 
                sescs = run.ieos[imat].sescs(self.rho[ioddimat],self.iev0[n,ioddimat]/self.rho0[ioddimat])
#                icheck = npwhere(self.alocal[ioddimat] < sescs)[0]
#                self.alocal[ioddimat[icheck]] = sescs[icheck]
                self.alocal[ioddimat] = sescs
            ########### Tillotson EOS #################
            if (run.ieosid[imat] == 'TIL'):
                for iii in ioddimat:
                    self.alocal[iii] = etab.Till_P(self.rho0[iii]/self.vr[n+2,iii], \
                                            self.iev0[n,iii]/self.rho0[iii],run.ieos[imat].params)[2]
        a = self.alocal[ioddj]
        # ----------------------------- END SOUND SPEED CALC ------------------
        igood = npwhere(np.square(a)+np.square(b) != 0.)[0]
        #print('len igood =',len(igood))
        if len(igood)>0:
            deltat_tmp = (2./3.)*(deltar[igood]/(npsqrt(np.square(a[igood])+np.square(b[igood]))))
            self.dtminj[ioddj[igood]]=deltat_tmp
        else:
            sys.exit('FATAL ERROR: zero time step') 
        # STS finds that the factor of 6 stabilizes planar impacts with Mie-Grueneisen EOS
        # and stops oscillations at the free surfaces
        deltat_tmp = deltat_tmp/run.tstepscale # lower the time step for stability
        # the Wilkins time step is OK for the sod test
        #deltat_tmp = deltat_tmp # arbitrarily lower the time step
        if np.min(deltat_tmp) > 1.1*dthalf:
            deltat_next = 1.1*dthalf
        else:
            deltat_next = np.min(deltat_tmp)
        #deltat_next = 0.001 # test code with constant time step
        #if verbose: print("Old, new time steps: ",self.deltat,deltat_next)
        self.deltat_past=self.deltat
        self.deltat_next=deltat_next
        #
        # ------------ write current n time step to file ---------------------
        # right now output every time step
        # later add condition for time steps from run object
        if True:
            if run.step_skip != -1:
                if np.mod(self.stepn,run.step_skip) == 0:
                    if run.binoutput:
                        self.binaryoutputpint(run)
                    else:
                        self.appendoutput(run)
            if run.time_skip > 0:
                if self.time[n+2] >= run.next_time_dump:
                    # only calculate phases when saving the time step
                    for imat in np.arange(run.nmat):
                        # select interior array points for this material
                        tmp = npwhere((self.matid == imat) & (self.ibc == 0))[0]
                        ioddimat = tmp[tmp %2 == 1] # select only odd indices for zone values
                        ########### SESAME EOS #################
                        if (run.ieosid[imat] == 'SES'):
                            sesphase = run.ieos[imat].sesphase(self.rho[ioddimat],self.iev0[n,ioddimat]/self.rho0[ioddimat])
                            self.phase[ioddimat] = sesphase
                            #print('sesphase = ',sesphase)
                    if run.binoutput:
                        self.binaryoutputpint(run)
                        self.binarydebugoutputpint(run)
                    else:
                        self.appendoutput(run)
                    run.next_time_dump += run.time_skip
        #
        # ------------ write current n time step to file ---------------------
        if verbose: 
            if np.mod(self.stepn, 1000) == 0:
                print("Count, time[n] sec, dt sec = ",self.stepn,self.time[n]/1.e6,self.deltat/1.e6)
        # ------------ check for very tiny time steps ---------------------
        if self.deltat_next < run.dtmin:
            print(f'FATAL ERROR: time step below minimum: dt, dtmin = {self.deltat_next}, {run.dtmin}') 
        #
        # END OF ONE TIME STEP 
    def shiftdata(self):
        # shift time(n) array t(0),t(1),t(2),t(3) half time steps
        self.time               = np.roll(self.time,2)
        self.time[2:4]          = 0.
        # ROLL ALL t(0),t(1),t(2),t(3) -> t(2),t(3),t(1),t(2) for (n,j) arrays
        # do not reset position values; it messes up the boundary conditions
        self.pos                = np.roll(self.pos,2,axis=0)
        # Then set future time step values to zero
        self.up                 = np.roll(self.up,2,axis=0)
        self.up[2:4,:]          = 0.
        self.phi                = np.roll(self.phi,2,axis=0)
        self.phi[2:4,:]         = 0.
        self.beta               = np.roll(self.beta,2,axis=0)
        self.beta[2:4,:]        = 0.
        self.vr                 = np.roll(self.vr,2,axis=0)
        self.vr[2:4,:]          = 0.
        self.eps1               = np.roll(self.eps1,2,axis=0)
        self.eps1[2:4,:]        = 0.
        self.eps2               = np.roll(self.eps2,2,axis=0)
        self.eps2[2:4,:]        = 0.
        self.s1                 = np.roll(self.s1,2,axis=0)
        self.s1[2:4,:]          = 0.
        self.s2                 = np.roll(self.s2,2,axis=0)
        self.s2[2:4,:]          = 0.
        self.s3                 = np.roll(self.s3,2,axis=0)
        self.s3[2:4,:]          = 0.
        self.pres               = np.roll(self.pres,2,axis=0)
        self.pres[2:4,:]        = 0.
        self.q                  = np.roll(self.q,2,axis=0)
        self.q[2:4,:]           = 0.
        self.iev0               = np.roll(self.iev0,2,axis=0)
        self.iev0[2:4,:]        = 0.
        self.temp               = np.roll(self.temp,2,axis=0)
        self.temp[2:4,:]        = 0.
        self.entropy           = np.roll(self.entropy,2,axis=0)
        self.entropy[2:4,:]    = 0.
        # sigmar[j] and sigmao[j] are recalculated at the beginning
        # of the time step. They are only calculated for present n, so
        # not updated in the time shifts above.
        #
    def createinteriorboundary(self,run,imat,ninteriorbc):
        jend = max(np.where(self.matid == imat)[0]) 
        print('CREATE interior boundary after imat,jend=',imat,jend,ninteriorbc)
        # jend is odd and interior
        # add even interior right boundary at jend+1 (pos bcid)
        # add odd void with -99 matid at jend+2
        # jend+3=convert even interior node to even left interior boundary node (neg bcid)
        # PREVIOUSLY FILLED VARIABLES that need to be preserved: pos,matid,ibc
        # initializing these values here; anything that needs the input file parameters in run class
        self.jtot = self.jtot+2
        # variables in time and space (n,j)
        n = self.ntot                 # running time step index
        j = self.jtot                 # spatial domain index
        # these variables are not yet filled in makegrid
        self.up       = np.zeros((n,j)) # particle velocity
        self.vr       = np.zeros((n,j)) # relative volume to initial state=rho0/rho=1/eta in Wilkins [dimless]
        self.phi      = np.zeros((n,j)) # rho0*(pos_j+1-pos_j)/vr_j+1/2 = (pos_j+1-pos_j)/rho_j+1/2=phi in Wilkins
        self.q        = np.zeros((n,j)) # artificial viscosity Mbar
        self.eps1     = np.zeros((n,j)) # velocity strain dup/dpos, per microsec
        self.eps2     = np.zeros((n,j)) # velocity strain up/pos, per microsec
        self.beta     = np.zeros((n,j)) # (sigmar-sigmao)/(0.5 dpos)*rho = beta in Wilkins
        self.iev0     = np.zeros((n,j)) # internal energy per initial volume 1e12 erg/cm3
        self.pres     = np.zeros((n,j)) # pressure Mbar
        self.s1       = np.zeros((n,j)) # s1 stress deviator deriv = 2G(deps1-dvr/vr/3)
        self.s2       = np.zeros((n,j)) # s2 stress deviator deriv
        self.s3       = np.zeros((n,j)) # s3 stress deviator deriv; not needed in 1D; kept for equation clarity
        self.temp     = np.zeros((n,j)) # temperature K
        self.entropy  = np.zeros((n,j)) # entropy [DEBUG check eu/cm3 units]
        self.yld      = np.zeros(j) # yield stress
        self.pfrac    = np.zeros(j) # fracture stress
        self.rho0     = np.zeros(j) # initial density rho0 g/cm3; DEBUG later can free up this memory by referring to imat definition
        self.rho      = np.zeros(j) # local density g/cm3
        self.mass     = np.zeros(j) # mass in node g
        self.deltaz   = np.zeros(j) # distortion energy Mbar
        self.dtminj   = np.zeros(j) # array used to find next time step and debug EOS problems
        self.sigmar   = np.zeros(j) # total radial stress sigma_r, only evaluated at n=1 current time
        self.sigmao   = np.zeros(j) # total tangential stress sigma_theta, only evaluated at n=1 current time
        self.alocal   = np.zeros(j) # intermediate variable for AV
        self.rholocal = np.zeros(j) # intermediate variable for AV
        self.phase    = np.zeros(j) # phase information
        # jend is odd and interior
        # jend+1=convert even interior node to even left interior boundary node (neg bcid)
        # first reclassify jend +1 (before adding two new nodes)
        self.ibc[jend+1] = -ninteriorbc
        #self.matid[jend+1] stays the same = imat+1
        #self.pos(n,jend+1) stays the same = xstart for imat+1
        # add even interior right boundary at j+1 (pos bcid)
        # add odd void with -99 matid at j+2
        # insert 2 new cells after jend
        self.ibc      = np.insert(self.ibc,jend+1,[ninteriorbc,-99])
        self.matid    = np.insert(self.matid,jend+1,[imat,-1])
        # KLUDGE DEBUG: hold the free surface boundary pressure in the entropy array
        self.bcpres[jend+2] = run.ipstart[imat]
        self.bcpres[jend+4] = run.ipstart[imat+1]
        posnew = np.zeros((n,j))
        # calculate position between the two boundaries
        posvoid = (run.ixstart[imat+1]-(run.ixstart[imat]+run.ilength[imat]))/2.+(run.ixstart[imat]+run.ilength[imat])
        for i in np.arange(4):
            posnew[i,:] = np.insert(self.pos[i,:],jend+1,[run.ixstart[imat]+run.ilength[imat],posvoid])
        self.pos = deepcopy(posnew)
        print('ipstart: ',run.ipstart[imat],run.ipstart[imat+1])
        print('check inner boundary pressures: ',jend+2,jend+4,self.phase[jend+2],self.phase[jend+4])
        return
    def closeinteriorboundary(self,jend):
        # this function is called at the beginning of a time step just after shifting time variables
        # so Wilkins n+1/2 is now t(0) variable
        # jend is even and interior and the last spatial element for imat
        print('CLOSE interior boundary for node jend=',jend,self.ibc[jend])
        # determine the new velocity at the closed interface using conservation of momentum
        # Wilkins equation in appendix B4(b)
        velnew = (self.mass[jend+3]*self.up[0,jend+2] + self.mass[jend-1]*self.up[0,jend]) \
                 / (self.mass[jend+3] + self.mass[jend-1])
        print('velnew, left m, v, right m, v: ', velnew,self.mass[jend-1],self.up[0,jend],self.mass[jend+3],self.up[0,jend+2],)
        #print('')
        # jend+3=convert from left interior boundary node to odd interior node 
        self.ibc[jend] = 0 # now interior even node for imat
        self.ibc[jend+3] = 0 # now interior odd node for imat+1
        # remove even interior right boundary at jend+1 (pos bcid)
        # remove odd void with -99 matid at jend+2
        self.jtot = self.jtot-2
        # variables in time and space (n,j)
        #n = self.ntot                 # running time step index
        #j = self.jtot                 # spatial domain index
        # these variables are all length j
        delarr = [jend+1,jend+2]
        tmp           = np.delete(self.matid,delarr)
        self.matid    = tmp
        tmp           = np.delete(self.ibc,delarr)
        self.ibc      = tmp
        tmp           = np.delete(self.yld,delarr)
        self.yld      = tmp
        tmp           = np.delete(self.pfrac,delarr)
        self.pfrac    = tmp 
        tmp           = np.delete(self.rho0,delarr)
        self.rho0     = tmp
        tmp           = np.delete(self.rho,delarr)
        self.rho      = tmp
        tmp           = np.delete(self.mass,delarr)
        self.mass     = tmp
        tmp           = np.delete(self.deltaz,delarr)        
        self.deltaz   = tmp
        tmp           = np.delete(self.dtminj,delarr)        
        self.dtminj   = tmp
        tmp           = np.delete(self.sigmar,delarr)        
        self.sigmar   = tmp
        tmp           = np.delete(self.sigmao,delarr)        
        self.sigmao   = tmp
        tmp           = np.delete(self.alocal,delarr)        
        self.alocal   = tmp
        tmp           = np.delete(self.rholocal,delarr)        
        self.rholocal = tmp
        tmp           = np.delete(self.phase,delarr)        
        self.phase    = tmp
        tmp           = np.delete(self.bcpres,delarr)        
        self.bcpres   = tmp
        # these are (n,j) arrays
        tmp           = np.delete(self.pos,delarr,axis=1)     
        self.pos = tmp
        tmp           = np.delete(self.up,delarr,axis=1)     
        self.up       = tmp
        self.up[:,jend] = velnew # set new velocity at jend even node
        tmp           = np.delete(self.vr,delarr,axis=1)                 
        self.vr       = tmp
        tmp           = np.delete(self.phi,delarr,axis=1)                 
        self.phi      = tmp
        tmp           = np.delete(self.q,delarr,axis=1)                 
        self.q        = tmp
        tmp           = np.delete(self.eps1,delarr,axis=1)                 
        self.eps1     = tmp
        tmp           = np.delete(self.eps2,delarr,axis=1)                 
        self.eps2     = tmp
        tmp           = np.delete(self.beta,delarr,axis=1)                 
        self.beta     = tmp
        tmp           = np.delete(self.iev0,delarr,axis=1)                 
        self.iev0     = tmp
        tmp           = np.delete(self.pres,delarr,axis=1)                 
        self.pres     = tmp
        tmp           = np.delete(self.s1,delarr,axis=1)                 
        self.s1       = tmp
        tmp           = np.delete(self.s2,delarr,axis=1)                 
        self.s2       = tmp
        tmp           = np.delete(self.s3,delarr,axis=1)                 
        self.s3       = tmp
        tmp           = np.delete(self.temp,delarr,axis=1)                 
        self.temp     = tmp
        tmp           = np.delete(self.entropy,delarr,axis=1)                 
        self.entropy  = tmp    
        return
    def checkforcontact(self,run):
        # following Wilkins Section B2(b) closing void
        # need to check for adjacent voids
        # DEBUG: need to sort out the half time step value delta-t^n
        n=1
        ivoid = np.where(self.matid == -1)[0]
        if len(ivoid) > 0:
            # there are interior boundaries
            for jvoid in ivoid:
                jend = jvoid-1 # should be even inner boundary of the void
                # if a fracture was just created there is no space
                if self.pos[n,jend] < self.pos[n,jend+2]:
                    # advance the even nodes on either side of the gap by 1.2*this time step
                    testposleft = self.pos[n,jend] + self.up[n-1,jend] * self.deltat*1.2
                    testposright = self.pos[n,jend+2] + self.up[n-1,jend+2] * self.deltat*1.2
                    #print(self.ibc[jend-2:jend+2])
                    #print(testposleft,testposright)
                    #print(self.pos[n,jend],self.up[n-1,jend])
                    #print(self.pos[n,jend+2],self.up[n-1,jend+2])
                    if testposleft >= testposright:
                        # very close; now solve for contact timestep
                        dxgap = self.pos[n,jend+2]-self.pos[n,jend]
                        dupgap = self.up[n-1,jend]-self.up[n-1,jend+2]
                        print('dxgap,dupgap=',dxgap,dupgap)
                        #A = self.sigmar[jend+3]/self.phi[n,jend+2] + self.sigmar[jend-1]/self.phi[n,jend] + \
                        #    (self.beta[n,jend+2]+self.beta[n,jend])*(run.dgeom-1)
                        #B = 2.*dupgap + A*self.deltat_past
                        # determine necessary deltat for contact
                        # DEBUG - THIS IS NOT THE WILKINS LOOP COME BACK TO THIS CALCULATION
                        # need to implement a convergence loop
                        dtgap = dxgap/dupgap # try to void overlap
                        print('WARNING: contact check does not include stress corrections! ibc=',self.ibc[jend])
                        print('old deltat, new dtgap [code time unit] = ',self.deltat,dtgap)
                        #print('')
                        # need to set the time step
                        self.deltat_past = self.deltat
                        self.deltat_next = dtgap
                        self.closeinteriorboundary(jend)
                        print('Contact! reduce the time step and remove the extra nodes ',self.time[1],self.deltat,dtgap)
        #print('')
        return
    def createinteriorfracture(self,run,imat,ispall,ninteriorbc):
        #jend = max(np.where(self.matid == imat)[0])
        jend = ispall
        print('CREATED new interior fracture after imat,jend,time=',imat,jend,ninteriorbc,self.time[3])
        #print('jend ibc',jend,self.ibc[jend-4:jend+4])
        #print('jend mat',jend,self.matid[jend-4:jend+4])
        # jend is odd and interior
        # add even interior right boundary at jend+1 (pos bcid)
        # add odd void with -99 matid at jend+2
        # jend+3=convert even interior node to even left interior boundary node (neg bcid)
        # PREVIOUSLY FILLED VARIABLES that need to be preserved: pos,matid,ibc
        # initializing these values here; anything that needs the input file parameters in run class
        self.jtot = self.jtot+2
        # variables in time and space (n,j)
        #n = self.ntot                 # running time step index
        j = self.jtot                 # spatial domain index
        # insert 2 cells in the j arrays
        n=1
        # arrays that need some thinking about values
        self.up       = np.insert(self.up,jend+1,self.up[:,jend+1],axis=1)
        self.up       = np.insert(self.up,jend+2,[0.,0.,0.,0.],axis=1)
        self.bcpres   = np.insert(self.bcpres,jend+1,[0.,0.])
        # arrays with only zeros added
        # these are (n,j) arrays
        self.vr       = np.insert(self.vr,jend+1,self.vr[:,jend+1],axis=1)
        self.vr       = np.insert(self.vr,jend+2,[0.,0.,0.,0.],axis=1)
        self.phi      = np.insert(self.phi,jend+1,self.phi[:,jend+1],axis=1)   
        self.phi      = np.insert(self.phi,jend+2,[0.,0.,0.,0.],axis=1)   
        self.q        = np.insert(self.q,jend+1,self.q[:,jend+1],axis=1)     
        self.q        = np.insert(self.q,jend+2,[0.,0.,0.,0.],axis=1)     
        self.eps1     = np.insert(self.eps1,jend+1,self.eps1[:,jend+1],axis=1)     
        self.eps1     = np.insert(self.eps1,jend+2,[0.,0.,0.,0.],axis=1)     
        self.eps2     = np.insert(self.eps2,jend+1,self.eps2[:,jend+1],axis=1)   
        self.eps2     = np.insert(self.eps2,jend+2,[0.,0.,0.,0.],axis=1)   
        self.beta     = np.insert(self.beta,jend+1,self.beta[:,jend+1],axis=1)    
        self.beta     = np.insert(self.beta,jend+2,[0.,0.,0.,0.],axis=1)    
        self.iev0     = np.insert(self.iev0,jend+1,self.iev0[:,jend+1],axis=1)   
        self.iev0     = np.insert(self.iev0,jend+2,[0.,0.,0.,0.],axis=1)   
        self.pres     = np.insert(self.pres,jend+1,self.pres[:,jend+1],axis=1)
        self.pres     = np.insert(self.pres,jend+2,[0.,0.,0.,0.],axis=1)
        self.s1       = np.insert(self.s1,jend+1,self.s1[:,jend+1],axis=1) 
        self.s1       = np.insert(self.s1,jend+2,[0.,0.,0.,0.],axis=1) 
        self.s2       = np.insert(self.s2,jend+1,self.s2[:,jend+1],axis=1) 
        self.s2       = np.insert(self.s2,jend+2,[0.,0.,0.,0.],axis=1) 
        self.s3       = np.insert(self.s3,jend+1,self.s3[:,jend+1],axis=1) 
        self.s3       = np.insert(self.s3,jend+2,[0.,0.,0.,0.],axis=1) 
        self.temp     = np.insert(self.temp,jend+1,self.temp[:,jend+1],axis=1)
        self.temp     = np.insert(self.temp,jend+2,[0.,0.,0.,0.],axis=1)
        self.entropy  = np.insert(self.entropy,jend+1,self.entropy[:,jend+1],axis=1)    
        self.entropy  = np.insert(self.entropy,jend+2,[0.,0.,0.,0.],axis=1)    
        # (j) arrays
        newarr = [0.,0.]
        self.yld      = np.insert(self.yld,jend+1,[self.yld[jend+1],0.])
        self.pfrac    = np.insert(self.pfrac,jend+1,[self.pfrac[jend+1],0.])
        self.rho0     = np.insert(self.rho0,jend+1,[self.rho0[jend+1],1.])
        self.rho      = np.insert(self.rho,jend+1,[self.rho[jend+1],0.])
        self.mass     = np.insert(self.mass,jend+1,[self.mass[jend+1],0.])
        self.deltaz   = np.insert(self.deltaz,jend+1,[self.deltaz[jend+1],0.])
        self.dtminj   = np.insert(self.dtminj,jend+1,[self.dtminj[jend+1],0.])
        self.sigmar   = np.insert(self.sigmar,jend+1,[self.sigmar[jend+1],0.])
        self.sigmao   = np.insert(self.sigmao,jend+1,[self.sigmao[jend+1],0.])
        self.alocal   = np.insert(self.alocal,jend+1,[self.alocal[jend+1],0.])
        self.rholocal = np.insert(self.rholocal,jend+1,[self.rholocal[jend+1],0.])
        self.phase    = np.insert(self.phase,jend+1,[self.phase[jend+1],0.])
        # more arrays that need thinking about values
        # jend is odd and interior
        # jend+1=convert even interior node to even left interior boundary node (neg bcid)
        # first reclassify jend +1 (before adding two new nodes)
        self.ibc[jend+1] = -ninteriorbc
        #self.matid[jend+1] stays the same = imat+1
        #self.pos(n,jend+1) stays the same = xstart for imat+1
        # add even interior right boundary at j+1 (pos bcid)
        # add odd void with -99 matid at j+2
        # insert 2 new cells after jend
        self.ibc      = np.insert(self.ibc,jend+1,[ninteriorbc,-99]) 
        self.matid    = np.insert(self.matid,jend+1,[imat,-1])
        # update positions
        # 3 array locations have the same position original pos[:,jend+1]
        # original self.pos[:,jend+1] is the outer fracture boundary = inner fracture boundary
        posbc = [self.pos[0,jend+1],self.pos[1,jend+1],self.pos[2,jend+1],self.pos[3,jend+1]]
        self.pos = np.insert(self.pos,jend+1,posbc,axis=1) # inner fracture boundary
        self.pos = np.insert(self.pos,jend+2,posbc,axis=1) # location of the void
        #print('check inner boundary pressures: ',jend,jend+3,self.bcpres[jend],self.bcpres[jend+3])
        #print('jend ibc',jend,self.ibc[jend-4:jend+4])
        #print('jend mat',jend,self.matid[jend-4:jend+4])

        return
    def writefirstoutput(self,run,checkexistsflag=True,verbose=False,overwriteflag=False):
        """ Tests for existing file and creates unique new file for output.
            Writes a header for the output columns.
            """
        iii=1
        testname = run.outputfilename
        testnamedebug = testname+'_debug'
        if checkexistsflag:
            while exists(testname):
                testname = run.outputfilename+str(iii)
                iii=iii+1
                if verbose: print('WARNING: Output file exists; creating new file: '+testname)
        else:
            if exists(testname):
                if verbose: print('WARNING: Overwriting output file: '+testname)
                if not(overwriteflag):
                    input(f"rm -f {testname}. Press Enter to continue...")
                system(f'rm -f {testname}')
            if exists(testnamedebug):
                system(f'rm -f {testnamedebug}')
        run.outputfilename=testname
        if verbose: print('Output filename = '+run.outputfilename)
        if run.debugflag:
            if verbose: print('Output debug filename = '+testnamedebug)
        # write the 0th time step initial condition: binary or ascii
        if run.binoutput:
            # pickle binary output
            self.binaryoutputpint(run,verbose=verbose)
            # optional debug output for full domain arrays
            #self.binarydebugoutputpint(run)
        else:
            # ascii output
            with open(run.outputfilename,"w") as f: 
                # headers for appendoutput
                f.write('step j ibc mat time r(0) r(1) r(2) pos '+ \
                        'up vr rho rho0 iev0 pres sigmar s1 s2 q dtminj phi beta eps1 eps2 deltaz aj rhoj sigmao temp etot\n')
            self.appendoutput(run)
        #
    def appendoutput(self,run,verbose=False):
        """ Appends data to ascii output file. NO UNITS CONVERSION.
            Primarily used for debugging and comparison to fortran KO."""
        n=1
        if verbose: print('Writing time ',self.time[n+2],' to file ',run.outputfilename)   
        with open(run.outputfilename,"a") as f: 
            for j in np.arange(0,self.jtot-1,2): # loop on even nodes
                f.write('%4i %4i %4i %4i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n' % \
                        (self.stepn,j,self.ibc[j],\
                         self.matid[j+1]+1,self.time[n+2],self.pos[n-1,j],self.pos[n,j],\
                         self.pos[n+1,j],self.pos[n+2,j],self.up[n+1,j],\
                         self.vr[n+2,j+1],self.rho[j+1],self.rho0[j+1],self.iev0[n+2,j+1],\
                         self.pres[n+2,j+1],self.sigmar[j+1],self.s1[n+2,j+1],\
                         self.s2[n+2,j+1],self.q[n+1,j+1],self.dtminj[j+1],\
                         self.phi[n,j],self.beta[n,j],self.eps1[n+1,j+1],self.eps2[n+1,j+1],\
                         self.deltaz[j+1],self.alocal[j+1],self.rholocal[j+1],\
                         self.sigmao[j+1],self.temp[n+1,j+1],self.ietotal+self.ketotal))
            # write the last boundary condition cells to be able to check it
            j=self.jtot-1
            f.write('%5i %4i %4i %4i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n' % \
                        (self.stepn,j,self.ibc[j],\
                         self.matid[j]+1,self.time[n+2],self.pos[n-1,j],self.pos[n,j],\
                         self.pos[n+1,j],self.pos[n+2,j],self.up[n+1,j],\
                         self.vr[n+2,j],self.rho[j],self.rho0[j],self.iev0[n+2,j],\
                         self.pres[n+2,j],self.sigmar[j],self.s1[n+2,j],\
                         self.s2[n+2,j],self.q[n+1,j],self.dtminj[j],\
                         self.phi[n,j],self.beta[n,j],self.eps1[n+1,j],self.eps2[n+1,j],\
                         self.deltaz[j],self.alocal[j],self.rholocal[j],\
                         self.sigmao[j],self.temp[n+1,j],self.ietotal+self.ketotal))
            #
        # update output log in the run object
        run.outputsteps = np.append(run.outputsteps,[self.stepn])
        run.outputtimes = np.append(run.outputtimes,[self.time[n+2]])
        # update the conservation values in the run object
        run.outputietot  = np.append(run.outputietot,[self.ietotal])
        run.outputketot  = np.append(run.outputketot,[self.ketotal])
        run.outputmvtot  = np.append(run.outputmvtot,[self.mvtotal])
        #
    def binaryoutputpint(self,run,verbose=False):
        """ Create time snapshot of cell centered values ready for pandas dataframe. 
            The output class is in the units defined by the user input in the yaml configuration file.
        """
        # set up special code units for pint
        ureg.define('eu = 1.0E12 ergs')   # energy unit
        #
        with open(run.inputfilename, 'r') as f:
            config = yaml.safe_load(f)
        # print the cell thermo parameters (ioddj between nodes)
        # print the corresponding averaged velocity and positions for these cell thermo parameters
        tmp = npwhere(self.ibc == 0)[0] # select all interior array indices
        ioddj = tmp[tmp %2 == 1] # select only odd indices for zone values
        # OUTPUT CLASS HOLDS
        dataout        = OutputClass()
        n=1
        dataout.j      = Q_(np.asarray(ioddj),'dimensionless')
        dataout.stepn  = Q_(np.full(len(ioddj),self.stepn),'dimensionless')
        timearr = Q_(np.full(len(ioddj),self.time[n+2]),config['codeunits']['time'])
        dataout.time   = timearr.to(config['units']['time'])
        dataout.mat    = Q_(np.asarray(self.matid[ioddj].astype(int)+1),'dimensionless')
        posarr = Q_(np.asarray(0.5*(self.pos[n+2,ioddj-1]+self.pos[n+2,ioddj+1])),config['codeunits']['length'])
        dataout.pos    = posarr.to(config['units']['length'])
        uparr = Q_(np.asarray(0.5*(self.up[n+1,ioddj-1]+self.up[n+1,ioddj+1])),config['codeunits']['velocity'])
        dataout.up     = uparr.to(config['units']['velocity'])
        rho0arr = Q_(np.asarray(self.rho0[ioddj]),config['codeunits']['density'])
        dataout.rho0   = rho0arr.to(config['units']['density'])
        rhoarr = Q_(np.asarray(self.rho[ioddj]),config['codeunits']['density'])
        dataout.rho    = rhoarr.to(config['units']['density'])
        ie = Q_(np.asarray(self.iev0[n+2,ioddj]/self.rho0[ioddj]),config['codeunits']['sp_energy'])
        dataout.ie     = ie.to(config['units']['sp_energy'])
        ent = Q_(np.asarray(self.entropy[n+2,ioddj]),config['codeunits']['sp_entropy'])
        dataout.entropy     = ent.to(config['units']['sp_entropy'])
        parr = Q_(np.asarray(self.pres[n+2,ioddj]),config['codeunits']['pressure'])
        dataout.pres   = parr.to(config['units']['pressure'])
        sigmararr = Q_(np.asarray(self.sigmar[ioddj]),config['codeunits']['pressure'])
        dataout.sigmar = sigmararr.to(config['units']['pressure'])
        sigmaoarr = Q_(np.asarray(self.sigmao[ioddj]),config['codeunits']['pressure'])
        dataout.sigmao = sigmaoarr.to(config['units']['pressure'])
        temparr = Q_(np.asarray(self.temp[n+2,ioddj]),config['codeunits']['temperature'])
        dataout.temp   = temparr.to(config['units']['temperature'])
        etotarr = Q_(np.full(len(ioddj),self.ietotal+self.ketotal),config['codeunits']['sp_energy'])
        dataout.etot   = etotarr.to(config['units']['sp_energy'])
        massarr = Q_(np.asarray(self.mass[ioddj]),config['codeunits']['mass'])
        dataout.mass   = massarr.to(config['units']['mass'])
        alocalarr = Q_(np.asarray(self.alocal[ioddj]),config['codeunits']['velocity'])
        dataout.alocal   = alocalarr.to(config['units']['velocity'])
        dataout.phase   = Q_(np.asarray(self.phase[ioddj]),'dimensionless')
        dtminjarr = Q_(np.asarray(self.dtminj[ioddj]),config['codeunits']['time'])
        dataout.dtminj  = dtminjarr.to(config['units']['time'])
        with open(run.outputfilename,"ab") as f: 
            #print('dumping dataout class with pickle ',self.stepn)
            pickle.dump(dataout,f)
        # update output log in the run object
        run.outputsteps = np.append(run.outputsteps,[self.stepn])
        run.outputtimes = np.append(run.outputtimes,[self.time[n+2]])
        # update the conservation values in the run object
        run.outputietot  = np.append(run.outputietot,[self.ietotal])
        run.outputketot  = np.append(run.outputketot,[self.ketotal])
        run.outputmvtot  = np.append(run.outputmvtot,[self.mvtotal])
        #
    def __str__(self):
       """ Print problem domain variables """
       return  'future scripts to print data object to screen and file'
    def binarydebugoutputpint(self,run,verbose=False):
        """ Create time snapshot of all array values ready for pandas dataframe. 
            The output class is in the units defined by the user input in the yaml configuration file.
        """
        # set up special code units for pint
        ureg.define('eu = 1.0E12 ergs')   # energy unit
        #
        with open(run.inputfilename, 'r') as f:
            config = yaml.safe_load(f)
        # print the cell thermo parameters (ioddj between nodes)
        # print the corresponding averaged velocity and positions for these cell thermo parameters
        #tmp = npwhere(self.ibc == 0)[0] # select all interior array indices
        #ioddj = tmp[tmp %2 == 1] # select only odd indices for zone values
        # DEBUG OUTPUT
        dataout        = DebugClass()
        nj = len(self.ibc)
        n=1
        dataout.j      = Q_(np.arange(nj),'dimensionless')
        dataout.stepn  = Q_(np.full(nj,self.stepn),'dimensionless')
        timearr = Q_(np.full(nj,self.time[n+2]),config['codeunits']['time'])
        dataout.time   = timearr.to(config['units']['time'])
        dataout.mat    = Q_(np.asarray(self.matid.astype(int)+1),'dimensionless')
        # position and velocity evaluated at n+1
        posarr = Q_(np.asarray(self.pos[n+1,:]),config['codeunits']['length'])
        dataout.pos    = posarr.to(config['units']['length'])
        uparr = Q_(np.asarray(self.up[n+1,:]),config['codeunits']['velocity'])
        dataout.up     = uparr.to(config['units']['velocity'])
        dataout.vr      = Q_(np.asarray(self.vr[n+1,:]),'dimensionless')
        # thermo variables evaluated at n+2
        # do not divide by zero if there is a void space
        imass = np.where(self.rho0 == 0.)[0]
        if len(imass) > 0:
            rhotmp = deepcopy(self.rho0)
            rhotmp[imass] = 1.
            ie = Q_(np.asarray(self.iev0[n+2,:]/rhotmp),config['codeunits']['sp_energy'])
        else:
            ie = Q_(np.asarray(self.iev0[n+2,:]/self.rho0),config['codeunits']['sp_energy'])
        dataout.ie     = ie.to(config['units']['sp_energy'])
        #
        parr = Q_(np.asarray(self.pres[n+2,:]),config['codeunits']['pressure'])
        dataout.pres   = parr.to(config['units']['pressure'])
        temparr = Q_(np.asarray(self.temp[n+2,:]),config['codeunits']['temperature'])
        dataout.temp   = temparr.to(config['units']['temperature'])
        tmparr = Q_(np.asarray(self.s1[n+2,:]),config['codeunits']['pressure'])
        dataout.s1   = tmparr.to(config['units']['pressure'])
        tmparr = Q_(np.asarray(self.s2[n+2,:]),config['codeunits']['pressure'])
        dataout.s2   = tmparr.to(config['units']['pressure'])
        tmparr = Q_(np.asarray(self.s3[n+2,:]),config['codeunits']['pressure'])
        dataout.s3   = tmparr.to(config['units']['pressure'])
        # variables in j only
        rho0arr = Q_(np.asarray(self.rho0),config['codeunits']['density'])
        dataout.rho0   = rho0arr.to(config['units']['density'])
        rhoarr = Q_(np.asarray(self.rho),config['codeunits']['density'])
        dataout.rho    = rhoarr.to(config['units']['density'])
        sigmararr = Q_(np.asarray(self.sigmar),config['codeunits']['pressure'])
        dataout.sigmar = sigmararr.to(config['units']['pressure'])
        sigmaoarr = Q_(np.asarray(self.sigmao),config['codeunits']['pressure'])
        dataout.sigmao = sigmaoarr.to(config['units']['pressure'])
        etotarr = Q_(np.full(nj,self.ietotal+self.ketotal),config['codeunits']['sp_energy'])
        dataout.etot   = etotarr.to(config['units']['sp_energy'])
        massarr = Q_(np.asarray(self.mass),config['codeunits']['mass'])
        dataout.mass   = massarr.to(config['units']['mass'])
        alocalarr = Q_(np.asarray(self.alocal),config['codeunits']['velocity'])
        dataout.alocal   = alocalarr.to(config['units']['velocity'])
        dataout.phase   = Q_(np.asarray(self.phase),'dimensionless')
        dataout.ibc     = Q_(np.asarray(self.ibc),'dimensionless')
        tmparr = Q_(np.asarray(self.bcpres),config['codeunits']['pressure'])
        dataout.bcpres   = tmparr.to(config['units']['pressure'])
        tmparr = Q_(np.asarray(self.pfrac),config['codeunits']['pressure'])
        dataout.pfrac   = tmparr.to(config['units']['pressure'])
        tmparr = Q_(np.asarray(self.yld),config['codeunits']['pressure'])
        dataout.yld   = tmparr.to(config['units']['pressure'])
        #
        if run.debugflag:
            debugfilename = run.outputfilename+'_debug'
            with open(debugfilename,"ab") as f: 
                #print('dumping dataout class with pickle ',self.stepn)
                pickle.dump(dataout,f)
        #
    def __str__(self):
       """ Print problem domain variables """
       return  'future scripts to print data object to screen and file'
##
#
############################################################## 
# INPUT FILE FUNCTIONS
# readinput_borg reads Borg-style I/O for fortran KO
# readinput_yaml reads new pyKO configuration file
############################################################## 
#
def readinput_borg(run,verbose=False):
    """############################################################## 
    # Read in John Borg KO style input file for KOv11.f
    # STS added input for geometry 
    # Format: 
    #   2 lines of comments for material layers
    #     one line for each material layer 
    #     required blank line to denote end of material layers
    #   2 lines of comments for boundary conditions
    #     one line for each boundary 
    #     required blank line to denote end of boundary conditions
    #   1 line of comment for stop time
    #     tstop in microseconds
    #   1 line of comment for geometry
    #     d = 1, 2, or 3 for planar, cylindrical, or spherical
    #   1 line of comment for aneos tables
    #     directory for aneos tables
    #   end of input line 'Do not erase this statement'
    #   any more lines are ignored
    # Borg EOS definitions:
    #   EOS provides the internal energy at V and P
    #   iEOS(2,imat)=1 - Mie Gruenisen
    #   iEOS(2,imat)=2 - Gamma law ideal gas
    #   iEOS(2,imat)=3 - Gamma law ideal gas with Newtonian Stress and heat transfer
    #   iEOS(2,imat)=4 - snow plow model with KO inputs for Hugoniot
    #   iEOS(2,imat)=5 - snow plow model with anamolous hugoniot
    #   iEOS(2,imat)=6 - P-alpha Model
    #
    # Here use EOS=7 for ANEOS SESAME tables for fortran compatibility 
    # if any EOS are ANEOS tables, the required files must be in 
    # single directory path for each EOS=7 in order
    #"""
    if exists(run.inputfilename):
        with open(run.inputfilename) as f: 
            inputlines = f.readlines() # read in all of input file at once
    else:
        sys.exit("FATAL ERROR: Cannot find input file: "+run.inputfilename)
    #
    # assign variables into run setup data class
    # count number of materials
    run.nmat=0 # start counter for the number of material layers in the input file
    for iline in inputlines[2::]:
        #print(run.nmat,bool(re.search(r'^\s*$', iline)))
        blankline = bool(re.search(r'^\s*$', iline)) # should match any white space in this blank line
        if blankline:
            break
        else:
            run.nmat=run.nmat+1
    #
    # read in descriptors for each layer
    for imat in range(run.nmat):
        matarr         = np.fromstring(inputlines[imat+2],sep=' ') # float numpy array
        run.inodes     = np.append(run.inodes,[int(matarr[1])])
        run.ilength    = np.append(run.ilength,[matarr[2]])
        run.ixstart    = np.append(run.ixstart,[matarr[3]])
        run.ipstart    = np.append(run.ipstart,[matarr[4]])
        run.iupstart   = np.append(run.iupstart,[matarr[5]])
        run.irhostart  = np.append(run.irhostart,[matarr[6]])
        run.iiev0start = np.append(run.iiev0start,[matarr[7]])
        run.itempstart = np.append(run.itempstart,[matarr[7]]/matarr[16])
        if matarr[0] == 1:
            run.ieosid.append('MGR')
            matnew = MGRClass()
            run.ieos.append(matnew)
            run.ieos[imat].name='mat'+str(imat+1)
            run.ieos[imat].p0=matarr[4] # Mbar
            run.ieos[imat].up0=matarr[5] # cm/us
            run.ieos[imat].rho0=matarr[6] # g/cm3
            run.ieos[imat].iev0=matarr[7] # internal energy per original volume, 1e12 ergs/cm3 
            run.ieos[imat].rhoref=matarr[8] # g/cm3
            run.ieos[imat].c0=matarr[9] # cm/us
            run.ieos[imat].s1=matarr[10] # dimless
            run.ieos[imat].s2=matarr[11] # us/cm
            run.ieos[imat].gamma0=matarr[12] # dimless
            run.ieos[imat].cv=matarr[16] # Mbar cm3 / (K g)
            run.ieos[imat].calccoefs() # calculate k1-k4
        if matarr[0] == 2:
            run.ieosid.append('IDG')
            matnew = IDGClass()
            run.ieos.append(matnew)
            run.ieos[imat].name='mat'+str(imat+1)
            run.ieos[imat].p0=matarr[4] # Mbar
            run.ieos[imat].up0=matarr[5] # cm/us
            run.ieos[imat].rho0=matarr[6] # g/cm3
            run.ieos[imat].iev0=matarr[7] # internal energy per original volume, 1e12 ergs/cm3 
            run.ieos[imat].gamma0=matarr[12] # dimless
            run.ieos[imat].cv=matarr[16] # Mbar cm3 / (K g)
            run.ieos[imat].t0 = matarr[7]/matarr[16] # K
        if (matarr[0] == 7):
            run.ieosid.append('SES')
            matnew = SESClass()
            run.ieos.append(matnew)
            run.ieos[imat].name='mat'+str(imat+1)
            run.ieos[imat].p0=matarr[4] # Mbar
            run.ieos[imat].up0=matarr[5] # cm/us
            run.ieos[imat].rho0=matarr[6] # g/cm3
            run.ieos[imat].iev0=matarr[7] # internal energy per original volume, 1e12 ergs/cm3 
            run.ieos[imat].gamma0=matarr[12] # dimless
            run.ieos[imat].cv=matarr[16] # Mbar cm3 / (K g)
            run.ieos[imat].t0 = matarr[7]/matarr[16] # K
            #run.ieos[imat].EOS.MODELNAME = 'test'        
        # read hydro parametrs
        if matarr[14] == 0.:
            run.istrid.append('HYDRO')
            strnew = HydroClass() # hydro with fracture
            run.istr.append(strnew)
            run.istr[imat].name='mat'+str(imat+1)
        else: # for now only strength model is VM; fracture not coded yet
            run.istrid.append('VM')
            strnew = VonMisesClass() # von Mises strength and fracture
            run.istr.append(strnew)
            run.istr[imat].name='mat'+str(imat+1)
            run.istr[imat].ys=matarr[13] # Mbar
            run.istr[imat].gmod=matarr[14] # Mbar
    #
    # read in boundary conditions
    # count number of boundary conditions
    run.nbc=0 # start counter for the number of boundary conditions in the input file
    istart = 2+run.nmat+3
    for iline in inputlines[istart::]:
        #print(run.nbc,bool(re.search(r'^\s*$', iline)))
        blankline = bool(re.search(r'^\s*$', iline)) # should match any white space in this blank line
        if blankline:
            break
        else:
            run.nbc=run.nbc+1
    #
    # read in descriptors for each boundary condition
    for ibc in range(run.nbc):
        bcarr = np.fromstring(inputlines[istart+ibc],sep=' ') # float numpy array
        if abs(bcarr[0]) == 1:
            # free surface boundary condition
            run.ibcid.append('FREE')
            bcnew = BCClass()
            run.ibc.append(bcnew)
            run.ibc[ibc].name = 'Free surface '+str(ibc+1)
            run.ibc[ibc].pbc = bcarr[1] # Mbar
            run.ibc[ibc].upbc = bcarr[2] # cm/us
            run.ibc[ibc].rrefbc = bcarr[3] # g/cm3
            run.ibc[ibc].iebc = bcarr[4] # 1e12 ergs/cm3
            run.ibc[ibc].vrbc = bcarr[5] # dimless relative volume to rhoref
        if abs(bcarr[0]) == 4:
            # free surface boundary condition
            run.ibcid.append('FIXED')
            bcnew = BCClass()
            run.ibc.append(bcnew)
            run.ibc[ibc].name = 'Fixed surface '+str(ibc+1)
            run.ibc[ibc].pbc = bcarr[1] # Mbar
            run.ibc[ibc].upbc = bcarr[2] # cm/us
            run.ibc[ibc].rrefbc = bcarr[3] # g/cm3
            run.ibc[ibc].iebc = bcarr[4] # 1e12 ergs/cm3
            run.ibc[ibc].vrbc = bcarr[5] # dimless relative volume to rhoref
    #
    # read in stop time
    run.tstop = float(inputlines[2+run.nmat+3+run.nbc+2])
    # read in geometry d=1., 2., or 3. Must be a float
    run.dgeom = float(inputlines[2+run.nmat+3+run.nbc+2+2])
    if run.dgeom == 1.0:
        run.dflag = 'PLA'
    elif run.dgeom == 2.0:
        run.dflag = 'CYL'
    elif run.dgeom == 3.0:
        run.dflag = 'SPH'
    else:
        sys.exit("FATAL ERROR: invalid geometry d="+str(run.dgeom))
    #
    if (run.ixstart[0] < 0) & (run.dflag != 'PLA'):
        # cannot have negative r values for CYL or SPH
        sys.exit('FATAL ERROR: invalid geometry, negative xstart with '+run.dflag)
    #            
    # if using SESAME tables, load the tables
    sescount = 0
    for imat in range(run.nmat):
        if run.ieosid[imat] == 'SES':
            run.ieos[imat].path = inputlines[2+run.nmat+3+run.nbc+2+2+2+sescount][0:-1]
            run.ieos[imat].EOS.R0REF = run.ieos[imat].rho0
            run.ieos[imat].readsesext(verbose=verbose) # Borg style input currently hard-codes unit conversion
            if verbose: print('imat = ',imat, ' reading SESAME STD+EXT tables in ',run.ieos[imat].path)
            sescount=sescount+1 # line number for next table path
    #
    if verbose:
        print(run) # print pyKO model run parameters
        for imat in range(run.nmat):
            print(run.ieos[imat]) # print material EOS parameters
            print(run.istr[imat]) # print material strength parameters
        for ibc in range(run.nbc):
            print(run.ibc[ibc]) # print boundary condition parameters
    return
    #
    # END read in John Borg KO style input file for KOv11.f
    ############################################################## 
    #
def readinput_yaml(run,verbose=False):
    """############################################################## 
       # Read in yaml configuration file for pyKO
       ##############################################################"""
#    import pint
#    ureg = pint.UnitRegistry()
#    Q_ = ureg.Quantity
    #
    if exists(run.inputfilename):
        with open(run.inputfilename, 'r') as f:
            config = yaml.safe_load(f)
    else:
        sys.exit("FATAL ERROR: Cannot find input file: "+run.inputfilename)
    #
    
    # set up special code units for pint
    ureg.define('eu = 1.0E12 ergs')   # energy unit
    #
    # assign variables into run setup data class
    #
    # count number of materials
    # count the number of materials
    keys = list(config)
    # start counter for the number of material layers in the input file
    run.nmat = 0
    matlist = []
    for key in keys:
        if key.startswith('mat'):
            matlist.append(key)
            run.nmat += 1
        if key.startswith('outputfilename'):
            run.outputfilename = config['outputfilename']
    #
    # read in descriptors for each material layer
    for imat in range(run.nmat):
        # required inputs for each layer
        # number of cells x 2 = number of nodes
        run.inodes     = np.append(run.inodes,int(2*config[matlist[imat]]['mesh']['cells']))
        # xstart
        xstart = Q_(np.asarray(config[matlist[imat]]['mesh']['xstart'],dtype='float'),config['units']['length'])
        run.ixstart    = np.append(run.ixstart,xstart.to(config['codeunits']['length']).magnitude)
        # length
        length = Q_(np.asarray(config[matlist[imat]]['mesh']['length'],dtype='float'),config['units']['length'])
        run.ilength    = np.append(run.ilength,length.to(config['codeunits']['length']).magnitude)
        # p0
        p0 = Q_(np.asarray(config[matlist[imat]]['init']['p0'],dtype='float'),config['units']['pressure'])
        run.ipstart    = np.append(run.ipstart,p0.to(config['codeunits']['pressure']).magnitude)
        # up0
        up0 = Q_(np.asarray(config[matlist[imat]]['init']['up0'],dtype='float'),config['units']['velocity'])
        run.iupstart   = np.append(run.iupstart,up0.to(config['codeunits']['velocity']).magnitude)
        # rho0
        rho0 = Q_(np.asarray(config[matlist[imat]]['init']['rho0'],dtype='float'),config['units']['density'])
        run.irhostart  = np.append(run.irhostart,rho0.to(config['codeunits']['density']).magnitude)
        # internal energy per initial volume
        # convert input specific energy to code units and then additional normalization by rho0 in code units
        # KO code units for iev0 =  eu/cm3
        spenergy = Q_(np.asarray(config[matlist[imat]]['init']['e0'],dtype='float'),config['units']['sp_energy'])
        iev0 = spenergy.to(config['codeunits']['sp_energy'])*rho0.to(config['codeunits']['density'])
        run.iiev0start = np.append(run.iiev0start,iev0.magnitude)
        # T0 --> move to material definition because population of this variable depends on the EOS model
        # run.itempstart = 0. #np.append(run.itempstart,[matarr[7]]/matarr[16])
        # -----------------------------------------------------------------------
        # Fill EOS parameters for each material
        if config[matlist[imat]]['eos']['type'] == 'MGR':
            run.ieosid.append('MGR')
            matnew = MGRClass()
            run.ieos.append(matnew)
            run.ieos[imat].name = config[matlist[imat]]['eos']['name'] # informal name
            # these variables extracted above
            run.ieos[imat].p0   = p0.to(config['codeunits']['pressure']).magnitude
            run.ieos[imat].up0  = up0.to(config['codeunits']['velocity']).magnitude
            run.ieos[imat].rho0 = rho0.to(config['codeunits']['density']).magnitude
            run.ieos[imat].iev0 = iev0.magnitude
            # these variables define the MGR model
            # reference density for the EOS
            rhoref = Q_(np.asarray(config[matlist[imat]]['eos']['rhoref'],dtype='float'),config['units']['density'])
            run.ieos[imat].rhoref = rhoref.to(config['codeunits']['density']).magnitude
            # c0 velocity
            c0 = Q_(np.asarray(config[matlist[imat]]['eos']['c0'],dtype='float'),config['units']['velocity'])
            run.ieos[imat].c0 = c0.to(config['codeunits']['velocity']).magnitude
            # s1 dimless
            run.ieos[imat].s1 = config[matlist[imat]]['eos']['s1'] # dimless
            # s2 units time/length
            s2 = Q_(np.asarray(config[matlist[imat]]['eos']['s2'],dtype='float'),config['units']['s2'])
            run.ieos[imat].s2 = s2.to(config['codeunits']['s2']).magnitude
            # gamma dimless
            run.ieos[imat].gamma0 = config[matlist[imat]]['eos']['gamma0'] 
            # specific heat capacity cv
            # convert input specific heat capacity to code units and then additional normalization by rho0 in code units
            # KO code units for cvv0 = eu/cm3/K
            spcv = Q_(np.asarray(config[matlist[imat]]['eos']['cv'],dtype='float'),config['units']['sp_heat_cap'])
            cv   = spcv.to(config['codeunits']['sp_heat_cap'])*rho0.to(config['codeunits']['density'])
            run.ieos[imat].cv = cv.magnitude
            # function to calculate the coefficients for the strain form for the Mie Grueneisen EOS
            # page 63 in Wilkins book
            # P = k1*x + k2*x^2 + k3*x^3 + k4*E [megabar]; x = 1 - relvolume; k2=0 for x<0; k4 = gamma0
            run.ieos[imat].calccoefs() # calculate k1-k4
            # set the initial temperature
            run.ieos[imat].t0 = iev0.magnitude/cv.magnitude # K
            # put temperature in the master intialization array
            run.itempstart = np.append(run.itempstart,iev0.magnitude/cv.magnitude)
        if config[matlist[imat]]['eos']['type'] == 'IDG':
            run.ieosid.append('IDG')
            matnew = IDGClass()
            run.ieos.append(matnew)
            run.ieos[imat].name   = config[matlist[imat]]['eos']['name'] # informal name
            # these variables extracted above
            run.ieos[imat].p0     = p0.to(config['codeunits']['pressure']).magnitude
            run.ieos[imat].up0    = up0.to(config['codeunits']['velocity']).magnitude
            run.ieos[imat].rho0   = rho0.to(config['codeunits']['density']).magnitude
            run.ieos[imat].iev0   = iev0.magnitude
            run.ieos[imat].gamma0 = config[matlist[imat]]['eos']['gamma0']
            # specific heat capacity cv
            # convert input specific heat capacity to code units and then additional normalization by rho0 in code units
            # KO code units for cvv0 = eu/cm3/K
            spcv = Q_(np.asarray(config[matlist[imat]]['eos']['cv'],dtype='float'),config['units']['sp_heat_cap'])
            cv   = spcv.to(config['codeunits']['sp_heat_cap'])*rho0.to(config['codeunits']['density'])
            run.ieos[imat].cv = cv.magnitude
            # set the initial temperature
            run.ieos[imat].t0 = iev0.magnitude/cv.magnitude # K
            # put temperature in the master intialization array
            run.itempstart = np.append(run.itempstart,iev0.magnitude/cv.magnitude)
        if config[matlist[imat]]['eos']['type'] == 'SES':
            run.ieosid.append('SES')
            matnew = SESClass()
            run.ieos.append(matnew)
            run.ieos[imat].name = config[matlist[imat]]['eos']['name'] # informal name
            # there variables are extracted above
            run.ieos[imat].p0     = p0.to(config['codeunits']['pressure']).magnitude
            run.ieos[imat].up0    = up0.to(config['codeunits']['velocity']).magnitude
            run.ieos[imat].rho0   = rho0.to(config['codeunits']['density']).magnitude
            run.ieos[imat].iev0   = iev0.magnitude
            # read in initial sound speed in case the table has problems
            cs0 = Q_(np.asarray(config[matlist[imat]]['init']['cs0'],dtype='float'),config['units']['velocity'])
            run.ieos[imat].cs0  = cs0.to(config['codeunits']['velocity']).magnitude
            #run.ieos[imat].gamma0 = config[matlist[imat]]['eos']['gamma0']
            # load the table
            run.ieos[imat].path   = config[matlist[imat]]['eos']['path']
            run.ieos[imat].EOS.R0REF = rho0.to(config['codeunits']['density']).magnitude
            run.ieos[imat].EOS.CS0REF = cs0.to(config['codeunits']['velocity']).magnitude
            run.ieos[imat].sesstd = config[matlist[imat]]['eos']['filestd']
            run.ieos[imat].sesext = config[matlist[imat]]['eos']['fileext']
            # the sesame reader file needs to know about units, so pass it config object
            run.ieos[imat].readtable(config,verbose=verbose)
            #if verbose: print('imat = ',imat, ' reading SESAME tables in ',run.ieos[imat].path)
            # t0 and p0 replaced based on rho0 and e0 to make consistent with SESAME table interpolation
            ptmp, ttmp, cstmp = run.ieos[imat].oneptc(run.ieos[imat].rho0,spenergy.to(config['codeunits']['sp_energy']).magnitude)
            run.ieos[imat].cv = 0. # not used with sesame
            run.ieos[imat].t0 = ttmp
            run.ieos[imat].p0 = ptmp
            run.ieos[imat].rhoref = run.ieos[imat].EOS.R0REF
            # put temperature in the master intialization array
            run.itempstart = np.append(run.itempstart,ttmp)
        if config[matlist[imat]]['eos']['type'] == 'TIL':
            run.ieosid.append('TIL')
            matnew = TILClass()
            run.ieos.append(matnew)
            run.ieos[imat].name = config[matlist[imat]]['eos']['name'] # informal name
            # these variables extracted above
            run.ieos[imat].p0   = p0.to(config['codeunits']['pressure']).magnitude
            run.ieos[imat].up0  = up0.to(config['codeunits']['velocity']).magnitude
            run.ieos[imat].rho0 = rho0.to(config['codeunits']['density']).magnitude
            run.ieos[imat].iev0 = iev0.magnitude
            # these variables define the Tillotson model
            # reference density for the EOS
            rhoref = Q_(np.asarray(config[matlist[imat]]['eos']['rhoref'],dtype='float'),config['units']['density'])
            run.ieos[imat].rhoref = rhoref.to(config['codeunits']['density']).magnitude
            # AA Bulk modulus K0 [pressure]
            AA = Q_(np.asarray(config[matlist[imat]]['eos']['AA'],dtype='float'),config['units']['pressure'])
            run.ieos[imat].AA = AA.to(config['codeunits']['pressure']).magnitude
            # BB [pressure]
            BB = Q_(np.asarray(config[matlist[imat]]['eos']['BB'],dtype='float'),config['units']['pressure'])
            run.ieos[imat].BB = BB.to(config['codeunits']['pressure']).magnitude
            # a dimless
            run.ieos[imat].a = config[matlist[imat]]['eos']['a'] # dimless
            # b dimless a+b = Grueneisen gamma at rhoref
            run.ieos[imat].b = config[matlist[imat]]['eos']['b'] # dimless
            # alpha dimless
            run.ieos[imat].alpha = config[matlist[imat]]['eos']['alpha'] # dimless
            # beta dimless
            run.ieos[imat].beta = config[matlist[imat]]['eos']['beta'] # dimless
            # E0; convert from specific energy input to energy per initial volume
            E0 = Q_(np.asarray(config[matlist[imat]]['eos']['E0'],dtype='float'),config['units']['sp_energy'])
            run.ieos[imat].E0 = E0.to(config['codeunits']['sp_energy']).magnitude*run.ieos[imat].rho0
            # Eiv; convert from specific energy input to energy per initial volume
            EIV = Q_(np.asarray(config[matlist[imat]]['eos']['EIV'],dtype='float'),config['units']['sp_energy'])
            run.ieos[imat].EIV = EIV.to(config['codeunits']['sp_energy']).magnitude*run.ieos[imat].rho0
            # Ecv; convert from specific energy input to energy per initial volume
            ECV = Q_(np.asarray(config[matlist[imat]]['eos']['ECV'],dtype='float'),config['units']['sp_energy'])
            run.ieos[imat].ECV = ECV.to(config['codeunits']['sp_energy']).magnitude*run.ieos[imat].rho0
            # set the initial temperature (K)
            # no temperature defined for original Tillotson with E=0 at P=0
            run.ieos[imat].t0 = 0. 
            # put temperature in the master intialization array (K)
            run.itempstart = np.append(run.itempstart,run.ieos[imat].t0)
            # create the parameter array for Tillotson EOS functions
            run.ieos[imat].fillparams()
        # -----------------------------------------------------------------------
        # Fill in strength model parameters for each material
        if config[matlist[imat]]['str']['type'] == 'HYDRO':
            run.istrid.append('HYDRO')
            strnew = HydroClass()
            run.istr.append(strnew)
            run.istr[imat].name='mat'+str(imat+1)
        if config[matlist[imat]]['str']['type'] == 'VM':
            run.istrid.append('VM')
            strnew = VonMisesClass() # von Mises strength and fracture
            run.istr.append(strnew)
            run.istr[imat].name = 'mat'+str(imat+1)
            # ys Yield strength
            ys = Q_(np.asarray(config[matlist[imat]]['str']['ys'],dtype='float'),config['units']['pressure'])
            run.istr[imat].ys   = ys.to(config['codeunits']['pressure']).magnitude
            # gmod Shear Modulus
            gmod = Q_(np.asarray(config[matlist[imat]]['str']['gmod'],dtype='float'),config['units']['pressure'])
            run.istr[imat].gmod = gmod.to(config['codeunits']['pressure']).magnitude
        # -----------------------------------------------------------------------
        # Dynamic fracture parameters
        fracnew = FractureClass()
        run.ifrac.append(fracnew)
        run.ifrac[imat].name='mat'+str(imat+1)
        pfrac = Q_(np.asarray(config[matlist[imat]]['frac']['pfrac'],dtype='float'),config['units']['pressure'])
        run.ifrac[imat].pfrac = pfrac.to(config['codeunits']['pressure']).magnitude
        # nrhomin is an optional input parameter
        if 'nrhomin' in config[matlist[imat]]['frac'].keys():
            #print('Key nrhomin exists')
            nrhomin = Q_(np.asarray(config[matlist[imat]]['frac']['nrhomin'],dtype='float'),'dimensionless')
            run.ifrac[imat].nrhomin = nrhomin.magnitude
        else:
            print('Key nrhomin does not exist. Using default 0.8')
            run.ifrac[imat].nrhomin = 0.8
        # -----------------------------------------------------------------------
        # if gravity, read in reference position for the initial state
        # gravity is an optional input parameter with each material
        if 'gravity' in config[matlist[imat]].keys():
            #print('Key gravity exists for mat ',imat)
            # flag yes gravity initialization for this layer
            run.grav.matflag = np.append(run.grav.matflag,[1])
            # position for gravity initialization
            refx = Q_(np.asarray(config[matlist[imat]]['gravity']['refpos'],dtype='float'),config['units']['length'])
            run.grav.refpos =  np.append(run.grav.refpos,[refx.to(config['codeunits']['length']).magnitude])
            #print('init position = ',run.grav.refpos)
            # pressure at refpos for gravity initialization
            refx = Q_(np.asarray(config[matlist[imat]]['gravity']['refpres'],dtype='float'),config['units']['pressure'])
            run.grav.refpres = np.append(run.grav.refpres,[refx.to(config['codeunits']['pressure']).magnitude])
            #print('init pressure = ',run.grav.refpres)
        else:
            # flag no gravity initialization for this layer
            run.grav.matflag = np.append(run.grav.matflag,[0])
            run.grav.refpos = np.append(run.grav.refpos,[0])
            #print('Key gravity does not exist for mat ',imat)
    # -----------------------------------------------------------------------
    # read in boundary conditions
    # at this time there is only an inner boundary condition and outer boundary condition
    # future boundaries will include interior gaps and fractures
    run.nbc=2 # number of boundary conditions
    run.ibcid.append(config['boundaries']['ibc'])
    run.ibcid.append(config['boundaries']['obc'])
    #
    # read in descriptors for each boundary condition
    # determine what is really needed as a user input
    # inner BC
    ibc=0
    bcnew = BCClass()
    run.ibc.append(bcnew)
    if config['boundaries']['ibc'] == 'FREE':
        # free surface boundary condition
        run.ibc[ibc].name = 'Free surface '+str(ibc+1)
        if 'ip0' in config['boundaries'].keys():
            pbc = Q_(np.asarray(config['boundaries']['ip0'],dtype='float'),config['units']['pressure'])
        else:
            # if not entered by user, set free surface pressure to zero
            pbc = Q_(np.asarray(0.0,dtype='float'),config['codeunits']['pressure'])
        run.ibc[ibc].pbc = pbc.to(config['codeunits']['pressure']).magnitude
        #run.ibc[ibc].upbc = bcarr[2] # cm/us
        #run.ibc[ibc].rrefbc = bcarr[3] # g/cm3
        #run.ibc[ibc].iebc = bcarr[4] # 1e12 ergs/cm3
        #run.ibc[ibc].vrbc = bcarr[5] # dimless relative volume to rhoref
    if config['boundaries']['ibc'] == 'FIXED':
        # fixed boundary condition
        run.ibc[ibc].name = 'Fixed surface '+str(ibc+1)
        #run.ibc[ibc].pbc = config['boundaries']['op0']
        #run.ibc[ibc].upbc = bcarr[2] # cm/us
        #run.ibc[ibc].rrefbc = bcarr[3] # g/cm3
        #run.ibc[ibc].iebc = bcarr[4] # 1e12 ergs/cm3
        #run.ibc[ibc].vrbc = bcarr[5] # dimless relative volume to rhoref
    # outer BC
    ibc=1
    bcnew = BCClass()
    run.ibc.append(bcnew)
    if config['boundaries']['obc'] == 'FREE':
        # free surface boundary condition
        run.ibc[ibc].name = 'Free surface '+str(ibc+1)
        if 'op0' in config['boundaries'].keys():
            pbc = Q_(np.asarray(config['boundaries']['op0'],dtype='float'),config['units']['pressure'])
        else:
            # if not entered by user, set free surface pressure to zero
            pbc = Q_(np.asarray(0.0,dtype='float'),config['codeunits']['pressure'])
        run.ibc[ibc].pbc = pbc.to(config['codeunits']['pressure']).magnitude
    if config['boundaries']['obc'] == 'FIXED':
        # fixed boundary condition
        run.ibc[ibc].name = 'Fixed surface '+str(ibc+1)    #
    # -----------------------------------------------------------------------
    # stop time
    tstop = Q_(np.asarray(config['tstop'],dtype='float'),config['units']['time'])
    run.tstop = tstop.to(config['codeunits']['time']).magnitude
    # first time step
    dtstart = Q_(np.asarray(config['dtstart'],dtype='float'),config['units']['time'])
    run.dtstart = dtstart.to(config['codeunits']['time']).magnitude
    # time between output snapshots
    dtoutput = Q_(np.asarray(config['dtoutput'],dtype='float'),config['units']['time'])
    run.time_skip = dtoutput.to(config['codeunits']['time']).magnitude
    run.next_time_dump = run.time_skip
    # output format
    if config['outputformat'] == 'BIN':
        run.binoutput = True
    else:
        run.binoutput = False
    # set geometry and set dgeom variable d=1., 2., or 3. d must be a float
    run.dflag = config['geometry']
    if run.dflag == 'PLA':
        run.dgeom = 1.0
    elif run.dflag == 'CYL':
        run.dgeom = 2.0
    elif run.dflag == 'SPH':
        run.dgeom = 3.0
    else:
        sys.exit("FATAL ERROR: invalid geometry "+str(run.dflag))
    if (run.ixstart[0] < 0) & (run.dflag != 'PLA'):
        # cannot have negative r values for CYL or SPH
        sys.exit('FATAL ERROR: invalid geometry, negative xstart with '+run.dflag)
    # gravitational acceleration
    gravity = Q_(np.asarray(config['gravity'],dtype='float'),config['units']['gravity'])
    run.grav.gravity = gravity.to(config['codeunits']['gravity']).magnitude
    # void pressure
    if 'pvoid' in config[matlist[imat]].keys():
        pvoid = Q_(np.asarray(config['pvoid'],dtype='float'),config['units']['pressure'])
        run.pvoid = pvoid.to(config['codeunits']['pressure']).magnitude
    else:
        # set a default value [in code units]
        run.pvoid = 1.E-9
    #
    # -----------------------------------------------------------------------
    # print parameters
    if verbose:
        print(run) # print pyKO model run parameters
        for imat in range(run.nmat):
            print(run.istr[imat])  # print material strength parameters
            print(run.ifrac[imat]) # print dynamic fracture parameters
            print(run.ieos[imat])  # print material EOS parameters
        for ibc in range(run.nbc):
            print(run.ibc[ibc]) # print boundary condition parameters
        #timestr = 'Time '.join(config['codeunits']['time'])+' \n'
        print('\nCODE UNITS:')
        print('  Time in ',config['codeunits']['time'])
        print('  Length in ',config['codeunits']['length'])
        print('  Pressure in ',config['codeunits']['pressure'])
        print('  Temperature in ',config['codeunits']['temperature'])
        print('  Density in ',config['codeunits']['density'])
        print('  Sp. energy in ',config['codeunits']['sp_energy'])
        print('  Sp. heat capacity in ',config['codeunits']['sp_heat_cap'])
        print('  Velocity in ',config['codeunits']['velocity'])
        print('  iev0 is specific internal energy * initial density')
    #print('ipstart = ',run.ipstart)
    #
    # END read in configuration yaml file
    ############################################################## 
#
############################################################## 
# RUN PROGRAM
############################################################## 
#
def run(fin='pyko.yml',fout='pyko-output.dat',ftype='YAML',verbose=True,userdtstart=0.,usertstepscale=0.,binoutput=True,debug=False):
    ############################################################## 
    # main program input arguments
    # initialize main run variables
    run                = RunClass(fin=fin,fout=fout,ftype=ftype)    # initialize run parameters class
    #
    ############################################################## 
    # read input file
    if run.inputfiletype == 'BORG':
        readinput_borg(run)
    # yaml configuration file will update outputfilename 
    if run.inputfiletype == 'YAML':
        readinput_yaml(run,verbose=verbose)
    #
    ##############################################################
    print('\n pyKO STARTING RUN \n')
    # initialize problem domain data structure
    # pass in the problem run parameters
    data = DomainClass() # create main data structure
    data.makegrid(run,verbose=verbose) # create initial grid; also calculates a guess for the initial time step
    if userdtstart > 0.:
        # for testing comparisons to fKO with fixed initial time step
        if verbose: print('Using user supplied initial time step (microsec): ',userdtstart)
        run.dtstart = userdtstart
        data.deltat   = run.dtstart      # current time step; while in time loop
        data.deltat_next = run.dtstart   # next time step; calculated at end of time loop
    if usertstepscale > 0.:
        if verbose: print('Using user supplied reduction factor for time steps: ',usertstepscale)
        run.tstepscale = usertstepscale
    run.debugflag = debug
    data.writefirstoutput(run,checkexistsflag=False,verbose=verbose,overwriteflag=True) # write initial grid to outputfile - right now just the header
    
    # testing: can call one time step only or a fixed number of time steps to check the simulation parameters
    # advance one time step
    #data.advancetime(run)
    #for steps in range(20):
    #    data.advancetime(run,verbose=verbose)
        
    # main time loop
    print('pyKO is running.....')
    t0 = time.time()
    while data.time[3] < run.tstop:
        data.advancetime(run,verbose=verbose)
    #
    if data.stepn > max(run.outputsteps): # write last time step if needed
        if run.binoutput:
            data.binaryoutputpint(run) # pickle output 
        else:
            data.appendoutput(run) # ascii columns
    t1 = time.time()
    
    # diagnostic energy checks and run time information
    print('Zeroth and final mvtotal: ',run.outputmvtot[0],run.outputmvtot[-1])
    print('Zeroth and final KE+IE: ', \
          run.outputietot[0]+run.outputketot[0],run.outputietot[-1]+run.outputketot[-1])
    print('KE+IE ratio (final/zeroth)):', \
          (run.outputietot[-1]+run.outputketot[-1])/(run.outputietot[0]+run.outputketot[0]))
    print('First output step and final output step KE+IE: ', \
          run.outputietot[1]+run.outputketot[1],run.outputietot[-1]+run.outputketot[-1])
    print('KE+IE ratio (final/first)):', \
          (run.outputietot[-1]+run.outputketot[-1])/(run.outputietot[1]+run.outputketot[1]))
    print('\n pyKO FINISHED RUN, simulation time = ',data.time[1+2],' us \n')
    print(' pyKO Wall Clock Main Loop Run Time=', t1-t0, ' s')
    return
#
############################################################## 
############################################################## 
# Default action for this python program
#
if __name__ == "__main__":
    # Useful habit to check the input parameters before running the code
    #filein = './test5/test5-water-free.yml'
    #test=RunClass(fin=filein)
    #test.checkinput()

    #filein = './test6/test6-gap.yml'
    #filein = './test7/test7-z-stagnation.yml'
    #filein = './test7/test7-gap-test.yml'
    #run(fin=filein,verbose=True)
    
    # The fortan code uses an initial time step of 0.001 microseconds; set as input here for code comparison
    # otherwise pyKO has an initialization section that estimates a good first time step
    filein = './test12/test12-gravity.yml'
    run(fin=filein,userdtstart=0.001,verbose=True)
#    
############################################################## 
############################################################## 
# END OF FILE
############################################################## 
############################################################## 
