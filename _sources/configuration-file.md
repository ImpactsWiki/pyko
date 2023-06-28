---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
   display_name: Python 3
   language: python
   name: python3
---

# Configuration file

pyKO uses a yaml format configuration file.

Yaml readers are not robust for exponential notion and numbers are sometimes interpreted as a string. Place a zero after the decimal point when using exponential notion.

The checkinput function prints the RunClass information for user inspection to confirm that the configuration file has been read in correctly. 

All the pyKO example configuration files use mks for the input parameters. Note that the checkinput function presents values in code units. See [](config:units).


```{code-block} python

# =======================
# pyKO configuration file
# =======================
# ===>>>> Exponential numbers: put a zero after the decimal point and before the e or E <<<<===
# yaml does not always process exponential notation as a number and can mistakenly load as a string
# check your input parameters, especially is using mixed units:
#       import pyko
#       import eos_table as etab      <<-- needed if using SESAME tables
#       run = RunClass()
#       run.checkinput('filename')
#
```

## Input-output parameters

Specify the problem name; output file name; output file format 'BIN' or 'ASC'

```{code-block} python
# --------------------------------------------------------------------
# I/O PARAMETERS
#
# problem name; output file name; output file format 'BIN' or 'ASC'
name           : 'Test 7 Z Stagnation'
outputfilename : './test7/pyko-test7-z-stagnation-bin2.dat'
outputformat   : 'BIN'
#
# tstop = stop time for the problem
# if dtstart is >0 then it overrides the internal calculation for the first time step; used to match fortran KO output
# dtoutput = output snapshot frequency in time
tstop      : 0.5E-6
dtstart    : 0.001E-6
dtoutput   : 0.01E-6
#
```

## Material properties and gridding parameters

Each material must be entered in order in the spatial domain. Negative position values are only permitted in planar geometry.

For each material, specify the inner edge, length and number of cells. The number of nodes is double the number of cells plus any necessary boundary node.

See [](page:eos) and [](page:strength) for explanation of material input parameters.

```{code-block} python
# --------------------------------------------------------------------
# MATERIAL PROPERTIES AND GRIDDING PARAMETERS
#
mat1:
    mesh:
        cells  : 200 
        xstart : -0.002
        length :  0.002
    init:
        up0    : 15000.0
        rho0   : 2700.0
        p0     : 1.e5
        e0     : 156808.87895946743
    eos: 
        name   : 'Aluminum'
        type   : 'SES'
        path    : '../aneos-Al/'
        filestd : 'NEW-SESAME-STD.TXT'
        fileext : 'NEW-SESAME-EXT.TXT'
    str:
        type   : 'VM'
        gmod   : 24.8E9
        ys     : 0.10E9
    frac:
        pfrac  : 1.0E20
#
mat2:
    mesh:
        cells  : 300
        xstart : 0.0
        length : 0.003
    init:
        up0    : 0.0
        rho0   : 3220.0
        p0     : 1.e5
        e0     : 73277.974035
    eos: 
        name   : 'Forsterite'
        type   : 'SES'
        path    : '../aneos-forsterite-2019/'
        filestd : 'NEW-SESAME-STD.TXT'
        fileext : 'NEW-SESAME-EXT.TXT'
    str:
        type   : 'VM'
        gmod   : 47.7E9
        ys     : 0.12E9
    frac:
        pfrac  : 1.0E20
#
mat3:
    mesh:
        cells  : 200
        xstart : 0.005
        length : 0.002
    init:
        up0    : 0.0
        rho0   : 2650.0
        p0     : 1.e5
        e0     : 29183.3658
    eos: 
        name   : 'Quartz'
        type   : 'SES'
        path    : '../aneos-silica-slvtv0.1/'
        filestd : 'NEW-SESAME-STD.TXT'
        fileext : 'NEW-SESAME-EXT.TXT'
    str:
        type   : 'VM'
        gmod   : 47.7E9
        ys     : 0.12E9
    frac:
        pfrac  : 1.0E20
#
```
## Domain parameters

Specify the inner (ibc) and outer (obc) boundary conditions as either FIXED or FREE surfaces.

Specify the problem geometry: planar, cylindrical or spherical. Cylindrical and spherical problems must have positive values for the spatial domain.

```{code-block} python
# --------------------------------------------------------------------
# DOMAIN PARAMETERS
#
# 'FIXED' or 'FREE'
# ibc = inner boundary condition (at minimum x); obc = outer boundary condition (at maximum x)
# fracture and interior gaps are not implemented yet
boundaries: 
    ibc    : 'FREE'
    obc    : 'FREE'
#
# 'PLA' = planar; 'CYL' = cylindrical; 'SPH' = spherical
geometry   : 'PLA'
#
# gravitational acceleration; negative means force toward negative position x
gravity    : 0.0
#
```

(config:units)=
## Units


```{code-block} python
# --------------------------------------------------------------------
# UNITS
#
# user input parameters units
# Provide the units for the user configuration parameters.
# This feature is provided for user convenience and allows for entering material parameters in more natural mixed units.
# The user parameter values will be converted to self-consistent code units during input processing.
# Output file units will be the same as the input units and included in the binary pickle files using the pint package.
# mesh definition units must be the same.
# s2 is the quadratic term for the Mie Grueneisen EOS
# here mks because life is easier in consistent units
units:
    time        : 'second'
    length      : 'meter'
    velocity    : 'meter/second'
    density     : 'kg/m^3'
    mass        : 'kg'
    pressure    : 'Pa'
    temperature : 'K'
    energy      : 'J'
    sp_energy   : 'J/kg'
    sp_entropy  : 'J/kg/K'
    sp_heat_cap : 'J/kg/K'
    gravity     : 'm/s^2'
    s2          : 'second/meter'
#
# EOS table units
# These are the units for the Stewart Group ANEOS tables.
# The table units will be converted to code units during input processing.
tableunits:
    # in NEW-SESAME-STD.TXT
    density      : 'g/cm^3'
    temperature  : 'K'
    pressure     : 'GPa'
    sp_energy    : 'MJ/kg'
    hfree_energy : 'MJ/kg'
    # in NEW-SESAME-EXT.TXT (also includes a phase identification flag)
    sp_entropy   : 'MJ/K/kg'
    sound_speed  : 'cm/s'
    sp_heat_cap  : 'MJ/kg/K'
#
# code units
# The routine that sets up the master run class object converts all inputs into code units,
# including converting any EOS table units to code units.
# These are the default set of code units for the original KO code.
# Energy: 1 eu = 10^12 ergs
# Pressure: 10^12 dynes/cm^2 = 100 kJ = 1 megabar     -> P=E/V : dynes/cm^2 = erg/cm^3
# Internal energy per original volume = 10^12 ergs/cm^3 = eu/cm^3 = 100 GJ/m^3
# Heat capacity per original volume = 10^12 ergs/cm^3/K = eu/cm^3/K = 100 GJ/m^3/K
# Velocity = cm/microseconds = 10 km/s
# Note that these units need to be declared with a coefficient of 1 (no scaling factor)
# Scaling factors are declared using the pint define command within function readinput_yaml:
#     ureg.define('eu = 1.0E12 ergs')   # energy unit
codeunits:
    time            : 'microseconds'
    length          : 'cm'
    mass            : 'g'
    density         : 'g/cm^3'
    relative_volume : 'dimensionless' 
    velocity        : 'cm/microsecond'
    pressure        : 'megabar'
    temperature     : 'K'
    energy          : 'eu'
    sp_energy       : 'eu/g'
    sp_heat_cap     : 'eu/K/g'
    sp_entropy      : 'eu/K/g'
    ie_perv0        : 'eu/cm^3'
    cv_perv0        : 'eu/cm^3/K'
    gravity         : 'cm/microseconds^2'
    s2              : 'microseconds/cm'
# It is possible to change the code units to a different self-consistent set of units.
# The code units presented here set up the initialization routine for the problem domain.
# Changing to a different set of units should work but has not been tested yet.
# Do not make modifications to the central code routines 
# that make the code dependent on a specific set of units.
# --------------------------------------------------------------------
#
# =========================
# end of configuration file
# =========================


```
