# README

AQUA WATER EOS TABLE 

AQUA EOS from https://github.com/mnijh/AQUA

Make-aneos-tables.ipynb converts AQUA to standard SESAME 201 + 301 table format.

NEW-SESAME-STD.TXT: Standard length Sandia-style SESAME file with 201 table and 301 table (density, temperature, pressure, sp. internal energy, Helmholtz free energy). 301 table units: g/cm3, K, GPa, MJ/kg, MJ/kg.
Note: No HFE tabulated.

NEW-SESAME-EXT.TXT: SESAME-style table with extra variables from ANEOS. Contains the standard 201 table and non-standard 301-extra-variables EOS table. The 301 table has: density grid values, temperature grid values, sp. entropy(T,rho), sound speed(T,rho), sp. heat capacity(T,rho), KPA flag(T,rho). 2-D arrays list all densities, looping over each temperature. 301 table units: g/cm3, K, MJ/K/kg, cm/s, MJ/K/kg, integer flag, integer flag. The KPA flag is an ANEOS output with phase information.

SESAME table format specification:
https://github.com/ststewart/aneos-forsterite-2019/blob/master/EOS-docs/Lyon-Johnson-1992-SESAME-database.pdf

