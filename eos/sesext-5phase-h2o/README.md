# README

5PHASE H2O EOS TABLE 

5PHASE H2O EOS by S. T. Stewart in  https://doi.org/10.1111/j.1945-5100.2008.tb00657.x

Make-aneos-tables.ipynb converts 5PHASE EOS to standard SESAME 201 + 301 table format.

NEW-SESAME-STD.TXT: Standard length Sandia-style SESAME file with 201 table and 301 table (density, temperature, pressure, sp. internal energy, Helmholtz free energy). 301 table units: g/cm3, K, GPa, MJ/kg, MJ/kg.
Note: No HFE tabulated.

NEW-SESAME-EXT.TXT: SESAME-style table with extra variables from ANEOS. Contains the standard 201 table and non-standard 301-extra-variables EOS table. The 301 table has: density grid values, temperature grid values, sp. entropy(T,rho), sound speed(T,rho), sp. heat capacity(T,rho), KPA flag(T,rho). 2-D arrays list all densities, looping over each temperature. 301 table units: g/cm3, K, MJ/K/kg, cm/s, MJ/K/kg, integer flag, integer flag. The KPA flag is an ANEOS output with phase information.

SESAME table format specification:
https://github.com/ststewart/aneos-forsterite-2019/blob/master/EOS-docs/Lyon-Johnson-1992-SESAME-database.pdf


## ORIGINAL 5PHASE EOS TABLE FORMAT
h2o_table_v8.3bNT.txt is a tabular EOS text file; the sav file is an IDL save file.


The IDL save file is in NIST units: K, MPa, kg/m3, kJ/kg

The 201-SESAME style table is in: K, GPa, g/cm3, MJ/kg, MJ/kg/K

In the h2o_table_v8.3bNT.txt file, the format is a 201 style table except that the Helmholtz free energy has been replaced with specific energy. 

The file 5PhaseEOS.zip contains both a no-tension table and a tension table in v8.0 gridding (no phase info).

