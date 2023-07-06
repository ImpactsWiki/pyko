# pyKO CHANGELOG
https://github.com/ImpactsWiki/pyko

## v0.6.x-dev 2023-07-05 in progress
* Beta version of Steinberg-Guinan model
* Fracture algorithm not used if the fracture input block is not used

## v0.6.1 2023-07-04 Update

* Added constant heat capacity temperature to Mie-Grueneisen EOS; initialize positive temperature with zero internal energy
* Added constant heat capacity temperature to Tillotson; initialize positive temperature with zero internal energy
* Now initial temperature and heat capacity are required inputs for TIL and MGR EOS
* Added functions to calculate the cold curves in MGR and TIL classes
* Fixed bugs in iSALE function for Tillotson EOS; now Hosono, iSALE, and Asphaug versions of Tillotson agree
* Added sound speed calculation to Mie-Grueneisen EOS; only used for time step calcs and output; artificial viscosity calculation is unchanged.
* Removed s2 input for MGR EOS; Wilkins implementation of MGR assumes linear Us=c+s*up and gamma/volume=constant
*   -- Wish list: option for non-linear Us, non-constant gamma/volume MGR in future
* See Test13 and Test14 jupyter notebooks for these updates and EOS function checks
* Updated some comments/units
* Changed IDG input key gamma0 to gamma.


## v0.6.0 2023-06-27 First Release

### Features

* Geometries
    * Planar, cylindrical, or spherical
* Boundary conditions
    * Fixed: symmetric, zero particle velocity
    * Free: zero pressure
* Elastic-plastic model
    * Von Mises
* Interfaces
    * Initial gaps
    * Dynamic fracture
    * Void closure on contact
* Gravitational acceleration
    * Initialization of isothermal layers only
* Equations of state
    * Ideal gas
    * Mie-Gr&uuml;neisen
    * Tillotson (needs more testing)
    * Tabular/SESAME, with tools for tabulated ANEOS; interpolation schemes are currently very simple
* Input and output
    * For testing and comparison to fortran KO:<br>
      Fixed formatting ascii input and output files
    * Preferred I/O:<br>
      Yaml configuration (input) files with unit conversion using pint<br>
      Binary output files using pickle with units
    * Customizable outputs: <br>
      User-defined cell data output in binary format<br>
      Debug output in binary format with all node variables

Fortran-style input files are deprecated for the yaml input file format. Not backwards compatible to earlier stages of development.

