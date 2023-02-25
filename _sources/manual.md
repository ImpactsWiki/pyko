# User Manual

pyKO is a 1-D elastic-plastic, 2nd order, Lagrangian hydrocode<br>
Adapted from John Borg's Fortran KO code v11 (https://www.eng.mu.edu/shockphysics/KO/)<br>
Based on Mark Wilkins, Computer Simulation of Dynamic Phenomena, Springer-Verlag, 1999

## Features
    * planar, cylindrical, spherical geometries
    * fixed and free boundary conditions
    * Von Mises shear strength
    * EOS
        * Ideal gas
        * Mie-Grueneisen
        * Tabular EOS
    * I/O
        1. For testing and comparison to fortran KO:
           Fixed formatting ascii input and output files
        2. Preferred I/O
           Yaml configuration (input) files with unit conversion using pint
           Binary output files using pickle with units

### Wish list
    * Fracture and gaps
    * Variable mesh spacing
    * Tillotson EOS
    * Gravity
    * More strength models

## pyKO files
    * pyko.py: primary code file
    * eos_table.py: module for using EOS tables developed for ANEOS model

## Required packages
    * Python 3.9+
    * Numpy
    * Pint
    * Yaml
    * Pickle

Required for test cases and visualization examples:
    * Pandas
    * Hvplot
    * Matplotlib
    * Scipy (for sodshock.py)
    


### Test cases in repository
Test cases are included as Jupyter Notebooks, for training and validation.<br>
The ANEOS models used here are for TESTING ONLY; they are not research-level EOS models.<br>
Each test has fortran-style ascii input file and new pyKO yaml configuration file.<p>

    * Test 1: Two Mie-Grueneisen plates planar impact, comparison to fortran KO
    * Test 1b: Three plates
    * Test 2: Sod shock tube test for ideal gases, comparison to analytic solution and fortran KO
    * Test 3: EOS table ideal gas Sod shock tube test, comparison to analytic solutoin
    * Test 4: EOS table planar plate impact, comparison to Mie-Grueneisen fortran KO
    * Test 5: EOS table expansion of shocked water into free space

