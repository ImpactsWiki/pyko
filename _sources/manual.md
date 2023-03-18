# User manual

pyKO is a one-dimensional Lagrangian elastic-plastic hydrocode written in python.

This code was developed from the book <a href="https://link.springer.com/book/10.1007/978-3-662-03885-7">Computer Simulation of Dynamic Phenomena by Mark Wilkins</a> (Springer-Verlag, 1999) and the <a href="https://www.eng.mu.edu/shockphysics/KO/">fortran KO code v11 by John Borg</a>.

## Features

* Planar, cylindrical, spherical geometries
* Fixed (symmetric, zero particle velocity) and free (zero pressure) boundary conditions
* Von Mises shear strength
* Initial gaps, dynamic fracture, void closure on contact
* Gravitational acceleration (needs more testing)
* EOS
    * Ideal gas
    * Mie-Grueneisen
    * Tillotson (needs more testing)
    * Tabular/SESAME, with features for tabulated ANEOS
* Input and output
    * For testing and comparison to fortran KO:<br>
      Fixed formatting ascii input and output files
    * Preferred I/O:<br>
      Yaml configuration (input) files with unit conversion using pint<br>
      Binary output files using pickle with units
    * Customizable outputs: <br>
      User-defined cell data output in binary format<br>
      Debug output in binary format with all node variables

### Wish list

* Full Wilkins convergence loop for time step to void closure
* Variable mesh spacing
* More strength models (will port commonly used Collins et al. MAPS 2004 with various additions)
* Material interface separation criteria
* Restart from binary
* Absorbing boundary condition
* Numba performance enhancements (this code is about 10x slower than fortran)

## pyKO source files

* `pyko.py`: primary hydrocode and I/O functions
* `eos_table.py`: module for EOS functions with ANEOS and Tillotson

## Required packages

* Python (developed in 3.9.16)
* Numpy (1.24.2)
* Pint (0.20.1)
* Yaml (6.0)
* Pickle (4.0)

Required python packages for test cases and visualization examples:<p>

* JupyterLab (3.5.0)
* Pandas (1.4.2)
* Hvplot (0.8.2)
* Matplotlib (3.7.0)
* Scipy (1.9.1) (needed for `sodshock.py`)

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

