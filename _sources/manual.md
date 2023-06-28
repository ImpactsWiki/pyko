# User manual

pyKO is a one-dimensional Lagrangian elastic-plastic hydrocode written in python. The code uses the Von Neumann finite difference scheme with second order accuracy. 

This code was developed from the book <a href="https://link.springer.com/book/10.1007/978-3-662-03885-7">Computer Simulation of Dynamic Phenomena by Mark Wilkins</a> (Springer-Verlag, 1999) and the <a href="https://www.eng.mu.edu/shockphysics/KO/">fortran KO code v11 by John Borg</a>.

## Features

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

### Wish list

* Full Wilkins convergence loop for time step to void closure
* Variable mesh spacing
* More strength models (will port commonly used Collins et al. MAPS 2004 with various additions)
* Material interface separation criteria
* Restart from binary
* Absorbing boundary condition
* Numba performance enhancements (this code is about 10x slower than fortran)
* Thermal diffusion

## pyKO source files

* `pyko.py`: primary hydrocode and I/O functions
* `eos_table.py`: module for EOS functions with ANEOS and Tillotson

## Required packages

See environment.yml and requirements.txt files in GitHub.

Tested with:

* Python (developed in 3.9; runs in 3.11)
* Numpy (1.24.2; 1.25.0)
* Pint (0.20.1;0.22)
* Yaml (6.0)
* Pickle (4.0)

Python packages for input examples and visualization in Jupyter notebooks:

* JupyterLab (3.5.0; 3.6.3) (hvplot not compatible with 4.0)
* Pandas (1.4.2; 2.0.2)
* Hvplot (0.8.2; 0.8.4)
* Matplotlib (3.7.0; 3.7.1)
* ipympl (0.9.3)
* Scipy (1.9.1) (needed for `sodshock.py`)
* ffmpeg (6.0.0)

### Example input files in repository
Example input files are included with Jupyter Notebooks, for training and visualization.

* Test 1: Two Mie-Grueneisen plates planar impact, comparison to fortran KO
* Test 1b: Three MGR plates
* Test 2: Sod shock tube test for ideal gases, comparison to analytic solution and fortran KO
* Test 3: EOS table ideal gas Sod shock tube test, comparison to analytic solutoin
* Test 4: EOS table planar plate impact, comparison to Mie-Grueneisen fortran KO
* Test 6: MGR plates with gaps, example of closing gaps
* Test 6b: SES plates with gaps, example of closing gaps
* Test 9: Spall test with MGR plates
* Test 9b: Spall test with TIL and MGR plates
* Test 12: Gravity Al target and impactor
* Test 12b: Gravity Forsterite target and impactor
