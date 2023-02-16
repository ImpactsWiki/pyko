# pyKO hydrocode

pyKO is a 1d Lagrangian Elastic-Plastic Hydrocode written in python.
This code was developed from the book <a href="https://link.springer.com/book/10.1007/978-3-662-03885-7">Computer Simulation of Dynamic Phenomena by Mark Wilkins</a> and the <a href="https://www.eng.mu.edu/shockphysics/KO/">fortran KO code by John Borg</a>.

This resource is available as part of the <a href="https://impacts.wiki">Impacts Community Wiki Project</a>.

This code is under development.

Current features include:
* Planar, cylindrical, and spherical geometry
* Elastic-plastic von Mises model
* EOS models: ideal gas, Mie-Grueneisen, and SESAME tables
* Fixed and free boundary conditions
* Ascii and binary (using pickle) output
* Jupyter notebook plotting and animation scripts
* Verification tests

Wish list:
* Fracture and gaps
* yaml input files
* Numba performance enhancements (the current python code is about 10x slower than fortran)
* More strength models
* Tillotson EOS
* Absorbing boundary condition

