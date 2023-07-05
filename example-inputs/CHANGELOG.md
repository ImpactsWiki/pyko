# CHANGELOG pyKO example inputs

## 2023-07-05
* Added Test8 notebook for a lab experiment on ice and two tabular water EOS that were not generated with ANEOS: the 5PHASE EOS and AQUA EOS. 
  See those eos directories for how to convert a table into a format compatible with pyKO.

## 2023-07-03
* Added Test13 notebooks to compare Tillotson, Mie Grueneisen, and tabular EOS. Examples for Al and Fe. Now including temperature calculations.
* Added Test14 notebook to compare Hosono, iSALE, and Asphaug implementations of Tillotson. Cross check of calculated sound speeds.

## 2023-06-28
* Added x-t diagram to test1b jupyter notebook

## v0.6.0 2023-06-27

Example input files are included with Jupyter Notebooks, for training and visualization.

* Test 1: Two Mie-Grueneisen plates planar impact, comparison to fortran KO
* Test 1b: Three MGR plates
* Test 2: Sod shock tube test for ideal gases, comparison to analytic solution and fortran KO
* Test 3: EOS table ideal gas Sod shock tube test, comparison to analytic solutoin
* Test 4: EOS table planar plate impact, comparison to Mie-Grueneisen fortran KO, includes movie visualization examples
* Test 6: MGR plates with gaps, example of closing gaps
* Test 6b: SES plates with gaps, example of closing gaps
* Test 9: Spall test with MGR plates
* Test 9b: Spall test with TIL and MGR plates
* Test 12: Gravity Al target and impactor
* Test 12b: Gravity Forsterite target and impactor
