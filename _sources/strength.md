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

(page:strength)=
# Material strength

## Von Mises strength

Von Mises strength requires two material parameters: the shear modulus and the yield strength.

Example configuration file entry (mks):
```
mat1:
    str:
        type   : 'VM'
        gmod   : 47.7E9
        ys     : 0.12E9
```

## Hydrodynamic material

A hydrodynamic material is indicated with this configuration file entry:
```
mat1:
    str:
        type   : 'HYDRO'
```

## Fracture

Fracture is not yet implemented but coming soon.

Fracture requires the Mie-Grueneisen EOS or a tabular EOS with a tension region.


