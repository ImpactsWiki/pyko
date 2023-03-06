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

A <a href="https://en.wikipedia.org/wiki/Von_Mises_yield_criterion">Von Mises material</a> yields when the second invariant of deviatoric stress 
$J_{2}$ reaches a critical value (input parameter ys=yield stress=$\sigma_{vm}$). Under uniaxial stress, the material yields when the principal stress exceeds the yield stress, $\sigma_1 > \sigma_{vm}$. Planar shock waves are uniaxial strain, and the Hugoniot Elastic Limit $\sigma_{HEL} = (1-\nu) \sigma_{vm} / (1-2\nu)$.


Von Mises strength requires two material parameters: the shear modulus (gmod) and the yield stress (ys).

Example configuration file entry (mks):
```
mat1:
    str:
        type   : 'VM'
        gmod   : 47.7E9
        ys     : 0.12E9
```



## Hydrodynamic material

A hydrodynamic material (shear modulus $G=0$) is indicated with this configuration file entry:
```
mat1:
    str:
        type   : 'HYDRO'
```

## Fracture

Fracture is not yet implemented but coming soon.

Fracture requires the Mie-Grueneisen EOS or a tabular EOS with a tension region.

Example configuration file entry (mks):
```
mat1:
    frac:
        pfrac  : 1.0E7
```