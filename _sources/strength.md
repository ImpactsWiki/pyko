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
# Material strength and fracture

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


## Steinberg-Guinan strength model

A very beta version of the Steinberg-Guinan strength model is in v0.6.1x-dev based on Wilkins book. It needs to be debugged. See Test15 notebook.

Example configuration file entry (mks):
```
mat1:
    str:
        name   : 'Al-Wilkins-SG'
        type   : 'SG'
        Y0     : 0.29E9
        Ymax   : 0.68E9
        beta   : 125.0
        n      : 0.1
        b      : 8.0
        h      : 6.2E-4
        Tm0    : 1220.0
        mu0    : 27.6E9
```


## Hydrodynamic material

A hydrodynamic material (shear modulus $G=0$) is indicated with this configuration file entry:
```
mat1:
    str:
        type   : 'HYDRO'
```

## Fracture and void space

Dynamic fracture is implemented as an optional feature. Fracture needs development to be more stable or ignored under extreme conditions.

Fracture requires an EOS with a tension region.

The fracture algorithm is not used if the frac input block is not included. Comment out or remove to prevent facture.

Input parameters: pfrac is the fracture stress. nrhomin is the maximum distension (rhomin/rhoref). The default value for nrhomin is 0.8.

Example configuration file entry (mks):
```
mat1:
    frac:
        pfrac   : 1.0E9
        nrhomin : 0.8
```

The code has the ability to close void spaces when surfaces come into contact. 
The current implementation of void closure uses a simplified calculation compared to Wilkins book for the time step for contact.
The full convergence loop for contact time should be implemented in the future.
