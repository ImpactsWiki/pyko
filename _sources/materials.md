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

(page:materials)=
# Example materials

Disclaimer: Most of these material models are appropriate for teaching/learning purposes. Research-level models must be vetted for the specific application. 


## Water Ice

[//]: # (https://link.springer.com/article/10.1023/A:1021134128038)

The properties of ice depend on the grain structure, temperature and strain rates. Nominal values for the Young's modulus $E=9$ GPa and poisson's ratio $\nu=0.32$ {cite}`Petrovic2003`. The corresponding shear modulus $G = \frac{E}{2(1+\nu)} = 3.38$ GPa. Cold ice under uniaxial strain has a Hugoniot Elastic Limit of 0.62 GPa {cite}`Stewart2005`. The average tensile strength for ice is 1.43 MPa at low strain rates {cite}`Petrovic2003`; the dynamic spall strength of ice is about 17 MPa {cite}`Lange1983`. The relationship between the Hugoniot Elastic Limit and the von Mises yield stress is:

$\sigma_{vm}  = \frac{(1-2\nu) \sigma_{HEL}}{(1-\nu)} = 0.328$ GPa.

Bulk modulus of ice is about 10 GPa at 150 K {cite}`You2018`.

There are multiple options for the equation of state for H<sub>2</sub>O, each with various strengths and weaknesses and should be chosen based on the desired application.

Example configuration file entry (mks):
```
mat1:
    str:
        type   : 'VM'
        gmod   : 3.38E9
        ys     : 0.328E9
        kmod   : 10.0E9
        cs0    : 3305.0
    frac:
        pfrac  : 17.0E6
```

## Aluminum 6061

<!--- (https://www.matweb.com/search/datasheet.aspx?matguid=b8d536e0b9b54bd7b69e4124d8f1d20a&ckck=1) --->
<!--- https://ascelibrary.org/doi/epdf/10.1061/%28ASCE%29EM.1943-7889.0000264 --->
<!--- https://www.researchgate.net/profile/John-Lewandowski-3/publication/229004194_Dynamic_Tensile_Deformation_of_Aluminum_Alloy_6061-T6_and_6061-OA/links/02e7e5240329c412df000000/Dynamic-Tensile-Deformation-of-Aluminum-Alloy-6061-T6-and-6061-OA.pdf --->


Al 6061 data from https://www.matweb.com/.

```
mat1:
    str:
        type   : 'VM'
        gmod   : 26.0E9
        ys     : 207.0E6
    frac:
        pfrac  : 276.0E6
```

## Copper

<!---  (https://www.matweb.com/search/datasheet.aspx?matguid=9aebe83845c04c1db5126fada6f76f7e) --->
<!--- https://aip.scitation.org/doi/full/10.1063/1.3607294 --->

Annealed copper shear modulus (46 GPa) and Poisson's ratio (0.343) from https://www.matweb.com/.

Hugoniot elastic limit and tensile strength are dependence on temperature, strain rate, propagation distance, etc. {cite}`Zaretsky2013`. These are nominal values for demonstration purposes only.
For an HEL of 0.2 GPa, $\sigma_{vm}=0.0956$ GPa.

Dynamic spall strength from https://aip.scitation.org/doi/full/10.1063/1.3607294

```
mat1:
    str:
        type   : 'VM'
        gmod   : 46.0E9
        ys     : 0.0956E9
    frac:
        pfrac  : 1.35E9
```

## Stainless steel 304

Most parameters from {cite}`Duffy1997`. $\nu_0=0.29$.

$c_v=502.416$ J/K/kg from https://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.htm

Need tensile strength; this is a guess.

All values mks.
```
mat1:
    eos:
        name   : 'Steel 304'
        type   : 'MGR'
        rhoref : 7870.0
        c0     : 4580.0
        s1     : 1.49
        s2     : 0.0
        gamma0 : 2.2
        cv     : 502.416
    str:
        type   : 'VM'
        gmod   : 78.0E9
        ys     : 0.2E9
    frac:
        pfrac  : 1.0E9
```

## Quartz

### alpha quartz

$\alpha$-quartz is not well represented by a simple material model. These values provide a coarse representation for teaching.

Wackerle (1962) HEL strongly dependent on orientation. Nominal $\sigma_{HEL}=10$ GPa. 

$\nu=0.08$ https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB014606

$\sigma_{vm}  = \frac{(1-2\nu) \sigma_{HEL}}{(1-\nu)} = 9.3$ GPa.

Shear modulus http://www-odp.tamu.edu/publications/204_SR/103/103_t1.htm
https://link.springer.com/article/10.1007/s00269-014-0711-z

Dynamic tensile strength approximately 1 GPa https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GL100468

<!--- #(1-2.*.085)/(1-.085)*10 --->

```
mat1:
    eos:
        name   : 'Quartz'
        type   : 'MGR'
        rhoref : 2648.0
        c0     : 620.0
        s1     : 1.74
        s2     : 0.0
        gamma0 : 0.84
        cv     : 740.0
    str:
        type   : 'VM'
        gmod   : 44.4E9
        ys     : 9.3E9
    frac:
        pfrac  : 1.0E9
        nrhomin: 0.9
```

### Fused quartz

The following US-up relationship for fused silica shock-compressed to 44–72 GPa was determined using the results from this study: US = 0.62(7) + 1.74(2)uP. {cite}`Berryman2019`

Poisson's ratio is 0.17
https://www.heraeus.com/en/hca/fused_silica_quartz_knowledge_base_1/properties_1/properties_hca.html#tabs-608478-6

Wackerle 1962, FS HEL about 9.5 GPa. Ys = (1-2*.17)/(1-.17)*9.5E9; assume spall is Ys/10.

Gamma0 from SESAME. https://sgp.fas.org/othergov/doe/lanl/lib-www/la-pubs/00296704.pdf

Shear modulus and heat capacity from www.accuratus.com

```
mat1:
    eos:
        name   : 'Fused Silica'
        type   : 'MGR'
        rhoref : 2201.0
        c0     : 620.0
        s1     : 1.74
        s2     : 0.0
        gamma0 : 0.65
        cv     : 740.0
    str:
        type   : 'VM'
        gmod   : 31.0E9
        ys     : 7.55E9
    frac:
        pfrac  : 0.755E9
```


## Forsterite

A tabular EOS is available from https://github.com/ststewart/aneos-forsterite-2019

## Olivine/Dunite

These olivine parameters were developed in {cite}`Marinova2011`.
```
mat1:
    eos:
        name   : 'Olivine'
        type   : 'TIL'
        rhoref : 3500.0
        E0     : 550.0E6
        EIV    : 4.5E6
        ECV    : 14.5E6
        AA     : 131.0E9
        BB     : 49.0E9
        a      : 0.5
        b      : 1.4
        alpha  : 5.0
        beta   : 5.0
```


## Lithium Flouride (LiF)

Hugoniot from {cite}`Hawreliak2023`. Grueneisen parameter from {cite}`Duffy1997`. Heat capacity from Wikipedia.

Poisson's ratio is 0.27.


```
mat1:
    eos:
        name   : 'LiF'
        type   : 'MGR'
        rhoref : 2640.0
        c0     : 5144.0
        s1     : 1.355
        s2     : 0.0
        gamma0 : 1.6
        cv     : 1507.0
    str:
        type   : 'VM'
        gmod   : 55.14E9
        ys     : 0.2E9
    frac:
        pfrac  : 0.02E9
```

https://www.crystran.co.uk/optical-materials/lithium-fluoride-lif

Specific Heat Capacity :	1562 J Kg-1 K-1
Dielectric Constant :	0.1
Youngs Modulus (E) :	64.97 GPa (2)
Shear Modulus (G) :	55.14 GPa (2)
Bulk Modulus (K) :	62.03 GPa (2)

