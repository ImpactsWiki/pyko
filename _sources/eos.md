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

(page:eos)=
# Equations of state

The following equation of state (EOS) models are supported:
* [](eos:idg)
* [](eos:mgr)
* [](eos:ses)

(eos:idg)=
## Ideal gas law

The ideal gas is defined by the initial conditions, the ratio of specific heat capacities, $\gamma$, and the specific heat capacity at constant volume $c_v$. The specific heat capacity is used to determine temperature assuming $T=E/cv$. $cv$ must be a non-zero value. The initial conditions must specify the initial density, pressure, energy and particle velocity.

Example configuration file entry (mks):
```
mat1:
    init:
        up0    : 0.0
        rho0   : 10.0E3
        p0     : 100.0E11
        e0     : 2.5E9
    eos: 
        name   : 'IDG 1'
        type   : 'IDG'
        gamma0 : 1.4
        cv     : 2.5E8
```

(eos:mgr)=
## Mie-Gr&uuml;neisen EOS

The Mie-Gr&uuml;neisen model uses a reference curve, usually the principal Hugoniot, and a thermal parameter to describe a pressure-volume-energy EOS surface (e.g., see Chapter 5 in {cite}`Forbes2012`). 

The material model requires a reference density, Hugoniot parameters, and the Gr&uuml;neisen gamma.
The principal Hugoniot is described by a quadratic equation: $U_S = c_0 + s_1 u_p + s_2 u_p^2$. The specific heat capacity is used to determine temperature assuming $T=E/cv$. $cv$ must be a non-zero value.

The initial conditions must specify the initial density, pressure, energy and particle velocity.

Example configuration file entry (mks):
```
mat1:
    init:
        up0    : 0.0
        rho0   : 8930.0
        p0     : 0.0
        e0     : 0.0
    eos:
        name   : 'Cu test'
        type   : 'MGR'
        rhoref : 8930.0
        c0     : 3900.0
        s1     : 1.49
        s2     : 0.0
        gamma0 : 1.99
        cv     : 1.0E6
```

The thermal Gr&uuml;neisen parameter $\gamma$ is given by

$$ 
\gamma(V) = V \left[\frac{dP}{dE}\right]_V, 
$$

where $\gamma$ is often assumed to only depend on volume. Note that parentheses are used to denote dependent variables.

In pyKO, it is assumed that $\gamma$ is solely dependent on V and

$$
\frac{\gamma(V)}{V} = \frac{\gamma_0}{V_0} = {\rm constant.}
$$

A word of caution: these simple equations for the shock Hugoniot and the Gr&uuml;neisen parameter are unlikely to be thermodynamically self-consistent over wide ranges of volume. Beware using the Mie-Gr&uuml;neisen model over wide ranges of pressure and density. 

The Mie-Gr&uuml;neisen equation of state defines a thermal pressure term, $P(V,E)$, at specific volume $V$ and specific internal energy $E$ referenced to a known state at $V$, which in this case is the principal shock Hugoniot with pressure $P_H(V)$ and specific internal energy $E_H(V)$. Then the thermal pressure term is given by 

$$ 
P(V,E) = P_H(V)+\frac{\gamma(V)}{V} \left[E(V)-E_H(V)\right]. 
$$

A single shock must satisfy the Rankine-Hugoniot conservation equations between an initial state, $s_0$, and the shock state, $s_H$. The initial state $s_0$ is defined by the state variables: specific volume $V_0$, presure $P_0$, and specific energy $E_0$ (which then specifies the initial specific entropy $S_0$ and temperature $T_0$). Note the tabulated initial value for the Gr&uuml;neisen parameter $\gamma_0$ must correspond to state $s_0$. The energy conservation equation is 

$$ 
E_H(V)=E_0+\left[\frac{P_0+P_H(V)}{2}\right]\left[V_0-V\right]. 
$$

Combining these two equations defines the Mie-Gr&uuml;neisen equation of state surface, which is limited to the volume range where $P_H(V)$ is defined,

$$
P(V,E) & =  P_H(V)+\frac{\gamma(V)}{V}\left[E(V)-E_0\right] - \left[P_0+P_H(V)\right]\frac{\gamma(V)}{2V}[V_0-V]  \\
       & =  \frac{\gamma(V)}{V}\left[E(V)-E_0\right] + P_H(V)\left\{ 1- \frac{\gamma(V)}{2V}\left[V_0-V\right] \right\}
        - \frac{P_0\gamma(V)}{2V}\left[V_0-V\right],\label{eq5} (5)
$$

and the last term is dropped when $P_0=0$. 

(eos:ses)=
## Tabular EOS

A tabular EOS contains the thermodynamic variables for the material over a rectangular grid of two independent variables, typically density and temperature. The pressure-volume-energy relationships are interpolated during problem execution to solve the energy conservation equation. pyKO uses bilinear interpolation on the density-temperature grid, but other tabular EOS schemes are easy to implement.

The SESClass.readses function assumes STD and EXT formats compatible with the eos_table.py module from the Stewart group.

* SESAME standard 201+301 format. This table format specficiation is described in {cite:t}`Lyon1992` (<a href="https://github.com/ststewart/aneos-forsterite-2019/blob/master/EOS-docs/Lyon-Johnson-1992-SESAME-database.pdf">PDF file</a>).

The units in the tabular EOS must be specified in the pyKO configuration file (see [](config:units)).

Stewart Group ANEOS tables (https://github.com/ststewart?tab=repositories) have the following format:

* NEW-SESAME-STD.TXT: Standard length SESAME file with 201 table and 301 table (density, temperature, pressure, sp. internal energy, Helmholtz free energy). 301 table units: g/cm<sup>3</sup>, K, GPa, MJ/kg, MJ/kg.

* NEW-SESAME-EXT.TXT: SESAME-style table with extra variables from ANEOS. Contains the standard 201 table and non-standard 301-extra-variables EOS table. The 301 table has: density grid values, temperature grid values, sp. entropy(T,rho), sound speed(T,rho), sp. heat capacity(T,rho), KPA flag(T,rho). 2-D arrays list all densities, looping over each temperature. 301 table units: g/cm<sup>3</sup>, K, MJ/K/kg, cm/s, MJ/K/kg, integer flag. The KPA flag is an ANEOS output with phase information.

To use a tabular EOS, specify the path, std and ext file names and the units for the file (in the units section). To use a different tabular format, customize the readses function to interface with the eos_tables.py module.
```
mat1:
    init:
        up0    : 0.0
        rho0   : 3220.0
        p0     : 1.e5
        e0     : 73277.974035
    eos: 
        name   : 'Forsterite'
        type   : 'SES'
        path    : '../aneos-forsterite-2019/'
        filestd : 'NEW-SESAME-STD.TXT'
        fileext : 'NEW-SESAME-EXT.TXT'
```

```{margin} Table interpolation
Future documentation will discuss pros and cons of different interpolation schemes.
```


(eos:til)=
## Tillotson EOS

```
# parameters: [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta]
# units:    [kg/m3, J/kg, J/kg, J/kg, Pa, Pa, [-]x4]
# these olivine parameters are from Marinova et al. 2011 Icarus 
# dunitetill = [3500.0, 550.0e+6, 4.500e+6, 14.50e+6, 131.00e+9,  49.00e+9, 0.5, 1.4, 5.0, 5.0]
# Basalt parameters from iSALE -- from where? Benz?
# basalttill = [2650.0, 4.87E8, 4.72E6, 18.2E6, 5.3E10, 5.3E10, 0.6, 0.6, 5., 5.]
```

Example configuration file entry (mks):
```
mat1:
    init:
        up0    : 0.0
        rho0   : 3500.0
        p0     : 0.0
        e0     : 0.0
    eos:
        name   : 'Dunite'
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

